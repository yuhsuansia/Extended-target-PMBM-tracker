function [PPP,MBM] = updatePMBM(PPP,MBM,W,time,model)

m = size(W,2);                      %number of measurements received

used_meas_u = false(m,1);           %measurement indices inside the gate of undetected objects
nu = length(PPP.w);                 %number of mixture components in PPP intensity
gating_matrix_u = false(m,nu);      %gating matrix for PPP components
for i = 1:nu
    %Perform gating for each mixture component in the PPP intensity
    gating_matrix_u(:,i) = ellipsoidalGating(W,PPP.GGIW(i),model);
    used_meas_u = used_meas_u | gating_matrix_u(:,i);
end

n_tt = length(MBM.track);           %number of pre-existing tracks
gating_matrix_d = cell(n_tt,1);     %gating matrix for each track
used_meas_d = false(m,1);           %measurement indices inside the gate of detected objects
for i = 1:n_tt
    %number of hypotheses in track i
    num_hypo = length(MBM.track{i});
    %construct gating matrix
    gating_matrix_d{i} = false(m,num_hypo);
    for j = 1:num_hypo
        %Perform gating for each single object hypothesis
        if MBM.track{i}(j).Bern.t_death(end) == time
            gating_matrix_d{i}(:,j) = ellipsoidalGating(W,MBM.track{i}(j).Bern.GGIW(end),model);
        end
    end
    used_meas_d = used_meas_d | sum(gating_matrix_d{i},2) >= 1;
end

%measurement indices inside the gate
used_meas = used_meas_d | used_meas_u;
%find indices of measurements inside the gate of undetected
%objects but not detected objects
used_meas_u_not_d = used_meas > used_meas_d;

%used measurements by detected targets
W1 = W(:,used_meas_d);
gating_matrix_u1 = gating_matrix_u(used_meas_d,:);
gating_matrix_d = cellfun(@(x) x(used_meas_d,:), gating_matrix_d, 'UniformOutput',false);

%Data association using stochastic optimisation
J = length(MBM.w);
wAssoc = [];
Nj = zeros(J,1);
P = cell(J,1);
for j = 1:J
    if isempty(W1)
        %Misdetection
        track_indices = find(MBM.table(j,:)>0);
        nj = length(track_indices);
        P{j}{1} = cell(nj,1);
        lik = 0;
        for i = 1:nj
            %Create new Bernoulli component due to missed detection
            track_miss = MBM.track{track_indices(i)}(MBM.table(j,track_indices(i)));
            %Compute the misdetection likelihood
            [~,lik_miss] = misdetectionBern(track_miss.Bern,model);
            lik = lik + lik_miss;
        end
    else
        %Find the most likely global hypotheses using Stochastic
        %Optimization
        [P{j},lik] = ObjectsSO(MBM,j,PPP,W1,gating_matrix_d,gating_matrix_u1,model);
    end
    Nj(j) = length(lik);
    wAssoc = [wAssoc;lik+MBM.w(j)];
end
%Normalise and prune components with low weights
[wAssoc,~] = normalizeLogWeights(wAssoc);
idx_keep = wAssoc > log(model.threshold_w);
wAssoc = wAssoc(idx_keep);
if length(wAssoc) == 1; wAssoc = 0; end

%Remove measurement partitions that correspond to low weight global hypotheses
true_meas_indices = find(used_meas_d==1);
m = size(W,2);
idx_0 = 0;
for j = 1:J
    idx_j = idx_keep(idx_0+1:idx_0+Nj(j));
    P{j} = P{j}(idx_j);
    %Convert set representation to boolean vector representation
    if ~isempty(P{j})
        for i = 1:length(P{j})
            for k = 1:length(P{j}{i})
                P{j}{i}{k} = ismember(1:m,true_meas_indices(P{j}{i}{k}));
            end
        end
    end
    idx_0 = idx_0 + Nj(j);
end

%For each single target hypothesis find its corresponding measurement cells
meas_track = cell(n_tt,1);
for i = 1:n_tt
    meas_track{i} = cell(length(MBM.track{i}),1);
end

meas_newtracks = false(1,m);
idx_new = 0;
for j = 1:J
    if ~isempty(P{j})
        track_indices = find(MBM.table(j,:)>0);
        nj = length(track_indices);
        for i = 1:length(P{j})
            for k = 1:nj
                if isempty(meas_track{track_indices(k)}{MBM.table(j,track_indices(k))})
                    meas_track{track_indices(k)}{MBM.table(j,track_indices(k))} = ...
                        [meas_track{track_indices(k)}{MBM.table(j,track_indices(k))};P{j}{i}{k}];
                elseif ~ismember(P{j}{i}{k},meas_track{track_indices(k)}{MBM.table(j,track_indices(k))},'rows')
                    meas_track{track_indices(k)}{MBM.table(j,track_indices(k))} = ...
                        [meas_track{track_indices(k)}{MBM.table(j,track_indices(k))};P{j}{i}{k}];
                end
            end
            num_mea_cell = length(P{j}{i});
            if num_mea_cell > nj
                for k = 1:num_mea_cell-nj
                    if ~ismember(P{j}{i}{k+nj},meas_newtracks,'rows')
                        idx_new = idx_new + 1;
                        meas_track{n_tt+idx_new,1}{1} = P{j}{i}{k+nj};
                        meas_newtracks = [meas_newtracks;P{j}{i}{k+nj}];
                    end
                end
            end
        end
    end
end

%Make sure boolean representation
for i = 1:n_tt
    for j = 1:length(meas_track{i})
        meas_track{i}{j} = logical(meas_track{i}{j});
    end
end
meas_newtracks = logical(meas_newtracks);

%Construct new global hypotheses look-up table
n_tt_upd = length(meas_track);
table = zeros(length(wAssoc),n_tt_upd);
idx = 0;
for j = 1:J
    if ~isempty(P{j})
        track_indices = find(MBM.table(j,:)>0);
        nj = length(track_indices);
        for i = 1:length(P{j})
            idx = idx + 1;
            for k = 1:nj
                [~,table(idx,track_indices(k))] = ...
                    ismember(P{j}{i}{k},meas_track{track_indices(k)}{MBM.table(j,track_indices(k))},'rows');
                if MBM.table(j,track_indices(k)) > 1
                    for p = 1:MBM.table(j,track_indices(k))-1
                        table(idx,track_indices(k)) = table(idx,track_indices(k)) + size(meas_track{track_indices(k)}{p},1);
                    end
                end
            end
            num_mea_cell = length(P{j}{i});
            if num_mea_cell > nj
                for k = 1:num_mea_cell-nj
                    [~,table_idx] = ismember(P{j}{i}{k+nj},meas_newtracks,'rows');
                    table(idx,table_idx-1+n_tt) = 1;
                end
            end
        end
    end
end

indices = 1:size(W,2);

%Update tracks
tracks = cell(n_tt_upd,1);
for i = 1:n_tt
    idx = 0;
    for j = 1:length(meas_track{i})
        for k = 1:size(meas_track{i}{j},1)
            idx = idx + 1;
            tracks{i}(idx,1) = MBM.track{i}(j); 
            if any(meas_track{i}{j}(k,:))
                [Bern,lik] = detectionBern(MBM.track{i}(j).Bern,W(:,meas_track{i}{j}(k,:)),model);
            else
                [Bern,lik] = misdetectionBern(MBM.track{i}(j).Bern,model);
            end
            tracks{i}(idx,1).Bern = Bern;
            tracks{i}(idx,1).lik = tracks{i}(idx,1).lik + lik;
            len = length(tracks{i}(idx,1).assocHistory);
            tracks{i}(idx,1).assocHistory(len+1,1).t = time;
            tracks{i}(idx,1).assocHistory(len+1,1).meas = indices(meas_track{i}{j}(k,:));
        end
    end
end
if n_tt_upd > n_tt
    for i = 1:n_tt_upd - n_tt
        in_gate = sum(gating_matrix_u(meas_track{n_tt+i}{1},:),1)>=1;
        if any(in_gate)
            [Bern,lik] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_track{n_tt+i}{1}),model);
            Bern.t_birth = time;
            Bern.t_death = time;
            Bern.w_death = 1;
        else
            Bern.r = 0; 
            Bern.GGIW = struct('a',0,'b',1,'v',0,'P',zeros(2,2),'V',zeros(2,2)); 
            lik = [];
        end
        tracks{n_tt+i,1} = struct('Bern',Bern,'lik',lik,'assocHistory',[]);
        tracks{n_tt+i,1}.assocHistory(1).t = time;
        tracks{n_tt+i,1}.assocHistory(1).meas = indices(meas_track{n_tt+i}{1});
    end
end

%Append new tracks
[tracks,table,wAssoc] = newObjectsSO(tracks,table,wAssoc,PPP,W,gating_matrix_u,used_meas_u_not_d,time,model);

%Remove Bernoulli components with low existence probability
n_tt = length(tracks);
for i = 1:n_tt
    %Find all Bernoulli components needed to be pruned
    idx = arrayfun(@(x) x.Bern.r < model.threshold_r, tracks{i});
%     idx = arrayfun(@(x) x.Bern.r < model.threshold_r | ...
%         (1-(x.Bern.GGIW(end).b/(x.Bern.GGIW(end).b+1))^x.Bern.GGIW(end).a) < 0.5 | ...
%         x.Bern.GGIW(end).v < 6 | x.Bern.GGIW(end).P(1,1) > 20^2 | x.Bern.GGIW(end).P(2,2) > 20^2 |...
%         x.Bern.GGIW(end).V(1,1)/(x.Bern.GGIW(end).v-6) > 20 | x.Bern.GGIW(end).V(2,2)/(x.Bern.GGIW(end).v-6) > 20, tracks{i});
    %Prune these Bernoulli components
    tracks{i} = tracks{i}(~idx);
    idx = find(idx);
    %Update hypothesis table, if a Bernoulli component is
    %pruned, set its corresponding entry to zero
    for j = 1:length(idx)
        temp = table(:,i);
        temp(temp==idx(j)) = 0;
        table(:,i) = temp;
    end
end

%Remove unused tracks
idx_empty = cellfun('isempty',tracks);
table = table(:,~idx_empty);
tracks = tracks(~idx_empty);

%Remove tracks that contains only null single object hypotheses
idx_keep = sum(table,1) > 0;
table = table(:,idx_keep);
tracks = tracks(idx_keep);
if isempty(table)
    wAssoc = [];
end

%Re-index hypothesis table
n_tt = length(tracks);
for i = 1:n_tt
    idx = table(:,i) > 0;
    [~,~,table(idx,i)] = unique(table(idx,i),'stable');
%     if any(idx)
%         table(idx,i) = 0;
%         table(~idx,i) = table(~idx,i) - 1;
%     end
end

%Merge duplicate hypothesis table rows
if ~isempty(table)
    [ht,~,ic] = unique(table,'rows','stable');
    if(size(ht,1)~=size(table,1))
        %There are duplicate entries
        w = zeros(size(ht,1),1);
        for i = 1:size(ht,1)
            indices_dupli = (ic==i);
            [~,w(i)] = normalizeLogWeights(wAssoc(indices_dupli));
        end
        table = ht;
        wAssoc = w;
    end
end

%%%%%%%%%%%%
%Merge similar Bernoulli components in the track
% 
% n_tt = length(tracks);
% for i = 1:n_tt
%     nb = length(tracks{i});
%     if nb > 1
%         %find the highest weight MB that contains this track
%         idx_mb = table(:,i) > 0;
%         idx_mb = find(idx_mb);
%         [~,idx] = max(wAssoc(idx_mb));
%         idx_b = table(idx_mb(idx),i);
%         I = idx_b;
%         for j = 1:nb
%             if j ~= idx_b && tracks{i}(idx_b).Bern.r == tracks{i}(j).Bern.r...
%                     && GGIW_KLdiv(tracks{i}(idx_b).Bern.GGIW,tracks{i}(j).Bern.GGIW) < model.merge
%                 I = [I j];
%             end
%         end
%         len = length(I);
%         if len > 1
%             Bern_temp = [tracks{i}(I).Bern];
%             w_temp = zeros(len,1);
%             for l = 1:len
%                 [~,w_temp(l)] = normalizeLogWeights(wAssoc(table(:,i)==I(l)));
%             end
%             %Assume the same weight
%             [~,GGIW_hat] = GGIW_merge_wrap(w_temp,[Bern_temp.GGIW]);
%             tracks{i}(idx_b).Bern.GGIW = GGIW_hat;
%             idx_remove = setdiff(I,idx_b);
%             tracks{i}(idx_remove) = [];
%             table(ismember(table(:,i),I),i) = idx_b;
%         end
%     end
% end
% 
% %Re-index hypothesis table
% n_tt = length(tracks);
% for i = 1:n_tt
%     idx = table(:,i)==0;
%     [~,~,table(:,i)] = unique(table(:,i),'stable');
%     if any(idx)
%         table(idx,i) = 0;
%         table(~idx,i) = table(~idx,i) - 1;
%     end
% end
% 
% %Merge duplicate hypothesis table rows
% if ~isempty(table)
%     [ht,~,ic] = unique(table,'rows','stable');
%     if(size(ht,1)~=size(table,1))
%         %There are duplicate entries
%         w = zeros(size(ht,1),1);
%         for i = 1:size(ht,1)
%             indices_dupli = (ic==i);
%             [~,w(i)] = normalizeLogWeights(wAssoc(indices_dupli));
%         end
%         table = ht;
%         wAssoc = w;
%     end
% end
% 
% if length(wAssoc) == 1; wAssoc = 0; end

%%%%%%%%%%%%

%Assign updated value
MBM.table = table;
MBM.track = tracks;
MBM.w = wAssoc;

%PPP misdetection update
PPP = misdetectionPPP(PPP,model);
%Prune PPP components
idx_keep = PPP.w > log(model.threshold_u);
PPP.w = PPP.w(idx_keep);
PPP.GGIW = PPP.GGIW(idx_keep);


end

