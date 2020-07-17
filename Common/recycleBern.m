% Code by Yuxuan Xia

function [PPP,MBM] = recycleBern(PPP,MBM,model)

%Recycle Bernoullis with low existence probability and add them to PPP

n_tt = length(MBM.track);
for i = 1:n_tt
    idx = arrayfun(@(x) x.Bern.r<model.recycle, MBM.track{i});
    if idx
        %Here, we should also consider the weights of different MBs
        idx_t = find(idx);
        n_h = length(idx_t);
        w_h = zeros(n_h,1);
        for j = 1:n_h
            idx_h = MBM.table(:,i) == idx_t(j);
            [~,w_h(j)] = normalizeLogWeights(MBM.w(idx_h));
        end
        %Recycle
        temp = [MBM.track{i}(idx).Bern];
        PPP.w = [PPP.w;log([temp.r]')+w_h];
        PPP.GGIW = [PPP.GGIW;[temp.GGIW]'];
        %Remove Bernoullis
        MBM.track{i} = MBM.track{i}(~idx);
        %Update hypothesis table, if a Bernoulli component is
        %pruned, set its corresponding entry to zero
        idxx = find(idx);
        for j = 1:length(idxx)
            temp = MBM.table(:,i);
            temp(temp==idxx(j)) = 0;
            MBM.table(:,i) = temp;
        end
    end
end

%Remove unused tracks
idx_empty = cellfun('isempty',MBM.track);
MBM.table = MBM.table(:,~idx_empty);
MBM.track = MBM.track(~idx_empty);

%Remove tracks that contains only null single object hypotheses
idx_keep = sum(MBM.table,1) > 0;
MBM.table = MBM.table(:,idx_keep);
MBM.track = MBM.track(idx_keep);
if isempty(MBM.table)
    MBM.w = [];
end

%Re-index hypothesis table
n_tt = length(MBM.track);
for i = 1:n_tt
    idx = MBM.table(:,i)==0;
    [~,~,MBM.table(:,i)] = unique(MBM.table(:,i),'rows');
    if any(idx)
        MBM.table(idx,i) = 0;
        MBM.table(~idx,i) = MBM.table(~idx,i) - 1;
    end
end

%Merge duplicate hypothesis table rows
if ~isempty(MBM.table)
    [ht,~,ic] = unique(MBM.table,'rows');
    if(size(ht,1)~=size(MBM.table,1))
        %There are duplicate entries
        w = zeros(size(ht,1),1);
        for i = 1:size(ht,1)
            indices_dupli = (ic==i);
            [~,w(i)] = normalizeLogWeights(MBM.w(indices_dupli));
        end
        MBM.table = ht;
        MBM.w = w;
    end
end

%Merge similar components in PPP
if length(PPP.w) > 1
    [PPP.w,PPP.GGIW] = mixtureReduction(PPP.w,PPP.GGIW,model);
end

end

