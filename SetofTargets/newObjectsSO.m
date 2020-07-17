% Code by Yuxuan Xia

function [tracks,table,wAssoc] = newObjectsSO(tracks,table,wAssoc,PPP,W,gating_matrix,used_meas_u_not_d,k,model)

%Perform stochastic optimisation to find measurement partitions with high
%likelihood for PPP components. Note that this can be directly applied to
%extended target PHD filter using stochastic optimisation

indices = 1:size(W,2);

W = W(:,used_meas_u_not_d);
gating_matrix = gating_matrix(used_meas_u_not_d,:);

m = size(W,2);                  %number of measurements
T = model.num_iterations*m;     %number of iterations

%Initialisation, each measurement corresponds to a new object
meas_cell = cell(m,1);
Lc = zeros(m,1);    %likelihood for each cell
for i = 1:m
    meas_cell{i} = i;
    [~,Lc(i)] = detectionPPP(PPP.w(gating_matrix(i,:)),...
        PPP.GGIW(gating_matrix(i,:)),W(:,i),model);
end
bestPartition = meas_cell;
lik = sum(Lc);
Lp = Lc;

max_repetition = max(20,ceil(model.max_repetition*m/2));
num_repetition = 0;
%Stochastic optimisation
for t = 1:T
    mea_selected = randi(m,1);                              %randomly select a measurement
    Loc = cellfun(@(x) x==mea_selected,meas_cell,'Un',0);
    cell_idx = find(cellfun(@(x) any(x(:)),Loc));           %find the corresponding cell index
    N = length(meas_cell);                                  %find the number of cells
    Wp = -inf(2*N+2,1);
    
    if length(meas_cell{cell_idx})==1 %selected cell has only one measurement
        L_move1 = -inf(N,1);
        for i = 1:N
            if i == cell_idx
                Wp(i) = 0;
            elseif any(sum(gating_matrix(meas_cell{i},:),1)>=1 & gating_matrix(mea_selected,:))
                in_gate = sum(gating_matrix([meas_cell{i} mea_selected],:),1)>=1;
                if any(in_gate)
                    [~,L_move1(i)] = detectionPPP(PPP.w(in_gate),...
                        PPP.GGIW(in_gate),W(:,[meas_cell{i} mea_selected]),model);
                    Wp(i) = L_move1(i) - Lp(i) - Lp(cell_idx);
                end
            end
        end
    else %selected cell has more than one measurement
        temp = meas_cell{cell_idx}(meas_cell{cell_idx}~=mea_selected);
        %selected measurement moved to an existing cell
        L_move1 = -inf(N,1);
        for i = 1:N
            if i == cell_idx
                Wp(i) = 0;
            elseif any(sum(gating_matrix(meas_cell{i},:),1)>=1 & gating_matrix(mea_selected,:))
                if length(temp) == 1
                    L_move2 = Lc(temp);
                else
                    in_gate = sum(gating_matrix(temp,:),1)>=1;
                    if any(in_gate)
                        [~,L_move2] = detectionPPP(PPP.w(in_gate),...
                            PPP.GGIW(in_gate),W(:,temp),model);
                    else
                        L_move2 = -inf;
                    end
                end
                in_gate = sum(gating_matrix([meas_cell{i} mea_selected],:),1)>=1;
                if any(in_gate)
                    [~,L_move1(i)] = detectionPPP(PPP.w(in_gate),...
                        PPP.GGIW(in_gate),W(:,[meas_cell{i} mea_selected]),model);
                    Wp(i) = L_move2 + L_move1(i) - Lp(i) - Lp(cell_idx);
                end
            end
        end
        %selected measurement moved to a new cell
        if length(temp) == 1
            Lp_move_new_cell = Lc(temp);
        else
            in_gate = sum(gating_matrix(temp,:),1)>=1;
            if any(in_gate)
                [~,Lp_move_new_cell] = detectionPPP(PPP.w(in_gate),...
                    PPP.GGIW(in_gate),W(:,temp),model);
            else
                Lp_move_new_cell = -inf;
            end
        end
        Wp(N+1) = Lp_move_new_cell + Lc(mea_selected) - Lp(cell_idx);
        %selected cell merged with an existing cell
        Lp_merge = -inf(N,1);
        for i = 1:N
            if i~=cell_idx
                in_gate = sum(gating_matrix([meas_cell{i} meas_cell{cell_idx}],:),1)>=1;
                if any(in_gate)
                    [~,Lp_merge(i)] = detectionPPP(PPP.w(in_gate),...
                        PPP.GGIW(in_gate),W(:,[meas_cell{i} meas_cell{cell_idx}]),model);
                    Wp(i+N+1) = Lp_merge(i) - Lp(i) - Lp(cell_idx);
                end
            end
        end
        %selected cell splited into two parts
        [IDX,~] = kmeanspp(W(:,meas_cell{cell_idx}),2);
        if length(meas_cell{cell_idx}(IDX==1)) == 1
            Lp_split1 = -inf;
        else
            in_gate = sum(gating_matrix(meas_cell{cell_idx}(IDX==1),:),1)>=1;
            if any(in_gate)
                [~,Lp_split1] = detectionPPP(PPP.w(in_gate),...
                    PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==1)),model);
            else
                Lp_split1 = -inf;
            end
        end
        if length(meas_cell{cell_idx}(IDX==2)) == 1
            Lp_split2 = -inf;
        else
            in_gate = sum(gating_matrix(meas_cell{cell_idx}(IDX==2),:),1)>=1;
            if any(in_gate)
                [~,Lp_split2] = detectionPPP(PPP.w(in_gate),...
                    PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==2)),model);
            else
                Lp_split2 = -inf;
            end
        end
        Wp(2*N+2) = Lp_split1 + Lp_split2 - Lp(cell_idx);
    end
    
    [Wp,~] = normalizeLogWeights(Wp);
    I = find(cumsum(exp(Wp))>rand(),1); %sampling
    if I == cell_idx %do nothing
        num_repetition = num_repetition + 1;
        if (num_repetition > max_repetition)
            break;
        end
    elseif I <= N %move selected measurement to an existing cell
        meas_cell{I} = [meas_cell{I} mea_selected];
        Lp(I) = L_move1(I);
        if length(meas_cell{cell_idx})==1
            meas_cell(cell_idx) = [];
            Lp(cell_idx) = [];
        else
            meas_cell{cell_idx} = temp;
            Lp(cell_idx) = L_move2;
        end
    elseif I == N+1 %move selected measurement to a new cell
        meas_cell{I} = mea_selected;
        Lp(I) = Lc(mea_selected);
        meas_cell{cell_idx} = temp;
        Lp(cell_idx) = Lp_move_new_cell;
    elseif I >= N+2 && I <= 2*N+1 %merge selected cell with an existing cell
        meas_cell{I-N-1} = [meas_cell{I-N-1} meas_cell{cell_idx}];
        Lp(I-N-1) = Lp_merge(I-N-1);
        meas_cell(cell_idx) = [];
        Lp(cell_idx) = [];
    elseif I == 2*N+2 %split selected cell
        meas_cell{N+1} = meas_cell{cell_idx}(IDX==2);
        Lp(N+1) = Lp_split2;
        meas_cell{cell_idx} = meas_cell{cell_idx}(IDX==1);
        Lp(cell_idx) = Lp_split1;
    end
    
    if I~=cell_idx
        num_repetition = 0;
    end
    
    if sum(Lp) > lik
        lik = sum(Lp); 
        bestPartition = meas_cell; 
    end
end

%Create new tracks and hypothesis look-up table

n_tt = length(bestPartition);
true_meas_indices = find(used_meas_u_not_d==1);

n = length(tracks);
for i = 1:n_tt
    in_gate = sum(gating_matrix(bestPartition{i},:),1)>=1;
    [Bern,lik] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,bestPartition{i}),model);
    tracks{n+i,1} = struct('Bern',Bern,'lik',lik,'assocHistory',[]);
    tracks{n+i,1}.assocHistory(1).t = k;
    tracks{n+i,1}.assocHistory(1).meas = indices(ismember(1:length(indices),true_meas_indices(sort(bestPartition{i}))));
end
if isempty(wAssoc)
    wAssoc = 0;
    table = ones(1,n_tt);
else
    table = [table ones(size(table,1),n_tt)];
end

end

