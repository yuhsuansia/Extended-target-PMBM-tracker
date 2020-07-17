% Code by Yuxuan Xia

function [tracks,table,wAssoc] = newObjectsCA(tracks,table,wAssoc,PPP,W,gating_matrix,used_meas_u_not_d,k,model)

%Perform data partitioning to find measurement partitions with high
%likelihood for PPP components. Note that this can be directly applied to
%extended target PHD filter using data partitioning

indices = 1:size(W,2);

W = W(:,used_meas_u_not_d);
gating_matrix = gating_matrix(used_meas_u_not_d,:);

%Perform DBSCAN with different hyperparameter settings to obtain multiple
%measurement partitionings
partitions = genPartitions_dbscan(W,model);

np = length(partitions);        %number of different measurement partitions
lik = zeros(np,1);

%Loop through each measurement partition
for j = 1:np
    nc = length(partitions{j}); %number of measurement cells
    Lc = zeros(nc,1);
    %Perform PPP update for each measurement cell in each partition
    for i = 1:nc
        in_gate = sum(gating_matrix(partitions{j}{i},:),1)>=1;
        [~,Lc(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,partitions{j}{i}),model);
    end
    lik(j) = sum(Lc); 
end

%Find the best partition with the highest likelihood
[~,idx] = max(lik);
bestPartition = partitions{idx};

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

