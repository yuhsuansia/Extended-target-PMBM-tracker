function P = genPartitions_dbscan(z,model)
% generate partitions of given measurement set

max_dist = model.max_dist;
min_dist = model.min_dist;
grid_dist = model.grid_dist;

if isempty(z)
    P = cell(0,1);
    return
elseif size(z,2) <= 1
    P = cell(1);
    P{1}{1} = 1:size(z,2);
    return
end

N_p = ceil((max_dist-min_dist)/grid_dist);
dist = linspace(min_dist,max_dist,N_p);
P = cell(N_p,1);

for i = 1:N_p
    p = dbscan(z, dist(i) ,1);
    % select unique partitions
    for j = 1:i
        if isequal(p,P{j})
            p = [];
            break;
        end
    end
    P{i} = p;
end

P = P(~cellfun(@isempty,P));

end