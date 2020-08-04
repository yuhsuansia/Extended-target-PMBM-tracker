function [tracks,table,wAssoc] = trackOrientedPMB(tracks,table,wAssoc)

%number of tracks
n_tt = length(tracks);
if length(wAssoc) == 1 || n_tt == 0
    return
else
    %merging local hypotheses in the same track
    for i = 1:n_tt
        %number of local hypotheses
        nb = length(tracks{i});
        %for each local hypothesis, find global hypotheses that contain it
        w_merge = zeros(nb,1);
        GGIW_to_merge = [];
        for j = 1:nb
            %add global hypotheses weights
            [~,w_merge(j)] = normalizeLogWeights(wAssoc(table(:,i) == j));
            %multiply with the existence probability (summation in logarithm)
            w_merge(j) = w_merge(j) + log(tracks{i}(j).Bern.r);
            GGIW_to_merge = [GGIW_to_merge;tracks{i}(j).Bern.GGIW];
        end
        %merge all GGIWs and existence probability
        [r_hat,GGIW_hat] = GGIW_merge_wrap(w_merge,GGIW_to_merge);
        tracks{i}(1).Bern.r = r_hat;
        tracks{i}(1).Bern.GGIW = GGIW_hat;
        tracks{i}(2:end) = [];
    end
    
    table = ones(1,n_tt);
    wAssoc = 0;
    
end

end

