% Code by Yuxuan Xia

function [estimates,trajectoryEstimates] = estimator(MBM,model)

%Initialise memory for trajectory estimates (both alive and dead)
trajectoryEstimates = [];

%Initialise memory for current state estimates
estimates.g = zeros(1,0);
estimates.x = zeros(4,0);
estimates.X = zeros(2,2,0);

%Target extent dimension
d = 2;

%Extract estimates from the highest weight global hypothesis (multi-Bernoulli)
[~,idx] = max(MBM.w);
table_entry = MBM.table(idx,:);
for i = 1:length(table_entry)
    %Only extract estimates from Bernoullis with existence probability
    %larger than a pre-defined threshold
    if table_entry(i) > 0 && MBM.track{i}(table_entry(i)).Bern.r >= model.exist_r
        %Select the most probable trajectory end time
        [~,ideath] = max(MBM.track{i}(table_entry(i)).Bern.w_death);
        %If the selected end time is equal to the maximum possible length,
        %this trajectory is assumed to be alive at current time step
        if ideath == length(MBM.track{i}(table_entry(i)).Bern.w_death)
            estimates.g = [estimates.g MBM.track{i}(table_entry(i)).Bern.GGIW(end).a/MBM.track{i}(table_entry(i)).Bern.GGIW(end).b];
            estimates.x = [estimates.x MBM.track{i}(table_entry(i)).Bern.GGIW(end).m];
            X = MBM.track{i}(table_entry(i)).Bern.GGIW(end).V/(MBM.track{i}(table_entry(i)).Bern.GGIW(end).v-2*d-2);
            estimates.X = cat(3,estimates.X,X);
        end
        
        %Trajectory death time
        t_death = MBM.track{i}(table_entry(i)).Bern.t_death(ideath);
        %Trajectory length
        tlen = t_death - MBM.track{i}(table_entry(i)).Bern.t_birth + 1;
        
        trajectoryEstimates(end+1,1).t_birth = [MBM.track{i}(table_entry(i)).assocHistory(1).t];
        trajectoryEstimates(end,1).t_death = t_death;
        
        trajectoryEstimates(end,1).g = [MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).a]./...
            [MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).b];
        
        trajectoryEstimates(end,1).x = [MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).m];
        
        temp = arrayfun(@(x) x.V/(x.v-2*d-2), MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen), 'un', false);
        trajectoryEstimates(end,1).X = zeros(d,d,length(temp));
        for j = 1:length(temp)
            trajectoryEstimates(end,1).X(:,:,j) = temp{j};
        end
    end
end  
    
    
end

