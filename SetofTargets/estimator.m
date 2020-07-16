% Code by Yuxuan Xia

function estimates = estimator(MBM,model)


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
        estimates.g = [estimates.g MBM.track{i}(table_entry(i)).Bern.GGIW.a/MBM.track{i}(table_entry(i)).Bern.GGIW.b];
        estimates.x = [estimates.x MBM.track{i}(table_entry(i)).Bern.GGIW.m];
        X = MBM.track{i}(table_entry(i)).Bern.GGIW.V/(MBM.track{i}(table_entry(i)).Bern.GGIW.v-2*d-2);
        estimates.X = cat(3,estimates.X,X);
    end
end


end

