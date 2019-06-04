function [PPP,MBM] = predictPMBM(PPP,MBM,model)

% Predict existing PPP
PPP.w = PPP.w + log(model.Ps);
PPP.GGIW = arrayfun(@(x) predictGGIWPPP(x,model), PPP.GGIW);

% Incorporate PPP birth
PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW];

% Predict MBM
n_track = length(MBM.track);
for i = 1:n_track
    nh = length(MBM.track{i});
    for h = 1:nh
        if MBM.track{i}(h).Bern.w_death(end) >= model.threshold_s
            MBM.track{i}(h).Bern.GGIW = predictGGIW(MBM.track{i}(h).Bern.GGIW,model);
            MBM.track{i}(h).Bern.t_death = [MBM.track{i}(h).Bern.t_death MBM.track{i}(h).Bern.t_death(end)+1];
            MBM.track{i}(h).Bern.w_death = [MBM.track{i}(h).Bern.w_death(1:end-1) MBM.track{i}(h).Bern.w_death(end)*(1-model.Ps) MBM.track{i}(h).Bern.w_death(end)*model.Ps];
        end
    end
end

end

