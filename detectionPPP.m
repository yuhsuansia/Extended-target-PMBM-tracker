function [Bern,lik] = detectionPPP(wp,GGIWp,C,model)

[GGIW_c,lik_c] = updateGGIWforPPP(GGIWp,C,model);

w_c = lik_c + wp + log(model.Pd);
[w_hat,Bern.GGIW] = GGIW_merge_wrap(w_c,GGIW_c);

num_meas = size(C,2);

if num_meas > 1
    Bern.r = 1;
    lik = log(w_hat);
else
    Bern.r = w_hat/(w_hat+model.lambda_fa);
    lik = log(w_hat+model.lambda_fa);
end

% lik = log(w_hat + model.lambda_fa^num_meas);


end

