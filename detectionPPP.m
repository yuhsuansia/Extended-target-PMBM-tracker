% Code by Yuxuan Xia

function [Bern,lik] = detectionPPP(wp,GGIWp,C,model)

%GGIW state update given measurement set
[GGIW_c,lik_c] = updateGGIWforPPP(GGIWp,C,model);

%Weight of new local hypothesis
w_c = lik_c + wp + log(model.Pd);

%Merge different components created by different PPP components
[w_hat,Bern.GGIW] = GGIW_merge_wrap(w_c,GGIW_c);

num_meas = size(C,2);

%For extented targets, if more than one measurement is assigned, the
%measurement set cannot correspond to clutter
if num_meas > 1
    Bern.r = 1;
    lik = log(w_hat);
else
    %For potential new target created using single measurement, take the
    %probability being clutter into account
    Bern.r = w_hat/(w_hat+model.lambda_fa);
    lik = log(w_hat+model.lambda_fa);
end


end

