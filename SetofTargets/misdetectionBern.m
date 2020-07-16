% Code by Yuxuan Xia

function [Bern,lik] = misdetectionBern(Bern,model)

% For extended target tracking, misdetection can either due to the fact
% that the target is not detected or the fact that the target is detected
% but generates 0 measurements. We should take both alternatives into
% account.

%Probability that the target is alive and detected but generates 0 measurement
temp = model.Pd*(Bern.GGIW.b/(Bern.GGIW.b+1))^Bern.GGIW.a;
%Porbability that the target is alive but not detected
% model.Qd;

%Total misdetection probability
qD = model.Qd + temp;

%Normalised weights of different events that cause misdetection
w1 = model.Qd/qD;
w2 = temp/qD;

%No change for trajectory is not detected
GGIW1 = Bern.GGIW;

%Update Gamma parameters for target generates 0 measurement
GGIW2 = GGIW1;
GGIW2.b = GGIW2.b + 1;

%Misdetection likelihood
lik = 1 - Bern.r + Bern.r*qD;

%Merge these two misdetection hypotheses
[~,Bern.GGIW.a,Bern.GGIW.b] = gammaMerge([w1;w2],[GGIW1.a;GGIW2.a],[GGIW1.b;GGIW2.b]);

%Updated Bernoulli existence probability
Bern.r = Bern.r*qD/lik;

lik = log(lik);

end
