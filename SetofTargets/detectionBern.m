% Code by Yuxuan Xia

function [Bern,lik] = detectionBern(Bern,C,model)

%local hypothesis state update
[Bern.GGIW,lik] = updateGGIW(Bern.GGIW,C,model);
%local hypothesis weight
lik = lik + log(Bern.r) + log(model.Pd);

Bern.r = 1;

end
