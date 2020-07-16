% Code by Yuxuan Xia

function [Bern,lik] = detectionBern(Bern,C,model)

%local hypothesis state update
[Bern.GGIW,lik] = updateGGIW(Bern.GGIW,C,model);
%local hypothesis weight
lik = lik + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));

Bern.r = 1;
Bern.t_death = Bern.t_death(end);
Bern.w_death = 1;

end
