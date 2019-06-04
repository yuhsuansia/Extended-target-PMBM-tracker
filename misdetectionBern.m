function [Bern,lik] = misdetectionBern(Bern,model)

temp = Bern.w_death(end)*model.Pd*(Bern.GGIW(end).b/(Bern.GGIW(end).b+1))^Bern.GGIW(end).a;
model.Qd = 1-(1-model.Qd)*Bern.w_death(end);

qD = model.Qd + temp;

w1 = model.Qd/qD;
w2 = temp/qD;

GGIW1 = Bern.GGIW(end);
GGIW2 = GGIW1;
GGIW2.b = GGIW2.b + 1;

lik = 1 - Bern.r + Bern.r*qD;

[~,Bern.GGIW(end).a,Bern.GGIW(end).b] = gammaMerge([w1;w2],[GGIW1.a;GGIW2.a],[GGIW1.b;GGIW2.b]);

Bern.r = Bern.r*qD/lik;

lik = log(lik);

%Updated time of death
Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)*(qD)]/(1-Bern.w_death(end)*(1-qD));

end
