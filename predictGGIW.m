function [GGIW] = predictGGIW(GGIW,model)

GGIW_.a = GGIW(end).a/model.eta;
GGIW_.b = GGIW(end).b/model.eta;

GGIW_.m = model.motionmodel.f(GGIW(end).m);
F = model.motionmodel.F(GGIW(end).m);
GGIW_.P = F*GGIW(end).P*F' + model.motionmodel.Q;

d = 2;
M = eye(2);
e = exp(-model.Ts/model.tao);
GGIW_.v = 2*d + 2 + e*(GGIW(end).v - 2*d - 2);
GGIW_.V = e*M*GGIW(end).V*M';

GGIW = [GGIW;GGIW_];

end

