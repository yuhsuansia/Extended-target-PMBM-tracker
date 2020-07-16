% Code by Yuxuan Xia

function [GGIW] = predictGGIW(GGIW,model)

% (heuristic) Gamma prediction
GGIW_.a = GGIW.a/model.eta;
GGIW_.b = GGIW.b/model.eta;

% (Kalman) Gaussian prediction
GGIW_.m = model.motionmodel.f(GGIW.m);
F = model.motionmodel.F(GGIW.m);
GGIW_.P = F*GGIW.P*F' + model.motionmodel.Q;

% (heuristic) Inverse-Wishart prediction
d = 2;
M = eye(2);
e = exp(-model.Ts/model.tao);
GGIW_.v = 2*d + 2 + e*(GGIW.v - 2*d - 2);
GGIW_.V = e*M*GGIW.V*M';

GGIW = GGIW_;

end

