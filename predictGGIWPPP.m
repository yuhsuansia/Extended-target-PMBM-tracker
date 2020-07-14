% Code by Yuxuan Xia

function [GGIW] = predictGGIWPPP(GGIW,model)

% (heuristic) Gamma prediction
GGIW.a = GGIW.a/model.eta;
GGIW.b = GGIW.b/model.eta;

% (Kalman) Gaussian prediction
GGIW.m = model.motionmodel.f(GGIW.m);
F = model.motionmodel.F(GGIW.m);
GGIW.P = F*GGIW.P*F' + model.motionmodel.Q;

% (heuristic) Inverse-Wishart prediction
d = 2;
M = eye(2);
e = exp(-model.Ts/model.tao);
GGIW.v = 2*d + 2 + e*(GGIW.v - 2*d - 2);
GGIW.V = e*M*GGIW.V*M';

end

