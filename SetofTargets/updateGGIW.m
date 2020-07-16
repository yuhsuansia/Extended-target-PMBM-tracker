% Code by Yuxuan Xia

function [GGIW,lik] = updateGGIW(GGIW,W,model)

%GGIW update for a trajectory sequence of GGIWs
%Only update the latest state

d = 2;

card_W = size(W,2);

GGIW_.a = GGIW.a + card_W;
GGIW_.b = GGIW.b + 1;

z_bar = mean(W,2);
epsilon = z_bar - model.measmodel.h(GGIW.m);
H = model.measmodel.H(GGIW.m);

X_hat = GGIW.V/(GGIW.v - 2*d - 2);
X_hat = (X_hat + X_hat')/2;

S = H*GGIW.P*H' + X_hat/card_W;
S = (S + S')/2;

Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

K = GGIW.P*H'*iS;

GGIW_.m = GGIW.m + K*epsilon;
GGIW_.P = GGIW.P - K*H*GGIW.P;

temp = (W - z_bar);
Z = temp*temp';

X_sqrt = sqrtm_2by2(X_hat);
S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.v = GGIW.v + card_W;
GGIW_.V = GGIW.V + N + Z;

%GGIW predicted likelihood

lik = (GGIW.v-d-1)/2*log(det2(GGIW.V)) - (GGIW_.v-d-1)/2*log(det2(GGIW_.V))...
        + gamma2ln((GGIW_.v-d-1)/2) - gamma2ln((GGIW.v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_.a)...
        - gammaln(GGIW.a) + GGIW.a*log(GGIW.b) - GGIW_.a*log(GGIW_.b)...
        - (card_W*log(pi)+log(card_W))*d/2;
    
GGIW = GGIW_;

end

