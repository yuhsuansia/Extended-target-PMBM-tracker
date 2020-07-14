% Code by Yuxuan Xia

function [GGIW_c,lik_c] = updateGGIWforPPP(GGIWp,C,model)

%Perform GGIW update for each mixture component

np = length(GGIWp);

d = 2;
card_W = size(C,2);
z_bar = mean(C,2);
temp = (C - z_bar);
Z = temp*temp';

GGIW_c = repmat(struct('a',[],'b',[],'m',[],'P',[],'v',[],'V',[]),[np,1]);
lik_c = zeros(np,1);

for i = 1:np
    
    GGIW_c(i,1).a = GGIWp(i,1).a + card_W;
    GGIW_c(i,1).b = GGIWp(i,1).b + 1;

    epsilon = z_bar - model.measmodel.h(GGIWp(i,1).m);
    H = model.measmodel.H(GGIWp(i,1).m);
    
    X_hat = GGIWp(i,1).V/(GGIWp(i,1).v - 2*d - 2);
    X_hat = (X_hat + X_hat')/2;
    
    S = H*GGIWp(i,1).P*H' + X_hat/card_W;
    S = (S + S')/2;
    
    Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
    
    K = GGIWp(i,1).P*H'*iS;
    
    GGIW_c(i,1).m = GGIWp(i,1).m + K*epsilon;
    GGIW_c(i,1).P = GGIWp(i,1).P - K*H*GGIWp(i,1).P;

    X_sqrt = sqrtm_2by2(X_hat);
    S_sqrt_inv = sqrtm(iS);
    N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';
    
    GGIW_c(i,1).v = GGIWp(i,1).v + card_W;
    GGIW_c(i,1).V = GGIWp(i,1).V + N + Z;
    
    %GGIW predicted likelihood
    lik_c(i) = (GGIWp(i,1).v-d-1)/2*log(det2(GGIWp(i,1).V)) - (GGIW_c(i,1).v-d-1)/2*log(det2(GGIW_c(i,1).V))...
        + gamma2ln((GGIW_c(i,1).v-d-1)/2) - gamma2ln((GGIWp(i,1).v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_c(i,1).a)...
        - gammaln(GGIWp(i,1).a) + GGIWp(i,1).a*log(GGIWp(i,1).b) - GGIW_c(i,1).a*log(GGIW_c(i,1).b)...
        - (card_W*log(pi)+log(card_W))*d/2;
    
end

end

