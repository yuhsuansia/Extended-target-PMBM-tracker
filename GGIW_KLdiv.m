function [KLdiv] = GGIW_KLdiv(GGIW1,GGIW2)

% Function that computes the KL-divergence between two
% Gamma-Gaussian-inverse Wishart (GGIW) distributions
% Code by Karl granstrom

% p(g,x,X) = Gam(g;a1,b1)N(x;m1,P1)IW(X;v1,V1)
% q(g,x,X) = Gam(g;a2,b2)N(x;m2,P2)IW(X;v2,V2)
%
% D_KL = D(p||q) = int p ln(p/q) dx

% Dimension of extension
d = size(GGIW1.V,1);
% Dimension of kinematical state
n_x = length(GGIW1.m);

KLdiv_G = GGIW1.a*log(GGIW1.b)-GGIW2.a*log(GGIW2.b)+gammaln(GGIW2.a)-gammaln(GGIW1.a)...
    +(GGIW1.a-GGIW2.a)*(psi(0,GGIW1.a)-log(GGIW1.b))+GGIW1.a*(GGIW2.b/GGIW1.b-1);

KLdiv_N = (...
    -0.5*log(det(GGIW1.P))+0.5*log(det(GGIW2.P))...
    -0.5*n_x+0.5*(GGIW1.m-GGIW2.m)'*(GGIW2.P\(GGIW1.m-GGIW2.m))...
    +0.5*trace(GGIW2.P\GGIW1.P)...
    );

KLdiv_IW = (...
    0.5*(GGIW1.v-d-1)*log(det(GGIW1.V))-0.5*(GGIW2.v-d-1)*log(det(GGIW2.V))...
    +sum(gammaln((GGIW2.v-d-(1:d))/2)-gammaln((GGIW1.v-d-(1:d))/2))...
    +0.5*(GGIW2.v-GGIW1.v)*(log(det(GGIW1.V))-sum(psi(0,(GGIW1.v-d-(1:d))/2)))...
    +trace(-0.5*(GGIW1.v-d-1)*(GGIW1.V\(GGIW1.V-GGIW2.V)))...
    );

KLdiv = KLdiv_G + KLdiv_N + KLdiv_IW;

end