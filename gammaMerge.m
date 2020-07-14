function [w2,a2,b2] = gammaMerge(w,a,b)

% Function used to merge a Gamma mixture, code by Karl Granstrom

% Number of components
N=length(a);

% Sum of weights
wb = sum(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct merging components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1 = 0;
c2 = 0;
for i = 1:N
    % Gammas
    c1 = c1+w(i)*(psi(0,a(i))-log(b(i)));
    c2 = c2+w(i)*a(i)/b(i);
end
c = c1/wb - log(c2/wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge the gammas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_k = sum(w(:).*a(:))/wb;

iter = 1;
while iter<100
    iter = iter+1;
    
    h_k = log(a_k)-psi(0,a_k)+c;
    hp_k = 1/a_k - psi(1,a_k);
    hb_k = -1/a_k^2 - psi(2,a_k);
    
    a_new = a_k-(2*h_k*hp_k)/(2*hp_k^2-h_k*hb_k); % Halley's
    
    if abs(a_new-a_k)<1e-5
        a_k = a_new;
        break
    else
        a_k = a_new;
    end
end
a2 = a_k;
b2 = a2/(sum(w(:).*a(:)./b(:))/wb);

w2 = wb;

end