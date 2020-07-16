function [w2,a2,b2,m2,P2,nu2,V2] = GGIW_merge(w,a,b,m,P,nu,V)

%Merge a mixture of GGIWs into a single GGIW, code by Karl Granstrom

nu_min = 7;
v_min = 1;

% Dimensions
n_x = size(m,1);
d = size(V,1);
% Number of components
N=length(nu);

% Sum of weights
wb = sum(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct merging components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1 = 0;
c2 = 0;
C1 = zeros(d,d);
C2 = 0;
C3 = 0;
for i = 1:N
    % Gammas
    c1 = c1+w(i)*(psi(0,a(i))-log(b(i)));
    c2 = c2+w(i)*a(i)/b(i);
    % inverse Wisharts
    C1 = C1+w(i)*(nu(i)-d-1)*(V(:,:,i)\eye(d));
    C2 = C2+w(i)*sum(psi(0,(nu(i)-d-(1:d))/2));
    C3 = C3 + w(i)*log(det(V(:,:,i)));
end
c = c1/wb - log(c2/wb);
C = d*wb*log(wb) -wb*log(det(C1)) +C2 -C3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge the gammas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
a_k = sum(w(:).*a(:))/wb;

iter = 1;
while iter<20
    iter = iter+1;
    
    h_k = log(a_k)-psi(0,a_k)+c;
    hp_k = 1/a_k - psi(1,a_k);
    hb_k = -1/a_k^2 - psi(2,a_k);
    
    a_new = a_k-(2*h_k*hp_k)/(2*hp_k^2-h_k*hb_k); % Halley's
    
    if abs(a_new-a_k)<1e-2
        a_k = a_new;
        break
    else
        a_k = a_new;
    end
    
    a_k=max(a_k,v_min);
end
a2 = a_k;
b2 = a2/(sum(w(:).*a(:)./b(:))/wb);
catch
    a2 = sum(w(:).*a(:))/wb;
    b2 = sum(w(:).*b(:))/wb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge the normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2 = sum(m.*repmat(w(:)',n_x,1),2)/wb;

P2 = zeros(n_x,n_x);
for i = 1:N
    P2 = P2 + w(i)*(P(:,:,i)+(m(:,i)-m2)*(m(:,i)-m2)');
end
P2 = P2/wb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge Inverse Wisharts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_k = mean(nu);

iter = 1;
while iter<100
    iter=iter+1;

    h_k  = d*wb*log(v_k-d-1)...
        -wb*sum(psi(0,(v_k-d-(1:d))/2))...
        +C;
    hp_k = d*wb/(v_k-d-1)...
        -0.5*wb*sum(psi(1,(v_k-d-(1:d))/2));
    hb_k = -d*wb/((v_k-d-1)^2)...
        -0.25*wb*sum(psi(2,(v_k-d-(1:d))/2));
   
%     v_new = v_k-h_k/hp_k; % Newtons
    v_new = v_k-(2*h_k*hp_k)/(2*hp_k^2-h_k*hb_k); % Halley's
    
    if abs(v_new-v_k)<1e-2
        v_k = v_new;
        break
    else
        v_k = v_new;
    end
    
    v_k=max(v_k,nu_min);
end
nu2=v_k;
V2 = zeros(d,d);
for i=1:N
    V2 = V2 + w(i)*(nu(i)-d-1)*(V(:,:,i)\eye(d));
end
V2 = (V2+V2')/2;
V2 = (nu2-d-1)*wb*(V2\eye(d));

w2 = wb;