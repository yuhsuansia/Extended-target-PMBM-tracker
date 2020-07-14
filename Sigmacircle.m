% Code by Karl Granstrom

function [x, y] = Sigmacircle(cx,cy,P,N,N2)

%Function used to draw an ellipse parameterised by matrix P, N is the sigma
%level

if nargin < 4
    N=1;
    N2 = 20;
elseif nargin < 5
    N2 = 20;
end

if size(P,1) == 2
    sqrtP = N*sqrtm(P);
else
    sqrtP = N*sqrtm(P);
end

phi = 0:pi/N2:2*pi;
xy = sqrtP*[cos(phi); sin(phi)];
xy = xy + [cx*ones(1,length(phi)) ; cy*ones(1,length(phi))];

x = xy(1,:)';
y = xy(2,:)';