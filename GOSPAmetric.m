function [ Error ] = GOSPAmetric(X,Y,c,p)

%
% Function that computes the Generalized Optimal Subpattern Assignment (OSPA) Metric
% for the random sets X and Y.
%
% Assumed gamma Gaussian inverse Wishart model
%
% Input X (and Y) is a struct with J estimates
% X.x -- Nx * J matrix with kinematic vectors
% X.X -- d * d * J tensor with random matrices

% Number of elements in X and Y.
Nx = size(X.x,2);
Ny = size(Y.x,2);

% if there are no elements at all
if Nx==0 && Ny==0
    gospa = 0;
    LocationError = 0;
    MissedError = 0;
    FalseError = 0;
    Error = [gospa,LocationError,MissedError,FalseError];
    return
elseif Nx==0
    gospa = (0.5*(c^p)*Ny)^(1/p);
    LocationError = 0;
    MissedError = Ny;
    FalseError = 0;
    Error = [gospa,LocationError,MissedError,FalseError];
    return
elseif Ny==0
    gospa = (0.5*(c^p)*Nx)^(1/p);
    LocationError = 0;
    MissedError = 0;
    FalseError = Nx;
    Error = [gospa,LocationError,MissedError,FalseError];
    return
end

% Distance matrix
D = repmat(c,[Nx Ny]);

for ix = 1:Nx
    for iy = 1:Ny
        % Gaussian Wasserstein Distance for kinematics and extent
        gwd = GaussianWassersteinDistance(X.x(1:2,ix),X.X(:,:,ix),Y.x(1:2,iy),Y.X(:,:,iy));
        % Poisson rates
        %prd = abs(X.g(ix)-Y.g(iy));
        
        % Apply threshold c
        D(ix,iy) = min(c,gwd);
    end
end

% Allocate memory
gospa = 0;
LocationError = 0;
absGamma = 0;

if Ny<=Nx
    % Compute assignment
    [Customer2Item,~] = auctionAlgorithm(-D');
    
    % Iterate over true targets
    for iy = 1:Ny
        % Check if distance is small enough
        if D(Customer2Item(iy),iy) < c
            % Location part of GOSPA
            gospa = gospa + D(Customer2Item(iy),iy)^p;
            % Location error
            LocationError = LocationError + D(Customer2Item(iy),iy);
            % Number of assignments
            absGamma = absGamma + 1;
        end
    end
    
else
    % Compute assignment
    [Customer2Item,~] = auctionAlgorithm(-D);
    
    % Iterate over estimates
    for ix = 1:Nx
        % Check if distance is small enough
        if D(ix,Customer2Item(ix)) < c
            % Location part of GOSPA
            gospa = gospa + D(ix,Customer2Item(ix))^p;
            % Location error
            LocationError = LocationError + D(ix,Customer2Item(ix));
            % Number of assignments
            absGamma = absGamma + 1;
        end
    end
    
end

gospa = (gospa + 0.5*(c^p)*(Nx+Ny-2*absGamma))^(1/p);

% Missed detection error
MissedError = Ny-absGamma;

% False alarm error
FalseError = Nx-absGamma;

Error = [gospa,LocationError,MissedError,FalseError];


% if Nx>Ny % assume X contains fewer elements than Y
%     tmp = X;
%     X = Y;
%     Y = tmp;
%     Nx = size(X,2);
%     Ny = size(Y,2);
% end
% 
% % Construct distance matrix D
% D = repmat(c,[Ny Ny]);
% for i = 1:Nx
%     xi = X(1:2,i);
%     for j = 1:Ny
%         xj = Y(1:2,j);
%         % Euclidean norm
%         D(i,j) = min(c,norm(xi-xj));
%     end
% end
% 
% % Use the auction algorithm to compute the best assignments
% % auction maximizes the reward. Maximizing the negative distance is equal
% % to minimizing the positive distance.
% [~,Item2Customer] = auctionAlgorithm(-D);
% 
% %D = dg+dx+dX;
% 
% % allocate memory
% alphas = repmat(c,[1 Ny]);
% % compute the alphas
% for j = 1:Ny
%     if Item2Customer(j) <= Nx
%         alphas(j) = D(Item2Customer(j),j);
%     end
% end
% 
% % Compute the OSPA metric
% gospa   = (mean(alphas.^p))^(1/p);


