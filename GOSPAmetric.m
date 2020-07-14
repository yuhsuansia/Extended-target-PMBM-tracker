function [ Error ] = GOSPAmetric(X,Y,c,p)

% Code by Karl Granstrom
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

% GOSPA metric 
gospa = (gospa + 0.5*(c^p)*(Nx+Ny-2*absGamma))^(1/p);

% Number of misdetection
MissedError = Ny-absGamma;

% Number of false detection
FalseError = Nx-absGamma;


Error = [gospa,LocationError,MissedError,FalseError];

end