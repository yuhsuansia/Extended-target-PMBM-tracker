function [dxy, wMat, loc_cost, miss_cost, fa_cost, switch_cost] ...
    = LPTrajMetric_sparse_extended(X, Y, c, p, gamma)
%% function [dxy, wMat, eVal] = LPTrajMetric(X, Y, c, p, gamma)
% This function computes the LP metric between sets of trajectories defined
% in https://arxiv.org/abs/1605.01177. 
% -------------------------------------------------------------------------
% Input:
% X, Y: sets of trajctories which are structs as follows:
%   X.tVec: 'nx x 1' dimensional vector that has start times of the 'nx'
%       trajectories in 'X'.
%   X.iVec: 'nx x 1' dimensional vector that has the duration of the 'nx'
%       trajectories in 'X'.
%   X.xState: 'stDim x T x nx' dimensional matrix, where 'stDim' is the 
%       state dimension, 'T' is the max length of the trajectories. The 
%       states of trajectory 'ind', 'X.xState(:, :, ind)' has '0' values 
%       outisde '[X.tVec(ind), X.tVec(ind)+X.iVec(ind)-1]'. Note that within the
%       window where X.xState is valid can have 'holes', with 'nan' values.
% c: >0, cut-off parameter
% p: >= 1, exponent parameter
% gamma: >0, track switch penalty
% -------------------------------------------------------------------------
% Output:
% dxy: Metric value
% wMat: Assignment matrix of dimension '(nx+1) x (ny+1) x T' has value
% between '0' and '1'
% loc_cost: localisation cost over time of dimension 'T x 1'
% miss_cost: cost for missed targets over time, dimension 'Tx1'
% fa_cost: cost for false targets over time, dimension 'Tx1'
% switch_cost: cost for switches over time, dimension '(T-1)x1'
% -------------------------------------------------------------------------
% Function modified by Angel F. Garcia-Fernandez: We use sparse matrices to compute the solution
% Also by Yuxuan Xia: Use Gaussian Wasserstein distance for base distance
%%%%%%%%%% Parameters of use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = size(X.xState, 3);
ny = size(Y.xState, 3);
nxny = nx*ny;
T = size(X.xState, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% localisation cost computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAB = locCostComp_v2(X, Y, c, p);
nxny2 = size(DAB, 1) * size(DAB, 2); % = (nx+1) * (ny+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  variables to be estimated in LP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [W_1,1(1) W_2,1(1) .. W_nx+1,1(1), .. W_1,ny+1(1) W_2,ny+1(1) ...
% W_nx+1,ny+1(1), W_1,1(T) W_2,1(T) .. W_nx+1,1(T), .. W_1,ny+1(T)
% W_2,ny+1(T) ... W_nx+1,ny+1(T) e(1) .. e(T-1) h_11(1) .. h_nx,ny(1) ...
% h_1,1(T-1) ... h_nx,ny(T-1)]'

%%% Length of the variable components in x
WtLen = nxny2*T; etLen = T-1;  htLen = nxny * (T-1);
nParam = WtLen + etLen + htLen; % total number of variables

%%% Position of the variable components in x
WtPos = (1:WtLen);  etPos = WtLen+(1:etLen);
htPos = WtLen + etLen + (1:htLen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  objective function f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(nParam, 1);
f(WtPos) = reshape(DAB, [WtLen,1]); % for vec(W(1)) to vec(W(T)), loc cost
f(etPos) = 0.5 * gamma^p * ones(T-1,1); %for e(1) to e(T-1), switch cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  equality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constraint 1 %%%%
%%% sum_i Wij(t) = 1 for j = 1 to ny, t = 1 to T
temp1 = []; 
tempOne = ones(1, nx+1);
for ind1 = 1:ny
    temp1 = sparse(blkdiag(temp1, tempOne));
end
temp1 = cat(2, temp1, zeros(ny, nx+1)); % sum_i Wij(t) for one t
Aeq1 = [];
for t  = 1:T
    Aeq1 = blkdiag(Aeq1, temp1);
end
Aeq1 = [Aeq1, sparse(T*ny, nParam-WtLen)]; % sum_i Wij(t) for all t
beq1 = ones(ny*T, 1);

%%%% Constraint 2 %%%%
%%% sum_j Wij(t) = 1 for i = 1 to nx, t = 1 to T

temp2= [speye(nx),sparse(nx,1)];

temp1 = [];
for ind1 = 1:ny+1
    temp1 = cat(2, temp1, temp2); % sum_j Wij(t) for one t
end
Aeq2 = [];
for t  = 1:T
    Aeq2 = blkdiag(Aeq2, temp1);
end
beq2 = ones(nx*T, 1);
Aeq2 = [Aeq2, sparse(T*nx, nParam-WtLen)]; % sum_j Wij(t) for all t

Aeq = [Aeq1; Aeq2]; beq = [beq1; beq2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  upper and lower bound constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 <= W < inf, 0 <= e < inf, 0 <= h < inf
lb = zeros(nParam, 1); 
ub = inf(nParam, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  inequality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constraint 1 %%%%
%%% sum_ij h_ij(t) - e(t) <= 0 for i = 1..nx, j = 1..ny, t = 1..T-1
%b1 = sparse((T-1),1); 
A1 = sparse((T-1),nParam);
for t = 1:T-1
    A1(t,WtLen+t) = -1;    % for -e(t)
    colInd = WtLen + etLen + (t-1)*nxny + (1:nxny); % for hij(t), all i,j
    A1(t,colInd) = 1;
end

%indeces_rows=t*ones(T-1,1);



%%%% Constraint 2 %%%%
%%% -W_ij(t+1)+W_ij(t)-h_ij(t) <= 0 for i = 1..nx, j = 1..ny, t = 1..T-1
temp1 = [speye(nx) sparse(nx,1)]; % W_ij(t) for i=1 to nx and one t,j
temp2 = [];
for ind = 1:ny % extending temp1 for all j
    temp2 = sparse(blkdiag(temp2, temp1));
end
temp2_t = [temp2 sparse(nxny,nParam-size(temp2,2))];
% W_ij(t) for i = 1..nx, j = 1..ny and one t
temp2_t1 = circshift(temp2_t, [0,nxny2]);
% W_ij(t+1) for i = 1..nx, j = 1..ny and one t

temp3_t = sparse(htLen, nParam);
temp3_t1 = sparse(htLen, nParam);
for t = 1:T-1 % extending temp2_t, temp2_t1 for t = 1..T-1
    temp3_t((t-1)*nxny+(1:nxny),:) = circshift(temp2_t, [0,(t-1)*nxny2]);
    temp3_t1((t-1)*nxny+(1:nxny),:) = circshift(temp2_t1, [0,(t-1)*nxny2]);
end
%b3 = sparse(nxny*(T-1), 1);
A3 = temp3_t - temp3_t1;
% W_ij(t) -W_ij(t+1) for i = 1..nx, j = 1..ny  and t = 1..T-1
A3(:, htPos) = -speye(htLen);
% -h_ij(t) for i = 1..nx, j = 1..ny  and t = 1..T-1

%%%% Constraint 3 %%%%
%%%  W_ij(t+1) - W_ij(t) - h_ij(t) <= 0 for i = 1 to nx, j = 1 to ny
%b4 = sparse(nxny*(T-1), 1);
A4 = - temp3_t + temp3_t1;
% - W_ij(t) + W_ij(t+1) for i = 1..nx, j = 1..ny  and t = 1..T-1
A4(:, htPos) = -speye(htLen);
% -h_ij(t) for i = 1..nx, j = 1..ny  and t = 1..T-1

A = [A1; A3; A4]; 
b = sparse((T-1)+nxny*(T-1)+nxny*(T-1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linProgOptions = optimoptions('linprog', 'Display','off');
[x, dxy] = linprog(f, A, b, Aeq, beq, lb, ub, [], linProgOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  Metric and the assignment values to be returned %%%%%%%%%%%%%%%
wMat = reshape(x(1:nxny2*T), [nx+1, ny+1, T]);
dxy = dxy.^(1/p);
[loc_cost, miss_cost, fa_cost, switch_cost] ...
    = computeLocFalseMissedSwitchCosts(wMat, DAB, X, Y, c, p, gamma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function locCostMat = locCostComp_v2(X, Y, c, p)
% function locCostMat = locCostComp(stMat, X, Y, c, p)
% computing the localisation cost at each time 't' for every (i,j)

tmpCost = c^p / 2; % cost for being unassigned
T = size(X.xState, 2); nx = size(X.xState, 3); ny = size(Y.xState, 3);
locCostMat  = zeros(nx+1, ny+1, T);
for t = 1:T
    for xind = 1:nx+1
        if (xind <= nx) % xind not dummy
            if (t>=X.tVec(xind)) && (t<=(X.tVec(xind)+X.iVec(xind)-1))
                % if X_xind exists at t
                for yind = 1:ny+1
                    if (yind <= ny) && (t >= Y.tVec(yind)) && ...
                            (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                        % if Y_yind exists at t
                        locCostMat(xind,yind,t) = computeLocCostPerTime( ...
                            X.xState(:,t,xind), Y.xState(:,t,yind), c, p);
                        
                    else % yind does not exist or yind is dummy
                        locCostMat(xind,yind,t) = tmpCost;
                    end
                end
            else        % if X_xind does not exist at t
                for yind = 1:ny
                    if (t >= Y.tVec(yind)) && ...
                            (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                        % if Y_yind exists at t
                        locCostMat(xind,yind,t) = tmpCost;
                    end
                end
            end
        else    % xind is dummy
            for yind = 1:ny
                if (t >= Y.tVec(yind)) && ...
                        (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                    % if Y_yind exists at t
                    locCostMat(xind,yind,t) = tmpCost;
                end
            end
        end
    end
end
end

function d = computeLocCostPerTime(x, y, c, p)
if all(~isnan(x)) && all(~isnan(y)) 
    % neither x nor y has hole
    m1 = x(1:2);
    P1 = reshape(x(3:6),2,2);
    m2 = y(1:2);
    P2 = reshape(y(3:6),2,2);
    gwd = GaussianWassersteinDistance(m1,P1,m2,P2);
	d = min(gwd^p,c^p);
elseif any(isnan(x) & ~isnan(y)) || any(~isnan(x) & isnan(y)) 
    % exactly one of x and y has hole
    d = c^p/2;
else
    d = 0;
end    
end

function [loc_cost, miss_cost, fa_cost, switch_cost] ...
    = computeLocFalseMissedSwitchCosts(w_mat, locCostMat, X, Y, c, p, gamma)
% computing the localisation cost, swtiching cost and cost for missed and false 
% targets.

tmp_cost = c^p / 2; % cost for being unassigned
T = size(X.xState, 2); nx = size(X.xState, 3); ny = size(Y.xState, 3);

switch_cost = 0.5 * gamma^p * ...
    squeeze(sum(sum(abs(diff(w_mat(1:nx, 1:ny, :), 1, 3)))));

loc_mask = zeros(size(w_mat));
miss_mask = zeros(size(w_mat));
fa_mask = zeros(size(w_mat));

fa_miss_mask = zeros(size(w_mat)); %Accounts for false and missed target costs that arise for a localisation cost of c^p

for t = 1:T
    % localisation and miss cost calculations
    for xind = 1:nx
        if (t>=X.tVec(xind)) && (t<=(X.tVec(xind)+X.iVec(xind)-1))
            % if X_xind exists at t
            for yind = 1:ny
                if (yind <= ny) && (t >= Y.tVec(yind)) && ...
                        (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                    % if Y_yind exists at t     
                    if all(~isnan(X.xState(:,t,xind))) && ...
                            all(~isnan(Y.xState(:,t,yind)))
                    % there is no hole in x or y at this time
                        %%% add to localisation cost at the time based on
                        %%% weight (unless the weight is c^p
                        
                        if(locCostMat(xind, yind, t)<2*tmp_cost)  
                            loc_mask(xind, yind, t) = 1;
                        else
                            fa_miss_mask(xind, yind, t) = 1;
                        end
                    elseif any(isnan(X.xState(:,t,xind))) && ...
                            all(~isnan(Y.xState(:,t,yind)))
                     % there is a hole in x but no hole in y
                        fa_mask(xind, yind, t) = 1;
                    elseif all(~isnan(X.xState(:,t,xind))) && ...
                            any(isnan(Y.xState(:,t,yind)))
                        % there is no hole in x but a hole in y
                        miss_mask(xind, yind, t) = 1;
                    end
                else % yind does not exist
                    
                    %%% add to miss cost at the time
                    miss_mask(xind, yind, t) = 1;
                end
            end
            %%% add to miss cost at the time for yind = ny+1
            yind = ny+1;
            miss_mask(xind, yind, t) = 1;
        end
    end
    
    for yind = 1:ny
        if (t>=Y.tVec(yind)) && (t<=(Y.tVec(yind)+Y.iVec(yind)-1))
            % if Y_yind exists at t
            if all(~isnan(Y.xState(:, t, yind))) % no hole in y
                for xind = 1:nx
                    if ~((xind <= nx) && (t >= X.tVec(xind)) && ...
                            (t <= (X.iVec(xind) + X.tVec(xind) - 1)))
                        % if X_xind does not exist at t
                        
                        %%% add to fa cost at the time
                        fa_mask(xind, yind, t) = 1;
                    end
                end
            end
            
            %%% add to fa cost at the time for xind = nx+1
            xind = nx+1;
            fa_mask(xind, yind, t) = 1;
        end
    end
end

loc_cost = squeeze(sum(sum(locCostMat .* w_mat .* loc_mask, 1), 2));
miss_cost = tmp_cost * squeeze(sum(sum(w_mat .* miss_mask, 1), 2))+ tmp_cost * squeeze(sum(sum(w_mat .* fa_miss_mask, 1), 2));
fa_cost = tmp_cost * squeeze(sum(sum(w_mat .* fa_mask, 1), 2))+ tmp_cost * squeeze(sum(sum(w_mat .* fa_miss_mask, 1), 2));
end
