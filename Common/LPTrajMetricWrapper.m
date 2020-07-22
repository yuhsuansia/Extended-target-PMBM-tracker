function [dxy, wMat, loc_cost, miss_cost, fa_cost, switch_cost] ...
    = LPTrajMetricWrapper(targetTracks, traj_est, c, p, gamma,K)

%reconstruct ground truth data
T = K;

nx = length(targetTracks);
X.tVec = [targetTracks.birthTime]';
X.iVec = [targetTracks.deathTime]' - X.tVec + 1;

X.xState = zeros(6,T,nx);

for i = 1:nx
    for t = X.tVec(i):X.tVec(i)+X.iVec(i)-1
        X.xState(1:2,t,i) = targetTracks(i).x(1:2,t-X.tVec(i)+1);
        X.xState(3:6,t,i) = reshape(targetTracks(i).X(:,:,t-X.tVec(i)+1),4,1);
    end
end

%reconstruct estimates

ny = length(traj_est);

Y.tVec = [traj_est.t_birth]';
Y.iVec = [traj_est.t_death]' - Y.tVec + 1;

Y.xState = zeros(6,T,ny);

for i = 1:ny
    for t = Y.tVec(i):Y.tVec(i)+Y.iVec(i)-1
        Y.xState(1:2,t,i) = traj_est(i).x(1:2,t-Y.tVec(i)+1);
        Y.xState(3:6,t,i) = reshape(traj_est(i).X(:,:,t-Y.tVec(i)+1),4,1);
    end
end

[dxy, wMat, loc_cost, miss_cost, fa_cost, switch_cost] ...
    = LPTrajMetric_sparse_extended(X, Y, c, p, gamma);

end