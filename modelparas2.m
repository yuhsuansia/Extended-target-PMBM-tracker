load('SplitMergeScenario.mat')
load('groundTruth_SplitMerge.mat')

K = length(Scenario.Z{1});
model.K = K;

% Effective window length for the gamma prediction
w_e_gamma = 20;
% Effective window length for the extent prediction
w_e_extent = 10;

% Parameters to control the agility of the prediction
model.tao = 1/(log(w_e_extent)-log(w_e_extent-1));
model.eta = 1/(1-1/w_e_gamma);

model.Ts = 1;   %sampling interval
sigma_v = 0.5;  %standard deviation of motion noise
sigma_r = 0.1;  %standard deviation of measurement noise

% Linear motion and measurement models
model.motionmodel = motionmodel.cvmodel(model.Ts,sigma_v);
model.measmodel = measmodel.cvmeasmodel(sigma_r);

% Reconstruct the data structure

% Generate tracks (ground truth)
X = cell(K,1);  %kinematic states
E = cell(K,1);  %extent states
N = zeros(K,1); %number of targets
for targetnum = 1:length(targetTracks)
    for k = targetTracks(targetnum).birthTime:targetTracks(targetnum).deathTime
        targetstate = targetTracks(targetnum).x(1:model.motionmodel.d,k-targetTracks(targetnum).birthTime+1);
        targetextent = targetTracks(targetnum).X(:,:,k-targetTracks(targetnum).birthTime+1);
        X{k} = [X{k} targetstate];
        E{k} = cat(3,E{k},targetextent);
        N(k) = N(k) + 1;
    end
end

for k = 1:K
    groundTruth{k}.x = X{k};
    groundTruth{k}.X = E{k};
end

% Target existence probability
model.Ps = 0.99;
% Target detection probability
model.Pd = Scenario.detection_prob;

% Range of the surveillance area
range_c = [-1 1;-1 1]*100;
% Poisson false alarm (clutter) rate
lambda_c = Scenario.false_alarm_rate;
% Poisson clutter intensity
model.lambda_fa= lambda_c/prod(range_c(:,2)-range_c(:,1));

% Target initial state (one broad prior)
nbirths = 1;
xstart = zeros(model.motionmodel.d,nbirths);

% Poisson birth model
model.birth.w = 0.01*ones(nbirths,1);   %weights of birth components
% Each GGIW component is described by the sufficient statistics of Gamma
% distribution (a,b), Gaussian distirbution (m,P) and inverse-Wishart
% distributino (v,V). For inverse-Wishart distribution, v is a scalar
% controls the distribution uncertainty. The larger v, the smaller
% uncertainty.
model.birth.GGIW = repmat(struct('a',5e3,'b',1e3,'m',[],'P',diag([ 50; 50; 5; 5 ])*diag([ 50; 50; 5; 5 ])','v',14,'V',10*eye(2)),[nbirths,1]);
for i = 1:nbirths
    model.birth.GGIW(i).m = xstart(:,i);
end

% Gating parameters
Pg = 0.999; %gating size in probability
model.gamma= chi2inv(Pg,model.measmodel.d);
% Effective missed detection probability after applying gating
model.Qd = 1 - model.Pd*Pg;

% Pruning thresholds
model.threshold_r = 1e-3;   %existence probability of Bernoulli component
model.threshold_u = 1e-3;   %weight of mixture component in PPP
model.threshold_w = 1e-3;   %weight of global hypothesis (multi-Bernoulli)
model.threshold_s = 1e-4;   %weight of trajectory is alive if exists

model.recycle = 1e-1;       %recycling threshold
model.merge = 4;            %merge threshold used to merge similar GGIWs
model.M = 20;              %cap of number of MBM components in PMBM
model.num_iterations = 10;  %controls the number of iterations used in SO
model.max_repetition = 2;   %controls the number of iterations used in SO

%extract target state from Bernoulli components with existence probability
%larger than this threshold
model.exist_r = 0.5;

% DBSCAN parameters, a grid search for hyperparameters
model.max_dist = 5;
model.min_dist = 0.1;
model.grid_dist = 0.1;

%% Plot ground truth

if ifplot
    
    screen_size = get(0, 'ScreenSize');
    f1 = figure(1);
    set(f1, 'Position', [0 0 screen_size(3) screen_size(4)]);
    grid on;box on;hold on
    
    cols = parula(length(targetTracks));
    for it = 1:length(targetTracks)
        xx = targetTracks(it).x(1,:);
        yy = targetTracks(it).x(2,:);
        plot(xx,yy,'linewidth',2)
        for ii = 1:size(xx,2)
            %illustrate the 3-sigma level of ellipse
            [cx,cy]=Sigmacircle(xx(ii),yy(ii),targetTracks(it).X(:,:,ii),3);
            plot(cx,cy,'-','color',cols(it,:),'linewidth',1);
        end
    end
    
    xlabel('x');ylabel('y')
    xlim([-100,100])
    ylim([-100,100])
    xlabel('x (m)','Interpreter','latex')
    ylabel('y (m)','Interpreter','latex')
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',16)
    
end
