clear;clc
dbstop if error

%Choose a scenario: Scenario 1: 27 targets born at four different
%locations; Scenario 2: targets move in proximity.
scenario = 1;

%Parameter setting
if scenario == 1
    modelparas1;
else
    modelparas2;
end

%Number of Monte Carlo Simulations
numMC = length(Scenario.Z);

%Parameters used in GOSPA metric
c = 20;
p = 1;

%Number of time steps
K = model.K;

GOSPA = zeros(K,4,numMC);
trajectoryEstimates = cell(numMC,1);
simulation_time = zeros(numMC,1);

for t = 1:numMC
    Z = Scenario.Z{t};
    % Initialisation
    PPP.w = log(model.birth.w);
    PPP.GGIW = model.birth.GGIW;
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Locl hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    
    estimates = cell(K,1);
    trajectoryEstimates{t} = cell(K,1);
    
    tic
    for k = 1:K
        %Print info
        [t,k]
        %Update step
        [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);
        %Evaluate filtering performance using GOSPA
        GOSPA(k,:,t) = GOSPAmetric(estimates{k},groundTruth{k},c,p);
        %Prediction Step
        if k < K
            [PPP,MBM] = predictPMBM(PPP,MBM,model);
        end
    end
    simulation_time(t) = toc;
end
