% Copyright (c) 2019, Yuxuan Xia E-mail: yuxuan.xia@chalmers.se
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% This code implements an extended target Poisson multi-Bernoulli mixture
% tracker based on sets of trajectories for tracking both alive and dead
% trajectories. One can also choose to implement an extended target Poisson
% multi-Bernoulli mixture filter based on sets of targets

clear;close all;clc
dbstop if error

addpath('Data','Third-party code','Common')

%Choose which multiple extended object tracking algorithm to implement:
%MEOT 1: PMBM filter; MEOT 2: PMBM tracker
MEOT = 1;

if MEOT == 1
    fprintf('You have chosen to implement the extended target PMBM filter.\n');
    addpath('SetofTargets')
    rmpath('SetofAllTrajectories')
elseif MEOT ==2
    fprintf('You have chosen to implement the extended target PMBM tracker.\n');
    addpath('SetofAllTrajectories')
    rmpath('SetofTargets')
else
    error('Please choose an existing algorithm!')
end

%Choose data assocation method: dataAssocMethod 1: conventional two-step
%approach clustering (DBSCAN) + assignment (MURTY); dataAssocMethod 2:
%likelihood-based one-step approach: stochastic optimisation based sampling
dataAssocMethod = 2;

%Choose a scenario: Scenario 1: 27 targets born at four different
%locations; Scenario 2: targets move in proximity (a broad birth prior).
scenario = 1;

%If plot
ifplot = false;

%Parameter setting
if scenario == 1
    fprintf('You have chosen to run the algorithm on Scenario 1.\n');
    modelparas1;
elseif scenario == 2
    fprintf('You have chosen to run the algorithm on Scenario 2.\n');
    modelparas2;
else
    error('Please choose a valid scenario!')
end

if dataAssocMethod == 1
    fprintf('You have chosen to use clustering and assignment for data association.\n');
    model.dataAssocMethod = 1;
elseif dataAssocMethod == 2
    fprintf('You have chosen to use stochastic optimisation for data association.\n');
    model.dataAssocMethod = 2;
else
    error('Please choose a valid data association method!')
end

%Choose whether to perform track-oriented PMB approximation
model.ifTOPMB = false;

%Parameters used in GOSPA (and LP trajectory) metric
alpha = 2; %for MTT
c = 20;
p = 1;
gamma = 2;


%Number of time steps
K = model.K;

%Initialise memory
GOSPA = zeros(K,4); %[GOSPA, Localisation, Missed, False]
estimates = cell(K,1);

if MEOT == 2
    trajectoryEstimates = cell(K,1);
end

%%%%%%
if ifplot
    % For illustration purposes
    screen_size = get(0, 'ScreenSize');
    f2 = figure(2);
    set(f2, 'Position', [0 0 screen_size(3) screen_size(4)]);
    grid on
    box on
    
    if scenario == 1
        xlim([-200,200])
        ylim([-200,200])
    elseif scenario == 2
        xlim([-100,100])
        ylim([-100,100])
    end
end
%%%%%%

%Load measurements
idx = 1; % range from 1 to 100
Z = Scenario.Z{idx};

%PPP initialisation
PPP.w = log(model.birth.w); %weights in logarithm
PPP.GGIW = model.birth.GGIW;

%MBM initialisation
MBM.w = [];     % Global hypotheses weights
MBM.track = {}; % Local hypotheses trees
MBM.table = []; % Global hypotheses look-up table

fprintf('Time step: ');

tic
for k = 1:K
    
    %Print info
    fprintf('%d ',k);
    
    %Update step
    [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
    
    %Extract estimates
    if MEOT == 1
        %estimate of the multi-target states
        estimates{k} = estimator(MBM,model);
    elseif MEOT == 2
        %both estimate of the multi-target states and the estimate of the
        %full trajectories (both alive and dead)
        [estimates{k},trajectoryEstimates{k}] = estimator(MBM,model);
    end
    
    %Recycling (NOT IMPLEMENTED IN PAPER)
    [PPP,MBM] = recycleBern(PPP,MBM,model);
    
    %Evaluate filtering performance using GOSPA
    [GOSPA(k,1), ~, decomposed_cost] = GOSPA_extended(groundTruth{k}, estimates{k}, p, c, alpha);
    GOSPA(k,2) = decomposed_cost.localisation;
    GOSPA(k,3) = decomposed_cost.missed;
    GOSPA(k,4) = decomposed_cost.false;
    
    %Prediction Step
    [PPP,MBM] = predictPMBM(PPP,MBM,model);
    
    %%%%%%
    if ifplot
        % For illustration purposes
        
        %true states
        xx = groundTruth{k}.x(1,:);
        yy = groundTruth{k}.x(2,:);
        
        %estimated states
        xx_est = estimates{k}.x(1,:);
        yy_est = estimates{k}.x(2,:);
        
        h1 = plot(Z{k}(1,:),Z{k}(2,:),'kx','linewidth',1);
        hold on
        for ii = 1:size(xx,2)
            [cx,cy]=Sigmacircle(xx(ii),yy(ii),groundTruth{k}.X(:,:,ii),3);
            h2 = plot(cx,cy,'r-','linewidth',2);
        end
        for ii = 1:size(xx_est,2)
            [cx_est,cy_est]=Sigmacircle(xx_est(ii),yy_est(ii),estimates{k}.X(:,:,ii),3);
            h3 = plot(cx_est,cy_est,'b-','linewidth',2);
        end
        
        if size(xx_est,2) > 0
            legend([h1 h2 h3], 'Measurements', 'True object positions and extents', ...
                'Estimated object positions and extents', 'Interpreter','latex');
        else
            legend([h1 h2], 'Measurements', 'True object positions and extents', 'Interpreter','latex');
        end
        
        xlabel('x (m)','Interpreter','latex')
        ylabel('y (m)','Interpreter','latex')
        set(gca,'TickLabelInterpreter', 'latex');
        set(gca,'FontSize',16)
        
        pause(0.1)
        
        hold off
    end
    %%%%%%
end

if MEOT == 2
    [dxy, wMat, loc_cost, miss_cost, fa_cost, switch_cost] = ...
        LPTrajMetricWrapper(targetTracks, trajectoryEstimates{end}, c, p, gamma,K);
end

simulation_time = toc;
fprintf('\nExecution time: %.2fs\n',simulation_time);

%%%%%
% Illustrate cardinality estimation
if ifplot
    card = cellfun(@(x) size(x.x,2), groundTruth);
    card_est = cellfun(@(x) size(x.x,2), estimates);
    figure(3)
    grid on;box on;hold on
    
    if scenario == 1
        plot(1:100,card,'Linewidth',2)
        plot(1:100,card_est,'Linewidth',2)
    elseif scenario == 2
        plot(1:40,card,'Linewidth',2)
        plot(1:40,card_est,'Linewidth',2)
    end
    
    xlabel('Time step','Interpreter','latex')
    ylabel('Number of targets','Interpreter','latex')
    
    legend('Ground truth', 'Estimates', 'Interpreter','latex');
    
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',16)
end
%%%%%
