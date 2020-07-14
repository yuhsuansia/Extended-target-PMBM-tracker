clear;close all;clc
dbstop if error

%Choose a scenario: Scenario 1: 27 targets born at four different
%locations; Scenario 2: targets move in proximity.
scenario = 1;

%If plot
ifplot = false;

%Parameter setting
if scenario == 1
    modelparas1;
elseif scenario == 2
    modelparas2;
else
    error('Please choose a valid scenario!')
end

%Number of Monte Carlo simulations
numMC = length(Scenario.Z);

%Parameters used in GOSPA metric
c = 20;
p = 1;

%Number of time steps
K = model.K;

%Initialise memory
GOSPA = zeros(K,4,numMC);               %GOSPA error
trajectoryEstimates = cell(numMC,1);    %Trajectory estimates
simulation_time = zeros(numMC,1);       %Simulation time for a complete run

for t = 1:1
    
    fprintf('\nMonte Carlo run: %d.\n',t);
    
    %%%%%%
    if ifplot
        % For illustration purposes
        screen_size = get(0, 'ScreenSize');
        f2 = figure;
        set(f2, 'Position', [0 0 screen_size(3) screen_size(4)]);
        grid on
        box on
        
        xlim([-200,200])
        ylim([-200,200])
    end
    %%%%%%
    
    %Load measurements
    Z = Scenario.Z{t};
    
    %PPP initialisation
    PPP.w = log(model.birth.w); %weights in logarithm
    PPP.GGIW = model.birth.GGIW;
    
    %MBM initialisation
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Local hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    
    estimates = cell(K,1);
    trajectoryEstimates{t} = cell(K,1);
    
    fprintf('Time step: ');
    
    tic
    for k = 1:K
        
        %Print info
        fprintf('%d ',k);
        
        %Update step
        [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
        
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory (both alive and dead))
        [estimates{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);
        
        %Recycling (NOT IMPLEMENTED IN PAPER)
        %         [PPP,MBM] = recycleBern(PPP,MBM,model);
        
        %Evaluate filtering performance using GOSPA
        GOSPA(k,:,t) = GOSPAmetric(estimates{k},groundTruth{k},c,p);
        
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
            
            legend([h1 h2 h3], 'Measurements', 'True object positions and extents', ...
                'Estimated object positions and extents', 'Interpreter','latex');
            
            xlabel('x (m)','Interpreter','latex')
            ylabel('y (m)','Interpreter','latex')
            axesH = gca;
            axesH.XAxis.TickLabelInterpreter = 'latex';
            axesH.YAxis.TickLabelInterpreter = 'latex';
            set(gca,'FontSize',16)
            
            pause(0.1)
            
            hold off
        end
        %%%%%%
    end
    
    simulation_time(t) = toc;
    fprintf('\nExecution time: %.2fs\n',simulation_time(t));
    
    %%%%%
    % Illustrate cardinality estimation
    if ifplot
        card = cellfun(@(x) size(x.x,2), groundTruth);
        card_est = cellfun(@(x) size(x.x,2), estimates);
        figure
        grid on
        box on
        hold on
        
        plot(1:100,card,'Linewidth',2)
        plot(1:100,card_est,'Linewidth',2)
        xlabel('Time step','Interpreter','latex')
        ylabel('Number of targets','Interpreter','latex')
        
        legend('Ground truth', 'Estimates', 'Interpreter','latex');
            
        axesH = gca;
        axesH.XAxis.TickLabelInterpreter = 'latex';
        axesH.YAxis.TickLabelInterpreter = 'latex';
        set(gca,'FontSize',16)
        
    end
    %%%%%
end
