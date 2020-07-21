%% Generate object tracks
%draw samples from RFS of objects (dynamic equations) to simulate object 
%trajectories. Assume that a 2D measurement model is used.
%Choose object detection probability
clc;clear
dbstop if error

%Choose object survival probability
P_S = 0.99;
%Choose range of surveillance area
sensor_model.range_c = [-1000 1000;-1000 1000];
%Choose gamma distribution parameter
alpha = 40;
beta = 4;
sensor_model.gamma = makedist('Gamma','a',alpha,'b',1/beta);

%Create linear motion model
T = 1;
sigma_q = 5;
motion_model = motionmodel.cvmodel(T,sigma_q);

%Specify inverse-Wishart parameters
motion_model.iwish.V = 200*eye(2);
motion_model.iwish.v = 20;

%Create linear measurement model
sigma_r = 10;
meas_model = measmodel.cvmeasmodel(sigma_r);

%Create ground truth model
nbirths = 3;
K = 100;
birth_model = repmat(struct('w',log(0.03),'x',[],'P',100*eye(motion_model.d)),[1,nbirths]);
birth_model(1).x  = [ 0; 0; 0; -10 ];
birth_model(2).x  = [ 400; -600; -10; 5 ];
birth_model(3).x  = [ -800; -200; 20; -5 ];

%Generate object tracks
[object_tracks] = trackgen(K,meas_model,motion_model,sensor_model,birth_model,P_S);

%Plot trajectories
% trajectory = [object_tracks.x];

figure
box on
grid on

hold on

cols = parula(length(object_tracks));
for it = 1:length(object_tracks)
    xx = object_tracks(it).x(1,:);
    yy = object_tracks(it).x(2,:);
    plot(xx,yy,'linewidth',3)
    for ii = 1:size(xx,2)
        %illustrate the 3-sigma level of ellipse
        [cx,cy]=Sigmacircle(xx(ii),yy(ii),object_tracks(it).X(:,:,ii),3);
        plot(cx,cy,'-','color',cols(it,:),'linewidth',2);
    end
end

% plot(trajectory(1,:), trajectory(2,:), 'bo', 'Linewidth', 1);
h1 = plot(birth_model(1).x(1),birth_model(1).x(2),'rx', 'Markersize',16, 'Linewidth',3);
plot(birth_model(2).x(1),birth_model(2).x(2),'rx', 'Markersize',16, 'Linewidth',3);
plot(birth_model(3).x(1),birth_model(3).x(2),'rx', 'Markersize',16, 'Linewidth',3);

xlabel('x (m)','Interpreter','latex'); ylabel('y (m)','Interpreter','latex')
legend(h1,'Potential birth locations','Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex','FontSize',16);

%%

function [object_tracks] = trackgen(K,measmodel,motionmodel,sensormodel,birthmodel,P_S)
%TRACKGEN draw samples from RFS of objects (birth model, dynamic equations) 
%to simulate object trajectories. Assume that a 2D measurement model is used. 
%INPUT: K: total tracking time --- scalar
%       sensormodel: a structure specifies the sensor parameters
%           P_D: object detection probability --- scalar
%           lambda_c: average number of clutter measurements per time scan, 
%                     Poisson distributed --- scalar
%           pdf_c: clutter (Poisson) intensity --- scalar
%       motionmodel: a structure specifies the motion model parameters
%           d: object state dimension --- scalar
%           F: function handle return transition/Jacobian matrix
%           f: function handle return predicted object state
%           Q: motion noise covariance matrix
%       measmodel: a structure specifies the measurement model parameters
%           d: measurement dimension --- scalar
%           H: function handle return transition/Jacobian matrix
%           h: function handle return the observation of the target state
%           R: measurement noise covariance matrix
%       birthmodel: a structure array specifies the birth model (Gaussian
%       mixture density) parameters --- (1 x number of birth components)
%           w: weights of mixture components (in logarithm)
%           x: mean of mixture components
%           P: covariance of mixture components
%       object survival probability --- scalar
%OUTPUT:object_tracks: a structure array specifies object tracks ---
%       (number of tracks x 1) 
%           tbirth: track initiate (object appear) time --- scalar
%           tdeath: track end time (object disappear) time --- scalar
%           x: object trajectory --- (object state dimension x time steps 
%              object exists in the scen)
%Note that if the number of tracks is zero, set the output to empty

n = 0;

%Convert birthmodel.w
[log_w,log_sum_w] = normalizeLogWeights([birthmodel.w]);

%surveillance area
range_c = sensormodel.range_c;

for k = 1:K
    %randomly sample number of births
    nb = poissrnd(exp(log_sum_w));
    for i = 1:nb
        %randomly sample birth component
        idx = find(rand<cumsum(exp(log_w)),1,'first');
        n = n+1;
        %randomly sample kinematic state
        object_tracks(n,1).x = mvnrnd(birthmodel(idx).x, birthmodel(idx).P)';
        %randomly sample Poisson rate
        object_tracks(n,1).g = random(sensormodel.gamma);
        %randomly sample extent matrix
        %Note that for coordinated turn motion model, it is common to
        %rotate the extent matrix using the rotation matrices R obtained
        %from the object heading, such that X = R'XR
        object_tracks(n,1).X = iwishrnd(motionmodel.iwish.V, motionmodel.iwish.v);
        %birth time
        object_tracks(n,1).birthTime = k;
        %death time
        object_tracks(n,1).deathTime = k;
        
        time = k;
        xk = object_tracks(n,1).x;
        
        xypos = mvnrnd(measmodel.h(xk)', measmodel.R)';
        xpos = xypos(1);
        ypos = xypos(2);
        
        ifalive = 1;
        while ifalive && (xpos>=range_c(1,1))&&(xpos<=range_c(1,2))...
                &&(ypos>=range_c(2,1))&&(ypos<=range_c(2,2))&&(time+1<=K)

            %simulate the kinematic state
            xk = mvnrnd(motionmodel.f(xk), motionmodel.Q)';
            object_tracks(n,1).x = [object_tracks(n,1).x xk];
            object_tracks(n,1).deathTime = object_tracks(n,1).deathTime + 1;
            
            %randomly sample Poisson rate
            object_tracks(n,1).g = random(sensormodel.gamma);
            
            %randomly sample extent matrix
            %Note that for coordinated turn motion model, it is common to
            %rotate the extent matrix using the rotation matrices R obtained
            %from the object heading, such that X = R'XR
            object_tracks(n,1).X = cat(3,object_tracks(n,1).X,...
                iwishrnd(motionmodel.iwish.V, motionmodel.iwish.v));
            
            %simulate measurements
            time = time + 1;
            xypos = mvnrnd(measmodel.h(xk)', measmodel.R)';
            xpos = xypos(1);
            ypos = xypos(2);
            
            if rand > P_S
                ifalive = 0;
            end
        end
    end
end

if n==0
    object_tracks = [];
end 

end