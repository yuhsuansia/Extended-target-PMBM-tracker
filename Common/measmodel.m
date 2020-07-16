% Code by Yuxuan Xia

classdef measmodel
    %MEASMODEL is a class containing different measurement models
    
    methods (Static)
        function obj = cvmeasmodel(sigma)
            %CVMEASMODEL creates the measurement model for a 2D nearly constant velocity motion model
            %INPUT:     sigma: standard deviation of measurement noise --- scalar
            %OUTPUT:    obj.d: measurement dimension --- scalar
            %           obj.H: function handle return an observation matrix --- 2 x 4 matrix
            %           obj.R: measurement noise covariance --- 2 x 2 matrix
            %           obj.h: function handle return a measurement --- 2 x 1 matrix
            obj.d = 2;
            obj.H = @(x) [1 0 0 0;0 1 0 0];
            obj.R = sigma^2*eye(2);
            obj.h = @(x) obj.H(x)*x;
        end
        
        function obj = ctmeasmodel(sigma)
            %CTMEASMODEL creates the measurement model for a 2D
            %coordinate turn motion model
            %INPUT:     sigma: standard deviation of measurement noise --- scalar
            %OUTPUT:    obj.d: measurement dimension --- scalar
            %           obj.H: function handle return an observation matrix --- 2 x 5 matrix
            %           obj.R: measurement noise covariance --- 2 x 2 matrix
            %           obj.h: function handle return a measurement --- 2 x 1 matrix
            % NOTE: the first two entries of the state vector represents
            % the X-position and Y-position, respectively.
            obj.d = 2;
            obj.H = @(x) [1 1 0 0 0;0 1 0 0 0];
            obj.R = sigma^2*eye(2);
            obj.h = @(x) obj.H(x)*x;
        end
        
        function obj = bearingmeasmodel(sigma, s)
            %BEARINGMEASUREMENT creats the bearing measurement model
            %INPUT: sigma: standard deviation of measurement noise --- scalar
            %       s: sensor position --- 2 x 1 vector
            %OUTPUT:obj.d: measurement dimension --- scalar
            %       obj.h: function handle to generate measurement --- scalar
            %       obj.H: function handle to call measurement model Jacobian --- 1 x (state dimension) vector
            %       obj.R: measurement noise covariance --- scalar
            %NOTE: the measurement model assumes that in the state vector the first two entries are the X-position and Y-position, respectively.
            obj.d = 1;
            %Range
            rng = @(x) norm(x(1:2)-s);
            %Bearing
            obj.h = @(x) atan2(x(2)-s(2),x(1)-s(1));
            %Measurement model Jacobian
            obj.H = @(x) [-(x(2)-s(2))/(rng(x)^2) (x(1)-s(1))/(rng(x)^2) zeros(1, length(x)-2)];
            %Measurement noise covariance
            obj.R = sigma^2;
        end
        
        function obj = dualbearingmeasmodel(sigma, s1, s2)
            %DUALBEARINGMEASUREMENT creats the dual bearing measurement model
            %INPUT: sigma: standard deviation of measurement noise --- scalar
            %       s1: sensor position for sensor 1 --- 2 x 1 vector
            %       s2: sensor position for sensor 2 --- 2 x 1 vector
            %OUTPUT:obj.d: measurement dimension --- scalar
            %       obj.h: function handle to generate measurement --- 2 x 1 vector
            %       obj.H: function handle to call measurement model Jacobian --- 2 x (state dimension) vector
            %       obj.R: measurement noise covariance --- 2 x 2 matrix
            %NOTE: the measurement model assumes that in the state vector the first two entries are the X-position and Y-position, respectively.
            obj.d = 2;
            %Range
            rng1 = @(x) norm(x(1:2)-s1);
            rng2 = @(x) norm(x(1:2)-s2);
            %Bearing
            obj.h = @(x) [atan2(x(2)-s1(2),x(1)-s1(1));
                          atan2(x(2)-s2(2),x(1)-s2(1))];
            %Measurement model Jacobian
            obj.H = @(x) [-(x(2)-s1(2))/(rng1(x)^2) (x(1)-s1(1))/(rng1(x)^2) zeros(1, length(x)-2);
                            -(x(2)-s2(2))/(rng2(x)^2) (x(1)-s2(1))/(rng2(x)^2) zeros(1, length(x)-2)];
            %Measurement noise covariance
            obj.R = sigma^2*eye(2);
        end
        
        function obj = rangebearingmeasmodel(sigma_r, sigma_b, s)
            %RANGEBEARINGMEASUREMENT creats the range/bearing measurement model
            %INPUT: sigma_r: standard deviation of measurement noise added to range --- scalar
            %       sigma_b: standard deviation of measurement noise added to bearing --- scalar
            %       s: sensor position --- 2 x 1 vector
            %OUTPUT:obj.d: measurement dimension --- scalar
            %       obj.h: function handle to generate measurement --- 2 x 1 vector
            %       obj.H: function handle to call measurement model Jacobian --- 2 x (state dimension) vector
            %       obj.R: measurement noise covariance --- 2 x 2 matrix
            %NOTE: the measurement model assumes that in the state vector the first two entries are the X-position and Y-position, respectively.
            obj.d = 2;
            %Range
            rng = @(x) norm(x(1:2)-s);
            %Bearing
            ber = @(x) atan2(x(2)-s(2),x(1)-s(1));
            %Measurement vector
            obj.h = @(x) [rng(x);ber(x)];
            %Measurement model Jacobian
            obj.H = @(x) [
                          (x(1)-s(1))/rng(x)      (x(2)-s(2))/rng(x)     zeros(1, length(x)-2);
                          -(x(2)-s(2))/(rng(x)^2) (x(1)-s(1))/(rng(x)^2) zeros(1, length(x)-2)
                          ];
            %Measurement noise covariance
            obj.R = [sigma_r^2 0;
                     0         sigma_b^2];
        end
    end
end