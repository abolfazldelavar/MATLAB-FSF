
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef chaosModel %Nonlinear Dynamic model

    % --------------------------------------------------------------------------
    % --- INSTRUCTION -------------------------------------------------------
    % INSTRUCTION
    % 1) Copy 'chaosModel.m' file in folder '+system' and rename it arbitrary,
    % 2) change the class name the same the file name,
    % 3) Edit properties according to your system detail,
    % 4) Inset dynamic fixed parameters in properties part in this file,
    % 5) Insert dynamic equations into 'dynamics' function,
    % 6) Write your output codes into 'external' function,
    % 7) If there is any state limitation, you can set them in 'limitations'
    % 8) Import the Jacobian into 'jacobian' function (Used in EKF and UKF),
    % 9) To use this model as an estimator, please refer to 'estimator' file
    %    Or if you have other ideas such as MPC controller, refer to related
    %    files in '+ru' directory.
    % --------------------------------------------------------------------------

    properties
        numStates   = 3         %Number of states
        numInputs   = 1         %Number of inputs
        numOutputs  = 2         %Number of outputs
        timeType    = 'c';      %c -> Continuous  ,  d -> Discrete
        solverType  = 'Euler';  % 'Euler', 'Runge'
        
        % EXTENTED KALMAN FILTER --------
        initialStates = ones(3,1);          % Initial value of states
        covariance    = 1e+3*eye(3);        % Covariance of states
        qMatrix       = eye(3)*1e0;         % Dynamic noise variance
        rMatrix       = eye(2)*1e0;         % Measurement noise variance
        
        % UNSKENTED KALMAN FILTER -------
        % Note that 'Extended KF' parameters is also useful for Unscented KF, 
        % so just put your values there.
        kappa         = 10;                 % A non-negative real number
        alpha         = 0.6;                % a \in (0, 1]
        
        % Other variables
        sigma = 10;
        ro    = 28;
        beta  = 8/3;
    end
    
    methods
        %% This part is internal dynamic functions that represents
        %  internal relations between states and inputs
        %  ~~> dx = f(x,u)
        function dx = dynamics(obj, x, u, k, ~, ~)
            % Parameters, States, Inputs, Current step, Sample-time, Current time
            dx(1) = obj.sigma*x(2, k) - obj.sigma*x(1, k) + u(1, k);
            dx(2) = obj.ro*x(1, k) - x(1, k)*x(3, k) - x(2, k);
            dx(3) = x(1, k)*x(2, k) - obj.beta*x(3, k);
        end
        
        %% Measurement functions 
        %  ~~> y = g(x,u)
        function y = external(~, x, ~, k, ~, ~)
            % Parameters, States, Inputs, Current step, Sample-time, Current time
            y(1) = x(1, k);
            y(2) = x(2, k);
        end
        
        %% All limitations before and after the state updating
        %  It can be useful for systems which have rules
        function x = limitations(~, x, mode)
            % Obj, States, Mode
            switch mode
                case 0  % before updating states
                    
                case 1  % After updating states
                    
            end
        end
        
        %% Jacobians
        function [A, L, H, M] = jacobians(obj, x, ~, k, ~, ~)
            % [A, L, H, M] <- (Parameters, States, Inputs, Current step, Sample-time, Current time)
            
            % INSTRUCTION:
            %   dx = Ax + Lw,     'x' is states and 'w' denotes process noise
            %   y  = Hx + Mv,     'x' is states and 'v' is measurement noise
            
            % A matrix, d(q(t))/dx(t)
            A      = zeros(3, 3);
            A(1,1) = -obj.sigma;
            A(1,2) = +obj.sigma;
            A(1,3) = 0;
            A(2,1) = obj.ro - x(3, k);
            A(2,2) = -1;
            A(2,3) = -x(1, k);
            A(3,1) = x(2, k);
            A(3,2) = x(1, k);
            A(3,3) = -obj.beta;
            % L matrix, d(q(t))/dw(t), Process noise effects
            L = eye(3);
            % H matrix, d(h(t))/dx(t)
            H = eye(2, 3);
            % M matrix, d(h(t))/dv(t), Measurement Noise effects
            M = eye(2);
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

