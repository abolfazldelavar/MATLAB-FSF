
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef estimator < handle
    
    %% --------------------------------------------------------------------------------
    % Author: A. Delavar, http://abolfazldelavar.com
    %  --------------------------------------------------------------------------
    % INSTRUCTION:
    %  1) Copy 'chaosModel.m' file in folder '+system' and rename it arbitrary,
    %  2) change the class name the same the file name,
    %  3) Edit properties according to your system detail,
    %  4) Insert dynamic fix parameters in properties part in the file,
    %  5) Insert dynamic equations into the 'dynamics' function,
    %  6) Write your output codes into the 'external' function,
    %  7) If there is any state limitation, you can set them in 'limitations'
    %  8) Import the Jacobian into 'jacobian' function (Used in EKF),
    %  9) Copy the below code into 'initialization.m':
    %         models.estimator = caller.estimator(systems.chaosModel,            ...
    %                                             params.step, signals.tLine,    ...
    %                                             'approach' , 'ekf',            ...
    %                                             'initialCondition', params.x0_real);
    %  9) Change the 'Approach' in above code if you want other approaches,
    % 10) Copy one piece of the below codes into 'simulation' to update your estimator,
    %       For EKF:
    %           models.estimator.nextstepEKF(u(:,k), y(:,k));
    %           nextstepEKF(models.estimator, u(:,k), y(:,k));
    %       For UKF:
    %           models.estimator.nextstepUKF(u(:,k), y(:,k));
    %           nextstepUKF(models.estimator, u(:,k), y(:,k));
    % ---------------------------------------------------------------------------------
    
    properties(SetAccess=private)
        block               % Get a copy of your model class
        sampleTime          % Simulation sample time
        numSteps            % The number of all time steps
        numStates           % The number of states
        numOutputs          % The number of measurements
        numInputs           % The number of inputs
        states              % The state signals over the time
        inputs              % The input signals over the time
        outputs             % The measurement signals over the time
        timeLine            % A copy of simulation time line
        currentStep         % The current step of simulation
        estAproach          % The estimation approach (EKF, UKF, ...)
        
        % EKF variables
        qMatrix             % Q mateix - Dynamic noise variance
        rMatrix             % R matrix - Measurement noise variance
        initialStates       % Initial belief of states
        covariance          % Initial belief of states covariance
        
        % UKF variables
        nUKF                % This is usually = Number of states
        kappa               % A non-negative real number
        alpha               % a \in (0, 1]
        lambda              % Unscented KF parameter
        betta               % Unscented KF parameter
        wm                  % Unscented KF weights. Used in state
        wc                  % Unscented KF weights. Used in covariance
    end
    
    methods
    %% Initial function run to set default properties and make initial signals
        function obj = estimator(inputModel, sampTime, timeline, varargin)
            % [Save internal] <- (input model class, Sample time, Time line, Options)
            % Options are: 'InitialCondition', 'Approach'
            
            obj.block           = inputModel;
            obj.sampleTime      = sampTime;
            obj.timeLine        = timeline(:)';
            obj.numSteps        = numel(obj.timeLine);
            obj.numStates       = obj.block.numStates;
            obj.numInputs       = obj.block.numInputs;
            obj.numOutputs      = obj.block.numOutputs;
            obj.inputs          = zeros(obj.numInputs,  obj.numSteps);
            obj.outputs         = zeros(obj.numOutputs, obj.numSteps);
            obj.currentStep     = 0;
            obj.initialStates   = obj.block.initialStates(:);
            obj.states          = [obj.initialStates ,zeros(obj.numStates, obj.numSteps)];
            obj.covariance      = obj.block.covariance;
            obj.estAproach      = 'ekf';    % Default approach
            
            % Extracting the arbitraty value of properties
            for i = 1:2:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        % Choosing the approach of estimation. e.g., EKF, UKF, ...
                        case 'approach'
                            obj.estAproach   = varargin{i+1};
                        case 'initialCondition'
                            obj.states(:, 1) = varargin{i+1};
                    end
                end
            end
            
            % This part initialize the estimator by setting parameters
            switch obj.estAproach
                case 'ekf'
                    obj.qMatrix = obj.block.qMatrix;
                    obj.rMatrix = obj.block.rMatrix;
                case 'ukf'
                    obj.qMatrix = obj.block.qMatrix;
                    obj.rMatrix = obj.block.rMatrix;
                    obj.kappa   = obj.block.kappa;
                    obj.alpha   = obj.block.alpha;
                    % Dependent variables
                    obj.nUKF    = obj.numStates;
                    obj.lambda  = obj.alpha^2*(obj.nUKF + obj.kappa) - obj.nUKF;
                    obj.betta   = 2;
                    % Making weights
                    obj.wm      = ones(2*obj.nUKF + 1, 1)/(2*(obj.nUKF + obj.lambda));
                    obj.wc      = obj.wm;
                    obj.wc(1)   = obj.wm(1) + (1 - obj.alpha^2 + obj.betta);                          
            end
        end
        
        %% Next step of Extended Kalman Filter (EKF)
        function obj = nextstepEKF(obj, u, y)
            % [Internal update] <- (Internal, Input, Output)
            
            %% Initialize parameters
            % The current time is calculated as below
            currentTime = obj.timeLine(obj.currentStep + 1);
            % Preparing and saving inputs and outputs to internal
            u = u(:); y = y(:);
            obj.inputs(:, obj.currentStep + 1)  = u;
            obj.outputs(:, obj.currentStep + 1) = y;
            % Using dynamics of system to calculate Jacobians
            [A, L, H, M] = obj.block.jacobians(obj.states,          ...
                                               obj.inputs,          ...
                                               obj.currentStep + 1, ...
                                               obj.sampleTime,      ...
                                               currentTime);
            
            %% Prediction step - Update xp
            %  This part tries to obtain a prediction estimate from dynamic
            %  model of your system directly from nonlinear equations
            xm = obj.states(:, obj.currentStep + 1);
            % Set before-state-limitations
            xv  = obj.block.limitations(obj.states, 0);
            handleDyn = @(xx) obj.block.dynamics(xx,                  ...
                                                 obj.inputs,          ...
                                                 obj.currentStep + 1, ...
                                                 obj.sampleTime,      ...
                                                 currentTime)';
            if obj.block.timeType == 'c'
                xp = libraries.functionLib.dynamicRunner(handleDyn, xv, xm, obj.sampleTime, obj.block.solverType);
            else
                xp  = handleDyn(xv);
            end
            % Set after-state-limitations
            xp = obj.block.limitations(xp, 1);
            
            % Prediction step - Update covariance matrix
            Pp = A*obj.covariance*A' + L*obj.qMatrix*L';        % Update the prediction covariance matrix
            
            %% Posterior step - Reciving measurements
            %  If there is not any measurement (y == NaN), posterior won't
            %  calculate and just prediction will be reported.
            if ~isnan(y)
                K  = Pp*H'/(H*Pp*H' + M*obj.rMatrix*M');        % Kalman Gain
                xm = xp + K*(y - H*xp);                         % Update the states
                Pm = (eye(size(obj.covariance, 1)) - K*H)*Pp;   % Update covariance matrix
            else
                xm = xp;
                Pm = Pp;
            end
            
            %% Update internal signals
            obj.states(:, obj.currentStep + 2) = xm;            % To save estimated states
            obj.covariance                     = Pm;            % To save covariance matrix
            obj.goAhead;                                        % Go to the next step
        end
        
        %% Next step of Unscented Kalman Filter (UKF)
        function obj = nextstepUKF(obj, u, y)
            % [Internal update] <- (Internal, Input, Output)
            % ------------------------------------------------------
            % To see how this algorithm works, refer to below source:
            %   A. Delavar and R. R. Baghbadorani, "Modeling, estimation, and
            %   model predictive control for Covid-19 pandemic with finite
            %   security duration vaccine," 2022 30th International Conference
            %   on Electrical Engineering (ICEE), 2022, pp. 78-83,
            %   doi: 10.1109/ICEE55646.2022.9827062.
            % ------------------------------------------------------
            
            %% Initialize parameters, STEP 0 & 1
            % The current time is calculated as below
            currentTime = obj.timeLine(obj.currentStep + 1);
            % Preparing and saving inputs and outputs to internal
            u = u(:); y = y(:);
            obj.inputs(:, obj.currentStep + 1)  = u;
            obj.outputs(:, obj.currentStep + 1) = y;
            % Using dynamics of system to calculate Jacobians
            [~, L, ~, M] = obj.block.jacobians(obj.states,          ...
                                               obj.inputs,          ...
                                               obj.currentStep + 1, ...
                                               obj.sampleTime,      ...
                                               currentTime);
            % Getting last states prior and its covariance
            xm = obj.states(:, obj.currentStep + 1);
            Pm = obj.covariance;
            
            %% Solving sigma points, STEP 2
            dSigma = sqrt(obj.nUKF + obj.lambda)*chol(Pm)';  % Calculating sqrt
            xmCopy = xm(:, ones(1, numel(xm)));              % Putting 'xm' is some column (copy)
            sp     = [xm, xmCopy + dSigma, xmCopy - dSigma]; % Obtaining sigma points
            
            %% Prediction states and their covariance, STEP 3
            %  This part tries to obtain a prediction estimate from dynamic
            %  model of your system directly from nonlinear equations
            nSpoints = size(sp, 2);
            xp       = zeros(obj.numStates, 1);
            Xp       = zeros(obj.numStates, nSpoints);
            for i = 1:nSpoints
                changedFullState = obj.states;
                changedFullState(:, obj.currentStep + 1) = sp(:, i);
                
                % Set before-state-limitations
                xv  = obj.block.limitations(changedFullState, 0);
                handleDyn = @(xx) obj.block.dynamics(xx,                  ...
                                                     obj.inputs,          ...
                                                     obj.currentStep + 1, ...
                                                     obj.sampleTime,      ...
                                                     currentTime)';
                if obj.block.timeType == 'c'
                    Xp(:,i) = libraries.functionLib.dynamicRunner(handleDyn,        ...
                                                                  xv,               ...
                                                                  xm,               ...
                                                                  obj.sampleTime,   ...
                                                                  obj.block.solverType);
                else
                    Xp(:,i)  = handleDyn(xv);
                end
                % Set after-state-limitations
                Xp(:,i) = obj.block.limitations(Xp(:,i), 1);
                
                xp = xp + obj.wm(i)*Xp(:, i);                % Prediction update
            end
            dPp = Xp - xp(:, ones(1, nSpoints));
            Pp  = dPp*diag(obj.wc)*dPp' + L*obj.qMatrix*L';  % Updating the covariance of states matrix
            
            %% Updating sigma points, STEP 4
            %dSigma = sqrt(obj.nUKF + obj.lambda)*chol(Pp)';  % Calculating sqrt
            %xmCopy = xp(:, ones(1, numel(xp)));              % Putting 'xp' is some column (copy)
            %sp     = [xp, xmCopy + dSigma, xmCopy - dSigma]; % Updating sigma points
            
            if ~isnan(y)
                %% Solving output estimation using predicted data, STEP 5
                %  This part tries to obtain a prediction output from sigma points
                zb = zeros(obj.numOutputs, 1);
                Zb = zeros(obj.numOutputs, nSpoints);
                for i = 1:nSpoints
                    changedFullState = obj.states;
                    changedFullState(:, obj.currentStep + 1) = Xp(:, i); % Or 'Xp(:, i)' instead of 'sp(:, i)'
                    Zb(:,i) = obj.block.external(changedFullState,      ...
                                                 obj.inputs,            ...
                                                 obj.currentStep + 1,   ...
                                                 obj.sampleTime,        ...
                                                 currentTime)';
                    zb = zb + obj.wm(i)*Zb(:, i);                   % Predicted output
                end
                dSt = Zb - zb(:, ones(1, nSpoints));
                St  = dSt*diag(obj.wc)*dSt' + M*obj.rMatrix*M';     % Updating the covariance of output matrix

                %% Solving Kalman gain, STEP 6
                SiG = dPp*diag(obj.wc)*dSt';
                K   = SiG*(St^-1);                                  % Kalman Gain
            end
            
            %% Solving posterior using measurement data, STEP 7
            %  If there is not any measurement (y == NaN), posterior won't
            %  calculate and just prediction will be reported.
            if ~isnan(y)
                xm = xp + K*(y - zb);                           % Update the states
                Pm = Pp - K*SiG';                               % Update covariance matrix
            else
                xm = xp;
                Pm = Pp;
            end
            
            %% Update internal signals
            obj.states(:, obj.currentStep + 2) = xm;            % To save estimated states
            obj.covariance                     = Pm;            % To save covariance matrix
            obj.goAhead;                                        % Go to the next step
        end
        
        % This function can make a jump in the step variable
        % If no arguments are available, jump 1 step
        function obj = goAhead(obj, i)
            if nargin < 1 + 1
                i = 1;
            end
            obj.currentStep = obj.currentStep + i;
        end
        
        % Reset Block by changing the current step to zero
        function obj = goFirst(obj)
            obj.currentStep = 0;
        end
        
        % The below function is used to plot the internal signals
        function show(obj, sel, varargin)
            % To illustrate states, inputs, or outputs, you might have to
            % use some varargins which are explained in 'scope' class
            
            % The default signal to draw is the states
            if ~exist('sel', 'var'); sel = 'x'; end
            switch sel
                case 'x'
                    signal   = obj.states(:,1:obj.numSteps);
                    nSignals = obj.numStates;
                case 'y'
                    signal   = obj.outputs(:,1:obj.numSteps);
                    nSignals = obj.numOutputs;
                case 'u'
                    signal   = obj.inputs(:,1:obj.numSteps);
                    nSignals = obj.numInputs;
            end
            % Make a scope
            scp = caller.scope(obj.timeLine, nSignals, signal);
            scp.show(varargin);
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

