
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef nonlinearSystem < handle
    
    %% --------------------------------------------------------------------------------
    % Author: A. Delavar, http://abolfazldelavar.com
    %  --------------------------------------------------------------------------
    % INSTRUCTION
    % 1) Copy 'chaos.m' file in folder '+system' and rename it arbitrary,
    % 2) change the class name the same the file name,
    % 3) Edit properties according to your system detail,
    % 4) Insert dynamic equations into 'dynamics' function,
    % 5) Write your output codes into 'external' function
    % 6) use below code in 'initialization.m' to set initial options
    %        models.chaos = caller.nonlinearSystem(systems.chaos, params.step, signals.tLine, initialStates);
    % 7) Use below code in 'simulation.m' to apply each step
    %        models.chaos.nextstep(inputSignal, xNoise, yNoise);
    %        nextstep(models.chaos, inputSignal, xNoise, yNoise);
    % --------------------------------------------------------------------------------- 
    
    properties (SetAccess=private)
        block               % Get a copy of your system class
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
        solverType          % The type of dynamic solver
    end
    
    methods
        function obj = nonlinearSystem(inputsystem, sampletime, timeline, initialcondition)
            % system, sample time, time line, initial states condition
            obj.block           = inputsystem;
            obj.sampleTime      = sampletime;
            obj.timeLine        = timeline(:)';
            obj.numSteps        = numel(obj.timeLine);
            obj.numStates       = obj.block.numStates;
            obj.numInputs       = obj.block.numInputs;
            obj.numOutputs      = obj.block.numOutputs;
            obj.solverType      = obj.block.solverType;
            obj.inputs          = zeros(obj.numInputs, obj.numSteps);
            obj.outputs         = zeros(obj.numOutputs, obj.numSteps);
            obj.currentStep     = 0;
            
            % If the initial input does not exist, set it zero
            % Else, put the initial condition in state matrix
            if nargin<4
                obj.states = zeros(obj.numStates, obj.numSteps + 1);
            else
                initialcondition = initialcondition(:);
                obj.states = [initialcondition ,zeros(obj.numStates, obj.numSteps)];
            end
        end
        
        % The 'nextstep' function can provide an easy way to call 
        % dydnamics of the system to calculate next sample states
        % To use this function, refer to the top INSTRUCTION part
        function obj = nextstep(obj, u, xNoise, yNoise)
            % input at time t, additive noise on states, additive noise on output
            
            % If the noises do not exist, consider them zero.
            if nargin < 4
                if nargin < 3
                    xNoise = 0;
                end
                yNoise = 0;
            end
            
            % The current time is calculated as below
            currentTime = obj.timeLine(obj.currentStep + 1);
            
            % Preparing the input signal and save to the internal array
            obj.inputs(:, obj.currentStep + 1) = u(:);
            
            % Set before-state-limitations:
            % This can be used if we want to process on states before
            % calculating the next states by dynamics.
            xv = obj.block.limitations(obj.states, 0);
            xo = obj.states(:, obj.currentStep + 1); % Getting the previous states
            
            % The below handle function is used in the following
            handleDyn = @(xx) obj.block.dynamics(xx,                  ...
                                                 obj.inputs,          ...
                                                 obj.currentStep + 1, ...
                                                 obj.sampleTime,      ...
                                                 currentTime)';
                                                 
            % This part calculates the states and outputs using the system dynamics
            if obj.block.timeType == 'c'
                % The type of solver can be under your control
                % To change your solver type, do not change any code here
                % Change the solver type in 'chaos.m' file or others
                x = libraries.functionLib.dynamicRunner(handleDyn, xv, xo, obj.sampleTime, obj.solverType);
            else
                % When the inserted system is discrete time, just the
                % dynamic must be solved as below
                x = handleDyn(xv);
            end
            
            % Set after-state-limitations
            x = obj.block.limitations(x, 1);
            
            % The output of the system is solved by the measurement
            % dynamics of the system which are available in 'chaos.m' file
            y = obj.block.external(obj.states,          ...
                                   obj.inputs,          ...
                                   obj.currentStep + 1, ...
                                   obj.sampleTime,      ...
                                   currentTime)';
            
            % Update internal signals which later can be used for plotting
            % and programming for other parts of the code
            obj.states(:, obj.currentStep + 2)  = x + xNoise(:);
            obj.outputs(:, obj.currentStep + 1) = y + yNoise(:);
            obj.goAhead;
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

