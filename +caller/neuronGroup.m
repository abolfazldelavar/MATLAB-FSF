
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef neuronGroup < handle
    
    %% --------------------------------------------------------------------------------
    % Author: A. Delavar, http://abolfazldelavar.com
    %  --------------------------------------------------------------------------
    % INSTRUCTION
    % 1) Copy 'Izhikevich.m' file in folder '+system' and rename it arbitrary,
    % 2) change the class name the same as the file name,
    % 3) Edit properties according to your model details,
    % 4) Insert dynamic equations into 'dynamics' function,
    % 5) Write your part of output codes into the 'external' function
    % 6) use the below code in 'initialization.m' to set initial options
    %        models.neuronNet = caller.neuronGroup(params.nNeurons, systems.Izhikevich, params.step, initialCondition);
    % 7) Use one of the below piece of code existed in 'simulation.m' to apply each step
    %        models.neuronNet.nextstep(inputSignal, xNoise, yNoise);
    %        nextstep(models.neuronNet, inputSignal, xNoise, yNoise);
    % --------------------------------------------------------------------------------- 
    
    properties(SetAccess=private)
        block               % Get a copy of your system class
        numberNeurons       % The number of neurons (Network size)
        sampleTime          % The simulation sample-time
        numStates           % The number of states
        numOutputs          % The number of measurements
        numInputs           % The number of inputs
        states              % State signals
        inputs              % Input signals
        outputs             % Measurement signals
        currentStep         % The current step of simulation
        solverType          % The type of dynamic solver
    end
    
    methods
        function obj = neuronGroup(nNeurons, inputsystem, sampletime, initialcondition)
            % the number of neurons, dynamic class, sample time, initial condition
            obj.block           = inputsystem;
            obj.sampleTime      = sampletime;
            obj.numberNeurons   = nNeurons;
            obj.numStates       = obj.block.numStates;
            obj.numInputs       = obj.block.numInputs;
            obj.numOutputs      = obj.block.numOutputs;
            obj.solverType      = obj.block.solverType;
            obj.inputs          = zeros(obj.numInputs, obj.numberNeurons);
            obj.outputs         = zeros(obj.numOutputs, obj.numberNeurons);
            obj.currentStep     = 0;
            
            % If the initial input does not exist, set it zero
            % Else, put the initial condition in the state matrix
            if nargin<4
                obj.states = zeros(obj.numStates, obj.numberNeurons);
            else
                if numel(initialcondition) == obj.numStates
                    initialcondition = initialcondition(:);
                    obj.states       = repmat(initialcondition, 1, obj.numberNeurons);
                else
                    if prod(size(initialcondition) == [obj.numStates obj.numberNeurons])
                        obj.states = initialcondition;
                    else
                        error('The dimential of initial value that inserted is wrong. Check it please.');
                    end
                end
            end
        end
        
        % The 'nextstep' function can provide an easy way to call 
        % dydnamics of the system to calculate next sample states
        % To use this function, refer to the top INSTRUCTION part
        function obj = nextstep(obj, u, xNoise, yNoise)
            % this object, input at time t, additive noise on states, additive noise on output
            
            % If the noise signals do not exist, set them zero.
            if nargin < 4
                if nargin < 3
                    xNoise = 0;
                end
                yNoise = 0;
            end
            
            % Preparing the input signal before do any further action
            if size(u,1) ~= obj.numInputs
                u = u';
            end
            
            % Set before-state-limitations:
            % This can be used if we want to process on states before
            % calculating the next states by dynamics.
            x = obj.block.limitations(obj.states, 0);
            
            % The below handle function is used in the following
            handleDyn = @(xx) obj.block.dynamics(xx, u);
            
            % This part calculates the states and outputs using the system dynamics
            if obj.block.timeType == 'c'
                % The type of solver can be under your control
                % To change your solver, do not change any code here
                % Change the solver type in 'Izhikevich.m' file or others
                x = libraries.functionLib.dynamicRunner(handleDyn, x, x, obj.sampleTime, obj.solverType);
            else
                % When the inserted system is discrete time, just the
                % dynamic must be solved as below
                x = handleDyn(x);
            end
            
            % Set after-state-limitations
            x = obj.block.limitations(x, 1);
            
            % The output of the system is solved by the measurement
            % dynamics of the system which are available in 'Izhikevich.m' file
            y = obj.block.external(x, u);
            
            % Updating internal signals
            obj.states  = x + xNoise;
            obj.outputs = y + yNoise;
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

