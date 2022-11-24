
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef Izhikevich %Nonlinear Dynamics of the system
    
    % --------------------------------------------------------------------------
    % --- INSTRUCTION -------------------------------------------------------
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
    % --------------------------------------------------------------------------
    
    properties
        numStates   = 2         % Number of states
        numInputs   = 1         % Number of inputs
        numOutputs  = 1         % Number of outputs
        timeType    = 'c';      % 'c' -> Continuous, 'd' -> Discrete
        solverType  = 'Euler';  % 'Euler', 'Runge'
        
        % Other variables
        a           = 0.1       % Time scale of the recovery variable
        b           = 0.2   	% Sensitivity of the recovery variable to the sub-threshold fluctuations of the membrane potential
        c           = -65    	% After-spike reset value of the membrane potential
        d           = 2       	% After-spike reset value of the recovery variable
        timescale   = 1e3;
    end
    
    methods
        %% This part is internal dynamic functions that represents
        %  internal relations between states and inputs
        %  ~~> dx = f(x,u)
        function dx = dynamics(obj, x, I)
            % Parameters, states, inputs
            dx      = zeros(2, size(x,2));
            dx(1,:) = obj.timescale*(0.04*x(1,:).^2 + 5*x(1,:) - x(2,:) + 140 + I);
            dx(2,:) = obj.timescale*(obj.a*(obj.b*x(1,:) - x(2,:)));
        end
        
        %% Measurement functions 
        %  ~~> y = g(x,u)
        function y = external(~, x, ~)
            % Parameters, states, inputs
            y(1,:) = x(1,:);
        end
        
        %% All limitations before and after the state updating
        %  It can be useful for systems which have rules
        function x = limitations(obj, x, mode)
            % Obj, States, Mode
            switch mode
                case 0  % before updating states
                    ind    = (x(1,:) == 30);
                    x(1,:) = x(1,:).*(~ind) + ind*obj.c;
                    x(2,:) = x(2,:) + ind*obj.d;
                case 1  % After updating states
                    x(1,:) = min(x(1,:), 30);
            end
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

