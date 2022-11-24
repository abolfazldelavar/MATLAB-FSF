
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef chaos %Nonlinear Dynamic System
    
    % --------------------------------------------------------------------------
    % --- INSTRUCTION -------------------------------------------------------
    % INSTRUCTION
    % 1) Copy 'chaos.m' file in folder '+system' and rename it arbitrary,
    % 2) change the class name the same the file name,
    % 3) Edit properties according to your system detail,
    % 4) Insert dynamic equations into 'dynamics' function,
    % 5) Write your output codes into 'external' function,
    % 6) If there is any state limitation, you can set them in 'limitations'
    % 7) use below code in 'initialization.m' to set initial options
    %        models.chaos = caller.nonlinearSystem(systems.chaos, Sample-time , Time Line, Initial condition);
    % 8) Use one of piece of below codes in 'simulation.m' to apply each step
    %        models.chaos.nextstep(Input Signal, xNoise, yNoise);
    %        nextstep(models.chaos, Input Signal, xNoise, yNoise);
    % --------------------------------------------------------------------------

    properties
        numStates   = 3         % Number of states
        numInputs   = 1         % Number of inputs
        numOutputs  = 2         % Number of outputs
        timeType    = 'c';      % 'c' -> Continuous, 'd' -> Discrete
        solverType  = 'Euler';  % 'Euler', 'Runge'
        
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
        
        % Here is the end of methods
    end
    % Here is the end of class
end

