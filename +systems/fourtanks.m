
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef fourtanks %Nonlinear Dynamic System

    % --------------------------------------------------------------------------
    % --- INSTRUCTION -------------------------------------------------------
    % INSTRUCTION
    % 1) Copy 'chaos.m' file in folder '+system' and rename it arbitrary,
    % 2) change the class name the same the file name,
    % 3) Edit properties according to your system detail,
    % 4) Insert dynamic equations into 'dynamics' function,
    % 5) Write your output codes into 'external' function
    % 6) If there is any state limitation, you can set them in 'limitations'
    % 7) use below code in 'initialization.m' to set initial options
    %        models.chaos = caller.nonlinearSystem(systems.chaos, Sample-time , Time Line, Initial condition);
    % 8) Use one of piece of below codes in 'simulation.m' to apply each step
    %        models.chaos.nextstep(Input Signal, xNoise, yNoise);
    %        nextstep(models.chaos, Input Signal, xNoise, yNoise);
    % --------------------------------------------------------------------------

    % This model is given from the below source ----------- >
    % K. H. Johansson, "The quadruple-tank process: a multivariable laboratory
    % process with an adjustable zero," in IEEE Transactions on Control Systems Technology,
    % vol. 8, no. 3, pp. 456-465, May 2000, doi: 10.1109/87.845876

    properties
        numStates   = 4         % Number of states
        numInputs   = 2         % Number of inputs
        numOutputs  = 2         % Number of outputs - the tanks 1 and 2 are measured
        timeType    = 'c';      % 'c' -> Continuous, 'd' -> Discrete
        solverType  = 'Euler';  % 'Euler', 'Runge'
        
        % Other variables
        a1  = 0.071;     % cm^2
        a2  = 0.057;     % cm^2
        a3  = 0.071;     % cm^2
        a4  = 0.057;     % cm^2
        A1  = 28;        % cm^2
        A2  = 32;        % cm^2
        A3  = 28;        % cm^2
        A4  = 32;        % cm^2
        g   = 981;       % cm/s^2
        k1  = 3.33;      % cm^3/Vs --- (3.14) is also possible
        k2  = 3.35;      % cm^3/Vs --- (3.29) is also possible
        ga1 = 0.7;       % (0.43) is also possible
        ga2 = 0.6;       % (0.34) is also possible
        kc  = 0.5;       % V/cm
    end
    
    methods
        %% This part is internal dynamic functions that represents
        %  internal relations between states and inputs
        %  ~~> dx = f(x,u)
        function dx = dynamics(obj, x, u, k, ~, ~)
            % Parameters, States, Inputs, Current step, Sample-time, Current time
            dx(1) = -obj.a1/obj.A1*sqrt(2*obj.g*x(1, k)) + obj.a3/obj.A1*sqrt(2*obj.g*x(3, k)) + obj.ga1*obj.k1/obj.A1*u(1, k);
            dx(2) = -obj.a2/obj.A2*sqrt(2*obj.g*x(2, k)) + obj.a4/obj.A2*sqrt(2*obj.g*x(4, k)) + obj.ga2*obj.k2/obj.A2*u(2, k);
            dx(3) = -obj.a3/obj.A3*sqrt(2*obj.g*x(3, k)) + (1 - obj.ga2)*obj.k2/obj.A3*u(2, k);
            dx(4) = -obj.a4/obj.A4*sqrt(2*obj.g*x(4, k)) + (1 - obj.ga1)*obj.k1/obj.A4*u(1, k);
        end
        
        %% Measurement functions 
        %  ~~> y = g(x,u)
        function y = external(obj, x, ~, k, ~, ~)
            % Parameters, States, Inputs, Current step, Sample-time, Current time
            y(1) = obj.kc*x(1, k);
            y(2) = obj.kc*x(2, k);
        end
        
        %% All limitations before and after the state updating
        %  It can be useful for systems which have rules
        function x = limitations(~, x, mode)
            % Obj, States, Mode
            switch mode
                case 0  % before updating states
                    
                case 1  % After updating states
                    x = max(x,0);
            end
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

