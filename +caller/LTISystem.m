
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef LTISystem < handle
    
    %% --------------------------------------------------------------------------------
    % Author: A. Delavar, http://abolfazldelavar.com
    %  --------------------------------------------------------------------------
    % INSTRUCTION
    % 1) Define your LTI system using 'tf', 'zpk', 'ss', ...
    % 2) use below code in 'initialization.m' to initial setting
    %        models.G = caller.LTISystem(The number of blocks, LTI_TF, Sample time, Initial condition);
    % 3) Use below order into 'simulation.m' to apply each step
    %        models.G.nextstep(Input Signal, xNoise, yNoise);
    %        nextstep(models.G, Input signal, xNoise, yNoise);
    % ---------------------------------------------------------------------------------
    
    properties (SetAccess=private)
        inputSystem     % Pure input system
        block           % Get a copy of your system class
        numberLTIs      % The number of LTI systems
        InitSampleTime  % Imported system sample time
        delay           % input delay of LTI system
        numStates       % The number of states
        numOutputs      % The number of measurements
        numInputs       % The number of inputs
        A               % Dynamic matrix A
        B               % Dynamic matrix B
        C               % Dynamic matrix C
        D               % Dynamic matrix D
        states          % The state signals over the time
        delayedInputs   % The input signals over the time
        outputs         % The measurement signals over the time
        sampleTime      % Simulation sample time
        currentStep     % The current step of simulation
    end
    
    methods
        function obj = LTISystem(nSystems, iSys, sampletime, initialcondition)
            % Number of systems, Block, Sample-time, Initial Condition
            obj.inputSystem = iSys;
            if iSys.Ts == sampletime && isStateSpace(iSys)
                %input is a discrete-time state-space with the same sample-time
                obj.block = iSys;
            elseif iSys.Ts == sampletime && ~isStateSpace(iSys)
                %input is a discrete-time transfer fuction with same sample time
                obj.block = minreal(ss(iSys));
            elseif isStateSpace(iSys)
                if iSys.Ts ~= 0
                    %input is a discrete-time state-space with different sample time   *
                    obj.block = d2d(iSys, sampletime);
                else
                    %input is a continuous-time state-space
                    obj.block = c2d(iSys, sampletime);
                end
            elseif iSys.Ts ~= 0
                %input is a discrete-time transfer fuction with different sample time  *
                obj.block = d2d(minreal(ss(iSys)), sampletime);
            else
                %input is a continuous-time transfer fuction
                obj.block = c2d(minreal(ss(iSys)), sampletime);
            end
            
            obj.delay           = obj.block.InputDelay;
            obj.numberLTIs      = nSystems;
            obj.sampleTime      = sampletime;
            obj.InitSampleTime  = iSys.Ts;
            obj.A               = obj.block.A;
            obj.B               = obj.block.B;
            obj.C               = obj.block.C;
            obj.D               = obj.block.D;
            obj.numStates       = size(obj.block.A, 1);
            obj.numInputs       = size(obj.block.B, 2);
            obj.numOutputs      = size(obj.block.C, 1);
            obj.delayedInputs   = zeros(obj.numInputs , obj.numberLTIs, obj.delay + 1);
            obj.outputs         = zeros(obj.numOutputs, obj.numberLTIs);
            obj.currentStep     = 0;
            
            % If the initial input does not exist, set it zero
            % Else, put the initial condition in the state matrix
            if nargin<4
                obj.states = zeros(obj.numStates, obj.numberLTIs);
            else
                if numel(initialcondition) == obj.numStates
                    initialcondition = initialcondition(:);
                    obj.states       = repmat(initialcondition, 1, obj.numberLTIs);
                else
                    if prod(size(initialcondition) == [obj.numStates obj.numberLTIs])
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
            % input at time t, additive noise on states, additive noise on output
            
            % If the noises do not exist, consider them zero.
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
            
            % Making delayed input signal
            obj.delayedInputs = circshift(obj.delayedInputs, -1, 3);
            obj.delayedInputs(:,:,end) = u;
            
            % Updating the states via dx = Ax + Bu
            x = obj.A*obj.states + obj.B*obj.delayedInputs(:,:,1);
            
            % Calculating outputs via y = Cx + Du
            y = obj.C*obj.states + obj.D*obj.delayedInputs(:,:,1);
            
            % Update internal signals which later can be used for plotting
            % and programming for other parts of the code
            obj.states  = x + xNoise;
            obj.outputs = y + yNoise;
            obj.goAhead;
        end
        
        % This function can make a jump in the step number variable
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
        
        % Here is the end of methods
    end
    % Here is the end of class
end

