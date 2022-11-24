
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef functionLib
    % Default functions, neccessary and essential functions that have been
    % provided to be used in projects. To use, you can call 'lib.func' in
    % anywhere you need.
    
    methods (Hidden=true)
        %% The below function print a comment every one second automatically
        function trig = disit(~, k, n, trig, HowSecond)
            % Trig: [Trigger, Previous k]
            % HowSecond is the span between each report (second)
            nowt      = toc;
            if floor(nowt) >= trig(1)
                trig(1) = trig(1) + HowSecond;     % To make a horizon
                mea     = k - trig(2);             % Calculatation of the mean
                trig(2) = k;
                trem    = (n - k)*HowSecond/mea;   % Remained time in second
                tremmin = min(floor(trem/60),1e5); % How minutes
                tremsec = round(rem(trem, 60));    % Remained seconds
                comPerc = round(k/n*100);          % Completed percentage
                
                % Making a graphical text
                txt     = repmat(' ', 1, 65);
                k       = num2str(k);
                txt(7 + (1:numel(k)) - floor(0.5*numel(k))) = k;
                txt(16) = '|';
                n       = num2str(n);
                txt(21 + (1:numel(n)) - floor(0.5*numel(n))) = n;
                txt(28) = '|';
                comPerc = num2str(comPerc);
                txt(34 + (1:numel(comPerc))) = comPerc;
                txt(43) = '|';
                tremmin = num2str(tremmin);
                txt(52 - (numel(tremmin):-1:1)) = tremmin;
                txt(52) = ':';
                tremsec = num2str(tremsec);
                txt(52 + (1:numel(tremsec))) = tremsec;
                disp(txt);
            end
        end
        
        %% Saying a hello at the start of the simulation
        function sayStart(obj)
            % Initialize
            params = valuation;
            % Check the directory and start a diary to save command window
            if ~isfolder(params.diaryDir); mkdir(params.diaryDir) ;end
            if params.makeDiary == true
                diary([params.diaryDir, '/', obj.getNow(1, '-'),'.txt']);
            end
            % Print the header
            disp('Faryadell Simulation Framework (FSF) - Version 1.0.0');
            disp(['The simulation has kicked off! (', obj.getNow(5,'/'), ', ', obj.getNow(6,':'), ')']);
            % Below codes make the header table
            disp(repmat('-', 1, 65));
            txt = repmat(' ', 1, 65);
            txt(3:63) = 'Current step | All steps | Progress (%) | Remained time (m:s)';
            disp(txt);
            disp(repmat('-', 1, 65));
            % The below order start a timer to obtain the time which this
            % simulation lasted.
            tic;
        end
        
        %% Report the simulation time when it finishes
        function sayEnd(~)
            % Stopping the timer which started before, and saving the time.
            ntime   = toc;
            tmin    = min(floor(ntime/60), 1e5); % How minutes
            tsec    = round(rem(ntime, 60), 4);  % Remained seconds
            disp(['_' char([7424, 7427, 7439, 7436, 4325, 7424, 4301, 7436]) ' ', ...
                      char([7429, 7431, 7436, 7424, 7456, 7424, 7450]), repmat('_', 1, 48)]);
            disp(['The simulation has been completed. (', num2str(tmin), ' minutes and ', num2str(tsec), ' seconds)']);
        end
        
        %% Returning a text contained date and time
        function txt = getNow(~, typeReport, splitchar)
            % [Output string] <- (Internal, They output style controller, Splitter)
            
            if ~exist('splitchar', 'var'); splitchar = '_'; end
            
            fullDate = clock; % Getting full information of time
            year     = num2str(fullDate(1));
            month    = num2str(fullDate(2));
            day      = num2str(fullDate(3));
            hour     = num2str(fullDate(4));
            minute   = num2str(fullDate(5));
            second   = num2str(round(fullDate(6)));
            
            % If the numbers are less than 10, add a zero before them
            if numel(month)<2 ; month = ['0', month]; end
            if numel(day)<2 ; day = ['0', day]; end
            if numel(hour)<2 ; hour = ['0', hour]; end
            if numel(minute)<2 ; minute = ['0', minute]; end
            if numel(second)<2 ; second = ['0', second]; end
            
            % If ipuut does not axist, set a value (=0)
            if nargin < 1 + 1; typeReport = 0; end
            
            % What style do you want? You can change arbitrary
            switch typeReport
                case 1
                    txt = [year, month, day, splitchar, hour, minute, second];
                case 2
                    txt = [year, month, day, splitchar, hour, minute, second];
                case 3
                    txt = [year, splitchar, month, splitchar, day, splitchar, hour, splitchar, minute, splitchar, second];
                case 4
                    txt = [year, month, day, splitchar, hour, minute];
                case 5
                    txt = [year, splitchar, month, splitchar, day];
                case 6
                    txt = [hour, splitchar, minute, splitchar, second];
                otherwise
                    txt = [year, month, day, hour, minute, second];
            end
        end
        
        % End of the methods
    end
    
    methods (Access=public)
        %% Delayed in a signal
        function y = delayed(~, u, k, pdelay)
            % (Signal, Current sample time, Delay number)
            if k - pdelay > 0
                y = u(k - pdelay);
            else
                y = 0.*u(k);
            end
        end
        
        %% Delete all files in the diary folder
        function clrLogs(~)
            params = valuation;
            delete([params.diaryDir, '/*.txt']);
        end
        
        %% Signal Generator
        function output = signalMaker(~, params, Tline)
            % SETUP -------------------------------------------------------------------------
            % 1) Insert the below code in 'setParameters.m' to use signal generator:
            %        %% Signal Generator parameters
            %        params.referAddType  = 'onoff';    % 'none', 'square', 'sin', 'onoff', ...                        
            %        params.referAddAmp   = 30;         % Amplitude of additive signal
            %        params.referAddFreq  = 2;          % Signal period at simulation time
            %
            % 2) Use the below code to call signal generation in 'modelMaker.m':
            %        signal = func.signalMaker(params, models.Tline);
            % -------------------------------------------------------------------------------
            
            Tline = Tline(1:params.n);
            switch params.referAddType
                case 'none'
                    output  = 0.*Tline;
                case 'square'
                    freq    = params.referAddFreq/params.Tout;
                    output  = square(2*pi*freq*Tline) * params.referAddAmp;
                case 'onoff'
                    freq    = params.referAddFreq/params.Tout;
                    output  = (2 - square(2*pi*freq*Tline) - 1)/2 * params.referAddAmp;
                case 'sin'
                    freq    = params.referAddFreq/params.Tout;
                    output  = sin(2*pi*freq*Tline) * params.referAddAmp;
            end
        end
        
        %% Making a signal of an exponential inverse
        function output = expinverse(~, Tline, bias, alph, areaa)
            % [Output signal] <- (Time line, The bias point, smoother, Signal domain)
            % Note that signal domain is a two component vector [a(1), a(2)]
            output = 1./(1 + exp(-alph*(Tline - bias)));
            output = (areaa(2) - areaa(1))*output + areaa(1);
        end
        
        %% Making an exponential signal
        function output = exponen(~, Tline, para)
            % [Output signal] <- (Time line, [Startpoint, Final point, Rate])
            Ss = para(1);        % Start point
            Sf = para(2);        % Final value
            Sr = para(3);        % Rate of unction
            output = (Ss - Sf)*exp(-Sr.*Tline) + Sf;            
        end
        
        %% Saturation a signal in a band area [band(1), band(2)]
        function y = satutation(~, u, band)
            if u > band(2)
                y = band(2);
            elseif u < band(1)
                y = band(1);
            else
                y = u;
            end
        end
        
        %% Linear mapping a number
        function output = lineMapping(~, x, from, to)
            % Mapping 'x' from band [w1, v1] to band [w2, v2]
            w1     = from(1);
            w2     = from(2);
            v1     = to(1);
            v2     = to(2);
            output = 2*((x - w1)/(w2 - w1)) - 1;
            output = (output + 1)*(v2 - v1)/2 + v1;
        end

        % End of the methods
    end
    
    
    methods (Static=true, Hidden=true)
        %% Dynamic solver
        function xtr = dynamicRunner(handleDyn, xv, xo, sTime, solverType)
            % [New states] <- Dynamic unction, Full-time states, Current states, Solver type
            switch solverType
                case 'Euler'
                    % Euler method properties is given below (T is sample time):
                    %   x(t+1) = x(t) + T*f(x(t))
                    xtr = xo + sTime*handleDyn(xv);

                case 'Runge'
                    % 4th oder of 'Runge Kutta' is described below (T is sample time):
                    %   K1     = T*f(x(t))
                    %   K2     = T*f(x(t) + K1/2)
                    %   K3     = T*f(x(t) + K2/2)
                    %   K4     = T*f(x(t) + K3)
                    %   x(t+1) = x(t) + 1/6*(K1 + 2*K2 + 2*K3 + K4)

                    K1 = sTime*handleDyn(xv);
                    K2 = sTime*handleDyn(xv + K1/2);
                    K3 = sTime*handleDyn(xv + K2/2);
                    K4 = sTime*handleDyn(xv + K3);
                    xtr  = xo + 1/6*(K1 + 2*K2 + 2*K3 + K4);

                otherwise
                    error(['The solver name is not correct, please change the word "', obj.solverType, '"']);
            end
        end
    end
    % End of class
end

