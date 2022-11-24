
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef scope < handle
    
    %% --------------------------------------------------------------------------------
    % Author: A. Delavar, http://abolfazldelavar.com
    %  --------------------------------------------------------------------------
    % This class is a kind of scope that can save your signals which can be
    % used for after running progress.
    % INSTRUCTION:
    % 1) use below code in 'initialization.m' to set initial options
    %        signals.V = caller.scope(Time line, The number of signals);
    % 2) Use one of the below piece of code existed in 'simulation.m' to save each step
    %        signals.V.getdata(Input signal at step k, Additive normal noise VARIANCE);
    %        getdata(signals.V, Input signal at step k, Additive normal noise VARIANCE);
    % --------------------------------------------------------------------------------- 
    
    properties(SetAccess=private)
        sampleTime          % Simulation sample-time
        timeLine            % Time line vector
        numSignals          % The number of signals
        signals             % Signal matrix
        n                   % The number of time steps
        currentStep         % The current step of simulation
    end
    
    methods
        function obj = scope(tLine, nSignals, initialcondition)
            % Time line, The number of signals, initial signal value
            obj.timeLine  	= tLine(:)';
            obj.sampleTime	= mean(tLine(2:end) - tLine(1:end-1));
            obj.numSignals  = nSignals;
            obj.currentStep = 0;
            obj.n         	= numel(obj.timeLine);
            
            % If the initial input does not exist, set it zero
            % Else, put the initial condition in state matrix
            if nargin<3
                obj.signals = zeros(obj.numSignals, obj.n);
            else
                if numel(initialcondition) == obj.numSignals
                    initialcondition = initialcondition(:);
                    obj.signals      = repmat(initialcondition, 1, obj.n);
                else
                    if prod(size(initialcondition) == [obj.numSignals obj.n])
                        obj.signals  = initialcondition;
                    else
                        error('The dimential of initial value that inserted is wrong. Check it please.');
                    end
                end
            end
        end
        
        % The 'getdata' function can receive value and save it.
        % To use this function, refer to the top INSTRUCTION part
        function obj = getdata(obj, insData, addNoiseVar)
            % Inserted data, additive noise Variance
            
            % If the noise signals do not exist, consider them zero.
            if nargin < 3
                noiseSig = 0;
            else
                noiseSig = addNoiseVar(:).*randn(obj.numSignals, 1);
            end
            
            % Preparing the imported data
            insData = insData(:);
            
            % Update internal signals which later can be used for plotting
            % and programming for other parts of the code
            obj.signals(:, obj.currentStep + 1) = insData + noiseSig;
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
        
        % The below function is used to plot the internal signals
        function show(obj, varargin)
            % To set your opyions you must use the below codes:
            %  'select' -> 1:3
            %  'ylabel' -> 'Voltage (v)'
            
            %Initialize
            plt = libraries.plotToolLib;
            if numel(varargin) == 1 && iscell(varargin{1}); varargin = varargin{1}; end
            
            % Default ylabel text and signals are set here. Indeed, all
            % signals are plotted in default mode.
            xtit = '$x$';
            ytit = '$y$';
            ztit = '$z$';
            sel = 1:obj.numSignals;
            
            % Extracting the arbitraty values of properties
            for i = 1:2:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'select'
                            sel      = varargin{i+1};
                        case 'derive'
                            derIndex = varargin{i+1};
                        case 'notime'
                            nosigs   = varargin{i+1};
                        case 'save'
                            fName    = varargin{i+1};
                        case 'ylabel'
                            ytit_c   = varargin{i+1};
                        case 'xlabel'
                            xtit_c   = varargin{i+1};
                        case 'zlabel'
                            ztit_c   = varargin{i+1};
                        case 'title'
                            ftitle   = varargin{i+1};
                        case 'legend'
                            flegend  = varargin{i+1};
                    end
                end
            end

            if ~exist('nosigs', 'var')
                % Pre-processing on data
                SIGNAL = obj.signals(sel,1:end-1);
                % If all derivatives are requested, make them here
                if exist('derIndex', 'var')
                    if derIndex == true
                        SIGNAL = SIGNAL - circshift(SIGNAL, +1, 2);
                    else
                        SIGNAL = SIGNAL - diag(derIndex>0)*circshift(SIGNAL, -1, 2);
                    end
                    SIGNAL(:,1)   = SIGNAL(:,2);
                    SIGNAL(:,end) = SIGNAL(:,end-1);
                end
                % General time line plot
                h = figure();
                    plot(obj.timeLine(1:end-1), SIGNAL);
                    xlabel('Time (s)');
                    if exist('xlabel', 'var'); xlabel(xtit); end
                    ylabel(ytit);
            else
                if size(nosigs, 2) == 2
                    % 2D signal line plot
                    h = figure();
                    hold on;
                    for i = 1:size(nosigs,1)
                        % Pre-processing on data
                        SIGNAL = obj.signals(nosigs(i,:),1:end-1);
                        % If all derivatives are requested, make them here
                        if exist('derIndex', 'var')
                            if derIndex == true
                                SIGNAL = SIGNAL - circshift(SIGNAL, +1, 2);
                                xtit = 'd$x/$d$t$';
                                ytit = 'd$y/$d$t$';
                            else
                                SIGNAL = SIGNAL - diag(derIndex>0)*circshift(SIGNAL, +1, 2);
                                if derIndex(1) ~= 0; xtit = 'd$x/$d$t$'; end
                                if derIndex(2) ~= 0; ytit = 'd$y/$d$t$'; end
                            end
                            SIGNAL(:,1)   = SIGNAL(:,2);
                            SIGNAL(:,end) = SIGNAL(:,end-1);
                        end
                        plot(SIGNAL(1,:), SIGNAL(2,:),  ...
                             'DisplayName', ['[', num2str(nosigs(i, 1)), ', ', num2str(nosigs(i, 2)), ']']);
                    end
                    xlabel(xtit)
                    ylabel(ytit);
                    legend;
                    
                elseif size(nosigs, 2) == 3
                    % 3D signal line plot
                    h = figure();
                    hold on;
                    for i = 1:size(nosigs,1)
                        % Pre-processing on data
                        SIGNAL = obj.signals(nosigs(i,:),1:end-1);
                        % If all derivatives are requested, make them here
                        if exist('derIndex', 'var')
                            if derIndex == true
                                SIGNAL = SIGNAL - circshift(SIGNAL, +1, 2);
                                xtit = 'd$x/$d$t$';
                                ytit = 'd$y/$d$t$';
                                ztit = 'd$z/$d$t$';
                            else
                                SIGNAL = SIGNAL - diag(derIndex>0)*circshift(SIGNAL, +1, 2);
                                if derIndex(1) ~= 0; xtit = 'd$x/$d$t$'; end
                                if derIndex(2) ~= 0; ytit = 'd$y/$d$t$'; end
                                if derIndex(3) ~= 0; ztit = 'd$z/$d$t$'; end
                            end
                            SIGNAL(:,1)   = SIGNAL(:,2);
                            SIGNAL(:,end) = SIGNAL(:,end-1);
                        end
                        plot3(SIGNAL(1,1:end-1),        ...
                              SIGNAL(2,1:end-1),        ...
                              SIGNAL(3,1:end-1),        ...
                             'DisplayName',             ...
                                            ['[', num2str(nosigs(i, 1)), ', ', ...
                                                  num2str(nosigs(i, 2)), ', ', ...
                                                  num2str(nosigs(i, 3)),']']);
                    end
                    xlabel(xtit);
                    ylabel(ytit);
                    zlabel(ztit);
                    legend;
                    grid on;
                    view(-22, 22);
                end
            end
            
            if exist('xtit_c', 'var'); xlabel(xtit_c); end
            if exist('ytit_c', 'var'); ylabel(ytit_c); end
            if exist('ztit_c', 'var'); zlabel(ztit_c); end
            if exist('ftitle', 'var'); title(ftitle);  end
            if exist('flegend', 'var')
                if iscell(flegend)
                    legend(flegend);
                elseif flegend == 1
                    legend;
                elseif flegend == 0
                    legend off;
                end
            end
            
            % Make it pretty
            plt.isi;
            % Save it, if it is requested.
            if exist('fName', 'var')
                plt.figureSaveCore(fName, h)
            end
        end
        
        % The below function is used to plot a raster plot
        function raster(obj, varargin)
            % To set your opyions you must use the below codes:
            %  'select'     -> 1:3
            %  'ylabel'     -> 'Voltage (v)'
            %  'colorLimit' -> [0, 1]
            %  'gradient'   -> Gradient variable that made by gradient func.
            
            plt = libraries.plotToolLib;
            
            % Default ylabel and signals are set here. Indeed, all
            % signals are plotted in default mode.
            xtit = 'Time (s)';
            ytit = 'Signal index';
            sel  = 1:obj.numSignals;
            
            % Making a gradient for default mode
            numcl = 256;    % The number of colors exist in gradient
            colrs = [255, 255, 255; 220, 140, 100; 29, 67, 80; 6, 30, 45]/255;
            locs  = [0, 0.333, 0.666, 1];
            cmap  = plt.gradient(numcl, colrs, locs, 0);
            colorLimit = [min(obj.signals(:)), max(obj.signals(:))];
            
            % Affecting the arbitraty values of properties
            for i = 1:2:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'select'
                            sel = varargin{i+1};
                        case 'xlabel'
                            xtit_c = varargin{i+1};
                        case 'ylabel'
                            ytit_c = varargin{i+1};
                        case 'colorlimit'
                            colorLimit = varargin{i+1};
                        case 'gradient'
                            cmap = varargin{i+1};
                    end
                end
            end
            
            % Plot part
            figure();
                % Plot the image
                imagesc(obj.signals(sel,:));
                % Affect the gradient
                colormap(cmap);
                % Show the gradient
                colorbar;
                % Limit the colors
                caxis(colorLimit);
                % Title
                title('Raster plot');
                % X label
                xlabel(xtit);
                if exist('xtit_c', 'var'); xlabel(xtit_c); end
                % Y label
                ylabel(ytit);
                if exist('ytit_c', 'var'); ylabel(ytit_c); end
                % Can be 'normal' or 'revese'
                set(gca, 'Ydir', 'normal');
                % Limit the X axis
                xlim([0, size(obj.signals, 2)] + 0.5);
                % Limit the Y axis
                ylim([0, sel(end)-sel(1)+1] + 0.5);
                % Hide all X axis ticks except **
                xticks([1, numel(obj.timeLine)]);
                % Rename the remind ticks whcih is show
                xticklabels({num2str(0), num2str(obj.timeLine(end))});
                % If there is just a signal, do:
                if numel(sel) == 1
                    yticks(1);
                    yticklabels(num2str(sel(1)));
                % If there are more that a signal, do:
                else
                    yticks([1, numel(sel)]);
                    yticklabels({num2str(sel(1)), num2str(sel(end))});
                end
                % Change the appearance to LaTeX form
                plt.isi;
        end
        
        % Here is the end of methods
    end
    % Here is the end of class
end

