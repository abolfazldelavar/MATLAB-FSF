
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

classdef plotToolLib
    % In this library, all functions are related to plotting and depiction are
    % provided which can be used in 'depiction.m' file.
    methods
        
        function isi(obj, varargin)
            % <- (Save or Not (string is the file name), Figure handle, Options)
            % The below function change the style similar to ISI journals.
            % It changes fonts and its sizes like LaTeX
            % OPTION ORDERS:
            %     opt.pictureWidth = 20;
            %     opt.hwRatio      = 0.65;
            %     opt.fontSize     = 17;
            %     isi(h, opt);
            % ---------------------------------------------------------
            %
            
            son          = 0;       % Save or not? What is its name?
            h            = gcf;     % The previous figure is considered
            fontSize     = 17;      % The default font size of the figure
            pictureWidth = 20;      % Set the default width of the considered figure
            hwRatio      = 0.65;    % The default ratio between height and width
            
            % Extracting the arbitraty values of properties
            for i = 1:2:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'save'
                            son          = varargin{i+1};
                        case 'figure'
                            h            = varargin{i+1};
                        case 'width'
                            pictureWidth = varargin{i+1};
                        case 'hwratio'
                            hwRatio      = varargin{i+1};
                        case 'fontsize'
                            fontSize     = varargin{i+1};
                    end
                end
            end

            % All font sizes contains labels and tick labels
            set(findall(h, '-property', 'FontSize'), 'FontSize', fontSize);
            % Hiding the top and the right borders
            set(findall(h, '-property', 'Box'), 'Box', 'off');
            % Changing the labels font and set  LaTeX format
            set(findall(h, '-property', 'Interpreter'), 'Interpreter', 'latex');
            % Changing the tick labels font and set  LaTeX format
            set(findall(h, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
            % Changing size of the figure as you ordered
            set(h, 'Units', 'centimeters', 'Position', [3 3 pictureWidth hwRatio*pictureWidth]);
            % Getting figure position and size
            pos = get(h, 'Position');
            % Below line changes the printed-version size
            set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)]);
            % Two below lines save the figure in formats PNG and PDF
            if ischar(son) %User must import the figure name as 'son' var
                obj.figureSaveCore(son, h);
            end
        end
        
        
        function figureSaveCore(~, son, h)
            % Loading requirements
            func     = libraries.functionLib;
            params   = valuation;
            % Set the focused figure if is not inserted
            if ~exist('h', 'var'); h = gcf; end
            % To get current PC time to use as a prefix in the name of file
            savePath   = [params.savePath, '/figs'];
            % Default saving format
            fFormat    = params.defaultImageFormat;
            AllFormats = {'jpg', 'png', 'pdf', 'fig'};
            isSetDirec = false;
            
            % Extracting the name and the directory which inported
            if ~exist('son', 'var') || ~ischar(son)
                % Set the time as its name, if there is no input in the arguments
                son = func.getNow(3,'-');
            else
                % Split folders
                tparts = split(son, '/');
                if numel(tparts) > 1
                    % Directory maker
                    savePath   = join(tparts(1:end-1), '/');
                    fparts     = split(tparts{end}, '.');
                    isSetDirec = 1;
                else
                    % If directory is not adjusted:
                    fparts = split(tparts, '.');
                end
                
                % Name and format
                if numel(fparts) > 1
                    % The name is also adjusted:
                    if numel(fparts{1}) == 0 || numel(fparts{end}) == 0
                        error('Please inter the correct notation for the file name.');
                    elseif ~(sum(contains(AllFormats, fparts{end})) && ...
                                 contains(fparts{end}, AllFormats))
                        error('You must input one of these formats: png/jpg/pdf/fig');
                    end
                    son     = join(fparts(1:end-1), '');
                    fFormat = fparts{2};
                else
                    % One of name or format just is inserted
                    % There is just name or format. It must be checked
                    if sum(contains(AllFormats, fparts{1})) && ...
                                 contains(fparts{1}, AllFormats)
                        fFormat = fparts{1};
                        % Set the time as its name, if there is no input in the arguments
                        son = func.getNow(3,'-');
                    else
                        % Just a name is imported, without directory
                        % and any formats
                        son = fparts{1};
                    end
                end
            end
            savePath = char(savePath);
            fFormat  = char(fFormat);
            son      = char(son);
            
            % Changing the file name 
            if params.uniqueSave == 1
                fName = [func.getNow, '_', son];
            else
                fName = son;
            end
            
            % Prepare direct
            if isSetDirec == 0
                fDir = [savePath, '/' fFormat];
            else
                fDir = savePath;
            end
            if ~isfolder(fDir); mkdir(fDir) ;end
            
            % Check the folders existance and make them if do not exist
            % Saving part
            switch fFormat
                case 'png'
                    print(h, [fDir, '/', fName], '-dpng', '-painters');
                case 'jpg'
                    print(h, [fDir, '/', fName], '-djpeg', '-painters');
                case 'pdf'
                    print(h, [fDir, '/', fName], '-dpdf', '-painters', '-fillpage');
                case 'fig'
                    savefig(h, [fDir, '/', fName]);
            end
            disp(['The graph named "' fName '.', fFormat, '" has been saved into "', fDir, '".']);
        end
        
        
        function resize(~, mode, h)
            % <- (Figure handle, Mode)
            if nargin < 2 + 1
                if nargin < 1 + 1
                    mode = 0;
                end
                h = gcf;
            end
            % Here the difult values were set.
            pictureWidth = 20;
            hwRatio      = 0.65;
            % This switch can provide some standard sizes
            % which you can change and personalize
            switch mode
                case 1
                    set(h, 'Units', 'centimeters', 'Position', [3 3 0.5*pictureWidth 0.5*pictureWidth]);
                case 2
                    set(h, 'Units', 'centimeters', 'Position', [3 3 1.5*pictureWidth hwRatio*pictureWidth]);
                otherwise
                    set(h, 'Units', 'centimeters', 'Position', [3 3 pictureWidth hwRatio*pictureWidth]);
            end
        end
        
        
        function dark(obj, varargin)
            % The function generates from a Matlab plot figure a version that can be
            % copied to a dark mode theme presentation or website.
            % The function replaces the default texts and box colors to
            % a user input color (default is white), and make the plot area transparent
            % to accept the dark background below it. The function also transform the
            % graphic colors that are not appropriate (low contrast) for a dark mode
            % theme to a version that is dark theme legible using a desaturation and
            % brightness approach.
            %
            % preparing this function I was inspired by https://material.io/design/color/dark-theme.html
            %
            % The function is a work in progess and may not support all figure objects


            %  Inputs:
            %  varargin(1)- The text color to modify (default is white)
            %  varargin(2)- The threshold from which to apply the cotrast correction (default is 4.5)
            %  varargin(3)- The dark background  (default is gray of value 0.16)
            %
            %
            %  How to the function:
            %  generate or import a Matlab figure and run the function:
            %
            %       plot(bsxfun(@times,[1:4],[2:5]'));xlabel('X');ylabel('Y');
            %       plot_darkmode
            %
            %  next copy the figure from the clipboard using Edit>Copy Figure and
            %  paste it on top of the dark background theme, for example in
            %  PowerPoint. Make sure that in the Copy Option, the  Transparent
            %  Background is enabled


            %   Ver 1.02 (2021-09-28)
            %   Adi Natan (natan@stanford.edu)

            %% defaults and initialize
            switch nargin-1
                case 3
                    textcolor           = varargin{1 + 1};
                    contrast_ratio      = varargin{2 + 1};
                    dark_bkg_assumption = varargin{3 + 1};
                case 2
                    textcolor           = varargin{1 + 1};
                    contrast_ratio      = varargin{2 + 1};
                    dark_bkg_assumption = ones(1,3)*39/255; % 0.16
                case 1
                    textcolor           = varargin{1 + 1};
                    contrast_ratio      = 4.5;
                    dark_bkg_assumption = ones(1,3)*39/255; % 0.16
                otherwise
                    textcolor           = [1,1,1]*0.95;
                    contrast_ratio      = 4.5;
                    dark_bkg_assumption = ones(1,3)*39/255; % 0.16
            end

            tcd = [{textcolor} , {contrast_ratio} , {dark_bkg_assumption}];


            % Getting the plot type
            g = get(get(gcf, 'children'), 'type');
            if ~strcmp(g, 'tiledlayout')
                % If it is not a kind of 'tiledlayout' do this

                % Getting figure children as a var named 'h'
                h             = get(gcf, 'children');
                % Providing a handle of each sections
                axes_ind      = findobj(h, 'type', 'Axes');
                legend_ind    = findobj(h, 'type', 'Legend');
                colorbar_ind  = findobj(h, 'type', 'Colorbar');
            else
                % If it is a kind of 'tiledlayout' do this

                % Getting figure children as a var named 'h'
                h0                = get(gcf, 'children');
                h0.Title.Color    = tcd{1};
                h0.Subtitle.Color = tcd{1};
                % Getting all children of main figure (like subsections of a subplot)
                h = get(get(gcf, 'children'), 'children');
                % Providing a handle of each sections
                axes_ind      = findobj(h, 'type', 'Axes');
                legend_ind    = findobj(h, 'type', 'Legend');
                colorbar_ind  = findobj(h, 'type', 'Colorbar');
            end

            %% modify Axes
            % Below loop repeats for each subplot of a mother figure
            for n = 1:numel(axes_ind)

                % % Edit x-ticks color
                % for m = 1:numel(axes_ind(n).XTickLabel)
                %     axes_ind(n).XTickLabel{m} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', textcolor), axes_ind(n).XTickLabel{m}];
                % end

                % % Edit y-ticks color
                % for m = 1:numel(axes_ind(n).YTickLabel)
                %     axes_ind(n).YTickLabel{m} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', textcolor), axes_ind(n).YTickLabel{m}];
                % end

                axes_ind(n).Color           = 'none';      % 'none' or tcd{3}; % make white area transparent
                axes_ind(n).XColor          = tcd{1};      % edit x axis color
                axes_ind(n).YColor          = tcd{1};      % edit y axis color
                axes_ind(n).ZColor          = tcd{1};      % edit z axis color

                axes_ind(n).XLabel.Color    = tcd{1};      % edit x label color
                axes_ind(n).YLabel.Color    = tcd{1};      % edit y label color
                axes_ind(n).ZLabel.Color    = tcd{1};      % edit z label color

                axes_ind(n).Title.Color     = tcd{1};      % edit title text color
                
                adjust_color                = @(x, y) obj.adjust_color(x, y);
                axes_ind(n).GridColor       = adjust_color(axes_ind(n).GridColor,      tcd);
                axes_ind(n).MinorGridColor  = adjust_color(axes_ind(n).MinorGridColor, tcd);
                % axes_ind(n).Subtitle.Color = textcolor;

                % take care of other axes children:
                h2              = get(axes_ind(n),'Children');
                g2              = get(axes_ind(n).Children,'type');
                text_ind        = find(strcmp(g2,'text'));
                patch_ind       = find(strcmp(g2,'patch'));
                line_ind        = find(strcmp(g2,'line'));
                errorbar_ind    = find(strcmp(g2,'errorbar'));
                area_ind        = find(strcmp(g2,'area'));
                bar_ind         = find(strcmp(g2,'bar'));
                hist_ind        = find(strcmp(g2,'histogram'));
                % contour_ind  = find(strcmp(g2,'contour'));
                % surface_ind = find(strcmp(g2,'surface'));

                % edit texts color
                for m = 1:numel(text_ind)
                    h2(text_ind(m)).Color = adjust_color( h2(text_ind(m)).Color ,tcd);
                    if ~strcmp( h2(text_ind(m)).BackgroundColor,'none')
                        %if text has some background color switch to dark bkg theme
                        h2(text_ind(m)).BackgroundColor = tcd{3};
                    end
                end

                % brighten patch colors if dim (use for the case of arrows etc)
                % this might not work well for all patch types so consider to comment
                for m = 1:numel(patch_ind)
                    h2(patch_ind(m)).FaceColor = adjust_color(h2(patch_ind(m)).FaceColor,tcd);
                    h2(patch_ind(m)).EdgeColor = adjust_color(h2(patch_ind(m)).EdgeColor,tcd);
                end

                for m = 1:numel(line_ind)
                    h2(line_ind(m)).Color = adjust_color(h2(line_ind(m)).Color,tcd);
                end


                for m = 1:numel(errorbar_ind)
                    h2(errorbar_ind(m)).Color = adjust_color(h2(errorbar_ind(m)).Color,tcd);
                    h2(errorbar_ind(m)).MarkerEdgeColor = adjust_color(h2(errorbar_ind(m)).MarkerEdgeColor,tcd);
                    h2(errorbar_ind(m)).MarkerFaceColor = adjust_color(h2(errorbar_ind(m)).MarkerFaceColor,tcd);
                end

                for m = 1:numel(area_ind)
                    h2(area_ind(m)).FaceColor = adjust_color(h2(area_ind(m)).FaceColor,tcd);
                    h2(area_ind(m)).EdgeColor = adjust_color(h2(area_ind(m)).EdgeColor,tcd);
                end

                for m = 1:numel(bar_ind)
                    h2(bar_ind(m)).FaceColor = adjust_color(h2(bar_ind(m)).FaceColor,tcd);
                    h2(bar_ind(m)).EdgeColor = adjust_color(h2(bar_ind(m)).EdgeColor,tcd);
                end


                for m = 1:numel(hist_ind)
                    h2(hist_ind(m)).FaceColor = adjust_color(h2(hist_ind(m)).FaceColor,tcd);
                    h2(hist_ind(m)).EdgeColor = adjust_color(h2(hist_ind(m)).EdgeColor,tcd);
                end

                %       for m=1:numel(contour_ind)
                %         h2(contour_ind(m)).FaceColor = adjust_color(h2(contour_ind(m)).FaceColor,tcd);
                %         h2(contour_ind(m)).EdgeColor = adjust_color(h2(contour_ind(m)).EdgeColor,tcd);
                %     end


            end
            %% modify Colorbars:
            for n = 1:numel(colorbar_ind)
                colorbar_ind(n).Color        =  textcolor;
                colorbar_ind(n).Label.Color  =  textcolor;
            end

            %% modify Legends:
            for n = 1:numel(legend_ind)
                legend_ind(n).Color     = 'none';     % make white area transparent
                legend_ind(n).TextColor = textcolor;  % edit text color
                legend_ind(n).Box       = 'off';      % delete box
            end


            %% modify annotations:
            ha = findall(gcf,'Tag','scribeOverlay');
            % get its children handles
            if ~isempty(ha)
                for n = 1:numel(ha)
                    hAnnotChildren = get(ha(n),'Children');
                    try
                        hAnnotChildrenType = get(hAnnotChildren, 'type');
                    catch
                        disp('annotation not available')
                        return
                    end

                    % edit lineType and shapeType colors
                    textboxshape_ind        = find(strcmp(hAnnotChildrenType,'textboxshape'));
                    ellipseshape_ind        = find(strcmp(hAnnotChildrenType,'ellipseshape'));
                    rectangleshape_ind      = find(strcmp(hAnnotChildrenType,'rectangleshape'));
                    textarrowshape_ind      = find(strcmp(hAnnotChildrenType,'textarrowshape'));
                    doubleendarrowshape_ind = find(strcmp(hAnnotChildrenType,'doubleendarrowshape'));
                    arrowshape_ind          = find(strcmp(hAnnotChildrenType,'arrowshape'));
                    arrow_ind               = find(strcmp(hAnnotChildrenType,'Arrow')); % older Matlab ver
                    lineshape_ind           = find(strcmp(hAnnotChildrenType,'lineshape'));


                    for m = 1:numel(textboxshape_ind)
                        hAnnotChildren(textboxshape_ind(m)).Color      =  textcolor;
                        hAnnotChildren(textboxshape_ind(m)).EdgeColor  =  adjust_color(hAnnotChildren(textboxshape_ind(m)).EdgeColor);
                    end

                    for m = 1:numel(ellipseshape_ind)
                        hAnnotChildren(ellipseshape_ind(m)).Color      =  adjust_color(hAnnotChildren(ellipseshape_ind(m)).Color,tcd);
                        hAnnotChildren(ellipseshape_ind(m)).FaceColor  =  adjust_color(hAnnotChildren(ellipseshape_ind(m)).FaceColor,tcd);
                    end

                    for m = 1:numel(rectangleshape_ind)
                        hAnnotChildren(rectangleshape_ind(m)).Color      =  adjust_color(hAnnotChildren(rectangleshape_ind(m)).Color,tcd);
                        hAnnotChildren(rectangleshape_ind(m)).FaceColor  =  adjust_color(hAnnotChildren(rectangleshape_ind(m)).FaceColor,tcd);
                    end

                    for m = 1:numel(textarrowshape_ind)
                        hAnnotChildren(textarrowshape_ind(m)).Color      =  adjust_color(hAnnotChildren(textarrowshape_ind(m)).Color,tcd);
                        hAnnotChildren(textarrowshape_ind(m)).TextColor  =  textcolor;
                        hAnnotChildren(textarrowshape_ind(m)).TextEdgeColor = adjust_color(hAnnotChildren(textarrowshape_ind(m)).TextEdgeColor,tcd);
                    end

                    for m = 1:numel(doubleendarrowshape_ind)
                        hAnnotChildren(doubleendarrowshape_ind(m)).Color = adjust_color(hAnnotChildren(doubleendarrowshape_ind(m)).Color,tcd);
                    end

                    for m = 1:numel(arrowshape_ind)
                        hAnnotChildren(arrowshape_ind(m)).Color = adjust_color(hAnnotChildren(arrowshape_ind(m)).Color,tcd);
                    end

                    for m = 1:numel(arrow_ind)
                        hAnnotChildren(arrow_ind(m)).Color = adjust_color(hAnnotChildren(arrow_ind(m)).Color,tcd);
                    end

                    for m = 1:numel(lineshape_ind)
                        hAnnotChildren(lineshape_ind(m)).Color = adjust_color(hAnnotChildren(lineshape_ind(m)).Color,tcd);
                    end
                end
            end
        end
        
        
        function outputColors = gradient(~, numofColors, mainColors, locs, showOrNot)
            % -----------------------------------------------------------------
            % 'numofColors' denotes the number of colors in output
            % 'mainColors' connotes the gradient colors from lower to upper
            % 'locs' provides an option to set position of colors
            % if 'showOrNot' is set true, the result of gradient will be shown.
            % the provided outputes is reported as 'outputColors' variable
            % Example of use:
            %       numcl = 256;
            %       colrs = [1 0 0; 1 1 1; 0 0 0];
            %       locs  = [0, 0.1, 1];
            %       cmap  = makeGradient(numcl, colrs, locs, 0);
            % -----------------------------------------------------------------

            numofGrads = size(mainColors, 1) - 1;
            cols       = zeros(numofGrads, 1);
            areaa      = [locs(1), locs(end)];
            
            func = libraries.functionLib;
            for i = 1:(numofGrads + 1)
                locs(i) = func.lineMapping(locs(i), areaa, [0, 1]);
            end

            for i = 1:numofGrads
                cols(i) = ceil(numofColors * (locs(i + 1) - locs(i)) );
            end

            initcolors = zeros(sum(cols), 3);
            shifft = 0;

            % Gradient maker
            for i = 1:numofGrads
                color1  = mainColors(i, :);
                color2  = mainColors(i + 1, :);
                gradian = interp1([0, 1], [color1; color2], linspace(0, 1, cols(i)));
                initcolors((1:cols(i)) + (i - 1) + shifft, :) = gradian;
                shifft  = shifft + cols(i) - 1;
            end
            outputColors = initcolors(1:numofColors, :);

            % Plot gradient
            if showOrNot == 1
                figure();
                img = repmat(1:numofColors, numofColors, 1);
                img = img(:,end:-1:1)';
                imshow(img, outputColors);
                colorbar;
            end
        end
        
        
        % End of methods
    end
    % End of class
    
    methods (Hidden=true)
        function out = adjust_color(~, in, tcd)
            % This function modifies an input color to fit a dark theme background.
            % For that a color needs to have sufficient contrast (WCAG's AA standard of at least 4.5:1)
            % The contrast ratio is calculate via :  cr = (L1 + 0.05) / (L2 + 0.05),
            % where L1 is the relative luminance of the input color and L2 is the
            % relative luminance of the dark mode background.
            % For this case we will assume a dark mode theme background of...
            % If a color is not passing this ratio, it will be modified to meet it
            % via desaturation and brightness to be more legible.
            % the function uses fminbnd, if you dont have the toolbox to use it you can
            % replace it with fmintx (avaiable in Matlab's file exchange)

            % if color is 'none' return as is
            if strcmp(in,'none')
                out=in;
                return
            end

            if isa(in,'char') % for inputs such as 'flat' etc...
                out=in;
                return
            end

            dark_bkg_assumption=tcd{3};

            % find the perceived lightness which is measured by some vision models
            % such as CIELAB to approximate the human vision non-linear response curve.
            % 1. linearize the RGB values (sRGB2Lin)
            % 2. find Luminance (Y)
            % 3. calc the perceived lightness (Lstar)
            % Lstar is in the range 0 to 1 where 0.5 is the perceptual "middle gray".
            % see https://en.wikipedia.org/wiki/SRGB ,

            sRGB2Lin=@(in) (in./12.92).*(in<= 0.04045) +  ( ((in+0.055)./1.055).^2.4 ).*(in> 0.04045);
            %Y = @(in) sum(sRGB2Lin(in).*[0.2126,  0.7152,  0.0722 ]);
            Y = @(in) sum(bsxfun(@times,sRGB2Lin( in ),[0.2126,  0.7152,  0.0722 ]),2 );
            Lstar = @(in)  0.01.*( (Y(in).*903.3).*(Y(in)<= 0.008856) + (Y(in).^(1/3).*116-16).*(Y(in)>0.008856));

            Ybkg = sum(sRGB2Lin(dark_bkg_assumption).*[0.2126,  0.7152,  0.0722 ]);

            cr = @(in)   (Y(in)' + 0.05) ./ (Ybkg + 0.05); % contrast ratio

            % rgb following desaturation of factor x
            ds=@(in,x) hsv2rgb( bsxfun(@times,rgb2hsv(in),[ones(numel(x),1) x(:) ones(numel(x),1)] ));

            % rgb following brightness change of factor x
            br=@(in,x) hsv2rgb( bsxfun(@power,rgb2hsv(in),[ones(numel(x),1) ones(numel(x),1) x(:)] ));


            if cr(in)<tcd{2} % default is 4.5

                %check if color is just black and replace with perceptual "middle gray"
                if ~sum(in)
                    fun0 = @(x) abs(Lstar( (ones(1,3)*x-dark_bkg_assumption ))-0.5);
                    L_factor=fminbnd(fun0,0.3,1);

                    out = ones(1,3)*L_factor;
                    return

                end


                % if saturation is what reduce contrast then desaturate
                in_hsv = rgb2hsv(in);
                if in_hsv(2) > 0.5
                    fun1 = @(x) abs(cr(ds(in,x)) - tcd{2});
                    [ds_factor, val] = fminbnd(fun1, 0, in_hsv(2));
                    if val < 1e-2
                        out = ds(in, ds_factor);
                        return
                    end
                end
                % desaturation alone didn't solve it, try to increase brightness
                fun2 = @(x) abs(cr(br(in,x)) - tcd{2});
                [br_factor, val] = fminbnd(fun2, 0, 1);

                if val<1e-2 && Lstar(br(in,br_factor))>0.5
                    out = br(in,br_factor);
                    return
                end
                % if niether worked then brightening + desaturation:
                fun3 = @(x) abs(cr(ds(br(in,br_factor),x))-tcd{2});
                [brds_factor, val]=fminbnd(fun3,0,1);

                if val<1e-2 && Lstar(ds(br(in,br_factor),brds_factor))>0.5
                    out = ds(br(in,br_factor),brds_factor);
                    return

                end
                % if all fails treat the color as black as above:
                fun0 = @(x) abs(Lstar( (ones(1,3)*x-dark_bkg_assumption ))-0.5);
                L_factor=fminbnd(fun0,0.3,1);
                out = ones(1,3)*L_factor;
            else
                out = in ;
            end
        end
        
        % End of the methods
    end
end

