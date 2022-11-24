
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

function params = valuation(extOpt) 
    % [Parameters] <- (External variables as a structure variable)
    % extOpt is used for external control on whole simulation by a
    % second m file. If you want to import data here, make a structure data and
    % put it into the 'validation()' function in 'main.m'.
    
    %% Simulation parameters
    % Tho below variables are time-step and simulation time, respectively.
    % You might need to modify them arbitrarily.
    params.step = 0.01; %(double)
    params.tOut = 10;   %(double)
    
    % The given line below is a dependent variable. You normally
    % should NOT change it.
    params.n = floor(params.tOut/params.step); %(DEPENDENT)
    
    % The two next variables carry the folders from which you can use to
    % call your fales and saving the results organizely.
    params.loadPath = 'data/inputs';  %(string)
    params.savePath = 'data/outputs'; %(string)
    
    % Do you want to save a diary after each simulation? So set the below logical
    % variable "true". The below directory is the place your logs are saved.
    params.makeDiary = true;   %(logical)
    params.diaryDir  = 'logs'; %(string)
    
    % The amount of time (in second) between each commands printed on CW.
    params.commandIntervalSpan = 5; %(double)
    
    % The below line is related to creating names when you are saving data
    % which include the time of saving, It is a logical variable which has
    % "false" value as its default
    params.uniqueSave = false; %(logical)
    
    % The below string can be set one of the following formats. It's the
    % default format whcih when you do not insert any formats, it will be
    % considered. These legall formats are "jpg", "png", "pdf", "fig"
    params.defaultImageFormat = 'png'; %(string)
    
    
    %% Your code ...
    
    
    %% External properties - This is used to control 'main.m from other files
    if nargin > 0
        params.extOpt = extOpt;
    end
end
