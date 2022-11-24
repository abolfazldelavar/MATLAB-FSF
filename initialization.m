
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

function [params, models, signals, lib] = initialization(params) 
    % [Parameters, Models, Signals, Libraries] <- (Parameters)
    % This function is created to support your signals, systems,
    % extra functions, and all initialization processes before 
    % the main simulation.
    
    %% Loading Libraries
    lib.func = libraries.functionLib();
    lib.draw = libraries.plotToolLib();
    lib.mfun = libraries.yourLib();
    
    %% Simulation Time Vectors
    signals.tLine = 0:params.step:params.tOut;
    
    
    %% Your Personalized Section
    %  Write your codes here ...
    
    
    %% Updating
    %  If you need to update or create some new parameters, finalize them here.
    models.updated = 1;
    params.updated = 1;
end