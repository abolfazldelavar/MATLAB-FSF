
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

function [params, models, signals, lib] = simulation(params, models, signals, lib)
    % [Parameters, Models, Signals, Libraries] <- (Parameters, Models, Signals, Libraries)
    % This file is your main core of simulation which there is a time-loop
    % Also, the model blocks and the signals must be updated here.
    % Before the main loop, you can initialize if you need, Also you can
    % finalize after that if you need, as well.
    
    %% Initial options
    func = lib.func;
    func.sayStart();
    % A trigger used to report steps in command
    trig = [2, 0];
    
    %% Main loop
    for k = 1:params.n
        % Displaying the iteration number on the command window
        trig = func.disit(k, params.n, trig, params.commandIntervalSpan);
        
        
        %% Write your codes here ...
        
        
    end
    
    %% Finalize options
    % To report the simulation time after running has finished
    func.disit(k, params.n, [0, trig(2)], params.commandIntervalSpan);
    func.sayEnd();
end
