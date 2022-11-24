
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%% // --------------------------------------------------------------

function [params, models, signals, lib] = depiction(params, models, signals, lib, plotOrNot)
    % [Parameters, Models, Signals, Libraries] <- (Parameters, Models, Signals, Libraries, Plot or Not (1 or 0))
    % All of your plots and depiction objects shoud be coded here
    % To use this, you have initialize, main, and finalize sections.
    
    if plotOrNot == true
        %% Initialize
        plt   = lib.draw;           % Import all plotting functions
        n     = params.n;
        nn    = 1:n;                % A vector from 1 to n
        tLine = signals.tLine(nn);  % A time-line vector
        
        
        %% Main part
        %  Insert your ideas here ...
        
        
        %% Finalize
        %  If you want to save images and graphs use this section to save
        %  them. All variables will be returned to the 'main' automatically,
        %  so you can change any parameters if you need.
        
    end
    % The end of 'depiction' function
end

