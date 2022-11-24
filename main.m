
%% // --------------------------------------------------------------
%    ***FARYADELL SIMULATION FRAMEWORK***
%    Creator:   Abolfazl Delavar
%    Web:       http://abolfazldelavar.com
%    Updated:   14 November 2022
%               2750 code lines
%
%    Description: This package has provided to simulate easily and quickly
%    Several useful tools have been written to simulate linear and nonlinear
%    dynamic systems, dependent or independent over time. To sum up, if
%    you are a student who wants to investigate and has decided to seek
%    dynamic and contrl systems, this is for you. Enjoy it :_)
% \\ ---------------------------------------------------------------

%% Closing all open windows, clearing variables and command
close all; clear; clc;

%% Setting parameters
params = valuation(); % extOpt

%% Making models and loading libraries
[params, models, signals, lib] = initialization(params);

%% Main simulation core
[params, models, signals, lib] = simulation(params, models, signals, lib);

%% Plotting and depiction
[params, models, signals, lib] = depiction(params, models, signals, lib, true);

%% Finalize, saving, and other after-simulated codes
finalize(params, models, signals, lib);

