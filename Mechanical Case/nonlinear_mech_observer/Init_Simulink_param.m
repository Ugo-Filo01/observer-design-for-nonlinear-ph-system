%% Mechanical-Case Simulink Parameters Code 
% This code must be launched before starting the Simulink scheme

% Clear 
clc; clear; close all;

%% Parameters 
% Mechanical parameters
params.m   = 1;
params.K   = 1e3;
params.eta = 50;

% Electrical parameters
params.epsllon = 2.8;   % epsilon
params.re  = 5;

% Initial condition
params.q0  = 1e-3;

%% Initial condition Mechanical system
% System initial condition
q_0 = 0;
p_0 = 0;

same_initial_cond = false;

if same_initial_cond

    qhat_0 = q_0;
    phat_0 = p_0;

else
    qhat_0 = 2e-4;
    phat_0 = 1e-3;

    % qhat_0 = 2e-5;
    % phat_0 = 1e-4;
end