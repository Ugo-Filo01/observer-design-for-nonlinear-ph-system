%% Electro-Mechanical system
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

params.u_op = 5e3; % for linearization

%%
% Electromechanical system
[H_sym, gradH_sym, HessH_sym, lambdas_sym] = H_electromech(params);

% Mechanical system (not used)
[H_sym_mech, gradH_sym_mech, HessH_sym_mech, lambdas_sym_mech, g_sym_mech, grad_g_mech] = H_mech(params);

% Start of the linearization of the system
%% Linearization around equilibrium for u_eq
% Electromechanical system
[A, B ,C ,Cnum,D ,x_eq ,u_eq ,y_eq, i_eq,Hxx] = linearize_electromech(params, HessH_sym);

% Mechanical system (not used)
[A_mech,B_mech,C_mech,D_mech,x_eq_mech,u_eq_mech,y_eq_mech, Hxx_mech, g_mech] = linearize_mech(params, gradH_sym_mech, ...
    HessH_sym_mech, g_sym_mech, grad_g_mech);

C_q = [1 0 0];


%% Initial condition Electromechanical system
% System initial condition
q_0 = 0;
p_0 = 0;
Q_0 = 0;

same_initial_cond = false; % Use boolean true instead of string 'true'
if same_initial_cond
    qhat_0 = 0;
    phat_0 = 0;
    Qhat_0 = 0;

    qtilde_0 = qhat_0 - x_eq(1); 
    ptilde_0 = phat_0 - x_eq(2); 
    Qtilde_0 = Qhat_0 - x_eq(3);

else
    qhat_0 = 2e-4;
    phat_0 = 1e-3;
    Qhat_0 = 5e-8;

    % qhat_0 = 2e-5;
    % phat_0 = 1e-4;
    % Qhat_0 = 5e-9;

    qtilde_0 = qhat_0 - x_eq(1); 
    ptilde_0 = phat_0 - x_eq(2); 
    Qtilde_0 = Qhat_0 - x_eq(3);

end







