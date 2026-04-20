clc; clear; close all;

%% ========================
% Physical parameters
% =========================
k   = 1000;
m   = 1;
eta = 50;
eps = 2.8;
q0  = 1e-3;

Umax = 5140;
ubar = Umax^2;

A0 = [0 1/m; -k -eta/m];

g2   = @(q) 2*eps*(q + q0).^3;
aFun = @(q) 6*eps*(q + q0).^2;

% Input applied at t=1s
U = @(t) (t>=1).*Umax;
u = @(t) U(t).^2;

fprintf('=== PAPER FIGURES GENERATION ===\n');
fprintf('Input applied at t=1s\n');

%% ========================
% Polytope construction
% =========================
% Value extract from Simulink simulation of the plant
qmin = -8.1257e-6;
qmax =  4.6754e-4;
pmin = -6.302e-3;
pmax =  2.228e-3;

cand_q = [qmin, qmax, -q0];
cand_q = cand_q(cand_q >= qmin & cand_q <= qmax);

amin = min(aFun(cand_q));
amax = max(aFun(cand_q));
gmin = min(g2([qmin qmax]));
gmax = max(g2([qmin qmax]));

beta_corners = [amin*pmin, amin*pmax, amax*pmin, amax*pmax];
betamin = min(beta_corners);
betamax = max(beta_corners);

a_set    = [amin, amax];
u_set    = [0, ubar];
beta_set = [betamin, betamax];
g_set    = [gmin, gmax];

Nv   = 16;
Abar = cell(Nv,1);
Cbar = cell(Nv,1);

idx = 0;
for ia = 1:2
    for iu = 1:2
        for ib = 1:2
            for ig = 1:2
                idx = idx + 1;
                Abar{idx} = A0 + [0 0; a_set(ia)*u_set(iu) 0];
                Cbar{idx} = [beta_set(ib)/m, g_set(ig)/m];
            end
        end
    end
end

betaP = 1e6;
rho   = 1e-6;

%% ==========================
% Find maximum lambda values
% ===========================
fprintf('\n=== FINDING MAXIMUM LAMBDA VALUES ===\n');

lambda_max_const = find_max_lambda_const(Abar, Cbar, Nv, betaP, rho);
fprintf('lambda_max (constant):  %.6f\n', lambda_max_const);

lambda_max_sched = find_max_lambda_sched(Abar, Cbar, Nv, betaP, rho);
fprintf('lambda_max (scheduled): %.6f\n', lambda_max_sched);

%% ===========================
% GROUP 1: Fair Comparison
% ============================
fprintf('\n=== GROUP 1: FAIR COMPARISON ===\n');

alpha = 0.1;
lambda_fair = alpha * lambda_max_const;

fprintf('Using lambda = %.6f for both observers\n', lambda_fair);

[ok1, P_const_fair, L_const_fair, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_fair, rho);
if ~ok1, error('Constant @ lambda_fair infeasible'); end

[ok2, P_sched_fair, K_cells_fair, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, lambda_fair, rho);
if ~ok2, error('Scheduled @ lambda_fair infeasible'); end

L_vertices_fair = cell(Nv,1);
for i = 1:Nv
    L_vertices_fair{i} = P_sched_fair \ K_cells_fair{i};
end

params_fair.L_vertices  = L_vertices_fair;
params_fair.K           = k;
params_fair.m           = m;
params_fair.eta         = eta;
params_fair.eps         = eps;
params_fair.q0          = q0;
params_fair.a_bounds    = [amin, amax];
params_fair.g_bounds    = [gmin, gmax];
params_fair.beta_bounds = [betamin, betamax];
params_fair.u_max       = ubar;

[T1_const, q1_const, p1_const, qh1_const, ph1_const, eq1_const, ep1_const, ...
    L1_const_history, T1_const_history] = ...
    simulate_observer('constant', L_const_fair, [], k, m, eta, eps, q0, u, g2);

[T1_sched, q1_sched, p1_sched, qh1_sched, ph1_sched, eq1_sched, ep1_sched, ...
    L1_sched_history, T1_sched_history] = ...
    simulate_observer('scheduled', [], params_fair, k, m, eta, eps, q0, u, g2);

fprintf('GROUP 1 simulations completed\n');

%% =============================
% GROUP 2: Scheduling Advantage
% ==============================
fprintf('\n=== GROUP 2: SCHEDULING ADVANTAGE ===\n');

lambda_const_max = lambda_max_const;
fprintf('Constant:  lambda = %.6f (100%% of max)\n', lambda_const_max);

[ok3, P_const_max, L_const_max, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_const_max, rho);
if ~ok3, error('Constant @ lambda_max infeasible'); end

lambda_sched_fair = lambda_max_const;
fprintf('Scheduled: lambda = %.6f (same as constant max)\n', lambda_sched_fair);

[ok4, P_sched_fair2, K_cells_fair2, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, lambda_sched_fair, rho);
if ~ok4, error('Scheduled @ lambda_max_const infeasible'); end

L_vertices_fair2 = cell(Nv,1);
for i = 1:Nv
    L_vertices_fair2{i} = P_sched_fair2 \ K_cells_fair2{i};
end

params_fair2 = params_fair;
params_fair2.L_vertices = L_vertices_fair2;

lambda_sched_max = lambda_max_sched;
fprintf('Scheduled: lambda = %.6f (100%% of max_sched)\n', lambda_sched_max);

[ok_fail, ~, ~, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_sched_max, rho);
if ok_fail
    warning('Constant still feasible at lambda=%.3f!', lambda_sched_max);
else
    fprintf(' Constant FAILS at lambda=%.3f\n', lambda_sched_max);
end

[ok5, P_sched_max, K_cells_max, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, lambda_sched_max, rho);
if ~ok5, error('Scheduled @ lambda_sched_max infeasible'); end

L_vertices_max = cell(Nv,1);
for i = 1:Nv
    L_vertices_max{i} = P_sched_max \ K_cells_max{i};
end

params_max = params_fair;
params_max.L_vertices = L_vertices_max;

[T2_const, q2_const, p2_const, qh2_const, ph2_const, eq2_const, ep2_const, ...
    L2_const_history, T2_const_history] = ...
    simulate_observer('constant', L_const_max, [], k, m, eta, eps, q0, u, g2);

[T2_sched_fair, q2_sched_fair, p2_sched_fair, qh2_sched_fair, ph2_sched_fair, ...
    eq2_sched_fair, ep2_sched_fair, L2_sched_fair_history, T2_sched_fair_history] = ...
    simulate_observer('scheduled', [], params_fair2, k, m, eta, eps, q0, u, g2);

[T2_sched_max, q2_sched_max, p2_sched_max, qh2_sched_max, ph2_sched_max, ...
    eq2_sched_max, ep2_sched_max, L2_sched_max_history, T2_sched_max_history] = ...
    simulate_observer('scheduled', [], params_max, k, m, eta, eps, q0, u, g2);

fprintf('GROUP 2 simulations completed\n');

%% ==============================================
% FIGURE GENERATION - GROUP 1
% Separate constant and scheduled gain subplots
% ===============================================
fprintf('\n=== GENERATING GROUP 1 FIGURES ===\n');
 
% Color scheme
col_real  = [0 0 0];           % Black
col_const = [0.85 0.33 0.10];  % Red
col_sched = [0 0.45 0.74];     % Blue
col_L2    = [0.93 0.69 0.13];  % Yellow/gold
col_L2c   = [0.50 0.50 0.50];  % Gray
 
% Font sizes
big_title = 22;
fs_title  = 20;
fs_label  = 19;
fs_legend = 17;
fs_tick   = 19;

% GROUP 1 - Main figure (2-column layout with 6 subplots)
fig1 = figure('Position', [50 50 1200 1200], 'Color', 'w');
sgtitle(sprintf('\\textbf{Group 1: Fair Comparison (Same $\\lambda$ = %.4f)}', lambda_fair), ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Row 1, Col 1: Position tracking
subplot(3,2,1);
plot(T1_const, q1_const*1e6, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T1_const, qh1_const*1e6, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/21):length(T1_const));
plot(T1_sched, qh1_sched*1e6, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$q$ [$\mu$m]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(a) Position and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', 'Const.', 'Sched.'}, 'Location', 'southeast', 'FontSize', fs_legend);
grid on; xlim([0 5]); ylim([-50 500]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 1, Col 2: Position error (zoom)
subplot(3,2,2);
plot(T1_const, zeros(size(T1_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T1_const, eq1_const*1e6, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/340):length(T1_const));
plot(T1_sched, eq1_sched*1e6, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$\tilde{q}$ [$\mu$m]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(b) Position Error (0--0.3 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', 'Const.', 'Sched.'}, 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.3]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 2, Col 1: Momentum tracking
subplot(3,2,3);
plot(T1_const, p1_const*1e3, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T1_const, ph1_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/21):length(T1_const));
plot(T1_sched, ph1_sched*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$p$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(c) Momentum and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', 'Const.', 'Sched.'}, 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 2, Col 2: Momentum error (zoom)
subplot(3,2,4);
plot(T1_const, zeros(size(T1_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T1_const, ep1_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/340):length(T1_const));
plot(T1_sched, ep1_sched*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$\tilde{p}$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(d) Momentum Error (0--0.3 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', 'Const.', 'Sched.'}, 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.3]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 3, Col 1: Constant observer gains
subplot(3,2,5);
plot(T1_const_history, L1_const_history(:,1), '--', 'LineWidth', 2.5, 'Color', col_const); hold on;
plot(T1_const_history, L1_const_history(:,2), '--', 'LineWidth', 2.5, 'Color', col_L2c);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(e) Observer Gains $L_1$ and $L_2$ Constant', ...
    'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.5, 1.1, 0]);
legend({'$L_1$ const.', '$L_2$ const.'}, 'Interpreter', 'latex', ...
    'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]); ylim([-1e8 17e8]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 3, Col 2: Scheduled observer gains
subplot(3,2,6);
plot(T1_sched_history, L1_sched_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched); hold on;
plot(T1_sched_history, L1_sched_history(:,2), '-', 'LineWidth', 2.5, 'Color', col_L2);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(f) Observer Gains $L_1(t)$ and $L_2(t)$ Scheduled', ...
    'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.5, 1.1, 0]);
legend({'$L_1$ sched.', '$L_2$ sched.'}, 'Interpreter', 'latex', ...
    'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Save Group 1 figure
savefig(fig1, 'Group1_Fair_Comparison_v2.fig');
exportgraphics(fig1, 'Group1_Fair_Comparison_v2.png', 'Resolution', 300);
exportgraphics(fig1, 'Group1_Fair_Comparison_v2.eps', 'ContentType', 'vector');
fprintf('✓ Saved: Group1_Fair_Comparison_v2.{fig,png,eps}\n');
 
%% ==============================
% FIGURE GENERATION - GROUP 2
% Tracking and error plots only
% ===============================
fprintf('\n=== GENERATING GROUP 2 FIGURES ===\n');
 
% GROUP 2 - Main figure (tracking and error plots)
fig2 = figure('Position', [50 50 1200 1000], 'Color', 'w');
sgtitle('\textbf{Group 2: Scheduling Advantage (Different $\lambda$)}', ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Colors for Group 2
col_const_2    = [0.85 0.33 0.10];  % Red (constant)
col_sched_fair = [0 0.45 0.74];     % Blue (scheduled @ same lambda)
col_sched_max  = [0.13 0.55 0.13];  % Green (scheduled @ max lambda)
col_L2_sf      = [0.53 0.81 0.92];  % Light blue
col_L2_sm      = [0.56 0.93 0.56];  % Light green
 
% Row 1, Col 1: Position tracking
subplot(2,2,1);
plot(T2_const, q2_const*1e6, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T2_const, qh2_const*1e6, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/21):length(T2_const));
plot(T2_sched_fair, qh2_sched_fair*1e6, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, qh2_sched_max*1e6, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$q$ [$\mu$m]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(a) Position and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', fs_legend);
grid on; xlim([0 5]); ylim([-50 500]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 1, Col 2: Position error (zoom)
subplot(2,2,2);
plot(T2_const, zeros(size(T2_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T2_const, eq2_const*1e6, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/130):length(T2_const)); 
plot(T2_sched_fair, eq2_sched_fair*1e6, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, eq2_sched_max*1e6, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$\tilde{q}$ [$\mu$m]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(b) Position Error (0--0.8 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.8]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 2, Col 1: Momentum tracking
subplot(2,2,3);
plot(T2_const, p2_const*1e3, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T2_const, ph2_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/21):length(T2_const));
plot(T2_sched_fair, ph2_sched_fair*1e3, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, ph2_sched_max*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$p$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(c) Momentum and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 2, Col 2: Momentum error (zoom)
subplot(2,2,4);
plot(T2_const, zeros(size(T2_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T2_const, ep2_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/130):length(T2_const)); 
plot(T2_sched_fair, ep2_sched_fair*1e3, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, ep2_sched_max*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$\tilde{p}$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(d) Momentum Error (0--0.8 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.8]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Save Group 2 main figure
savefig(fig2, 'Group2_Scheduling_Advantage_v2.fig');
exportgraphics(fig2, 'Group2_Scheduling_Advantage_v2.png', 'Resolution', 300);
exportgraphics(fig2, 'Group2_Scheduling_Advantage_v2.eps', 'ContentType', 'vector');
fprintf('✓ Saved: Group2_Scheduling_Advantage_v2.{fig,png,eps}\n');
 
%% =================================
% GROUP 2 - Gains comparison figure
% ==================================
fprintf('\n=== GENERATING GROUP 2 GAINS COMPARISON FIGURE ===\n');
 
fig3 = figure('Position', [200 100 700 1200], 'Color', 'w');
sgtitle('\textbf{Group 2: Observer Gains Comparison}', ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Subplot 1: Constant gains (top)
subplot(3,1,1);
plot(T2_const_history, L2_const_history(:,1), '--', 'LineWidth', 2.5, 'Color', col_const_2); hold on;
plot(T2_const_history, L2_const_history(:,2), '--', 'LineWidth', 2.5, 'Color', col_L2c);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(a) Constant Observer Gains', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ const.($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ const.($\\lambda$=%.3f)', lambda_sched_fair)'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]); ylim([-20e5 1e5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Subplot 2: Scheduled gains (both cases) (middle)
subplot(3,1,2);
plot(T2_sched_fair_history, L2_sched_fair_history(:,1), ':', 'LineWidth', 2.5, 'Color', col_sched_fair); hold on;
plot(T2_sched_fair_history, L2_sched_fair_history(:,2), ':', 'LineWidth', 2.5, 'Color', col_L2_sf);
plot(T2_sched_max_history, L2_sched_max_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched_max);
plot(T2_sched_max_history, L2_sched_max_history(:,2), '-', 'LineWidth', 2.5, 'Color', col_L2_sm);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(b) Scheduled Observer Gains', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_1$ sched. ($\\lambda$=%.2f)', lambda_sched_max), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
 
% Subplot 3: Zoom of scheduled gains (excluding L2 max component) (bottom)
subplot(3,1,3);
plot(T2_sched_fair_history, L2_sched_fair_history(:,1), ':', 'LineWidth', 2.5, 'Color', col_sched_fair); hold on;
plot(T2_sched_fair_history, L2_sched_fair_history(:,2), ':', 'LineWidth', 2.5, 'Color', col_L2_sf);
plot(T2_sched_max_history, L2_sched_max_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched_max);
% NOTE: Deliberately NOT plotting L2_sched_max_history(:,2)
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(c) Scheduled Gains Zoom', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_1$ sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
 
% Improve vertical spacing
drawnow;
axesHandles = findobj(fig3, 'Type', 'axes');
for ax = axesHandles'
    pos = get(ax, 'Position');
    pos(2) = pos(2) + 0.02;
    pos(4) = pos(4) - 0.02;
    set(ax, 'Position', pos);
end
 
 
% Save Group 2 gains figure
savefig(fig3, 'Group2_Gains_Comparison_v2.fig');
exportgraphics(fig3, 'Group2_Gains_Comparison_v2.png', 'Resolution', 300);
exportgraphics(fig3, 'Group2_Gains_Comparison_v2.eps', 'ContentType', 'vector');
fprintf('Saved: Group2_Gains_Comparison_v2.{fig,png,eps}\n');
 
%% =========================
% Debug print
% ==========================
fprintf('Figure Generation completed\n\n');
fprintf('LAMBDA VALUES:\n');
fprintf('  λ_max (constant):  %.6f\n', lambda_max_const);
fprintf('  λ_max (scheduled): %.6f\n', lambda_max_sched);
fprintf('  Improvement ratio: %.2fx\n', lambda_max_sched/lambda_max_const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lambda_max = find_max_lambda_const(Abar, Cbar, Nv, betaP, rho)
    % Binary search to find maximum feasible lambda for constant gain observer
    lambda_lo = 0.0;
    lambda_hi = 1.0;

    [ok_lo, ~, ~, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_lo, rho);
    if ~ok_lo, error('Even lambda=0 infeasible'); end

    % Find upper bound
    for k = 1:30
        [ok_hi, ~, ~, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_hi, rho);
        if ~ok_hi, break; end
        lambda_hi = 2*lambda_hi;
    end

    % Binary search
    lambda_max = 0;
    for it = 1:50
        lambda_mid = 0.5*(lambda_lo + lambda_hi);
        [ok_mid, ~, ~, ~] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda_mid, rho);

        if ok_mid
            lambda_max = lambda_mid;
            lambda_lo = lambda_mid;
        else
            lambda_hi = lambda_mid;
        end

        if abs(lambda_hi - lambda_lo) < 1e-3, break; end
    end
end

function lambda_max = find_max_lambda_sched(Abar, Cbar, Nv, betaP, rho)
    % Binary search to find maximum feasible lambda for scheduled gain observer
    lambda_lo = 0.0;
    lambda_hi = 1.0;

    [ok0, ~, ~, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, 0.0, rho);
    if ~ok0, error('Even lambda=0 infeasible'); end

    % Find upper bound
    for k = 1:25
        [okhi, ~, ~, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, lambda_hi, rho);
        if ~okhi, break; end
        lambda_hi = 2*lambda_hi;
    end

    % Binary search
    lambda_max = 0;
    for iter = 1:60
        lambda_mid = 0.5*(lambda_lo + lambda_hi);
        [okmid, ~, ~, ~] = solve_scheduled_LMI(Abar, Cbar, Nv, betaP, lambda_mid, rho);
        if okmid
            lambda_max = lambda_mid;
            lambda_lo = lambda_mid;
        else
            lambda_hi = lambda_mid;
        end
        if abs(lambda_hi - lambda_lo) < 1e-3, break; end
    end
end

function [T, q, p, qh, ph, eq, ep, L_history, T_history] = ...
    simulate_observer(mode, L_const, params, k, m, eta, eps, q0, u, g2)
    % Simulate plant and observer dynamics
    
    % Plant dynamics
    fPlant = @(t,x) [x(2)/m; -k*x(1) - (eta/m)*x(2) + g2(x(1))*u(t)];

    % Initial conditions
    x0    = [0; 0];       % Plant initial state
    xhat0 = [2e-4; -2e-3]; % Observer initial state (with error)
    xx0   = [x0; xhat0];

    tspan = [0 5];
    opts  = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',1e-3);

    if strcmp(mode, 'constant')
        % Constant gain observer
        A0 = [0 1/m; -k -eta/m];
        fObs = @(t,xx) obs_const(t, xx, fPlant, A0, u, g2, m, L_const);
        [T, XX] = ode15s(fObs, tspan, xx0, opts);
        L_history = repmat(L_const(:).', length(T), 1);
        T_history = T;
    else
        % Scheduled gain observer
        global L_sched_global T_sched_global
        L_sched_global = [];
        T_sched_global = [];

        fObs = @(t,xx) obs_sched(t, xx, fPlant, k, m, eta, eps, q0, u, params);
        [T, XX] = ode15s(fObs, tspan, xx0, opts);

        L_history = L_sched_global;
        T_history = T_sched_global;
    end

    % Extract state components
    q  = XX(:,1);
    p  = XX(:,2);
    qh = XX(:,3);
    ph = XX(:,4);
    eq = q - qh;
    ep = p - ph;
end

function dxx = obs_const(t, xx, fPlant, A0, u, g2, m, L)
    % Constant gain observer dynamics
    x    = xx(1:2);
    xhat = xx(3:4);
    
    % Plant dynamics
    dx = fPlant(t, x);
    
    % Output and its estimate
    y  = g2(x(1))    * x(2)    / m;
    yh = g2(xhat(1)) * xhat(2) / m;
    
    % Observer dynamics
    dxhat = A0*xhat + [0; g2(xhat(1))*u(t)] + L*(y - yh);
    dxx = [dx; dxhat];
end

function dxx = obs_sched(t, xx, fPlant, K, m, eta, eps, q0, u, params)
    % Scheduled gain observer dynamics
    global L_sched_global T_sched_global

    x    = xx(1:2);
    xhat = xx(3:4);
    
    % Plant dynamics
    dx = fPlant(t, x);

    % Output and its estimate
    y  = g2_local(x(1),    x(2),    m, eps, q0);
    yh = g2_local(xhat(1), xhat(2), m, eps, q0);

    % Compute scheduled gain
    L = compute_scheduled_gain(xhat, u(t), params);

    % Store gain history
    L_sched_global = [L_sched_global; L(:).'];
    T_sched_global = [T_sched_global; t];

    % Observer dynamics
    qh  = xhat(1);
    ph  = xhat(2);
    g2h = 2*eps*(qh + q0)^3;

    dxhat = [ph/m; -K*qh - eta*(ph/m)] + [0; g2h*u(t)] + L*(y - yh);
    dxx = [dx; dxhat];
end

function y = g2_local(q, p, m, eps, q0)
    % Local output function
    g2v = 2*eps*(q + q0)^3;
    y = g2v * (p/m);
end

function L = compute_scheduled_gain(xhat, uval, params)
    % Compute scheduled observer gain using polytopic interpolation
    qh = xhat(1);
    ph = xhat(2);

    % Compute scheduling variables
    a_hat    = 6*params.eps*(qh + params.q0)^2;
    g_hat    = 2*params.eps*(qh + params.q0)^3;
    beta_hat = a_hat*ph;

    % Normalization factors
    da = params.a_bounds(2)    - params.a_bounds(1);
    dg = params.g_bounds(2)    - params.g_bounds(1);
    db = params.beta_bounds(2) - params.beta_bounds(1);

    if da < eps, da = 1; end
    if dg < eps, dg = 1; end
    if db < eps, db = 1; end

    % Normalized scheduling parameters (0 to 1)
    mu_a = (a_hat    - params.a_bounds(1))    / da;
    mu_g = (g_hat    - params.g_bounds(1))    / dg;
    mu_b = (beta_hat - params.beta_bounds(1)) / db;
    mu_u = uval / params.u_max;

    % Saturation
    mu_a = max(0,min(1,mu_a));
    mu_g = max(0,min(1,mu_g));
    mu_b = max(0,min(1,mu_b));
    mu_u = max(0,min(1,mu_u));

    % Interpolation weights
    w_aL = 1-mu_a; w_aH = mu_a;
    w_uL = 1-mu_u; w_uH = mu_u;
    w_bL = 1-mu_b; w_bH = mu_b;
    w_gL = 1-mu_g; w_gH = mu_g;

    % Compute weights for all 16 vertices
    h = zeros(16,1);
    idx = 0;
    for ia = 1:2
        wa = (ia==1)*w_aL + (ia==2)*w_aH;
        for iu = 1:2
            wu = (iu==1)*w_uL + (iu==2)*w_uH;
            for ib = 1:2
                wb = (ib==1)*w_bL + (ib==2)*w_bH;
                for ig = 1:2
                    wg = (ig==1)*w_gL + (ig==2)*w_gH;
                    idx = idx + 1;
                    h(idx) = wa*wu*wb*wg;
                end
            end
        end
    end

    % Normalize weights
    s = sum(h);
    if s <= 0
        h = ones(16,1)/16;
    else
        h = h / s;
    end

    % Compute interpolated gain
    L = zeros(2,1);
    for i = 1:16
        L = L + h(i)*params.L_vertices{i};
    end
end

function [ok, P, L, tmin] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda, rho)
    % Design constant gain observer using LMI
    setlmis([]);
    n = 2;

    % Define LMI variables
    [Pvar, ~] = lmivar(1, [n 1]);
    [Kvar, ~] = lmivar(2, [n 1]);

    % LMI constraints for each vertex
    for i = 1:Nv
        lmi = newlmi;
        lmiterm([lmi 1 1 Pvar], 1, Abar{i}, 's');
        lmiterm([lmi 1 1 Kvar], -1, Cbar{i}, 's');
        lmiterm([lmi 1 1 Pvar], 2*lambda, 1);
        lmiterm([lmi 1 1 0], -rho*eye(n));
    end

    % Lower bound on P
    lmiPmin = newlmi;
    lmiterm([lmiPmin 1 1 Pvar], -1, 1);
    lmiterm([lmiPmin 1 1 0], eye(n));

    % Upper bound on P
    lmiPmax = newlmi;
    lmiterm([lmiPmax 1 1 Pvar], 1, 1);
    lmiterm([lmiPmax 1 1 0], -betaP*eye(n));

    lmis = getlmis;
    options = [1e-7, 800, -1, 0, 1];
    [tmin, xfeas] = feasp(lmis, options);

    if tmin > 0
        ok = false; P = []; L = []; return;
    end

    P = dec2mat(lmis, xfeas, Pvar);
    K = dec2mat(lmis, xfeas, Kvar);

    if min(eig(P)) <= 0
        ok = false; P = []; L = []; return;
    end

    L = P \ K;
    ok = true;
end

function [ok, P, K_cells, tmin] = solve_scheduled_LMI(Abar_cells, Cbar_cells, Nv, betaP, lambda, rho)
    % Design scheduled gain observer using LMI (separate K for each vertex)
    setlmis([]);
    n = size(Abar_cells{1},1);

    % Define LMI variables
    [Pvar, ~] = lmivar(1, [n 1]);

    Kvar_cells = cell(Nv,1);
    for i = 1:Nv
        [Kvar_cells{i}, ~] = lmivar(2, [n 1]);
    end

    % LMI constraints for each vertex
    for i = 1:Nv
        lmi = newlmi;
        lmiterm([lmi 1 1 Pvar], 1, Abar_cells{i}, 's');
        lmiterm([lmi 1 1 Kvar_cells{i}], -1, Cbar_cells{i}, 's');
        lmiterm([lmi 1 1 Pvar], 2*lambda, 1);
        lmiterm([lmi 1 1 0], -rho*eye(n));
    end

    % Lower bound on P
    lmiPmin = newlmi;
    lmiterm([lmiPmin 1 1 Pvar], -1, 1);
    lmiterm([lmiPmin 1 1 0], eye(n));

    % Upper bound on P
    lmiPmax = newlmi;
    lmiterm([lmiPmax 1 1 Pvar], 1, 1);
    lmiterm([lmiPmax 1 1 0], -betaP*eye(n));

    lmis = getlmis;
    options = [1e-7, 4000, -1, 0, 1];
    [tmin, xfeas] = feasp(lmis, options);

    if tmin > 0
        ok = false; P = []; K_cells = []; return;
    end

    P = dec2mat(lmis, xfeas, Pvar);

    if min(eig(P)) <= 0
        ok = false; P = []; K_cells = []; return;
    end

    K_cells = cell(Nv,1);
    for i = 1:Nv
        K_cells{i} = dec2mat(lmis, xfeas, Kvar_cells{i});
    end

    ok = true;
end