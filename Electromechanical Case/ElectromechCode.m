clc; clear; close all;

%% =========================
% Physical parameters
% =========================
k   = 1000;
m   = 1;
eta = 50;
eps = 2.8;
q0  = 1e-3;
re  = 5;

J = [ 0  1  0;
     -1  0  0;
      0  0  0];

R = [ 0    0     0;
      0   eta    0;
      0    0   1/re];

g    = [0; 0; 1/re];
A_JR = (J - R);

%% Input configuration
Umax = 5140;  % Maximum voltage [V]

% Input type: 'step' or 'sine'
input_type = 'step';
f_sine = 0.5;         % Frequency for sine [Hz]
phi_sine = 0;         % Phase for sine [rad]
t_step_start = 1.0;   % Step start time [s]

switch lower(input_type)
    case 'step'
        U = @(t) Umax * (t >= t_step_start);
    case 'sine'
        U = @(t) Umax * sin(2*pi*f_sine.*t + phi_sine);
    otherwise
        error('Unknown input_type. Use ''step'' or ''sine''.');
end

u = @(t) U(t);

fprintf('═══════════════════════════════════════════════\n');
fprintf('  ELECTROMECHANICAL OBSERVER COMPARISON\n');
fprintf('  Constant vs Scheduled Gain\n');
fprintf('═══════════════════════════════════════════════\n\n');
fprintf('Input type: %s\n', input_type);
if strcmp(input_type, 'step')
    fprintf('Step applied at t = %.1f s\n\n', t_step_start);
else
    fprintf('Sine wave: f = %.1f Hz, φ = %.2f rad\n\n', f_sine, phi_sine);
end

%% ============================================================
% Equilibrium-based parameter bounds
% ============================================================
fprintf('Computing equilibrium-based parameter bounds...\n');

Nu = 200;
uGrid = linspace(0, Umax, Nu);

s_of_u = nan(size(uGrid));
q_of_u = nan(size(uGrid));
Q_of_u = nan(size(uGrid));

for iu = 1:numel(uGrid)
    uu = uGrid(iu);
    Fs = @(s) 2*eps*(uu^2).*s.^3 - k*s + k*q0;

    if uu == 0
        s = q0;
    else
        sL = 0.05*q0;
        sR = 20*q0;
        fL = Fs(sL); fR = Fs(sR);

        expandCount = 0;
        while sign(fL) == sign(fR) && expandCount < 80
            sL = 0.5*sL;
            sR = 2.0*sR;
            fL = Fs(sL); fR = Fs(sR);
            expandCount = expandCount + 1;
        end

        if sign(fL) == sign(fR), continue; end
        s = fzero(Fs, [sL sR]);
    end

    if ~isfinite(s) || s <= 0, continue; end

    q = s - q0;
    Q = eps*s^4 * uu;

    s_of_u(iu) = s;
    q_of_u(iu) = q;
    Q_of_u(iu) = Q;
end

ok = isfinite(s_of_u) & isfinite(q_of_u) & isfinite(Q_of_u);
s_of_u = s_of_u(ok);
q_of_u = q_of_u(ok);
Q_of_u = Q_of_u(ok);

if isempty(s_of_u)
    error('No feasible equilibrium.');
end

% Compute TS parameter ranges from equilibrium
t1_vals = 10*(Q_of_u.^2)./(eps*(s_of_u.^6));
t2_vals = -4*(Q_of_u)./(eps*(s_of_u.^5));
t3_vals = 1./(eps*(s_of_u.^4));

inflate = 1.10;

t1_min = 0;
t1_max = inflate * max(t1_vals);

t2_abs = inflate * max(abs(t2_vals));
t2_min = -t2_abs;
t2_max =  t2_abs;

t3_min = (1/inflate) * min(t3_vals);
t3_max = inflate * max(t3_vals);

fprintf('\n═══ EQUILIBRIUM-BASED BOUNDS ═══\n');
fprintf('s ∈ [%.3e, %.3e] m\n', min(s_of_u), max(s_of_u));
fprintf('q ∈ [%.3e, %.3e] m\n', min(q_of_u), max(q_of_u));
fprintf('Q ∈ [%.3e, %.3e] C\n\n', min(Q_of_u), max(Q_of_u));

fprintf('═══ TS PARAMETER BOUNDS ═══\n');
fprintf('t1 ∈ [%.6e, %.6e]\n', t1_min, t1_max);
fprintf('t2 ∈ [%.6e, %.6e]\n', t2_min, t2_max);
fprintf('t3 ∈ [%.6e, %.6e]\n\n', t3_min, t3_max);

%% =========================
% Build 8 vertices (t1, t2, t3)
% =========================
Nv = 8;
Abar_vertices = cell(Nv,1);
Cbar_vertices = cell(Nv,1);

t1_grid = [t1_min, t1_max];
t2_grid = [t2_min, t2_max];
t3_grid = [t3_min, t3_max];

idx = 0;
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            idx = idx + 1;

            t1v = t1_grid(i1);
            t2v = t2_grid(i2);
            t3v = t3_grid(i3);

            % Hessian at vertex
            AH_v = [k + t1v,  0,    t2v;
                    0,        1/m,  0;
                    t2v,      0,    t3v];

            % System matrices
            Abar_vertices{idx} = A_JR * AH_v;
            Cbar_vertices{idx} = g' * AH_v;
        end
    end
end

fprintf('Created %d vertices\n\n', Nv);

kap0 = cellfun(@(A) cond(A), Abar_vertices);
fprintf('═══ FAMILY CONDITIONING ═══\n');
fprintf('  Nv = %d\n', Nv);
fprintf('  max κ = %.3e, median κ = %.3e\n\n', max(kap0), median(kap0));

%% =========================
% Common LMI parameters
% =========================
betaP = 1e6;
rho   = 1e-6;

%% =========================
% Find maximum lambda values
% =========================
fprintf('\n=== FINDING MAXIMUM LAMBDA VALUES ===\n');

lambda_max_const = find_max_lambda_const(Abar_vertices, Cbar_vertices, Nv, betaP, rho);
fprintf('lambda_max (constant):  %.6f\n', lambda_max_const);

lambda_max_sched = find_max_lambda_sched(Abar_vertices, Cbar_vertices, Nv, betaP, rho);
fprintf('lambda_max (scheduled): %.6f\n', lambda_max_sched);

%% =========================
% GROUP 1: Fair Comparison
% =========================
fprintf('\n=== GROUP 1: FAIR COMPARISON ===\n');

alpha = 0.1;
lambda_fair = alpha * lambda_max_const;

fprintf('Using lambda = %.6f for both observers\n', lambda_fair);

[ok1, P_const_fair, L_const_fair, ~] = design_observer_LMI_constant(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_fair, rho);
if ~ok1, error('Constant @ lambda_fair infeasible'); end

[ok2, P_sched_fair, K_cells_fair, ~] = solve_scheduled_LMI(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_fair, rho);
if ~ok2, error('Scheduled @ lambda_fair infeasible'); end

L_vertices_fair = cell(Nv,1);
for i = 1:Nv
    L_vertices_fair{i} = P_sched_fair \ K_cells_fair{i};
end

params_fair.L_vertices  = L_vertices_fair;
params_fair.k           = k;
params_fair.m           = m;
params_fair.eta         = eta;
params_fair.eps         = eps;
params_fair.q0          = q0;
params_fair.re          = re;
params_fair.t1_bounds   = [t1_min, t1_max];
params_fair.t2_bounds   = [t2_min, t2_max];
params_fair.t3_bounds   = [t3_min, t3_max];

[T1_const, q1_const, p1_const, Q1_const, qh1_const, ph1_const, Qh1_const, ...
    eq1_const, ep1_const, eQ1_const, L1_const_history, T1_const_history] = ...
    simulate_observer('constant', L_const_fair, [], k, m, eta, eps, q0, re, u, A_JR, g);

[T1_sched, q1_sched, p1_sched, Q1_sched, qh1_sched, ph1_sched, Qh1_sched, ...
    eq1_sched, ep1_sched, eQ1_sched, L1_sched_history, T1_sched_history] = ...
    simulate_observer('scheduled', [], params_fair, k, m, eta, eps, q0, re, u, A_JR, g);

fprintf('GROUP 1 simulations completed\n');

%% =========================
% GROUP 2: Scheduling Advantage
% =========================
fprintf('\n=== GROUP 2: SCHEDULING ADVANTAGE ===\n');

lambda_const_max = lambda_max_const;
fprintf('Constant:  lambda = %.6f (100%% of max)\n', lambda_const_max);

[ok3, P_const_max, L_const_max, ~] = design_observer_LMI_constant(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_const_max, rho);
if ~ok3, error('Constant @ lambda_max infeasible'); end

lambda_sched_fair = lambda_max_const;
fprintf('Scheduled: lambda = %.6f (same as constant max)\n', lambda_sched_fair);

[ok4, P_sched_fair2, K_cells_fair2, ~] = solve_scheduled_LMI(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_sched_fair, rho);
if ~ok4, error('Scheduled @ lambda_max_const infeasible'); end

L_vertices_fair2 = cell(Nv,1);
for i = 1:Nv
    L_vertices_fair2{i} = P_sched_fair2 \ K_cells_fair2{i};
end

params_fair2 = params_fair;
params_fair2.L_vertices = L_vertices_fair2;

lambda_sched_max = lambda_max_sched;
fprintf('Scheduled: lambda = %.6f (100%% of max_sched)\n', lambda_sched_max);

[ok_fail, ~, ~, ~] = design_observer_LMI_constant(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_sched_max, rho);
if ok_fail
    warning('Constant still feasible at lambda=%.3f!', lambda_sched_max);
else
    fprintf(' Constant FAILS at lambda=%.3f\n', lambda_sched_max);
end

[ok5, P_sched_max, K_cells_max, ~] = solve_scheduled_LMI(...
    Abar_vertices, Cbar_vertices, Nv, betaP, lambda_sched_max, rho);
if ~ok5, error('Scheduled @ lambda_sched_max infeasible'); end

L_vertices_max = cell(Nv,1);
for i = 1:Nv
    L_vertices_max{i} = P_sched_max \ K_cells_max{i};
end

params_max = params_fair;
params_max.L_vertices = L_vertices_max;

[T2_const, q2_const, p2_const, Q2_const, qh2_const, ph2_const, Qh2_const, ...
    eq2_const, ep2_const, eQ2_const, L2_const_history, T2_const_history] = ...
    simulate_observer('constant', L_const_max, [], k, m, eta, eps, q0, re, u, A_JR, g);

[T2_sched_fair, q2_sched_fair, p2_sched_fair, Q2_sched_fair, ...
    qh2_sched_fair, ph2_sched_fair, Qh2_sched_fair, ...
    eq2_sched_fair, ep2_sched_fair, eQ2_sched_fair, ...
    L2_sched_fair_history, T2_sched_fair_history] = ...
    simulate_observer('scheduled', [], params_fair2, k, m, eta, eps, q0, re, u, A_JR, g);

[T2_sched_max, q2_sched_max, p2_sched_max, Q2_sched_max, ...
    qh2_sched_max, ph2_sched_max, Qh2_sched_max, ...
    eq2_sched_max, ep2_sched_max, eQ2_sched_max, ...
    L2_sched_max_history, T2_sched_max_history] = ...
    simulate_observer('scheduled', [], params_max, k, m, eta, eps, q0, re, u, A_JR, g);

fprintf('GROUP 2 simulations completed\n');

%% =========================
% FIGURE GENERATION - GROUP 1
% Separate constant and scheduled gain subplots
% =========================
fprintf('\n=== GENERATING GROUP 1 FIGURES ===\n');
 
% Color scheme
col_real  = [0 0 0];           % Black
col_const = [0.85 0.33 0.10];  % Red
col_sched = [0 0.45 0.74];     % Blue
col_L1    = [0.93 0.69 0.13];  % Yellow/gold (for L2)
col_L2    = [0.47 0.67 0.19];  % Green (for L3)
col_L1c   = [0.50 0.50 0.50];  % Gray (for const L2)
col_L2c   = [0.30 0.30 0.30];  % Dark gray (for const L3)
 
% Font sizes
big_title = 22;
fs_title  = 20;
fs_label  = 19;
fs_legend = 17;
fs_tick   = 19;

% GROUP 1 - Main figure (2-column layout with 8 subplots)
fig1 = figure('Position', [50 50 1200 1400], 'Color', 'w');
sgtitle(sprintf('\\textbf{Group 1: Fair Comparison (Same $\\lambda$ = %.4f)}', lambda_fair), ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Row 1, Col 1: Position tracking
subplot(4,2,1);
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
subplot(4,2,2);
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
subplot(4,2,3);
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
subplot(4,2,4);
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

% Row 3, Col 1: Charge tracking
subplot(4,2,5);
plot(T1_const, Q1_const*1e9, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T1_const, Qh1_const*1e9, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/21):length(T1_const));
plot(T1_sched, Qh1_sched*1e9, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$Q$ [nC]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(e) Charge and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', 'Const.', 'Sched.'}, 'Location', 'southeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 3, Col 2: Charge error (zoom)
subplot(4,2,6);
plot(T1_const, zeros(size(T1_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T1_const, eQ1_const*1e9, '-.', 'LineWidth', 2.5, 'Color', col_const, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T1_const)/340):length(T1_const));
plot(T1_sched, eQ1_sched*1e9, ':', 'LineWidth', 2.5, 'Color', col_sched);
ylabel('$\tilde{Q}$ [nC]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(f) Charge Error (0--0.3 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', 'Const.', 'Sched.'}, 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.3]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 4, Col 1: Constant observer gains
subplot(4,2,7);
plot(T1_const_history, L1_const_history(:,1), '--', 'LineWidth', 2.5, 'Color', col_const); hold on;
plot(T1_const_history, L1_const_history(:,2), '--', 'LineWidth', 2.5, 'Color', col_L1c);
plot(T1_const_history, L1_const_history(:,3), '--', 'LineWidth', 2.5, 'Color', col_L2c);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(g) Observer Gains $L_1$, $L_2$, $L_3$ Constant', ...
    'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'$L_1$ const.', '$L_2$ const.', '$L_3$ const.'}, 'Interpreter', 'latex', ...
    'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 4, Col 2: Scheduled observer gains
subplot(4,2,8);
plot(T1_sched_history, L1_sched_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched); hold on;
plot(T1_sched_history, L1_sched_history(:,2), '-', 'LineWidth', 2.5, 'Color', col_L1);
plot(T1_sched_history, L1_sched_history(:,3), '-', 'LineWidth', 2.5, 'Color', col_L2);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(h) Observer Gains $L_1(t)$, $L_2(t)$, $L_3(t)$ Scheduled', ...
    'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'$L_1$ sched.', '$L_2$ sched.', '$L_3$ sched.'}, 'Interpreter', 'latex', ...
    'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Save Group 1 figure
savefig(fig1, 'Group1_Electromech_Fair_Comparison.fig');
exportgraphics(fig1, 'Group1_Electromech_Fair_Comparison.png', 'Resolution', 300);
exportgraphics(fig1, 'Group1_Electromech_Fair_Comparison.eps', 'ContentType', 'vector');
fprintf('✓ Saved: Group1_Electromech_Fair_Comparison.{fig,png,eps}\n');
 
%% =========================
% FIGURE GENERATION - GROUP 2
% Tracking and error plots only
% =========================
fprintf('\n=== GENERATING GROUP 2 FIGURES ===\n');
 
% GROUP 2 - Main figure (tracking and error plots)
fig2 = figure('Position', [50 50 1200 1200], 'Color', 'w');
sgtitle('\textbf{Group 2: Scheduling Advantage (Different $\lambda$)}', ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Colors for Group 2
col_const_2    = [0.85 0.33 0.10];  % Red (constant)
col_sched_fair = [0 0.45 0.74];     % Blue (scheduled @ same lambda)
col_sched_max  = [0.13 0.55 0.13];  % Green (scheduled @ max lambda)
 
% Row 1, Col 1: Position tracking
subplot(3,2,1);
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
subplot(3,2,2);
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
subplot(3,2,3);
plot(T2_const, p2_const*1e3, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T2_const, ph2_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/21):length(T2_const));
plot(T2_sched_fair, ph2_sched_fair*1e3, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, ph2_sched_max*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$p$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(c) Momentum and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 2, Col 2: Momentum error (zoom)
subplot(3,2,4);
plot(T2_const, zeros(size(T2_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T2_const, ep2_const*1e3, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/130):length(T2_const)); 
plot(T2_sched_fair, ep2_sched_fair*1e3, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, ep2_sched_max*1e3, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$\tilde{p}$ [g$\cdot$m/s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(d) Momentum Error (0--0.8 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.8]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);

% Row 3, Col 1: Charge tracking
subplot(3,2,5);
plot(T2_const, Q2_const*1e9, '-', 'LineWidth', 3.0, 'Color', col_real); hold on;
plot(T2_const, Qh2_const*1e9, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/21):length(T2_const));
plot(T2_sched_fair, Qh2_sched_fair*1e9, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, Qh2_sched_max*1e9, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$Q$ [nC]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(e) Charge and its estimation', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'Real', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Row 3, Col 2: Charge error (zoom)
subplot(3,2,6);
plot(T2_const, zeros(size(T2_const)), '-', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
plot(T2_const, eQ2_const*1e9, '-.', 'LineWidth', 2.5, 'Color', col_const_2, ...
      'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', ...
      'MarkerIndices', 1:round(length(T2_const)/130):length(T2_const)); 
plot(T2_sched_fair, eQ2_sched_fair*1e9, '-.', 'LineWidth', 2.5, 'Color', col_sched_fair);
plot(T2_sched_max, eQ2_sched_max*1e9, ':', 'LineWidth', 2.5, 'Color', col_sched_max);
ylabel('$\tilde{Q}$ [nC]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(f) Charge Error (0--0.8 s zoom)', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({'0 ref.', sprintf('Const. ($\\lambda$=%.3f)', lambda_const_max), ...
    sprintf('Sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
    sprintf('Sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 0.8]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Save Group 2 main figure
savefig(fig2, 'Group2_Electromech_Scheduling_Advantage.fig');
exportgraphics(fig2, 'Group2_Electromech_Scheduling_Advantage.png', 'Resolution', 300);
exportgraphics(fig2, 'Group2_Electromech_Scheduling_Advantage.eps', 'ContentType', 'vector');
fprintf('✓ Saved: Group2_Electromech_Scheduling_Advantage.{fig,png,eps}\n');
 
%% =========================
% GROUP 2 - Gains comparison figure
% =========================
fprintf('\n=== GENERATING GROUP 2 GAINS COMPARISON FIGURE ===\n');

% Colors for gains  
col_L2_sf = [0.53 0.81 0.92];  % Light blue
col_L2_sm = [0.56 0.93 0.56];  % Light green
col_L3_sf = [0.93 0.69 0.13];  % Gold
col_L3_sm = [0.80 0.60 0.10];  % Dark gold
 
fig3 = figure('Position', [200 100 700 1200], 'Color', 'w');
sgtitle('\textbf{Group 2: Observer Gains Comparison}', ...
    'Interpreter', 'latex', 'FontSize', fs_title+1, 'FontWeight', 'bold');
 
% Subplot 1: Constant gains (top)
subplot(3,1,1);
plot(T2_const_history, L2_const_history(:,1), '--', 'LineWidth', 2.5, 'Color', col_const_2); hold on;
plot(T2_const_history, L2_const_history(:,2), '--', 'LineWidth', 2.5, 'Color', col_L1c);
plot(T2_const_history, L2_const_history(:,3), '--', 'LineWidth', 2.5, 'Color', col_L2c);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(a) Constant Observer Gains', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ const.($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ const.($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_3$ const.($\\lambda$=%.3f)', lambda_sched_fair)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
 
% Subplot 2: Scheduled gains (both cases) (middle)
subplot(3,1,2);
plot(T2_sched_fair_history, L2_sched_fair_history(:,1), ':', 'LineWidth', 2.5, 'Color', col_sched_fair); hold on;
plot(T2_sched_fair_history, L2_sched_fair_history(:,2), ':', 'LineWidth', 2.5, 'Color', col_L2_sf);
plot(T2_sched_fair_history, L2_sched_fair_history(:,3), ':', 'LineWidth', 2.5, 'Color', col_L3_sf);
plot(T2_sched_max_history, L2_sched_max_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched_max);
plot(T2_sched_max_history, L2_sched_max_history(:,2), '-', 'LineWidth', 2.5, 'Color', col_L2_sm);
plot(T2_sched_max_history, L2_sched_max_history(:,3), '-', 'LineWidth', 2.5, 'Color', col_L3_sm);
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(b) Scheduled Observer Gains', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_3$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_1$ sched. ($\\lambda$=%.2f)', lambda_sched_max), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.2f)', lambda_sched_max), ...
        sprintf('$L_3$ sched. ($\\lambda$=%.2f)', lambda_sched_max)}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs_legend-1);
grid on; xlim([0 5]);
set(gca, 'FontSize', fs_tick, 'LineWidth', 1.2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
 
% Subplot 3: Zoom of scheduled gains (excluding L2 and L3 max components) (bottom)
subplot(3,1,3);
plot(T2_sched_fair_history, L2_sched_fair_history(:,1), ':', 'LineWidth', 2.5, 'Color', col_sched_fair); hold on;
plot(T2_sched_fair_history, L2_sched_fair_history(:,2), ':', 'LineWidth', 2.5, 'Color', col_L2_sf);
plot(T2_sched_fair_history, L2_sched_fair_history(:,3), ':', 'LineWidth', 2.5, 'Color', col_L3_sf);
plot(T2_sched_max_history, L2_sched_max_history(:,1), '-', 'LineWidth', 2.5, 'Color', col_sched_max);
% NOTE: Deliberately NOT plotting L2 and L3 from the max scheduled case
ylabel('Gain $L$', 'Interpreter', 'latex', 'FontSize', fs_label, 'FontWeight', 'bold');
title('(c) Scheduled Gains Zoom', 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
legend({sprintf('$L_1$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_2$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
        sprintf('$L_3$ sched. ($\\lambda$=%.3f)', lambda_sched_fair), ...
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
savefig(fig3, 'Group2_Electromech_Gains_Comparison.fig');
exportgraphics(fig3, 'Group2_Electromech_Gains_Comparison.png', 'Resolution', 300);
exportgraphics(fig3, 'Group2_Electromech_Gains_Comparison.eps', 'ContentType', 'vector');
fprintf('✓ Saved: Group2_Electromech_Gains_Comparison.{fig,png,eps}\n');
 
%% =========================
% Summary
% =========================
fprintf('Figure generation completed\n');
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

function [T, q, p, Q, qh, ph, Qh, eq, ep, eQ, L_history, T_history] = ...
    simulate_observer(mode, L_const, params, k, m, eta, eps, q0, re, u, A_JR, g)
    % Simulate electromechanical plant and observer dynamics

    % Plant dynamics
    fPlant = @(t,x) electromech_plant_rhs(t, x, A_JR, g, u, k, m, eps, q0);

    % Initial conditions
    x0    = [0; 0; 0];             % Plant initial state
    xhat0 = [2e-4; -2e-3; 5e-9];   % Observer initial state (with error)
    xx0   = [x0; xhat0];

    tspan = [0 5];
    opts  = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',1e-3);

    if strcmp(mode, 'constant')
        % Constant gain observer
        fObs = @(t,xx) obs_const(t, xx, fPlant, A_JR, g, u, k, m, eps, q0, re, L_const);
        [T, XX] = ode15s(fObs, tspan, xx0, opts);
        L_history = repmat(L_const(:).', length(T), 1);
        T_history = T;
    else
        % Scheduled gain observer
        global L_sched_global T_sched_global
        L_sched_global = [];
        T_sched_global = [];

        fObs = @(t,xx) obs_sched(t, xx, fPlant, A_JR, g, u, k, m, eps, q0, re, params);
        [T, XX] = ode15s(fObs, tspan, xx0, opts);

        L_history = L_sched_global;
        T_history = T_sched_global;
    end

    % Extract state components
    q  = XX(:,1);
    p  = XX(:,2);
    Q  = XX(:,3);
    qh = XX(:,4);
    ph = XX(:,5);
    Qh = XX(:,6);
    eq = q - qh;
    ep = p - ph;
    eQ = Q - Qh;
end

function dxx = obs_const(t, xx, fPlant, A_JR, g, u, k, m, eps, q0, re, L)
    % Constant gain observer dynamics
    x    = xx(1:3);
    xhat = xx(4:6);
    
    % Plant dynamics
    dx = fPlant(t, x);
    
    % Output and its estimate
    y  = y_electromech(x,    re, eps, q0);
    yh = y_electromech(xhat, re, eps, q0);
    
    % Observer dynamics
    qh = xhat(1); ph = xhat(2); Qh = xhat(3);
    sh = qh + q0;
    
    gHhat = [k*qh - 2*Qh^2/(eps*sh^5);
             ph/m;
             Qh/(eps*sh^4)];
    
    innov = y - yh;
    dxhat = A_JR * gHhat + g * u(t) + L * innov;
    
    dxx = [dx; dxhat];
end

function dxx = obs_sched(t, xx, fPlant, A_JR, g, u, k, m, eps, q0, re, params)
    % Scheduled gain observer dynamics
    global L_sched_global T_sched_global

    x    = xx(1:3);
    xhat = xx(4:6);
    
    % Plant dynamics
    dx = fPlant(t, x);

    % Output and its estimate
    y  = y_electromech(x,    re, eps, q0);
    yh = y_electromech(xhat, re, eps, q0);

    % Compute scheduled gain
    L = compute_scheduled_gain(xhat, params);

    % Store gain history
    L_sched_global = [L_sched_global; L(:).'];
    T_sched_global = [T_sched_global; t];

    % Observer dynamics
    qh = xhat(1); ph = xhat(2); Qh = xhat(3);
    sh = qh + q0;
    
    gHhat = [k*qh - 2*Qh^2/(eps*sh^5);
             ph/m;
             Qh/(eps*sh^4)];

    innov = y - yh;
    dxhat = A_JR * gHhat + g * u(t) + L * innov;
    
    dxx = [dx; dxhat];
end

function dx = electromech_plant_rhs(t, x, A_JR, g, u, k, m, eps, q0)
    % Electromechanical plant dynamics
    q = x(1); p = x(2); Q = x(3);
    s = q + q0;

    gH = [k*q - 2*Q^2/(eps*s^5);
          p/m;
          Q/(eps*s^4)];

    dx = A_JR * gH + g * u(t);
end

function y = y_electromech(x, re, eps, q0)
    % Output function: y = Q²/(re*eps*s^4)
    q = x(1); Q = x(3);
    s = q + q0;
    y = Q^2 / (re * eps * s^4);
end

function L = compute_scheduled_gain(xhat, params)
    % Compute scheduled observer gain using polytopic interpolation
    qh = xhat(1); Qh = xhat(3);
    sh = qh + params.q0;

    % Compute scheduling variables
    t1_hat = 10*Qh^2 / (params.eps*sh^6);
    t2_hat = -4*Qh   / (params.eps*sh^5);
    t3_hat = 1       / (params.eps*sh^4);

    % Normalized scheduling parameters (0 to 1)
    mu_t1 = (t1_hat - params.t1_bounds(1)) / (params.t1_bounds(2) - params.t1_bounds(1));
    mu_t2 = (t2_hat - params.t2_bounds(1)) / (params.t2_bounds(2) - params.t2_bounds(1));
    mu_t3 = (t3_hat - params.t3_bounds(1)) / (params.t3_bounds(2) - params.t3_bounds(1));

    % Saturation
    mu_t1 = max(0, min(1, mu_t1));
    mu_t2 = max(0, min(1, mu_t2));
    mu_t3 = max(0, min(1, mu_t3));

    % Interpolation weights
    w_t1L = 1 - mu_t1; w_t1H = mu_t1;
    w_t2L = 1 - mu_t2; w_t2H = mu_t2;
    w_t3L = 1 - mu_t3; w_t3H = mu_t3;

    % Compute weights for all 8 vertices
    h = zeros(8,1);
    idx = 0;
    for i1 = 1:2
        w1 = (i1==1)*w_t1L + (i1==2)*w_t1H;
        for i2 = 1:2
            w2 = (i2==1)*w_t2L + (i2==2)*w_t2H;
            for i3 = 1:2
                w3 = (i3==1)*w_t3L + (i3==2)*w_t3H;
                idx = idx + 1;
                h(idx) = w1 * w2 * w3;
            end
        end
    end

    % Normalize weights
    s = sum(h);
    if s <= 0
        h = ones(8,1)/8;
    else
        h = h / s;
    end

    % Compute interpolated gain
    L = zeros(3,1);
    for i = 1:8
        L = L + h(i) * params.L_vertices{i};
    end
end

function [ok, P, L, tmin] = design_observer_LMI_constant(Abar, Cbar, Nv, betaP, lambda, rho)
    % Design constant gain observer using LMI
    setlmis([]);
    n = 3;

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
    P = 0.5*(P+P');

    if min(eig(P)) <= 1e-12
        ok = false; P = []; K_cells = []; return;
    end

    K_cells = cell(Nv,1);
    for i = 1:Nv
        K_cells{i} = dec2mat(lmis, xfeas, Kvar_cells{i});
    end

    ok = true;
end