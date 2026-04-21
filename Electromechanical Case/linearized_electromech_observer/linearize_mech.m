function [A,B,C,D,x_eq_mech,u_eq_mech,y_eq_mech, Hxx, g] = linearize_mech(params, gradH_sym, ...
    HessH_sym, g_sym, grad_g)
% LINEARIZE_ELECTROMECH  Linearization of the pH electromech plant
%   around the equilibrium associated to u_eq.

    syms q p real

    % Equilibrium
    [x_eq_mech, u_eq_mech, y_eq_mech] = equilibrium_mech(params);
    q_eq = x_eq_mech(1);
    p_eq = x_eq_mech(2);

    % Evaluate Hessian and gradient at equilibrium
    gradH = double(subs(gradH_sym, {q,p}, {q_eq, p_eq}));
    Hxx   = double(subs(HessH_sym, {q,p}, {q_eq, p_eq}));

    eta  = params.eta;

    J      = [ 0  1;
              -1  0 ];
    R      = diag([0, eta]);
    g      = double(subs(g_sym, {q,p}, {q_eq, p_eq}));
    grad_g = double(subs(grad_g, {q,p}, {q_eq, p_eq}));

    % % Linearized plant matrices
    A = (J - R) * Hxx + grad_g * u_eq_mech;      % 2x2
    B = g;                                  % 2x1
    % C = g' * Hxx + grad_g' * gradH;         % 1x2
    C = g.' * Hxx + gradH.' * grad_g;       % 1x2
    D = 0;                                  % scalar

    % the previous structure is not pH since it has spurious terms
    % A = (J - R) * Hxx;    % 2x2
    % B = g;                % 2x1
    % C = g' * Hxx;         % 1x2
    % D = 0;                % scalar

    % with the notation I have 0 error so the linear system is tracking
    % prefectly
    % A = (J - R) * Hxx;    % 2x2
    % B = g;                % 2x1
    % C = g' * Hxx;         % 1x2
    % D = 0;                % scalar
end
