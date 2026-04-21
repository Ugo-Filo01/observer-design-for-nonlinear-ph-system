function [A,B,C, C_num,D,x_eq,u_eq,y_eq, i_eq, Hxx] = linearize_electromech(params, HessH_sym)
% LINEARIZE_ELECTROMECH  Linearization of the pH electromech plant
%   around the equilibrium associated to u_eq.

    syms q p Q u real

    % 1) Equilibrium
    [x_eq, u_eq, y_eq, i_eq] = equilibrium_electromech(params);
    q_eq = x_eq(1);
    p_eq = x_eq(2);
    Q_eq = x_eq(3);

    % 2) Evaluate Hessian at equilibrium
    Hxx = double(subs(HessH_sym, {q,p,Q}, {q_eq, p_eq, Q_eq}));

    % 3) Constant structure matrices
    eta     = params.eta;
    re      = params.re;
    epsllon = params.epsllon;
    q0 = params.q0;

    Cq   = epsllon * (q + q0)^4;
    dHdQ = Q / Cq;
    i = - (1/re) * dHdQ + (1/re) * u;

    J = [ 0  1  0;
         -1  0  0;
          0  0  0 ];

    R = diag([0, eta, 1/re]);
    g = [0; 0; 1/re];

    % 4) Linearized plant matrices
    A = (J - R) * Hxx;      % 3x3
    B = g;                  % 3x1
    C = g' * Hxx;         % 1x3
    D = 0;                  % scalar

    C_jac = jacobian(i, [q p Q]);   % oppure [q; p; Q]

    C_num = double(subs(C_jac, {q,p,Q,u},{q_eq, p_eq, Q_eq, u_eq}));
end
