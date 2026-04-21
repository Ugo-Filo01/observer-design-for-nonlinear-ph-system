function [x_eq, u_op, y_eq, i_eq] = equilibrium_electromech(params)
% EQUILIBRIUM_ELECTROMECH  Find equilibrium (q,p,Q) for constant input u_eq.

    K       = params.K;
    epsllon = params.epsllon;
    q0      = params.q0;
    u_op    = params.u_op;
    re      = params.re;

    % Solve for s = q + q0
    fun_s = @(s) 2*epsllon*u_op^2 .* s.^3 - K*s + K*q0;

    s0  = q0;             % initial guess near nominal gap
    s_eq = fzero(fun_s, s0);

    q_eq = s_eq - q0;
    p_eq = 0;
    Q_eq = epsllon * s_eq^4 * u_op;

    x_eq = [q_eq; p_eq; Q_eq];
    y_eq = (Q_eq) / (re * epsllon * (q_eq + q0)^4);

    Cq_eq   = epsllon * s_eq^4;
    dHdQ_eq = Q_eq / Cq_eq;
    i_eq = - (1/re) * dHdQ_eq + (1/re) * u_op;

end
