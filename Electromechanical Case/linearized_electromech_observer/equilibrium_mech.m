function [x_eq_mech, u_op_mech, y_eq_mch] = equilibrium_mech(params)
% EQUILIBRIUM_ELECTROMECH  Find equilibrium (q,p,Q) for constant input u_eq.

    K       = params.K;
    epsllon = params.epsllon;
    q0      = params.q0;
    u_op_mech    = params.u_op;
    m            = params.m;

    fun_s = @(s) 2*epsllon*u_op_mech^2 .* s.^3 - K*s + K*q0;

    s0   = q0;
    s_eq = fzero(fun_s, s0);

    q_eq = s_eq - q0;
    p_eq = 0;

    x_eq_mech = [q_eq; p_eq];
    y_eq_mch = (p_eq / m) * 2 * epsllon * (q_eq + q0)^3;
    u_op_mech = u_op_mech^2;

end
