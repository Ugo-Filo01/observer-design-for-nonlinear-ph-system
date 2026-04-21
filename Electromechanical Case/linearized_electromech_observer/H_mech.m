function [H_sym, gradH_sym, HessH_sym, lambdas_sym, g_sym, grad_g] = H_mech(params)
% H_ELECTROMECH  Symbolic Hamiltonian, gradient and Hessian for the
% electro-mechanical system, with parameters substituted.
%
%   [H_sym, gradH_sym, HessH_sym, lambdas_sym] = H_electromech(params)
%
% params must contain fields: m, K, eps, q0.

    % declare symbolic state variables
    syms q p real
    % syms m K epsllon q0 real

    % read parameters
    m       = params.m;
    K       = params.K;
    epsllon = params.epsllon;
    q0      = params.q0;

    g_sym = [0; 2 * epsllon * (q + q0)^3];

    % symbolic Hamiltonian with numeric parameters
    H_sym   = 0.5*K*q^2 + 0.5*(p^2)/m;

    disp('Mechanical Hamiltonian H(q,p,):');
    disp(H_sym);

    % symbolic gradient
    gradH_sym = jacobian(H_sym, [q p]).';   % 3x1

    % symbolic Hessian (Jacobian of the gradient)
    HessH_sym = jacobian(gradH_sym, [q p]); % 3x3

    % symbolic eigenvalues of the Hessian
    lambdas_sym = eig(HessH_sym);
    % disp('Electro-mechanical system eig:');
    % disp(lambdas_sym);
    detH = simplify(det(HessH_sym));
    disp('Hessian determinant: ')
    disp(detH);

    grad_g = jacobian(g_sym, [q, p]);
end