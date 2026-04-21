function [H_sym, gradH_sym, HessH_sym, lambdas_sym] = H_electromech(params)
% H_ELECTROMECH  Symbolic Hamiltonian, gradient and Hessian for the
% electro-mechanical system, with parameters substituted.
%
%   [H_sym, gradH_sym, HessH_sym, lambdas_sym] = H_electromech(params)
%
% params must contain fields: m, K, eps, q0.

    % declare symbolic state variables
    syms q p Q real
    % syms m K epsllon q0 real

    % read parameters
    m       = params.m;
    K       = params.K;
    epsllon = params.epsllon;
    q0      = params.q0;

    % capacitance
    Cq      = epsllon*(q + q0)^4;

    % symbolic Hamiltonian with numeric parameters
    H_sym   = 0.5*K*q^2 + 0.5*(p^2)/m + 0.5*(Q^2)/Cq;

    disp('Electro-mechanical Hamiltonian H(q,p,Q):');
    disp(H_sym);

    % symbolic gradient
    gradH_sym = jacobian(H_sym, [q p Q]).';   % 3x1

    % symbolic Hessian (Jacobian of the gradient)
    HessH_sym = jacobian(gradH_sym, [q p Q]); % 3x3

    % symbolic eigenvalues of the Hessian
    lambdas_sym = eig(HessH_sym);
    % disp('Electro-mechanical system eig:');
    % disp(lambdas_sym);
    detH = simplify(det(HessH_sym));
    disp('Hessian determinant: ')
    disp(detH);
end
