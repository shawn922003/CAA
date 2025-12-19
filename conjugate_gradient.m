function [x, iter, res_norm] = conjugate_gradient(A, b, x0, tol, max_iter)
% ---------------------------------------------------------
% Conjugate Gradient Method for solving A*x = b
%
% Input:
%   A         - Symmetric Positive Definite (SPD) matrix
%   b         - Right-hand side vector
%   x0        - Initial guess (can be zeros(size(b)))
%   tol       - Convergence tolerance (e.g. 1e-6)
%   max_iter  - Maximum number of iterations
%
% Output:
%   x         - Solution vector
%   iter      - Number of iterations performed
%   res_norm  - Final residual norm
%
% Example:
%   A = [4,1;1,3];
%   b = [1;2];
%   [x, iter, res] = conjugate_gradient(A, b, zeros(2,1), 1e-6, 100)
% ---------------------------------------------------------

    % Initialization
    x = x0;
    r = b - A * x;          % initial residual
    p = r;                  % initial direction
    rs_old = r' * r;        % squared residual norm

    for iter = 1:max_iter
        Ap = A * p;
        alpha = rs_old / (p' * Ap);   % step size

        % Update x and r
        x = x + alpha * p;
        r = r - alpha * Ap;

        rs_new = r' * r;
        res_norm = sqrt(rs_new);

        disp(["iteration" iter "residual norm: " res_norm]);

        % Check convergence
        if res_norm < tol
            fprintf('Converged at iteration %d, residual = %.3e\n', iter, res_norm);
            break;
        end

        % Update direction
        beta = rs_new / rs_old;
        p = r + beta * p;

        rs_old = rs_new;
    end

    % If not converged, report
    if iter == max_iter && res_norm >= tol
        fprintf('Max iterations reached, residual = %.3e\n', res_norm);
    end
end
