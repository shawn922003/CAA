clc; clear; close all;

% =============================================
% 1. 建立一個 SPD (Symmetric Positive Definite) 系統
% =============================================
n = 50;
A = gallery('poisson', n);   % 產生 SPD 矩陣 (size = n^2 x n^2)
b = ones(size(A,1), 1);
x0 = zeros(size(b));

% =============================================
% 2. 用共軛梯度法求解
% =============================================
tol = 1e-6;
max_iter = 1000;

[x_cg, iter, res] = conjugate_gradient(A, b, x0, tol, max_iter);

% =============================================
% 3. 用 MATLAB 內建直接法求解 (標準答案)
% =============================================
x_direct = A \ b;

% =============================================
% 4. 結果比較
% =============================================
fprintf('\n===== Conjugate Gradient vs Direct Solve =====\n');
fprintf('Iterations: %d\n', iter);
fprintf('Residual norm: %.3e\n', res);
fprintf('||x_cg - x_direct|| = %.3e\n', norm(x_cg - x_direct));
fprintf('||A*x_cg - b|| = %.3e\n', norm(A*x_cg - b));
fprintf('||A*x_direct - b|| = %.3e\n', norm(A*x_direct - b));

% =============================================
% 5. 顯示結果
% =============================================
disp('x_cg (Conjugate Gradient):');
disp(x_cg);

disp('x_direct (A\\b):');
disp(x_direct);

max_iter = 1000;
[x_pcg, flag, relres, iter, resvec] = pcg(A, b, tol, max_iter);

% 比較標準解
x_direct = A \ b;

fprintf('PCG finished in %d iterations, relative residual = %.3e\n', iter, relres);
fprintf('||x_pcg - x_direct|| = %.3e\n', norm(x_pcg - x_direct));





