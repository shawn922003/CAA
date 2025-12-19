n = 50;
A = randn(n);
q1 = randn(n,1);
m = 10;

[Q,H] = arnoldi(A, q1, m);

% 檢查正交性
disp(norm(Q'*Q - eye(m), 'fro'));

% 檢查 Arnoldi 關係：A*Q(:,1:m-1) ≈ Q*H
disp(norm(A*Q(:,1:m-1) - Q*H, 'fro'));
