A = rand(100, 1);
B = A * transpose(A);

E = eig(B);

disp("A: ");
disp(A);
disp("B: ");
disp(B);
disp("Eig: ");
disp(E);

% ===== Histogram analysis =====
figure;
histogram(E, 10); % 分成 10 個桶 (bins)
xlabel('Eigenvalue');
ylabel('Frequency');
title('Histogram of Eigenvalues of B = A*A^T');
grid on;
