A = [0 2; 3 4];
b = [1 1];
x = guassian_elimination(A, b);
disp(x);
disp(A \ b')