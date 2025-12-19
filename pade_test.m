syms x
G = [1/3+1/4 -1/4; -1/4 1/4];

C = [1 0; 0 2];

Is = [1/3 0]';

m0=G\Is;
m1=G\(C*m0);
m2=G\(C*m1);

sys = m0(1) + m1(1) * x + m2(1) * x^2;

psys = pade(sys, x, 0, 'Order', [ 2]);

disp(psys);