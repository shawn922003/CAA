syms r1 r2 c1 c2 s t real

G = [1/r1+1/r2, -1/r2;
     -1/r2,      1/r2];
C = [c1 0; 0 c2];
I = [1/r1; 0];

K = 6;
M = size(G,1);

taylor_I = taylor_coeffs(I, s, 0, K-1);   % 假設回傳 K×M，每列是長度 M 的係數

m = sym(zeros(K, M));

% 關鍵：把 1×2 轉成 2×1 再做左除；解完再轉回 1×2
m(1,:) = G \ taylor_I(:, 1);

for i = 2:K
    m(i, :) = G\(taylor_I(:, i) - C * m(i - 1, :)');
end

disp("================ moment ================");
disp(m);

V = sym(zeros(M, 1)); 
for i = 1:K
    V = V + m(i, :)' * (s ^ (i-1));
end

disp("================ moment V ================")
disp(V);


disp("================ real V ================")
disp((G + s * C)\I);


% ---- 用 moments 做 Padé：每個分量各自處理
% 把 KxM 的 m 轉成 MxK，mom(j,:) = 該分量的 [m0 m1 ... m_{K-1}]
mom = m.';   % M x K

% 你現在 K=3 -> 只能做 [0/2]（需要 L+M+1 = 3 個係數）
L1 = 0; M1 = 2;
L2 = 0; M2 = 2;

% 第1個分量
[num1, den1] = pade_from_moments_sym(mom(1,1:(L1+M1+1)), L1, M1);
V1_pade = simplify(poly2sym(fliplr(num1), s) / poly2sym(fliplr(den1), s));

% 第2個分量
[num2, den2] = pade_from_moments_sym(mom(2,1:(L2+M2+1)), L2, M2);
V2_pade = simplify(poly2sym(fliplr(num2), s) / poly2sym(fliplr(den2), s));

disp("=============== Padé from moments ===============");
disp("V1_pade(s) = "); disp(V1_pade);
disp("V2_pade(s) = "); disp(V2_pade);



% === 反拉普拉斯 + 代入 r1=r2=c1=c2=1 ===
assumeAlso(t>=0);  % 幫助 ilaplace/Heaviside 符號推理

% 精確解（2x1）
V_exact = simplify((G + s*C)\I);

% 先把所有參數代 1
subs_map = struct('r1',1,'r2',1,'c1',1,'c2',1);
V1p_sub = simplify(subs(V1_pade, subs_map));
V2p_sub = simplify(subs(V2_pade, subs_map));
Vex_sub = simplify(subs(V_exact,   subs_map));   % 2x1

% 做反拉普拉斯
v1_pade_t = simplify(ilaplace(V1p_sub, s, t));
v2_pade_t = simplify(ilaplace(V2p_sub, s, t));
v1_exact_t = simplify(ilaplace(Vex_sub(1), s, t));
v2_exact_t = simplify(ilaplace(Vex_sub(2), s, t));

disp("=== time-domain (Padé) with all params = 1 ===")
disp(v1_pade_t); disp(v2_pade_t);
disp("=== time-domain (Exact) with all params = 1 ===")
disp(v1_exact_t); disp(v2_exact_t);

% === 繪圖（Padé vs Exact）===
% 轉成可數值評估的 function handle
v1p_fun = matlabFunction(v1_pade_t, 'Vars', t);
v2p_fun = matlabFunction(v2_pade_t, 'Vars', t);
v1e_fun = matlabFunction(v1_exact_t, 'Vars', t);
v2e_fun = matlabFunction(v2_exact_t, 'Vars', t);

tgrid = linspace(0, 10, 1000);

figure;  % 第一個分量
plot(tgrid, v1p_fun(tgrid), 'LineWidth', 1.5); hold on;
plot(tgrid, v1e_fun(tgrid), 'LineWidth', 1.5, 'LineStyle','--');
grid on; xlabel('t (s)'); ylabel('v_1(t)');
title('Response of V_1(t): Padé vs Exact (r1=r2=c1=c2=1)');
legend('Padé [1/2]', 'Exact');

figure;  % 第二個分量
plot(tgrid, v2p_fun(tgrid), 'LineWidth', 1.5); hold on;
plot(tgrid, v2e_fun(tgrid), 'LineWidth', 1.5, 'LineStyle','--');
grid on; xlabel('t (s)'); ylabel('v_2(t)');
title('Response of V_2(t): Padé vs Exact (r1=r2=c1=c2=1)');
legend('Padé [0/2]', 'Exact');