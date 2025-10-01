Vd = 1;
C = 1;
R = 1;
delta_t = 0.1;

% 定義迭代函式 (用 nested function 或 anonymous function 皆可)
forward_next_step = @(v_curr) (1 - delta_t / R / C) * v_curr + delta_t / R / C * Vd;
backward_next_step = @(v_curr) (R * C * v_curr + Vd * delta_t) / (delta_t + R * C);
trapezoid_next_step = @(v_curr) ((1 - delta_t / 2 / R /C) * v_curr + delta_t / R / C * Vd) / (1 + delta_t / 2 / R / C);

step = 100;
forward_vol_t = zeros(step,1);   % 初始化
for i = 2:step
    forward_vol_t(i) = forward_next_step(forward_vol_t(i - 1));
end

backward_vol_t = zeros(step, 1);
for i = 2:step
    backward_vol_t(i) = backward_next_step(backward_vol_t(i - 1));
end

trapezoid_vol_t = zeros(step, 1);
for i = 2: step
    trapezoid_vol_t(i) = trapezoid_next_step(trapezoid_vol_t(i - 1));
end

% 繪圖
time = (0:step-1) * delta_t;   % 時間軸
plot(time, forward_vol_t, 'LineWidth', 1.5);
plot(time, backward_vol_t, 'LineWidth', 1.5);
plot(time, trapezoid_vol_t, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('RC Charging Curve');
grid on;


