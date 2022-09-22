close all;

h = importdata('local_h.csv');
local_e = importdata('local.csv');
global_e = importdata('global.csv');
num_sol = importdata('sol.csv');
per8 = importdata('per.csv');
per16 = importdata('per2.csv');
per32 = importdata('per3.csv');
per_i = importdata('per_i.csv');
vol = importdata('vol.csv');
lab5_local = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_5\NumMethods2_5\local.csv');
lab5_global = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_5\NumMethods2_5\global.csv');
lab5_vol = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_5\NumMethods2_5\vol.csv');

%y = @(x)(x.*(x.*x+1));
y = @(x)((2 .* x + 1).*log(2 .* x + 1) + 1);
a = 0;
b = 4;
x = (a:(b-a)/32:b);

figure
xlabel('x');
hold on;
grid on;
plot(x, y(x), 'LineWidth', 2);
plot(x, num_sol, '--', 'LineWidth', 2);
ylabel('y');
legend('Точное решение', 'Численное решение', 'Location', 'SouthEast');
title('Точное и численное решение при h = 0.25');

figure
xlabel('x');
hold on;
grid on;
semilogy(x, abs(y(x)-num_sol), 'LineWidth', 2);
ylabel('Абсолютная ошибка');
title('График ошибки');

figure
loglog(h, local_e, 'LineWidth', 2);
hold on;
grid on;
loglog(h, global_e, 'LineWidth', 2);
loglog(h, lab5_local, 'LineWidth', 2);
loglog(h, lab5_global, 'LineWidth', 2);
xlabel('h');
ylabel('Ошибка');
title('График ошибки');
legend('local', 'global', 'lab5 local', 'lab5 global', 'Location', 'SouthEast');

figure
xlabel('x');
hold on;
grid on;
ylabel('y');
plot(x, y(x), 'LineWidth', 1);
for i = 1:1:11
    plot(x, per32(i,:), 'LineWidth', 1)
end
legend('10^{-1}', '10^{-2}', '10^{-3}', '10^{-4}','10^{-5}','10^{-6}','10^{-7}', '10^{-8}', '10^{-9}', '10^{-10}', '10^{-11}', 'Location', 'NorthWest');

dots1 = zeros(9);
dots1(1,1) = a;
h1 = (b-a)/8;
for i = 2:1:9
    dots1(1,i) = dots1(1,i-1) + h1;
end
dots2 = zeros(17);
dots2(1,1) = a;
h2 = (b-a)/16;
for i = 2:1:17
    dots2(1,i) = dots2(1,i-1) + h2;
end
dots3 = zeros(33);
dots3(1,1) = a;
h3 = (b-a)/32;
for i = 2:1:33
    dots3(1,i) = dots3(1,i-1) + h3;
end

figure
for i = 1:1:11
    error = max(abs(y(dots1(1,:))-per(i,:)));
    loglog(10^-i, error, 'b*', 'LineWidth', 1);
    hold on;
    
    error = max(abs(y(dots2(1,:))-per16(i,:)));
    loglog(10^-i, error, 'r*', 'LineWidth', 1);
    
    error = max(abs(y(dots3(1,:))-per32(i,:)));
    loglog(10^-i, error, 'g*', 'LineWidth', 1);
end

grid on;
xlabel('Возмущение');
ylabel('Ошибка');

step = [8,16,32,64,128];
vol = vol./4;
figure
plot(step, vol, 'b', 'LineWidth', 1);
hold on;
grid on;
plot(step, lab5_vol, 'g', 'LineWidth', 1);
xlabel('h');
ylabel('Количество вызовов функции');
legend('Adams', 'RK', 'location', 'northwest');