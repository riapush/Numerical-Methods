%close all;

h = importdata('local_h.csv');
local_e = importdata('local.csv');
global_e = importdata('global.csv');
num_sol = importdata('sol.csv');

y = @(x)(x.*(x.*x+1));
a = 0;
b = 2;
x = (a:(b-a)/32:b);

figure
xlabel('x');
hold on;
grid on;
plot(x, y(x), 'LineWidth', 2);
plot(x, num_sol, '--', 'LineWidth', 2);
ylabel('y');
legend('Точное решение', 'Численное решение', 'Location', 'SouthEast');
title('Точное и численное решение при h = 0.0625');

figure
xlabel('x');
hold on;
grid on;
semilogy(x, abs(y(x)-num_sol), 'LineWidth', 2);
ylabel('Абсолютная ошибка');
title('График ошибки');

figure
xlabel('h');
loglog(h, local_e, 'LineWidth', 2);
hold on;
grid on;
loglog(h, global_e, 'LineWidth', 2);
ylabel('Локальная и глобальная ошибка');
title('График локальной ошибки');
legend('local', 'global')