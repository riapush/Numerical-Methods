h = importdata('local_h.csv');
local_e = importdata('local.csv');
global_e = importdata('global.csv');
num_sol = importdata('sol.csv');
pert = importdata('pert.csv');
pert_small = importdata('pert_small.csv');
pert_glob = importdata('pert_glob.csv');

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


figure
xlabel('x');
hold on;
grid on;
ylabel('y');
title('Влияние ошибки в исходных данных на решение');
for i = 1:1:4
    plot(x, pert(i,:), 'LineWidth', 1);
end
plot(x, y(x), 'LineWidth', 1);
legend('8%', '16%', '32%', '64%', 'Точное решение', 'Location', 'NorthWest');

figure
xlabel('Возмущение');
hold on;
grid on;
ylabel('Относительная ошибка');
per = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 20.48, 40.96, 81.92];
k = b - (b-a)/32;
title({'Влияние ошибки в исходных данных на решение' '(все данные в процентах)'});
for i = 1:1:14
    semilogx(per(1,i), abs(pert_glob(1,i)-y(k))*100/abs(y(k)), '*');
end
