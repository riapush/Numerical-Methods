iter = importdata('iter.csv');
I_num = importdata('integral.csv');
iter_non = importdata('iter_non.csv');
I_num_non = importdata('integral_non.csv');
lab3_iter = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\splits.csv');
lab3_err = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\err.csv');

F = @(x)(x .* x .* x ./ 3 + 0.1.*cos(10 .* x));
F_non = @(x)(0.1.*cos(10 .* x) + 0.75*x^(4/3));
a = 0;
b = 2;
I_acc = F(b)-F(a);
f = @(x)(x.^2 - sin(10.*x));
I_acc_non = F_non(b)-F_non(a);

eps = [10^1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12];

figure
loglog(eps, iter, 'LineWidth', 1);
hold on
grid on;
loglog(eps, iter_non, 'LineWidth', 1);
loglog(eps, lab3_iter, 'LineWidth', 1);
xlabel('Заданная точность');
ylabel('Количество разбиений');
legend('smooth', 'non-smooth', 'lab 3');
title('Зависимость кол-ва разбиений от заданной точности');

figure
loglog(eps, abs(I_acc-I_num), 'LineWidth', 1);
hold on
grid on;
loglog(eps, abs(I_acc_non-I_num_non), 'LineWidth', 1);
loglog(eps, abs(lab3_err - I_acc), 'LineWidth', 1);
loglog(eps,eps);
xlabel('Заданная точность');
ylabel('Погрешность');
title('Зависимость погрешности от заданной точности');
legend('smooth', 'non-smooth', 'lab 3', 'Location', 'SouthEast');

figure
loglog((b-a)./iter, abs(I_acc-I_num), 'LineWidth', 1);
hold on;
grid on;
loglog((b-a)./iter_non, abs(I_acc_non-I_num_non), 'LineWidth', 1);
legend('smooth', 'non-smooth', 'Location', 'SouthEast');
ylabel('Погрешность');
xlabel('h');
