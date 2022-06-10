dots_trapezoid = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\dots.csv');
dots_gauss = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_4\NumMethods2_4\dots.csv');

splits_trapezoid = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\splits.csv');
splits_trapezoid_non = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\splits_non.csv');
err_trapezoid = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\err.csv');
err_trapezoid_non = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\err_non.csv');

err_gauss = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_4\NumMethods2_4\err.csv');
err_gauss_non = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_4\NumMethods2_4\err_non.csv');
splits_gauss = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_4\NumMethods2_4\splits.csv');
splits_gauss_non = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_4\NumMethods2_4\splits_non.csv');

eps = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12];

h = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\h.csv');
err_check = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2_3\NumMethods2_3\err_check.csv');
figure
loglog(h,err_check);
hold on
grid on
xlabel('Длина отрезка');
ylabel('Погрешность');


figure
loglog(eps, dots_trapezoid, 'LineWidth', 2);
hold on
grid on
loglog(eps, dots_gauss, 'LineWidth', 2);
legend('Метод трапеции', 'Квадратурная формула Гаусса');
title('Зависимость объема вычислений от заданной точности');
xlabel('Заданная точность');
ylabel('Количество вызовов функции f');

figure
loglog(eps, splits_gauss, 'LineWidth', 1);
hold on
grid on;
loglog(eps, splits_gauss_non, 'LineWidth', 1);
loglog(eps, splits_trapezoid, 'LineWidth', 1);
loglog(eps, splits_trapezoid_non, 'LineWidth', 1);
xlabel('Заданная точность');
ylabel('Количество разбиений');
legend('Гладкая Гаусс', 'Негладкая Гаусс', 'Гладкая метод трапеции', 'Негладкая метод трапеции');
title('Зависимость кол-ва разбиений от заданной точности');

figure
loglog(eps, err_gauss, 'g', 'LineWidth', 2);
hold on
grid on;
loglog(eps, err_gauss_non, 'g--', 'LineWidth', 2);
loglog(eps, err_trapezoid, 'b', 'LineWidth', 2);
loglog(eps, err_trapezoid_non, 'b--', 'LineWidth', 2);
xlabel('Заданная точность');
ylabel('Погрешность');
loglog(eps, eps, 'r',  'LineWidth', 1);
legend('Гладкая Гаусс', 'Негладкая Гаусс', 'Гладкая метод трапеции', 'Негладкая метод трапеции', 'Биссектриса', 'Location', 'SouthEast');
title('Зависимость погрешности от заданной точности');