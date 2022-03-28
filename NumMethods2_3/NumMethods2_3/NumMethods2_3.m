err = importdata("err.csv");
time = importdata("iter.csv");
err_opt = importdata("err_opt.csv");
time_opt = importdata("iter_opt.csv");
splits = importdata("splits.csv");


degree = linspace(1,13,13);
eps = 10.^-degree;

I = @(x)(x.^3/3+cos(10.*x)/10);
acc_I = I(2)-I(0);
err = abs(err-acc_I);
err_opt = abs(err_opt-acc_I);

figure
loglog(eps, err, 'b', 'LineWidth', 2);
hold on;
grid on;
loglog(eps, eps, 'LineWidth', 2);
%loglog(eps, err_opt, 'b--', 'LineWidth', 2);
title('Зависимость ошибки от заданной точности');
xlabel('Заданная точность');
ylabel('Ошибка приближения');
legend('Метод трапеций', 'Биссектриса', 'Location', 'SouthEast');

figure
loglog(eps, time, 'b--', 'LineWidth', 2);
hold on;
grid on;
loglog(eps, time_opt, 'LineWidth', 2);
xlabel('Заданная точность');
ylabel('Количество вычислений');
title('Зависимость количества вычислений от заданной точности');
legend( 'Метод трапеций', 'Оптимизация', 'Location', 'SouthWest');

figure
loglog(eps, splits);
hold on;
grid on;
xlabel('Заданная точность');
ylabel('Количество разбиений');
title('Зависимость количества разбиений от заданной точности');
