err = importdata("err.csv");
time = importdata("iter.csv");
err_opt = importdata("err_opt.csv");
time_opt = importdata("iter_opt.csv");


degree = linspace(1,13,13);
eps = 10.^-degree;

I = @(x)(x.^3/3+cos(10.*x)/10);
acc_I = I(2)-I(0);
err = abs(err-acc_I);
err_opt = abs(err_opt-acc_I);

figure
loglog(eps, err);
hold on;
grid on;
loglog(eps,eps);
loglog(eps, err_opt);
title('Зависимость ошибки от заданной точности');
xlabel('Заданная точность');
ylabel('Ошибка приближения');

figure
loglog(eps, time, 'b--');
hold on;
grid on;
loglog(eps, time_opt, 'r');
xlabel('Заданная точность');
ylabel('Количество итераций');
title('Зависимость количества итераций от заданной точности');
