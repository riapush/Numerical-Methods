f = @(x)exp(x);
x = linspace(0,1,1000);

csv1 = readmatrix("solution1.csv");
csv2 = readmatrix("solution2.csv");

x1 = csv1(:,1);
y1 = csv1(:,2);
x2 = csv2(:,1);
y2 = csv2(:,2);

figure
plot(x, f(x), 'LineWidth', 1)
hold on
plot(x1, y1, '--*', 'LineWidth', 1)
plot(x2, y2, '--o', 'LineWidth', 1)
title("Некоторые решения")
legend("Точное решение", "6 точек (h=0.2)", "11 точек (h=0.1)",'location', 'northwest')
grid on

figure 
plot(x1, abs(f(x1)-y1), 'LineWidth', 1)
hold on
plot(x2, abs(f(x2)-y2), 'LineWidth', 1)
title("Ошибки")
legend("6 точек (h=0.2)", "11 точек (h=0.1)", 'location', 'northwest')
grid on



csv5 = readmatrix("per_A.csv");
da = csv5(:,1);
err = csv5(:,2);

figure
subplot(2,1,1)
loglog(da, err)
title("Зависимость фактической ошибки от величины возмущения A")
xlabel("\delta A")
ylabel("Фактическая ошибка")
grid on

csv3 = readmatrix("per_B.csv");
db = csv3(:,1);
err = csv3(:,2);

subplot(2,1,2)
loglog(db, err)
title("Зависимость фактической ошибки от величины возмущения B")
xlabel("\delta B")
ylabel("Фактическая ошибка")
grid on

csv4 = readmatrix("error_h.csv");
h = csv4(:,1);
err = csv4(:,2);

figure
loglog(abs(h), err, abs(h), h.^4, 'LineWidth', 1)
title("Зависимость фактической ошибки от h")
xlabel("h")
ylabel("Фактическая ошибка")
legend("Зависимость", "y=h^4", 'Location', 'northwest')
grid on