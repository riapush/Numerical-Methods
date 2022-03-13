x = [-2:0.001:2];
f = @(x)(x.^2 - sin(10.*x));
cheb_grid = @(x)(0.183.*x.^3+1.001.*x.^2-1.074.*x-0.003);
uni_grid = @(x)((-7.25.*x.^3+16.*x.^2+21.72.*x)/16);

figure
plot(x, f(x), 'LineWidth', 2);
hold on
grid on
plot(x, cheb_grid(x), 'LineWidth', 2);
plot(x, uni_grid(x), 'LineWidth', 2);
plot(1.414,0.999, 'r*', 'LineWidth', 2);
plot(-1.414,2.999, 'r*', 'LineWidth', 2);
plot(-2, 4.91, 'go', 'LineWidth', 2);
plot(2, 3.09, 'go', 'LineWidth', 2);
legend('График функции', 'Полином Эрмита (Чебышевская сетка)', 'Полином Эрмита (Равномерная сетка)')
uni_max = abs(f(x)-uni_grid(x));
cheb_max = abs(f(x)-cheb_grid(x));

figure
semilogy(x, uni_max, 'LineWidth', 2);
hold on
grid on
semilogy(x, cheb_max, 'LineWidth', 2);
legend('Равномерная сетка', 'Чебышевская сетка');