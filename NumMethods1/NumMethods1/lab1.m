p = @(x)(x.^6 + x.^5 - 13 .* x.^3 - 9 .* x + 2);
p_d1 = @(x)(6.*x.^5+5.*x.^4-(13*3).*x.^2-9);
p_d2 = @(x)((6*5).*x.^4+20.*x.^3-(13*6).*x);
p_root = fzero(p,[1,3])
t = @(x)(5.^x-6.*x-7);
t_d1 = @(x)((5.^x)*log(5)-6)
t_d2 = @(x)((5.^x)*(log(5))^2)
t_root = fzero(t, [1,2])

a = [1:0.000005:3];
figure
hold on
grid on
title('polynom')
plot(a,p(a), 'r', 'LineWidth', 2)
plot(a,p_d1(a), 'g', 'LineWidth', 2)
plot(a,p_d2(a), 'm', 'LineWidth', 2)
legend('polynom', 'first deriative', 'second deriative', 'Location', 'NorthWest')

a = [0:0.000005:3];
figure
hold on
grid on
title('transcendental')
plot(a, t(a), 'r', 'LineWidth', 2)
plot(a, t_d1(a), 'g', 'LineWidth', 2)
plot(a, t_d2(a), 'm', 'LineWidth', 2)
line([1,1],[-50,350], 'Color', 'black', 'LineWidth', 1)
line([2,2],[-50,350], 'Color', 'black', 'LineWidth', 1)
legend('transcendental', 'first deriative', 'second deriative', 'Location', 'NorthWest')





% построим зависимость абс.погрешности от кол-ва итераций.
% для этого импортируем данные
p_bisection_converge = importdata('p_bisection_converge.csv')
iter_p_b_con = p_bisection_converge(:,1)';
val_p_b_con = abs(p_bisection_converge(:,2)' - p_root);

t_bisection_converge = importdata('t_bisection_converge.csv');
iter_t_b_con = t_bisection_converge(:,1)'
val_t_b_con = abs(t_bisection_converge(:,2)' - t_root)

p_chord_converge = importdata('p_chord_converge.csv');
iter_p_c_con = p_chord_converge(:,1)';
val_p_c_con = abs(p_chord_converge(:,2)' - p_root);

t_chord_converge = importdata('t_chord_converge.csv');
iter_t_c_con = t_chord_converge(:,1)';
val_t_c_con = abs(t_chord_converge(:,2)' - t_root);

figure
semilogy(iter_p_b_con, val_p_b_con, 'Color','r','LineWidth',2)
hold on
grid on
xlabel ('k');
ylabel (' |x* - x_k| ');
title ('Зависимость абсолютной погрешности от кол-ва итераций');
semilogy(iter_t_b_con, val_t_b_con, 'Color','g','LineWidth',2);
semilogy(iter_p_c_con, val_p_c_con, 'Color','b','LineWidth',2);
semilogy(iter_t_c_con, val_t_c_con, 'Color','m','LineWidth',2);
legend('bisection polynom','bisection transcendental', 'chord polynom', 'chord transcendental', 'Location', 'NorthEast');

% построим зависимость точности от кол-ва итераций
% для этого импортируем данные

p_bisection_acc = importdata('p_bisection_accuracy.csv');
iter_p_b_acc = p_bisection_acc(:,1)'; % необходимое число итераций
iter_p_b_acc(1) = [];
p_b_acc = p_bisection_acc(:,2)'; % степени 10
p_b_acc(1) = [];

p_b_val = abs(p_bisection_acc(:,3)' - p_root);
p_b_val(1) = [];

t_bisection_acc = importdata('t_bisection_accuracy.csv');
iter_t_b_acc = t_bisection_acc(:,1)';
iter_t_b_acc(1) = [];
t_b_acc = t_bisection_acc(:,2)';
t_b_acc(1) = [];

t_b_val = abs(t_bisection_acc(:,3)' - t_root);
t_b_val(1) = [];

p_chord_acc = importdata('p_chord_accuracy.csv');
iter_p_c_acc = p_chord_acc(:,1)'
iter_p_c_acc(1) = [];
p_c_acc = p_chord_acc(:,2)'
p_c_acc(1) = [];

p_c_val = abs(p_chord_acc(:,3)' - p_root);
p_c_val(1) = [];

t_chord_acc = importdata('t_chord_accuracy.csv');
iter_t_c_acc = t_chord_acc(:,1)';
iter_t_c_acc(1) = [];
t_c_acc = t_chord_acc(:,2)';
t_c_acc(1) = [];

t_c_val = abs(t_chord_acc(:,3)' - t_root);
t_c_val(1) = [];


figure
semilogx(p_b_acc, iter_p_b_acc, 'Color','g','LineWidth',2);
hold on
semilogx(t_b_acc, iter_t_b_acc, 'Color','r','LineWidth',1);
semilogx(p_c_acc, iter_p_c_acc, 'Color','y','LineWidth',3);
semilogx(t_c_acc, iter_t_c_acc, 'Color','m','LineWidth',2);
grid on
title ('Зависимость кол-ва итераций от необходимой точности')
xlabel ('epsilon');
ylabel ('k');
legend('bisection polynom','bisection transcendental', 'chord polynom', 'chord transcendental', 'Location', 'NorthWest');

% Достигается ли точность?

figure
loglog(p_b_acc, p_b_val, 'Color','g','LineWidth',2);
hold on
grid on
xlabel ('epsilon');
ylabel ('|x^{k^*} - x^*|');
title ('Достигается ли точность?')
loglog(p_c_acc, p_c_val, 'Color','r','LineWidth',2);
loglog(t_b_acc, t_b_val, 'Color','y','LineWidth',3);
loglog(t_c_acc, t_c_val, 'Color','m','LineWidth',2);
syms x;
y=x;
fplot(y);
legend('bisection polynom', 'chord polynom', 'bisection transcendental', 'chord transcendental', 'y = x', 'Location', 'SouthEast');

% зависимость от x_0

i = [3,5,7,9,11,13];
p_b_x0 = importdata('p_bisection_root.csv');
p_b_x0_iter = p_b_x0(:,1);

t_b_x0 = importdata('t_bisection_root.csv');
t_b_x0_iter = t_b_x0(:,1);

p_c_x0 = importdata('p_chord_root.csv');
p_c_x0_iter = p_c_x0(:,1);

t_c_x0 = importdata('t_chord_root.csv');
t_c_x0_iter = t_c_x0(:,1);

figure
subplot(2,2,1)
plot(i-p_root, p_b_x0_iter, 'Color','g','LineWidth',2);
hold on
grid on
ylabel ('i');
xlabel ('|x^0 - x^*|');
title('polynom: bisection method')
subplot(2,2,2)
plot(i-t_root, t_b_x0_iter, 'Color','g','LineWidth',2);
grid on
ylabel ('i');
xlabel ('|x^0 - x^*|');
title('transcendental: bisection method')
subplot(2,2,3)
semilogy(i-p_root, p_c_x0_iter, 'Color','g','LineWidth',2);
grid on
ylabel ('i');
xlabel ('|x^0 - x^*|');
title('polynom: chord method')
subplot(2,2,4)
semilogy(i-t_root, t_c_x0_iter, 'Color','g','LineWidth',2);
grid on
ylabel ('i');
xlabel ('|x^0 - x^*|');
title('transcendental: chord method')

