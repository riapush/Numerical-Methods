p = @(x)(x.^6 + x.^5 - 13 .* x.^3 - 9 .* x + 2);
p_root = fzero(p,[1,3])
t = @(x)(5.^x-6.*x-7);
t_root = fzero(t, [1,2])

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
legend('bisection polynom', 'chord polynom', 'bisection transcendental', 'chord transcendental', 'Location', 'SouthEast');
syms x;
y=x;
fplot(y);

% зависимость от x_0

%p_b_x0 = importdata('p_bisection_root.csv');
%p_b_x0_iter = p_x0(:,1);
%p_b_x0_root = p_x0(:,3);