fixed_degree = importdata("fixed_degree.csv");
cond_val = fopen("for_cond_val.csv");

f = @(x)(x.^2 - sin(10.*x));
x = (-2:0.00001:2);
error_f_d = zeros(90);
error_f_n = zeros(30);
num_of_points = zeros(90);
num_of_points(1) = 10;

for i = 1:1:90
    coef = fixed_degree(i,:);
    polynom = @(x)(coef(1)*x.^0 + coef(2)*x.^1 + coef(3)*x.^2 + coef(4)*x.^3 + coef(5)*x.^4+coef(6)*x.^5);
    error_f_d(i) = max(abs(polynom(x)-f(x)));
    if i ~= 1
        num_of_points(i) = num_of_points(i-1)+1;
    end
end

figure
loglog(num_of_points,error_f_d);
hold on
grid on;
title({'Зависимость макс. отклонения от числа узлов' 'при фиксированной степени полинома (5)'});
xlabel('Количество узлов');
ylabel('Модуль максимального отклонения');

x = linspace(-2,2,50);
fixed_n = fopen('fixed_n.csv', 'r');
fixed_n_array = fscanf(fixed_n, "%f;")';
k = 1;
degree = zeros(30);
for i = 1:1:30
    coef = zeros(i+1);
    for j = 1:1:i+1
        coef(j) = fixed_n_array(k);
        k = k + 1;
    end
    coef = coef(end:-1:1);
    p = polyval(coef,x);
    error_f_n(i) = max(abs(p-f(x)));
    degree(i) = i;
end

figure
loglog(degree,error_f_n);
hold on
grid on;
title({'Зависимость макс. отклонения от степени полинома' 'при фиксированном количестве узлов (50)'});
xlabel('Cтепень полинома');
ylabel('Модуль максимального отклонения');

conds = zeros(30);
all_vals = fscanf(cond_val, "%lf; ");
k = 1;
for i = 2:1:31
    matrix = zeros(i,i);
    for m = 1:1:i
        for n = 1:1:i
            matrix(m,n) = all_vals(k);
            k = k + 1;
        end
    end
    conds(i-1) = cond(matrix);
end

figure
semilogy(degree, conds);
hold on;
grid on;
title({'Зависимость числа обусловленности матрицы от степени полинома' 'при фиксированном количестве узлов (50)'});
xlabel('Cтепень полинома');
ylabel('Число обусловленности');

x = linspace(-2,2,1000);
p2_import = fopen("comparison_2.csv", "r");
p_1 = fscanf(p2_import, "%lf; ");
p_1 = p_1(end:-1:1);
p_2 = polyval(p_1,x);
p1_import = importdata("D:/Git/GitHub/Numerical-Methods/NumMethods2_1/NumMethods2_1/comparison_1.csv");
f = @(x)(x.^2 - sin(10.*x));

nodes = linspace(-2,2,15);
figure
plot(x, p_2, 'LineWidth', 2);
hold on
grid on
plot(x, p1_import, 'LineWidth', 2);
plot(x,f(x), '--', 'LineWidth', 2);
title('Сравнение с интерполяцией');
for i = 1:1:15
    plot(nodes(i), f(nodes(i)), 'bo', 'LineWidth', 2);
end
legend('Интерполяция', 'МНК', 'f(x)', 'Location', 'South');

fclose('all');