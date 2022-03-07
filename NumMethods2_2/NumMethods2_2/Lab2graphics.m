fixed_degree = importdata("fixed_degree.csv");

f = @(x)(x.^2 - sin(10.*x));
x = (-2:0.00001:2);
error_f_d = zeros(40);
error_f_n = zeros(48);
num_of_points = zeros(40);
num_of_points(1) = 10;

for i = 1:1:40
    coef = fixed_degree(i,:);
    polynom = @(x)(coef(1)*x.^0 + coef(2)*x.^1 + coef(3)*x.^2 + coef(4)*x.^3 + coef(5)*x.^4+coef(6)*x.^5);
    error_f_d(i) = max(abs(polynom(x)-f(x)));
    if i ~= 1
        num_of_points(i) = num_of_points(i-1)+5;
    end
end

figure
semilogy(num_of_points,error_f_d);
hold on
grid on;
title({'Зависимость макс. отклонения от числа узлов' 'при фиксированной степени полинома (5)'});
xlabel('Количество узлов');
ylabel('Модуль максимального отклонения');

fixed_n = importdata("fixed_n.csv");
degree = zeros(48);
for i = 1:1:34
    coef = fixed_n(i,:);
    p = polyval(coef,x);
    error_f_n(i) = max(abs(p-f(x)));
    degree(i) = i+1;
end

figure
semilogy(degree,error_f_n);
hold on
grid on;
title({'Зависимость макс. отклонения от степени полинома' 'при фиксированном количестве узлов (50)'});
xlabel('Cтепень полинома');
ylabel('Модуль максимального отклонения');

fclose('all');