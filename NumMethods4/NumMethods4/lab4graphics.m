% matrix = importdata('matrix_accuracy.csv');
% solutions = importdata('accuracy.csv')
% eig_val = sort(eig(matrix))'
% eps = [1,2,3,4,5,6,7,8,9,10,11,12,13, 14];
% err = zeros(10,1);
% max_err = zeros(14,1)
% for i=1:1:14
%     eps(i) = 10^-i;
% end
% for i=1:1:14
%     err = sort(abs(solutions(i,:) - eig_val))
%     max_err(i) = err(10)
% end
% figure
% loglog(eps, max_err, 'LineWidth', 2)
% hold on
% grid on
% syms x;
% y=x;
% fplot(y);
% title('Достигается ли точность?');
% xlabel('Заданная точность');
% ylabel('Полученная точность');

n = 10;
%sep_values = importdata('sep_values.csv');
sep_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];
iter = importdata('iterations.csv');
iter_hes = importdata('iterations_hes.csv');
iter_shift = importdata('iterations_shift.csv');
iter_shift(:,11:end)=[];
iter_hes(:,11:end)=[];

figure
semilogy(sep_values, iter, 'LineWidth', 2);
hold on
semilogy(sep_values, iter_hes, 'LineWidth', 2);
semilogy(sep_values, iter_shift, 'LineWidth', 2);
grid on
title({'График зависимости количества итераций';'от числа отделимости при заданной точности'});
xlabel('Число отделимости');
ylabel('Количество итераций');
legend('QR алгоритм', 'QR с приведением к форме Хессенберга', 'QR со сдвигом', 'Location', 'NorthWest')

matrix = importdata('matrix_per.csv');
eigen_vec = importdata('eigen_vectors.csv')
[~, ~, eigen_vec_matlab] = eig(matrix)
sol = importdata('solution_delta.csv');
delta = [0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10];
delta_len = length(delta);
error = zeros(delta_len, 1)';
accurate_sol = sort(eig(matrix)');
medium_err = 0;
for i =1:1:delta_len
        medium_err = norm(accurate_sol - sort(sol(i,:)));
        error(i) = medium_err/norm(abs(accurate_sol));
end

figure
loglog(delta, error, 'LineWidth', 2)
hold on
grid on
title('График влияния внесенного возмущения на решение')
xlabel('Возмущение');
ylabel('Относительная погрешность нахождения соб.значений');

for j = 1:1:delta_len
    matrix1 = zeros(n,n);
    for i=1:1:n
       matrix1(i,:) = eigen_vec((i-1)*10+i,:);
    end
    for i = 1:1:n
        matrix1(i,:) = matrix1(i,:)/norm(matrix1(i,:));
    end
    vec_err = abs(matrix1-V);
end