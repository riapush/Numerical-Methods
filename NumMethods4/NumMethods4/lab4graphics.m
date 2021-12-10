matrix = importdata('matrix_accuracy.csv');
solutions = importdata('accuracy.csv')
eig_val = sort(eig(matrix))'
eps = [1,2,3,4,5,6,7,8,9,10,11,12,13];
err = zeros(10,1);
max_err = zeros(13,1)
for i=1:1:13
    eps(i) = 10^-i;
end
for i=1:1:13
    err = sort(abs(solutions(i,:) - eig_val))
    max_err(i) = err(10)
end
figure
loglog(eps, max_err, 'LineWidth', 2)
hold on
grid on
syms x;
y=x;
fplot(y);
title('Достигается ли точность?');
xlabel('Заданная точность');
ylabel('Полученная точность');

% clc
% clear on
% sep_values = importdata('sep_values.csv');
% iter = importdata('iterations.csv');
% iter(:,17:end)=[];
% 
% figure
% plot(sep_values,iter);
% hold on
% grid on
% title({'График зависимости количества итераций';'от числа отделимости при заданной точности'});
% xlabel('Число отделимости');
% ylabel('Количество итераций');
% 
% matrix = importdata('matrix1.csv');
% sol = importdata('solution.csv')
% sol1 = zeros(5,1)';
% delta = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 0];
% delta_len = length(delta);
% error = zeros(delta_len, 1)';
% figure
% eig_num = sort(eig(matrix)')
% eig_num(3)
% medium_err = 0;
% for i =1:1:delta_len
%     for j=1:1:5
%         for k = 1:1:5
%             eig_num(k)
%             sol(j,k)
%             medium_err = medium_err + abs(eig_num(k) - sol(j,k))
%         end
%     end
% end
