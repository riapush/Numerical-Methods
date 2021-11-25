% найдем точные корни 
n = 10;
num = 2;
rang = 10;
A = importdata('matrix.csv');
b = importdata('rp.csv');
root = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2\householders_method\householders_method\x.csv');
cond_val = [10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10];
tmpA = zeros(rang);
tmpb = zeros(rang, 1);
tmproot = zeros(rang, 1);
err = zeros(n, 1);
err1 = zeros(n, 1);

for i = 1:1:n
    for j = 1:1:rang
            tmpA(j, :) = A(i, (j-1)*rang + 1:1:j*rang);
    end
    tmpb = b(i, :)';
    if (i == 1)
        Achoosen = tmpA;
        b_wo_pertrebution = tmpb;
    end
    tmproot = root(i, :)';
    matlab_root = tmpA \ tmpb;
    err1(i) = norm(tmpA * tmproot - tmpb);
    err(i) = norm(tmproot - matlab_root);
end
figure
loglog(cond_val, err, 'LineWidth', 2);
grid on
hold all
loglog(cond_val, err1, 'LineWidth', 2);
title({'График зависимости норм невязки и фактической ошибки';'от числа обусловленности матрицы'});
xlabel('Cond(A)');
ylabel('$||x - x^*||$ or $||Ax - b||$', 'interpreter', 'LaTex');
legend('Фактическая ошибка', 'Невязка', 'Location', 'southeast');


y1 = @(x) 10 * x;
delta = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 0];
delta_len = length(delta);
root = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2\householders_method\householders_method\x1.csv');
%Achoosen = importdata('matrix1.csv');
Bchoosen = importdata('rp1.csv');
matlab_root = Achoosen \ b_wo_pertrebution; %root(1) is root with delta = 0
err2 = zeros(n, 1);
errb = zeros(n, 1);
for i = 1:delta_len
err2(i) = norm(root(i,:)' - matlab_root)/norm(matlab_root);
errb(i) = norm(delta(i) * Bchoosen(i,:))/norm(Bchoosen(i,:));

end

for i = 1:delta_len-1
if (err2(i) <= errb(i)*10)
    display('TRUE')
else
    display('FALSE')
end
end

figure('name', 'зависимость отн. погрешности от возмущения');
loglog(errb, err2, 'Color', 'green', 'LineWidth', 2);
hold on
grid on
loglog(errb, y1(errb), 'Color', 'red', 'LineWidth', 1);
title({'График зависимости относительной погрешности'; 'от возмущения'});
xlabel('Возмущение');
ylabel('Отн. погрешность норм');
legend('Хорошо обусловленная матрица', 'y = cond\_well*x', 'Location', 'southeast')


Rang = importdata('rang3.csv');
Time = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2\householders_method\householders_method\time.csv');
figure
loglog(Rang, Time, 'LineWidth', 2);
hold on
title({'График зависимости времени выполнения метода'; 'от ранга матрицы'});
xlabel('Ранг матрицы');
ylabel('Время выполнения');
p1 =   1.088e-06;
p2 =  -0.0002761;
p3 =     0.02851;
p4 =     -0.6696;
f = @(x) (p1*x.^3 + p2*x.^2 + p3.*x + p4);
loglog(Rang,f(Rang), 'LineWidth', 2)
legend('time', 'fit', 'Location', 'NorthWest')