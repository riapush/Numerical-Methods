graph1 = importdata('iter.csv');
eps = graph1(:,1);
iter = graph1(:,2);

figure
semilogx(eps,iter, 'LineWidth', 2);
hold on
grid on
title({'График зависимости заданной точности';'от количества итераций'});
xlabel('Заданная точность');
ylabel('Количество итераций');

Rang = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2\rang3.csv');
hh = importdata('D:\Git\GitHub\Numerical-Methods\NumMethods2\householders_method\householders_method\time.csv');
gs = importdata('time.csv')
figure
subplot(2,1,1)
plot(Rang, hh);
hold on
title({'График зависимости времени выполнения метода'; 'от ранга матрицы (метод Хаусхолдера)'});
xlabel('Ранг матрицы');
ylabel('Время выполнения');
grid on
subplot(2,1,2)
plot(Rang, gs);
grid on
title({'График зависимости времени выполнения метода'; 'от ранга матрицы (метод Зейделя)'});
xlabel('Ранг матрицы');
ylabel('Время выполнения');