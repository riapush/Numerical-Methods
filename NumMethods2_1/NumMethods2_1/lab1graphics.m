all_cheb = importdata('cheb.csv');
all_uni = importdata('uni.csv');
all_cheb_nodes = fopen('cheb_points.csv', "r");
all_uni_nodes = fopen('uni_points.csv', "r");
cheb_x1 = fscanf(all_cheb_nodes, "%f;")'
uni_x1 = fscanf(all_uni_nodes, "%f;")'

f = @(x)(x.^2 - sin(10.*x));
x = importdata('x.csv');
error_cheb = zeros(12);
error_uni = zeros(12);

figure
k = 1;
m = 1;
for i = 1:1:22
    cheb1 = all_cheb(i,:);
    uni1 = all_uni(i,:);
    error_cheb(i) = max(abs(cheb1 - f(x)));
    error_uni(i) = max(abs(uni1 - f(x)));
    
    subplot(3,2,1)
    hold off
    plot(x,f(x), 'b', 'LineWidth', 2);
    hold on
    grid on
    plot(x, cheb1, 'r--', 'LineWidth', 2);
    for j = 1:1:i+2
        plot(cheb_x1(k), f(cheb_x1(k)), 'g*');
        k = k + 1;
    end
    title(['Чебышевская сетка при ', num2str(i+2), ' узлах'])
    legend('function', 'Hermit', 'Location', 'SouthWest')
    hold off
    
    subplot(3,2,3)
    plot(x,f(x), 'r', 'LineWidth', 2);
    hold on
    grid on
    plot(x, uni1, 'b--.', 'LineWidth', 2);
    for j = 1:1:i+2
        plot(uni_x1(m), f(uni_x1(m)), 'g*');
        m = m + 1;
    end
    title(['Равномерная сетка при ', num2str(i+2), ' узлах'])
    legend('function', 'Hermit', 'Location', 'SouthWest');
    hold off
    
    subplot(3,2,5)
    semilogy(x, abs(cheb1 - f(x)), 'LineWidth', 2, 'Color', 'blue')
    hold on
    grid on
    semilogy(x, abs(uni1 - f(x)), 'LineWidth', 2, 'Color', 'red')
    legend('Чебышевская сетка', 'Равномерная сетка', 'Location', 'SouthWest');
    hold off
    
    subplot(3,2,[2,6])
    semilogy(i+2, error_cheb(i), 'or');
    hold on
    grid on
    semilogy(i+2, error_uni(i), '*g');
    title('Зависимость максимальной ошибки от количества узлов')
    xlabel('Количество узлов');
    ylabel('Модуль ошибки');
    legend('Чебышевская сетка', 'Равномерная сетка', 'Location', 'SouthWest');
    
    if i == 1
        pause (2)
    else
        pause(0.5)
    end
end
fclose('all');