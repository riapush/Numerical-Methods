f = @(x)(x.^2 - sin(10.*x));
x = (-2:0.01:2);
p = @(x)(x.^2-0.256.*x);

figure
plot(x,f(x));
hold on;
grid on;
plot(x,p(x));
% coef = [0.000000 -0.225104 1.000000 0.286324 0.000000 -0.064950];
% coef = coef(end:-1:1);
% 
% p = polyval(coef, x);
% 
% figure
% plot(x,f(x));
% hold on
% grid on;
% plot(x,p);
% for j = 1:1:100
%     plot(uni_x1(j), f(uni_x1(j)), 'g*');
% end