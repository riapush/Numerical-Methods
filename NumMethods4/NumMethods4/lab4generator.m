%  clc
%  clear on
%  file = fopen('matrix_accuracy.csv', 'w');
% 
% n = 10;
% [Q,~] = qr(rand(n));
% Achoosen = zeros(n,1);
% sep_val = i/20;
% lambda(10) = randi([1 20]);
% for j=9:-1:1
%     lambda(j) = lambda(j+1)/sep_val;
% end
% D = diag([lambda(1), lambda(2), lambda(3), lambda(4), lambda(5), lambda(6), lambda(7), lambda(8), lambda(9), lambda(10)]);
% A = Q'*D*Q;
% if i == 1
%     Achoosen = A;
% end
% for i=1:1:n
%     fprintf(file, '%.15f;', A(:,i));
%     fprintf(file,'\n');
% end
% fprintf(file,'\n');
% fclose(file);

n = 10;
count = 10;
Achoosen = zeros(n,1);
fileN = fopen('rang.csv', 'w');
fprintf(fileN, '%i', n);
fclose(fileN);
file1 = fopen('matrix.csv', 'w');
file2 = fopen('sep_values.csv', 'w');
fclose(file1);
sep_value = 1 - 10^-10;
file1 = fopen('matrix.csv', 'a');
new_sep_value = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
D = zeros(n,1);
for i=1:count
    fprintf(file2, '%.10f;', new_sep_value(i));
    for j=1:n-1
        if j == 1
            D(j,j) = 0.3;
            D(j+1, j+1) = D(j,j)*(new_sep_value(i));
        else
            D(j+1, j+1) = D(j,j)*(new_sep_value(i));
        end
    end
    [Q, R] = qr(rand(n)*20);
    A = Q*D*Q';
    for j=1:n
        AR = A(j,:);
        fprintf(file1, '%.17f;', AR);
        fprintf(file1, '\n');
    end
    fprintf(file1, '\n');
    %new_sep_value = sep_value/(i+1);
    if (i == 1)
        Achoosen = A;
    end
end
fclose(file1);
fclose(file2);


delta = [0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10];
delta_len = length(delta);
file1 = fopen('matrix_per.csv', 'w');
for i=1:1:n
    for j=1:1:n
        fprintf(file1, '%.17f;', Achoosen(i,j));
    end
    fprintf(file1, '\n');
end
fileA = fopen('perturbation.csv', 'w');
for i = 1:delta_len
    A = Achoosen + delta(i)*rand(n, 1);
    fprintf(fileA, '%.17f;', A);
    fprintf(fileA, '\n');
end
fclose(fileA);
fclose(file1);

