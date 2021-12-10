 clc
 clear on
 file = fopen('matrix_accuracy.csv', 'w');

n = 10;
[Q,~] = qr(rand(n));
Achoosen = zeros(n,1);
sep_val = i/20;
lambda(10) = randi([1 20]);
for j=9:-1:1
    lambda(j) = lambda(j+1)/sep_val;
end
D = diag([lambda(1), lambda(2), lambda(3), lambda(4), lambda(5), lambda(6), lambda(7), lambda(8), lambda(9), lambda(10)]);
display(D);
A = Q'*D*Q;
if i == 1
    Achoosen = A;
end
for i=1:1:n
    fprintf(file, '%.15f;', A(:,i));
    fprintf(file,'\n');
end
fprintf(file,'\n');
fclose(file);

%  clc
%  clear on
%  file = fopen('matrix.csv', 'w');
%  file_sep = fopen('sep_values.csv','w');
% sep_vals = zeros(16,1);
% 
% n = 10;
% [Q,~] = qr(rand(n));
% Achoosen = zeros(n,1);
% for i = 1:16
%     sep_val = i/20;
%     lambda(10) = randi([1 20]);
%     for j=9:-1:1
%         lambda(j) = lambda(j+1)/sep_val;
%     end
%     D = diag([lambda(1), lambda(2), lambda(3), lambda(4), lambda(5), lambda(6), lambda(7), lambda(8), lambda(9), lambda(10)]);
%     display(D);
%     A = Q'*D*Q;
%     if i == 1
%         Achoosen = A;
%     end
%     sep_vals(i) = 1/sep_val;
%     fprintf(file, '%.15f;', A);
% fprintf(file,'\n');
% end
% fprintf(file_sep, '%.3f;', sep_vals);
% fclose(file_sep);
% fclose(file);
% 
% delta = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 0];
% delta_len = length(delta);
% file1 = fopen('matrix1.csv', 'w');
% for i=1:1:5
%     for j=1:1:5
%         fprintf(file1, '%.20f;', Achoosen(i,j));
%     end
%     fprintf(file1, '\n');
% end
% fileA = fopen('perturbation.csv', 'w');
% for i = 1:delta_len
% A = Achoosen + delta(i)*rand(n, 1);
% fprintf(fileA, '%.20f;', A);
% fprintf(fileA, '\n');
% end
% fclose(fileA);
% fclose(file1);


