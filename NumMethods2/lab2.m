n = 10;
File1 = fopen('rp.csv', 'w');
File2 = fopen('matrix.csv', 'w');
File3 = fopen('rang.csv', 'w');
for cond_val = [10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10]
    [Q, R] = qr(rand(n));                       
    D = diag(linspace(1, cond_val, n));
    A = Q * D * Q';
    root = rand(n,1);
    b = A * root;
    if cond_val == 10
        Achoosen = A;
        Bchoosen = b;
    end
    for i = 1:1:n
        fprintf(File1, '%0.20d ', b(i, 1));
    end
    fprintf(File1, '\n');    
    for i = 1:1:n
        for j = 1:1:n
            fprintf(File2, '%0.20d ', A(i, j));
        end
    end
    fprintf(File2, '\n'); 
end
fprintf(File3, '%d', n);
fclose(File1);
fclose(File2);
fclose(File3);



File1 = fopen('matrix1.csv', 'w');
File2 = fopen('rp1.csv', 'w');
File3 = fopen('delta.csv', 'w');
cond_well = 10;
delta = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 0];
delta_len = length(delta);

A = Achoosen; % + delta(i)* Achoosen;
for j=1:n
AR = A(j,:);
fprintf(File1, '%.20f;', AR);
fprintf(File1, '\n');
end
for i = 1:delta_len
B = Bchoosen + delta(i)*rand(n, 1); % Bchoosen;
fprintf(File2, '%.20f;', B);
fprintf(File2, '\n');
fprintf(File2, '\n');
end
fclose(File1);
fclose(File2);
fclose(File3);


File1 = fopen('matrix2.csv', 'w');
File2 = fopen('rp2.csv', 'w');
File3 = fopen('rang3.csv', 'w');
File4 = fopen('matrixnum3.csv', 'w');
cond = 10;
k = 0;
for n = 10:10:400
    fprintf(File3, '%d ', n);
    [Q, R] = qr(rand(n));                       
    D = diag(linspace(1, cond, n));
    A = Q * D * Q';
    root = rand(n, 1);
    b = A * root;
    for i = 1:1:n
        for j = 1:1:n
            fprintf(File1, '%0.20d ', A(i, j));
        end
    end
    fprintf(File1, '\n');
    for i = 1:1:n
        fprintf(File2, '%0.20d ', b(i, 1));
    end
    fprintf(File2, '\n');
    k = k + 1;
end
fprintf(File4, '%d', k);
fclose(File1);
fclose(File2);
fclose(File3);
fclose(File4);