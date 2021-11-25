n = 10;
count = 1;
cond_val = 10;
fileB = fopen('rp.csv', 'w');
fileA = fopen('matrix.csv', 'w');
fileR = fopen('rang.csv', 'w');
fprintf(fileR, '%i', 10);
sum = 0;
X = rand(n,1);
for i = 1:count
    D = diag(linspace(1, cond_val, n));
    [Q, R] = qr(rand(n)*20);
    A = Q*D*Q';
    for j = 1:n
        for k = 1:n
            if j ~= k
                sum = sum + A(j, k);
            end
        end
        for k = 1:n
            if j ~= k
                A(j, k) = A(j, k)/(sum);
            end
        end
        sum = 0;
    end
    B = A * X;
    fprintf(fileB, '%.17f;', B);
    for j=1:n
        AR = A(j,:);
        fprintf(fileA, '%.17f;', AR);
        fprintf(fileA, '\n');
    end
    fprintf(fileA, '\n');
    fprintf(fileB, '\n');
    fprintf(fileB, '\n');
end
fclose(fileA);
fclose(fileB);
fclose(fileR);
