clear all
close all
path(path,'../matlabUtils')

for i = 1:10

N = 10;
A = rand(N);
A = A + A';

NRIG = 4;
P = rand(N,NRIG);
L = chol(P'*P);
P = (L' \ P')';

P = eye(size(A)) - P*P';
A = P' * A * P;

PREC = 1.e-6;
[Ap,R,B] = GeneralizedInverse(A,PREC);
if (size(R,2) ~= NRIG)
  save test_GeneralizedInverse A PREC;
  error('A tester en detail... voir test_GeneralizedInverse.mat')
end
err1 = norm(A*R)/norm(A);
err2 = norm(A * Ap * A - A)/norm(A);
if (err1 > 1.e-6) | (err2 > 1.e-6)
  i
  err1
  err2
  error('PB')
end

end

disp('TEST PASSE AVEC SUCCES')
quit
