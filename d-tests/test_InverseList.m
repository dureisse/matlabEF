clear all
close all
path(path,'../matlabUtils');

numer1 = [1:10];
numer_inv1 = InverseList(numer1,max(numer1));
err1 = norm(numer1 - numer_inv1)/norm(numer1);
if (err1 > 1.e-10)
  err1
  error('problem 1')
end

numer1 = [1:5 10:-1:7 30];
numer_inv1 = InverseList(numer1,max(numer1));
numer2 = InverseList(numer_inv1,max(numer_inv1));
err2 = norm(numer1 - numer2)/norm(numer1);
if (err2 > 1.e-10)
  err2
  error('problem 2')
end

disp('TEST PASSE AVEC SUCCES')
quit
