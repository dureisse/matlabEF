clear all
close all
path(path,'../matlabUtils');

x1 = [1 2 3 4 5];
y1 = [1 2 3 4 5];
x2 = x1;
y2 = [5 4 3 2 1];
[x,y] = EvolIntersect(x1,y1,x2,y2);
if ~all([x,y] == [3,3])
  error('stop 1')
end

x1 = [1   2   3 4 5 6   7 8]; x2 = x1;
y1 = [1   1   3 4 3 2   1 0];
y2 = [0.9 1.5 2 2 2 1.9 2 2];
[x,y] = EvolIntersect(x1,y1,x2,y2);
%plot(x1,y1,x2,y2,x,y,'o')
err1 = max(abs(x-[1.166667 2.333333 6.090909]), ...
           abs(y-[1.0000 1.666667 1.9090909]));
if (err1 > 1.e-6)
  error('stop2')
end

x1 = [1.:1.:7.];
x2 = x1;
y1 = [300. 17.5  1.78 0.192 0.0187 0.001576 0.000113];
y2 = [8.62 1.55 0.313 0.0535 0.00816 0.00107 0.000122];
[x,y] = EvolIntersect(x1,y1,x2,y2);
%plot(x1,y1,x2,y2,x,y,'o')
%axis([6. 7 0 0.0002]) 
err1 = max(abs([x y] - [6.982524272 1.3857e-04]));
if (err1 > 1.e-8)
  error('stop3a')
end
[x,y] = EvolIntersect(x1,y1,x2,y2,'semilogy');
%semilogy(x1,y1,x2,y2,x,y,'o')
%axis([6. 7 0 0.0002]) 
err1 = max(abs([x y] - [6.834793965 1.7464e-04]));
if (err1 > 1.e-8)
  error('stop3b')
end
[x,y] = EvolIntersect(x1,y1,x2,y2,'loglog');
%loglog(x1,y1,x2,y2,x,y,'o')
%axis([6. 7 0 0.0002]) 
err1 = max(abs([x y] - [6.823984416 1.7464e-04]));
if (err1 > 1.e-8)
  error('stop3c')
end

disp('TEST PASSE AVEC SUCCES')
quit
