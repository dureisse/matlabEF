clear all
close all
path(path,'../matlabUtils')
path(path,'../matlabEF')

X1 = 1; Y1 = 2; 
X2 = 3; Y2 = 4; 
X3 = 5; Y3 = 6; 

x1 = 10; y1 = 11; z1 = 100;
x2 = 12; y2 = 13; z2 = 101;
x3 = 14; y3 = 15; 
x4 = 16; y4 = 17; 


numer1 = [1:3];
list1 = [{'UX'} {'UY'}];
map1 = [1 2
        3 4
        5 6];
vec1 = [X1 Y1 X2 Y2 X3 Y3]';
numer2 = [2:-1:1];
list2 = [{'FY'} {'FZ'} {'FX'}];
vec2 = [y2 z2 x2 y1 z1 x1]';
map2 = [1 2 3
        4 5 6];

listps1 = [{'UX'} {'UY'}];
listps2 = [{'FX'} {'FY'}];
numer3 = [1 2 4];
X4 = 0; Y4 = 0; 
vec3 = PscalVect(vec1,numer1,map1,list1, ...
                 vec2,numer2,map2,list2, ...
                 listps1,listps2,numer3);
vec3a = [X1*x1+Y1*y1 X2*x2+Y2*y2 X4*x4+Y4*y4]';

err1 = norm(vec3-vec3a)/norm(vec3a);
if (err1 > 1.e-10)
  err1
  error('pb')
end


% Cas d'un champ scalaire
T1 = 0.1; T2 = 0.2; T3 = 0.3;
numer3 = [1:3];
list3 = [{'T'}];
vec3 = [T1 T2 T3]';
map3 = numer3';

listps1 = [{'UX'} {'UY'}];
listps3 = [{'T'}];
vec4 = PscalVect(vec3,numer3,map3,list3, ...
                 vec1,numer1,map1,list1, ...
                 listps3,listps1,numer3);
vec4a = [T1*X1 T1*Y1 T2*X2 T2*Y2 T3*X3 T3*Y3]';

err2 = norm(vec4-vec4a)/norm(vec4a);
if (err2 > 1.e-10)
  err2
  error('pb 2')
end

disp('TEST PASSE AVEC SUCCES')
quit
