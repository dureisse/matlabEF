function [LL,MM,CC,err,varargout] = GeometricQuadraticFit(xcoor,varargin)
% Fit a tangent quadratic surface to a cloud of nodes
%
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 15 / 12 / 2006
%
% Inputs
%   xcoor(nbno,3)	coordinates of the nodes
% Optional inputs
%   w(nbno)		weights (default: 1)
% Ouputs
%   LL(1,3)		
%   MM(3,3)		
%   CC			
%   err			minimum square error
%
% The fitted surface equation is:
% LL [x y z]' = [x y z] MM [x y z]' + CC

  nin1 = nargin - 1;
  switch nin1
    case 0,
      w = 1;
    case 1,
      w = varargin{1};
      error('Weight not yet done')
    otherwise,
      nin1
      error('Bad number of optional input arguments')
  end
  nou1 = nargout-4;
  if nou1 > 1
    error('Too any output arguments')
  end

  idim = size(xcoor,2);
  if (idim ~= 3)
    idim
    error('Only in 3D')
  end

nbno0 = size(xcoor,1);
if nbno0 <=6
  nbno0
  error('Not enough nodes')
end

% Average plane
%[N,I,G] = GeometricForm(xcoor,w);
[N,I,G] = GeometricForm(xcoor);

% Local problem
NN = inv(N');
xcoorl = xcoor * NN;

% Quadratic fit, with linear and constant term
M = zeros(3+2+1,3+2+1);
C = zeros(3+2+1,1);
for i = 1:nbno0
  Pi = [xcoorl(i,1)^2 ; xcoorl(i,2)^2 ; sqrt(2)*xcoorl(i,1)*xcoorl(i,2)];
  Qi = [xcoorl(i,1) ; xcoorl(i,2)];
  Vi = [Pi ; Qi ; 1];
  M = M + Vi*Vi';
  C = C + xcoorl(i,3)*Vi;
end
s = svd(M);
if abs(s(end)/s(1)) < 1.e-8
  error('Bad conditioning')
end
B = M \ C;
a = B(1,1); b = B(2,1); c = B(3,1)/sqrt(2);
d = B(4,1); e = B(5,1);
f = B(6,1);

% Final surface equation: LL [x y z]' = [x y z] MM [x y z]' + CC
L = ([0 0 1] - [d e 0])*NN';
nn = [1 0 0 ; 0 1 0]*NN';
MM = nn' * [a c ; c b] * nn;

LL = L;
CC = f;

% Erreur
num = 0.; den = 0.;
for i = 1:nbno0
  x = xcoor(i,1); y = xcoor(i,2); z = xcoor(i,3);
  err1 = (LL * [x y z]') - ([x y z] * MM * [x y z]') - CC;
  err2 = (LL * [x y z]') + ([x y z] * MM * [x y z]') + CC;
  num = num + err1^2;
  den = den + err2^2;
end
err = sqrt(num/den);

if nou1 == 1
  varargout{1} = N;
end
