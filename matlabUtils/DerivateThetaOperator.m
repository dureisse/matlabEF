function [D1,V1] = DerivateThetaOperator(lt1,theta1)
% Derivation operator for continuous time functions and theta method
%
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 03 / 04 / 2006
%
% Inputs
%   lt1(nt)	Time instant values
%   theta1	Parameter of the integration/derivation scheme
% Outputs
%   D1(nt,nt)	Derivate operator
%   V1(1,nt)	For taking into account initial value
%
%   If a function is given by a vector line a(1,nt)
%   its derivate will be da = a.D1 + v0 * V1
%   v0 is the derivative value at lt1(1)

nt = length(lt1);
I1 = eye(nt);
%%h = wshift('1D',lt1,1) - lt1; h = h(1:end-1);
h = lt1(2:end) - lt1(1:end-1);
uth = 1./(theta1*h);
A = spdiags([[-uth 0]'  [0 uth]'],[-1 0],nt,nt);
B = (1. - 1./theta1) * spdiags(ones(nt,1),-1,nt,nt);

% Solves X' = A*X + B*X' + B*[v0 0. ... 0.]'
% i.e.   (1-B) X' = A*X + B(:,1)*v0
% i.e. gives X' = [(1-B)^-1 * A] * X + [(1-B)^-1 * B(:,1)]*v0

M = I1 - B;
D1 = (M \ A)';
V1 = (M \ B(:,1))';


