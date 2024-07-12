function [vec1] = StdToCowin(tens1)
% Cowin form of a symmetric tensor
%
% The tensor is symmetrized before being transformed
%
% Input
%   tens1(3,3)
% Output
%   vec1(6,1)

r2 = sqrt(2.);
vec1 = zeros(6,1);
vec1(1,1) = tens1(1,1);
vec1(2,1) = tens1(2,2);
vec1(3,1) = tens1(3,3);
vec1(4,1) = r2*0.5*(tens1(2,3) + tens1(3,2));
vec1(5,1) = r2*0.5*(tens1(3,1) + tens1(1,3));
vec1(6,1) = r2*0.5*(tens1(1,2) + tens1(2,1));

