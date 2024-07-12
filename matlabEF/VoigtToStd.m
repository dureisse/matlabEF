function [tens1] = StdTovoigt(vec1)
% symmetric tensor from its Voigt form
% Actually Cowin, not Voigt!
%
% The tensor is symmetric!
%
% Input
%   vec1(6,1)
% Output
%   tens1(3,3)

r22 = sqrt(2.)/2.;
tens1 = [vec1(1,1)     r22*vec1(6,1)  r22*vec1(5,1)      
         r22*vec1(6,1) vec1(2,1)      r22*vec1(4,1)
         r22*vec1(5,1) r22*vec1(4,1)  vec1(3,1)];
