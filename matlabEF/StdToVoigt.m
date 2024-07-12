function [vec1] = StdTovoigt(tens1)
% Voigt form of a symmetric tensor
% Actually Cowin, not Voigt!

vec1 = StdToCowin(tens1);
