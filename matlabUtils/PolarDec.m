function [V,R] = PolarDec(F);
% Polar decomposition of F = V * R
% F should be square and regular

V = sqrtm(F * F');
R = V \ F;
