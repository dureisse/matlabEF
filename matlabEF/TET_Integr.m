function [COOR,WEIGHT] = TET_Integr(nbptg);
% Location and weights of integration points in a tetrahedron
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 01 / 12 / 2007
%
% This routine has been factorized to get consistent integrations
% on each kind of triangle (TET4, TE10, TE10d)
% Inputs
%   nbptg		Number of integration points
% Outputs
%   COOR(nbptg,2)	Reference coordinates of the points
%   WEIGHT(nbptg)	Weights


switch nbptg

  case 1,
%   centroid
    COOR   = [0.25 0.25 0.25];
    WEIGHT = [1./6.];

  case 4,
%%       Cools
%        XX = 0.36180339887498948;
%        COOR   = [XX         XX         XX
%                  (1.-3.*XX) XX         XX
%                  XX         (1.-3.*XX) XX
%                  XX         XX         (1.-3.*XX)];
%        WEIGHT = 1./24.*ones(1,4);
%       Felippa AFEM et Cast3M et Zienkievicz p. 202
        BETA  = 0.13819660;
        COOR   = [BETA         BETA         BETA
                  (1.-3.*BETA) BETA         BETA
                  BETA         (1.-3.*BETA) BETA
                  BETA         BETA         (1.-3.*BETA)];
        WEIGHT = 1./24.*ones(1,4);

%  case 5,
%%       5 points d'integration
%%       1 poids negatif !
%        nbptg  = 5; % Nombre de points d'integration
%        UNQUA = 0.25; UNDEMI = 0.5; UNSIX = 1./6.;
%        COOR = [UNQUA UNQUA UNQUA
%                UNSIX UNSIX UNSIX
%                UNSIX UNSIX UNDEMI
%                UNSIX UNDEMI UNSIX
%                UNDEMI UNSIX UNSIX];
%        WEIGHT = [-.80D0*UNSIX .45D0*UNSIX .45D0*UNSIX .45D0*UNSIX .45D0*UNSIX];

  otherwise,
    nbptg
    error('Number of integration points not implemented')
end
