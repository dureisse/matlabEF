function [COOR,WEIGHT] = TRI_Integr(nbptg);
% Location and weights of integration points in a triangle
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 06 / 11 / 2007
%
% This routine has been factorized to get consistent integrations
% on each kind of triangle (TRI3, TRI6, TRI6d)
% Inputs
%   nbptg		Number of integration points
% Outputs
%   COOR(nbptg,2)	Reference coordinates of the points
%   WEIGHT(nbptg)	Weights

%     INTEGRATION 7 POINTS TRIANGLE - Cast3M
%      X101=.1012865073235D0;
%      X797=.7974269853531D0;
%      X470=.4701420641051D0;
%      X059=.0597158717898D0;
%      P125=.1259391805448D0;
%      P132=.1323941527885D0;
%    COOR   = [X101 X101
%              X797 X101
%              X101 X797
%              X470 X059
%              X470 X470
%              X059 X470
%              1./3. 1./3.];
%    WEIGHT = [P125*0.5 P125*0.5 P125*0.5 P132*0.5 P132*0.5 P132*0.5 .1125D0];


switch nbptg

  case 1,
%   centroid
    COOR   = [1./3. 1./3.];
    WEIGHT = [0.5];

  case 3,
%   cf Zienkiewicz
    COOR   = [1./6. 1./6.
              2./3. 1./6.
              1./6. 2./3.];
    WEIGHT = [1./6. 1./6. 1./6.];
%   other choice: Hammer
%    COOR   = [0.5 0.
%              0.  0.5
%              0.5 0.5];
%    WEIGHT = [1./6. 1./6. 1./6.];

  case 4,
%   problem: negative weight
%   Version poids central negatif Zienckiewick p. 201
%    COOR   = [1./3. 1./3.
%              1./5. 1./5.
%              3./5. 1./5.
%              1./5. 3./5.];
%    WEIGHT = [-27./96.  25./96.  25./96. 25./96.];
%    dans Hillion, Num. Integr. on a Triangle, IJNME 11:797-815, 1977
     disp('Integration 4 points Hillion en test...')
     COOR = [0.666390246 0.178558728
             0.178558728 0.666390246
             0.280019915 0.075031109
             0.075031109 0.280019915];
     WEIGHT = [0.159020691 0.159020691 0.090979309 0.090979309];

  case 6,
%   dans BATOZ tome 1 p200 : integration de poly d'ordre 4
    A = 0.445948490915965;
    B = 0.091576213509771;
    C = 0.111690794839005;
    D = 0.054975871827661;
    COOR = [A          A
            (1. -2.*A) A
            A          (1.-2.*A)
            B          B
            (1.-2.*B)  B
            B          (1.-2.*B)];
    WEIGHT = [C C C D D D];

  case 7,
%   cf Castem et Bathe p. 280, mais numerotation differente d'ici
%   dans BATOZ tome 1 p200 : integration de poly d'ordre 5
    COOR = [...
      1./3.                 1./3.
      (9.-2.*sqrt(15.))/21. (6.+sqrt(15.))/21.
      (6.+sqrt(15.))/21.    (9.-2.*sqrt(15.))/21.
      (6.+sqrt(15.))/21.    (6.+sqrt(15.))/21.
      (9.+2.*sqrt(15.))/21. (6.-sqrt(15.))/21.
      (6.-sqrt(15.))/21.    (9.+2.*sqrt(15.))/21.
      (6.-sqrt(15.))/21.    (6.-sqrt(15.))/21.];
    WEIGHT = [ 9./80. ...
      (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. ...
      (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400.];

  otherwise,
    nbptg
    error('Number of integration points not implemented')
end
