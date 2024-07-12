function [Ap,R,A] = GeneralizedInverse(A,PREC);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 08 / 05 / 00
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 01 / 2003
%
% Inversion d'une matrice carree symetrique A par la methode de Crout,
% sans pivotage, mais avec recuperation des modes rigides.
% A contient en sortie sa factorisation LDLT generalisee 
% (si mode rigide detecte en position I, D(I) est mis a 1)
% L et D ecrasent le triangle inferieur, le triangle superieur
% reste indemne.
% (d'apres la routine invgn2.eso de style Fortran)
% 
% Input
%   A(N,N)      Matrice a inverser
%   PREC        Precision pour la detection des modes rigides
% Output
%   Ap(N,N)     Inverse generalisee
%   R(N,NRIG)   Modes rigides
%   A(N,N)      Contient L,D de la factorisation generalisee LDLT de A

N = size(A,1);
B = zeros(N,2*N);
V = zeros(N,1);

% -   Factorisation LDLT de A
      NRIG  = 0;
listrig = [];
      for J = 1:N
        LRIG = 0;
        for I = 1:J
          V(I,1) = A(J,I) * A(I,I);
        end
% -	On calcule le pivot
	ILON = J - 1;
	V(J,1) = A(J,J) - A(J,1:ILON) * V(1:ILON,1);
	if (abs(V(J,1)) <= PREC)
% -	   Detection d'un pivot nul...
	   LRIG = 1;
	   V(J,1) = 1.;
	   NRIG = NRIG + 1;
	   B(J,N+NRIG) = 1.;
listrig(NRIG) = J;
        end
% -	On modifie la colonne (on y stocke L)
	A(J,J) = V(J,1);
	for K = J + 1:N
	  ILON = J - 1;
          XX = A(K,1:ILON) * V(1:ILON,1);
	  A(K,J) = (A(K,J) - XX) / A(J,J);
	end
	if (LRIG)
% -	  ... pouvait-on pivoter ?
	  for I = J+1:N
	    if (abs(A(I,J)) >= PREC)
              J
              I
              error('il fallait pivoter')
	    end
	  end
	end
      end
% -
% -   Pour inverser A, on resoud les N systemes lineaires LDLT B = 1
      for J = 1:N
	B(J,J) = 1.;
      end
% -   L
      for K = 1:N+NRIG
        for J = 2:N
	  ILON = J - 1;
          XX = A(J,1:ILON) * B(1:ILON,K);
	  B(J,K) = B(J,K) - XX;
        end
      end
% -   D
      for K = 1:N
	for J = 1:N
	  B(J,K) = B(J,K) / A(J,J);
	end
      end
% -   LT
      for K = 1:N+NRIG
	for J = N-1:-1:1
	  ILON = N - J;
          XX = A(J+1:J+ILON,J)' * B(J+1:J+ILON,K);
	  B(J,K) = B(J,K) - XX;
	end
      end


Ap = B(:,1:N);
%for i = 1:NRIG
%  J = listrig(i);
%  Ap(i,i) = 0.;
%end
R  = B(:,N+1:N+NRIG);


for i = 1:NRIG
  J = listrig(i);
  A(J,J) = 0.;
end
clear B;
