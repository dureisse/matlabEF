function [B] = LocalB_Pout(DPHI,NNIP,NDDP,PHI,idim,DDPHI);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2003
%
% Local operator derivative of kinematic torsor of the basis functions
% suitable for Euler-Bernoulli beam
%
% Entrees
%   DPHI(i,nni) = valeur de la derivee par rapport a la coordonne i dans
%                 l'espace de reference, de la fonction de forme nni
%   NNIP = liste des numeros  des fonctions de formes associees aux ddl
%          locaux de l'element fini
%   NDDP = liste des degres de liberte primaux associes
%   PHI(nni) = valeur de la fonction de forme nni
%   idim = Dimension of the physical space
%   DDPHI(i,j,nni) = valeur de la derivee seconde par rapport aux
%                    coordonnees i et j dans l'espace reel,
%                    de la fonction de forme nni
% Sorties
%   B(nbcopB,nbddl)	Operateur gradient pour le deplacement axial,
%                       derivee seconde pour la fleche
%
% Assumption: dofs are in right order (1 to idim)
%             'Voight' notation for kinematic torsor is used for B


  nbddl = length(NNIP); % Number of local ddl (within the element)
  idimr = size(DPHI,1); % Dimension of the reference space
  if (idimr ~= 1)
    idimr
    error('Bad idimr')
  end
  if (idim == 2)
    nbcop  = 3; % Number of primal components for the kinematic field
%                 U1,U2,R3
    nbcopB = 2; % Number of components for the 'strain' field
%                 EPS1,C3
    Q      = [-1.];        % Operator e1 vectoriel .
  elseif (idim == 3)
    nbcop  = 6; % Number of primal components for the kinematic field
%                 U1,U2,U3,R1,R2,R3
    nbcopB = 4; % Number of components on the derivated field
%                 EPS1,C1,C2,C3
    Q      = [1. 0. 0.
              0. 0. -1.
              0. 1. 0.]; % Operator e1 vectoriel .
  else
    idim
    error('Bad idim')
  end

  B = zeros(nbcopB,nbddl);

if ((idim == 2) && (nbddl == 6))
%   1    2    1  2  3  4  5  6
% B(EPS1 C3 , U1 U2 R3 U1 U2 R3)
  B(1,1) = DPHI(1,1);
  B(2,2) = DDPHI(1,1,2);
  B(2,3) = DDPHI(1,1,3);
  B(1,4) = DPHI(1,4);
  B(2,5) = DDPHI(1,1,5);
  B(2,6) = DDPHI(1,1,6);
  B(2,:) = -Q' * B(2,:);
elseif ((idim == 3) && (nbddl == 12))
%   1    2  3  4    1  2  3  4  5  6  7  8  9  10 11 12
% B(EPS1 C1 C2 C3 , U1 U2 U3 R1 R2 R3 U1 U2 U3 R1 R2 R3)
% traction
  B(1,1) = DPHI(1,1);
  B(1,7) = DPHI(1,7);
% torsion
  B(2,4)  = DPHI(1,4);
  B(2,10) = DPHI(1,10);
% flexion / 2
  B(3,3) = DDPHI(1,1,3);
  B(3,9) = DDPHI(1,1,9);
  B(3,5) = -DDPHI(1,1,5);
  B(3,11) = -DDPHI(1,1,11);
% flexion / 3
  B(4,2) = DDPHI(1,1,2);
  B(4,8) = DDPHI(1,1,8);
  B(4,6) = DDPHI(1,1,6);
  B(4,12) = DDPHI(1,1,12);

%  B(2:4,:) = -Q' * B(2:4,:);
else
  error('Pas prevu dans le cas general...')
end
