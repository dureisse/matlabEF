function [B] = LocalB_Dkir(DPHI,NNIP,NDDP,PHI,idim,DDPHI);
% Local derivative of shape functions for discrete Kirchhoff plates
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 07 / 2004
%
% Entrees
%   DPHI(i,nni) = valeur de la derivee par rapport a la coordonne i dans
%                 l'espace reel, de la fonction de forme nni
%   NNIP = liste des numeros  des fonctions de formes associees aux ddl
%          locaux de l'element fini
%   NDDP = liste des degres de liberte primaux associes
%   PHI(nni) = valeur de la fonction de forme nni
%   idim = Dimension of the physical space
%   DDPHI(i,j,nni) = valeur de la derivee seconde par rapport aux
%                    coordonnees i et j dans l'espace reel,
%                    de la fonction de forme nni
% Sorties
%   B(nbcopB,nbddl)	Operateur de deformation
%
% Assumption: dofs are in right order
%             'Voight' notation is used
% Beware that (unlike Batoz) the dofs are true rotations:
% RX is rotation around X axis.


r2 = sqrt(2.)/2.; % Due to Voigt notation for membrane

nbddl = length(NNIP); % Number of local ddl (within the element)
idimr = size(DPHI,1); % Dimension of the reference space
if (idimr ~= 2)
  idimr
  error('Bad idimr')
end

if (idim == 3)
% Number of primal dofs for the kinematic field
  nbcop  = 5; % membrane: U1,U2 , bending: U3,R1,R2
% Number of components on the derivated field
  nbcopB = 6; % membrane: EP11,EP22,EP12 , bending: N11,N22,N12
else
  idim
  error('Bad idim')
end

if ((idim == 3) && (nbddl == 28))
% DKQ
% 28 ddl : aux 4 coins (U1,U2,U3,R1,R2),
%          aux noeuds virtuels milieux des 4 aretes (R1,R2).
% Les ddl des noeuds virtuels seront ensuite a eliminer
% par la condition discrete de Kirchhoff.

% membrane : 8 ddl : aux 4 coins (U1,U2)
  Bm = zeros(3,28);
  nbno = 4;
  for ino = 1:nbno
    ddl1 = (ino-1)*5+1; % U1 du noeud ino
    Bm(1,ddl1)   = DPHI(1,ddl1);
    Bm(2,ddl1+1) = DPHI(2,ddl1+1);
    Bm(3,ddl1)   = r2 * DPHI(2,ddl1);
    Bm(3,ddl1+1) = r2 * DPHI(1,ddl1+1);
  end

% flexion : 20 ddl : aux 4 coins (U3,R1,R2) et aux 4 aretes (R1,R2)
% mais sans la condition discrete de Kirchhoff, on se sert uniquement
% de 16 ddl : aux 4 coins (R1,R2) et aux 4 aretes (R1,R2).
  Bf = zeros(3,28);
  for ino = 1:nbno
    ddl1 = (ino-1)*5+4; % R1 du noeud ino
    ddl2 = ddl1 + 1;    % R2 du noeud ino
    Bf(1,ddl2) = DPHI(1,ddl2);
    Bf(2,ddl1) = -DPHI(2,ddl1);
    Bf(3,ddl1) = -r2 * DPHI(1,ddl1);
    Bf(3,ddl2) = r2 * DPHI(2,ddl2);
  end
  for ino = 1:nbno
    ddl1 = (ino-1)*2+5*nbno+1; % R1 du noeud virtuel arete ino
    ddl2 = ddl1 + 1;           % R2 du noeud virtuel arete ino
    Bf(1,ddl2) = DPHI(1,ddl2);
    Bf(2,ddl1) = -DPHI(2,ddl1);
    Bf(3,ddl1) = -r2 * DPHI(1,ddl1);
    Bf(3,ddl2) = r2 * DPHI(2,ddl2);
  end

  B = [Bm ; Bf];
  clear Bm Bf;

elseif ((idim == 3) && (nbddl == 21))
% DKT
% 21 ddl : aux 3 coins (U1,U2,U3,R1,R2),
%          aux noeuds virtuels milieux des 3 aretes (R1,R2).
% Les ddl des noeuds virtuels seront ensuite a eliminer
% par la condition discrete de Kirchhoff.

% membrane : 6 ddl : aux 3 coins (U1,U2)
  Bm = zeros(3,21);
  nbno = 3;
  for ino = 1:nbno
    ddl1 = (ino-1)*5+1; % U1 du noeud ino
    Bm(1,ddl1)   = DPHI(1,ddl1);
    Bm(2,ddl1+1) = DPHI(2,ddl1+1);
    Bm(3,ddl1)   = r2 * DPHI(2,ddl1);
    Bm(3,ddl1+1) = r2 * DPHI(1,ddl1+1);
  end

% flexion : 15 ddl : aux 3 coins (U3,R1,R2) et aux 3 aretes (R1,R2)
% mais sans la condition discrete de Kirchhoff, on se sert uniquement
% de 12 ddl : aux 3 coins (R1,R2) et aux 3 aretes (R1,R2).
  Bf = zeros(3,21);
  for ino = 1:nbno
    ddl1 = (ino-1)*5+4; % R1 du noeud reel coin ino
    ddl2 = ddl1 + 1;    % R2 du noeud reel coin ino
    Bf(1,ddl2) = DPHI(1,ddl2);
    Bf(2,ddl1) = -DPHI(2,ddl1);
    Bf(3,ddl1) = -r2 * DPHI(1,ddl1);
    Bf(3,ddl2) = r2 * DPHI(2,ddl2);
%    Bf(1,ddl1) = DPHI(1,ddl1);
%    Bf(2,ddl2) = DPHI(2,ddl2);
%    Bf(3,ddl1) = r2 * DPHI(2,ddl1);
%    Bf(3,ddl2) = r2 * DPHI(1,ddl2);
  end
  for ino = 1:nbno
    ddl1 = (ino-1)*2+5*nbno+1; % R1 du noeud virtuel arete ino
    ddl2 = ddl1 + 1;           % R2 du noeud virtuel arete ino
    Bf(1,ddl2) = DPHI(1,ddl2);
    Bf(2,ddl1) = -DPHI(2,ddl1);
    Bf(3,ddl1) = -r2 * DPHI(1,ddl1);
    Bf(3,ddl2) = r2 * DPHI(2,ddl2);
%    Bf(1,ddl1) = DPHI(1,ddl1);
%    Bf(2,ddl2) = DPHI(2,ddl2);
%    Bf(3,ddl1) = r2 * DPHI(2,ddl1);
%    Bf(3,ddl2) = r2 * DPHI(1,ddl2);
  end

  B = [Bm ; Bf];
  clear Bm Bf;

else
  error('Pas prevu dans le cas general autre que DKQ...')
end
