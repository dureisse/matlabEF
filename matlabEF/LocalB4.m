function [B] = LocalB4(DPHI,NNIP,NDDP);
% Local operator symmetric gradient of the basis functions
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 09 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 12 / 2004
%   Ajout au 3D
%
% DPHIX(i,nni) = valeur de la derivee par rapport a la coordonne i dans
%                l'espace de reference, de la fonction de forme nni
% NNIP = liste des numeros  des fonctions de formes associees aux ddl
%        locaux de l'element fini
% NDDP = liste des degres de liberte primaux associes
%
% Assumption: dofs are in right order (1 to idim)
%             Voigt notation is used for B
% Beware that for 2D plane stress, the last component EPZZ is set to 0
%   in fact, it could not be derived from displacement, because the
%   constitutive relation is needed to prescribe SMZZ = 0.

  r2 = sqrt(2.)/2.; % Due to Voigt notation

%DD 10/2/03  nbddl = length(NDDP); % Number of local ddl (within the element)
  nbddl = length(NNIP); % Number of local ddl (within the element)
  idimr = size(DPHI,1); % Dimension of the reference space
  idim = idimr;         % Dimension of the physical space
  nbcop = idim;         % Number of primal components for the field U
% Number of components on the derivated field (Voigt notation)
  if (idim == 2)
%% DD 03/11/25    nbcopB = 3;
    nbcopB = 4;
  elseif (idim == 3)
    nbcopB = 6;
  else
    idim
    error('bad idim')
  end

  B = zeros(nbcopB,nbddl);

% Loop on dof
  for ddl = 1:nbddl
    name = NDDP(ddl);
    i    = name;  % Assuming right order for dof
    B(i,ddl) = DPHI(i,ddl);
    if (idim == 2)
      B(3,ddl) = r2 * DPHI(3-i,ddl);
      B(4,ddl) = 0;
    elseif (idim == 3)
      switch i
        case 1,
          B(5,ddl) = r2 * DPHI(3,ddl);
          B(6,ddl) = r2 * DPHI(2,ddl);
        case 2,
          B(4,ddl) = r2 * DPHI(3,ddl);
          B(6,ddl) = r2 * DPHI(1,ddl);
        case 3,
          B(4,ddl) = r2 * DPHI(2,ddl);
          B(5,ddl) = r2 * DPHI(1,ddl);
        otherwise,
          error('Bad dof order...')
      end
    end
  end
