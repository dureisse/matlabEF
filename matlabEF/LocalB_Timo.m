function [B] = LocalB_Timo(DPHI,NNIP,NDDP,PHI,idim);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 12 / 2003
%
% Local operator derivative of kinematic torsor of the basis functions
% suitable for Timoshenko beam
%
% DPHIX(i,nni) = valeur de la derivee par rapport a la coordonne i dans
%                l'espace de reference, de la fonction de forme nni
% NNIP = liste des numeros  des fonctions de formes associees aux ddl
%        locaux de l'element fini
% NDDP = liste des degres de liberte primaux associes
% PHI(nni) = valeur de la fonction de forme nni
% idim = Dimension of the physical space
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
    nbcop  = 3;            % Number of primal components for the field U
    nbcopB = 3;            % Number of components on the derivated field
    Q      = [0. ; -1.];   % Operator e1 vectoriel .
  elseif (idim == 3)
    nbcop  = 6;         % Number of primal components for the field U
    nbcopB = 6;         % Number of components on the derivated field
    Q      = [0. 0. 0.
              0. 0. -1.
              0. 1. 0.]; % Operator e1 vectoriel .
  else
    idim
    error('Bad idim')
  end

  B = zeros(nbcopB,nbddl);

% Loop on node
  ddl = 0;
  nbno = nbddl / nbcop;
  for ino = 1:nbno
    for i = 1:nbcop
      ddl = ddl + 1;
      B(i,ddl) = DPHI(1,ddl); % Assuming right order for dof
    end
    ddltheta = (ino-1)*nbcop+idim+1:(ino-1)*nbcop+nbcop;
    B(1:idim,ddltheta) = Q * PHI(NNIP(ddltheta))';
  end
