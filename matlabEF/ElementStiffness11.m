function [KE,BE,bbE,varargout] = ElementStiffness11( ...
                                          modle1,D1,intge1, ...
                                          xcoorel,mode,varargin)
% Stiffness matrix and related operators of one element
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 05 / 08 / 2002
%   Ajout des B
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 11 / 2002
%   Bug corrige au niveau des numerotations
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%   Ajout du mode d'analyse BARR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   Modification composantes 2D ZZ
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 12 / 2003
%   Ajout du mode d'analyse TIMO
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2003
%   Ajout du mode d'analyse POUT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 04 / 2004
%   Matrice de rotation locale pour les elements de poutre
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 04 / 2004
%   Argument optionnel de geometrie locale pour elements de structure
%   (pour la rotation locale)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 04 / 2004
%   Passage au 3D pour les barres
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 07 / 2004
%   Ajout du mode d'analyse DKIR (Kirchhoff discret, plaques planes)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 05 / 2005
%   Ajout du mode d'analyse JOIN
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 12 / 2005
%   Ajout de la rigidite de spin pour les plaques DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 17 / 01 / 2006
%   Passage au 3D pour les POUT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2006
%   Ajout de la reconstruction d'une rotation quadratique pour DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  LaMCoS                           le 08 / 10 / 2010
%   Ajout du mode d'analyse DSHE (cisaillement discret, plaques planes)
%
% Integrale du produit des gradients symetriques des fonctions de forme
% Et matrice B associee BE
% Et pour la champ de deformation, matrice bbE
%
% Inputs
%   modle1		Elemental model
%   D1(:,:,nptg)	Hooke's operator (material behavior)
%   intge1		Elemental integration description
%   xcoorel(nbno,idim)	Coordinates of nodes of current element
%   mode		Mode of analysis
% Optional inputs
%   T(ino*nvec,idim)	Elemental local geometry information at nodes,
%                       mandatory for structural elements
%   spin1		Percentage of artificial spin stiffness for plates
% Outputs
%   KE(nbddl,nbddl)	Stiffness matrix of the element
%   BE(:,nbddl)		Dual generalized operator
%   bbE(:,nbddl)	Primal operator
% Optional outputs
%   CCE(nbddlQuadRot,nbddl)	For DKIR elements, reconstruction of
%                               quadratic rotations
%
% Comments
%   For beam elements (BARR,POUT,TIMO), T(ino,idim) contains the tangential
%   vector at each node.
%   For plate elements (DKIR), T(nvec,idim) contains the local basis,
%   constant on all the element. T(3,:) is the normal.
%   For joint elements (JOIN) it depends if one is in 2D or 3D.
%
%   For structural elements,
%   the rotation matrix is not local (i.e. at material point)
%   due to the interpolation of local dofs.
%   This is strange.. to be checked for SEG3 TIMO or POUT elements!
%   The Hooke operator is expressed in the geometric local basis
%   i.e. the global basis for massive elements, and the local
%   basis T for structural elements

% POUR LA VERSION SUIVANTE : nbnni NE SERT PAS

nin1 = nargin-5;
clear T;
spin1 = 0.;
for in1 = 1:nin1
  if ~all(size(varargin{in1})==[1 1])
    T = varargin{in1};
  else
    spin1 =  varargin{in1};
    if ((spin1 > 0.1) || (spin1 < 0.))
      spin1
      error('Check the artifical added spin stiffness')
    end
    if ~strcmp(mode,'DKIR')
      mode
      error('Artifical spin stiffness only for DKIR')
    end
  end
end

nou1 = nargout-3;
if (nou1 && ~strcmp(mode,'DKIR'))
  nou1
  mode
  error('Optional outputs for DKIR only')
end

% nbno : Number of primal or dual nodes
% idim : Physical dimension (real space)
[nbno idim] = size(xcoorel);

% idimr : Reference dimension (reference space)
% nbnni : Number of primal and dual basis functions
% nbptg : Number of integration points
[idimr nbnni nbptg] = size(intge1.DPHI);

if (length(modle1.DDLP) ~= length(modle1.DDLD))
  modle1.DDLP
  modle1.DDLD
  error('Different primal and dual physical components not implemented')
end

%% % nbcop : Number of primal and dual components for physical field
%% nbcop = length(modle1.DDLP);

if (length(modle1.NDDP) ~= length(modle1.NDDD))
  modle1.NDDP
  modle1.NDDD
  error('Different primal and dual ddl not implemented')
end

if ~all(modle1.NNOP == modle1.NNOD)
  modle1.NNOP
  modle1.NNOD
  error('Different primal and dual nodes not implemented')
end
if ~all(modle1.NNIP == modle1.NNID)
  modle1.NNIP
  modle1.NNID
  error('Different primal and dual basis functions not implemented')
end
% nbddl : Number of primal and dual ddl for local physical field
nbddl = length(modle1.NNIP);


% Size within local (physical) space
KE  = zeros(nbddl,nbddl);
BE  = zeros(0,nbddl);
bbE = zeros(0,nbddl);

% Coordinate of nodes that participate to the transformation
% of the element
xcoort1 = xcoorel(modle1.NNOT,:);

% list of node numbers
nnop = modle1.NNOP;
% list of dof names
nddp  = modle1.NDDP;
% list of basis functions
nnip  = modle1.NNIP;

% list of basis functions for element transformation
nnit = modle1.NNIT;


% Case of structural element with local basis coordinates
% """""""""""""""""""""""""""""""""""""""""""""""""""""""
% Rotation operator Q to get from local dof to global dof
if strcmp(mode,'BARR')

% Local basis rotation at nodes of the element
% Assumption of correct ordering of dof
  nbddl_loc = 1 * nbno; % 1 local dof (U1)
  Q = sparse(nbddl_loc,nbddl);
  iddl = 0;
  for ino = 1:nbno
    Q(ino,iddl+1:iddl+idim) = T(ino,:);
    iddl = iddl + idim;
  end

elseif (strcmp(mode,'TIMO') || strcmp(mode,'POUT'))
  if (idim == 2)

%   Local basis rotation at nodes of the element
%   Assumption of correct ordering of dof
    nbddl_loc = nbddl; % local dofs
    Q = sparse(nbddl_loc,nbddl);
    for ino = 1:nbno
      N = ProdVect(T(ino,:)'); % normal vector (2D: directly rotated)
      lddl = [(ino-1)*3+1:(ino-1)*3+3];
      Q(lddl,lddl) = [T(ino,:) 0.
                      N'       0.
                      0.  0.   1.];
    end

  else
    nbddl_loc = nbddl; % local dofs
    Q = sparse(nbddl_loc,nbddl);
    for ino = 1:nbno
%      N1 = T((ino-1)*3+1,:);
%      N2 = T((ino-1)*3+2,:);
%      N3 = T((ino-1)*3+3,:);
      N1 = T(ino,:);
      N2 = ProdVect([[0. 0. 1.]' N1']); % if N1=x then N2=y
      n2 = norm(N2);
      if (n2 < 1.e-5);
        N2 = ProdVect([[0. 1. 0.]' N1']); % if N1=z then N2=x
        n2 = norm(N2);
        if (n2 < 1.e-5);
	  N1
	  error('pb de norme de vecteurs locaux')
	end
      end
      N2 = N2' / n2;
      N3 = ProdVect([N1' N2']); N3 = N3' / norm(N3);
      lddl = [(ino-1)*6+1:(ino-1)*6+6];
      Q(lddl,lddl) = [N1 0. 0. 0.
                      N2 0. 0. 0.
                      N3 0. 0. 0.
                      0. 0. 0. N1
                      0. 0. 0. N2
                      0. 0. 0. N3];
    end
  end

elseif strcmp(mode,'JOIN')
  if (idim == 2)

%   Local basis rotation at nodes of the element
%   Assumption of correct ordering of dof
    nbddl_loc = nbddl; % local dofs
    Q = sparse(nbddl_loc,nbddl);
    for ino = 1:nbno
      N = ProdVect(T(ino,:)'); % normal vector (2D: directly rotated)
      lddl = [(ino-1)*2+1:(ino-1)*2+2];
      Q(lddl,lddl) = [T(ino,:)
                      N'      ];
    end

  else
    error('more than 2D not implemented for JOIN')
  end

elseif strcmp(mode,'DKIR')
  if (idim == 3)

%   Local basis rotation at nodes of the element
%   Assumption of correct ordering of dof
    if (nbddl ~= 7*nbno)
      nbddl
      nbno
      error('test failed....')
    end
%    nbddl = 6*nbno;     % 6 global dofs (UX,UY,UZ,RX,RY,RZ)
    nbddl_loc = 5*nbno; % 5 local dofs (2 for membrane, 3 for bending)
                        % (U1,U2,U3,R1,R2)
    N1 = T(1,:);
    N2 = T(2,:);
    N3 = T(3,:);
    Q = sparse(nbddl_loc,6*nbno);
    for ino = 1:nbno
      lddl = [(ino-1)*6+1:(ino-1)*6+6];
      lddl_loc = [(ino-1)*5+1:(ino-1)*5+5];
      if (1 == 0)
      Q(lddl_loc,lddl) = [N1       0. 0. 0.
                          N2       0. 0. 0.
                          N3       0. 0. 0.
                          0. 0. 0. N1+spin1*N3
                          0. 0. 0. N2+spin1*N3];
      end
      Q(lddl_loc,lddl) = [N1       0. 0. 0.
                          N2       0. 0. 0.
                          N3       0. 0. 0.
                          0. 0. 0. N1
                          0. 0. 0. N2];
    end
    KE  = zeros(nbddl_loc,nbddl_loc);
    BE  = zeros(0,nbddl_loc);
    bbE = zeros(0,nbddl_loc);

  else
    error('less than 3D not implemented for DKIR')
  end

elseif strcmp(mode,'AXIS')
  if (idim ~= 2)
    idim
    error('Dimension should be 2 for AXIS')
  end
% Real coordinates of integration points
  xcoorptg = intge1.PHI' * xcoorel;

end


% Stiffness in local basis
% """"""""""""""""""""""""
% Loop on integration points
for ptg = 1:nbptg
  wptg  = intge1.WEIGHT(ptg);

% Transformation for isoparametric elements
% """""""""""""""""""""""""""""""""""""""""
  dphixt1 = intge1.DPHI(:,nnit,ptg);
  [Mjaco,Jaco] = LocalJaco2(dphixt1,xcoort1);
  if (ptg == 1)
    Jaco0 = Jaco;
  else
    if (Jaco0 * Jaco < 0.)
      Jaco0
      ptg
      Jaco
      error('Change of sign in Jacobian')
    end
  end

% B-matrix (operator to be integrated)
% """"""""
% (here: same for dual and primal)
  if (idim == idimr)
    Ijaco = inv(Mjaco);
    dphix = intge1.DPHI(:,nnip,ptg);
    dphiX = Ijaco * dphix;
    B = LocalB4(dphiX,nnip,nddp);
    if strcmp(mode,'AXIS')
      if idim ~= 2
        idim
        error('Only dimension 2 for AXIS mode')
      end
      junk = intge1.PHI(:,ptg) / xcoorptg(ptg,1);
      B(4,1:2:end) = junk';
    end

  elseif ((idim == 2) && (idimr == 1))
%   Element 1D dans un espace 2D
    Ijaco = 1./Jaco;
    if strcmp(mode,'BARR')
      dphix = intge1.DPHI(:,nnip,ptg);
      dphiX = Ijaco * dphix;
      B = dphiX;
    elseif strcmp(mode,'TIMO')
      dphix = intge1.DPHI(:,nnip,ptg);
      dphiX = Ijaco * dphix;
      phi   = intge1.PHI(nnip,ptg);
      B = LocalB_Timo(dphiX,nnip,nddp,phi,idim);
    elseif strcmp(mode,'POUT')
      phi    = intge1.PHI(nnip,ptg);
      dphix  = intge1.DPHI(:,nnip,ptg);
      dphiX  = Ijaco * dphix;
      ddphix = intge1.DDPHI(:,:,nnip,ptg);
      ddphiX = Ijaco * Ijaco * ddphix;
      B = LocalB_Pout(dphiX,nnip,nddp,phi,idim,ddphiX);
%     Problem of rotational dofs: has to be multiply by transformation
%     gradient to be interpreted as derivates in real space
      for ino = 1:nbno
        lddl = [(ino-1)*3+1:(ino-1)*3+3];
        B(:,lddl(3)) = Jaco * B(:,lddl(3));
      end
    elseif strcmp(mode,'JOIN')
      if idim~=2
        error('2D only...')
      end
      phi    = intge1.PHI(nnip,ptg);
      for i=1:2
        j=find(nddp==i);
        B(i,j) = phi(j,:)';
      end
    else
      mode
      error('1D mode not known...')
    end

  elseif ((idim == 3) && (idimr == 1))
%   Element 1D dans un espace 3D
    Ijaco = 1./Jaco;
    if strcmp(mode,'BARR')
      dphix = intge1.DPHI(:,nnip,ptg);
      dphiX = Ijaco * dphix;
      B = dphiX;
    elseif strcmp(mode,'POUT')
%%      A VERIFIER
      phi    = intge1.PHI(nnip,ptg);
      dphix  = intge1.DPHI(:,nnip,ptg);
      dphiX  = Ijaco * dphix;
      ddphix = intge1.DDPHI(:,:,nnip,ptg);
      ddphiX = Ijaco * Ijaco * ddphix;
      B = LocalB_Pout(dphiX,nnip,nddp,phi,idim,ddphiX);
%     Problem of rotational dofs: has to be multiply by transformation
%     gradient to be interpreted as derivates in real space
%     but not for torsional dof
      for ino = 1:nbno
        lddl = [(ino-1)*6+1:(ino-1)*6+6];
        B(:,lddl(5:6)) = Jaco * B(:,lddl(5:6));
      end
    else
      mode
      error('1D mode not known... bis')
    end

  elseif ((idim == 3) && (idimr == 2))
%   Element 2D dans un espace 3D
    N1 = T(1,:)';
    N2 = T(2,:)';
    N3 = T(3,:)';
    F = Mjaco';
    Fnew = [N1 N2]' * F;
    Ijaco = inv(Fnew');
    if strcmp(mode,'DKIR')
      phi    = intge1.PHI(nnip,ptg);
      dphix  = intge1.DPHI(:,nnip,ptg);
      dphiX  = Ijaco * dphix;
      ddphix = intge1.DDPHI(:,:,nnip,ptg);
      ddphiX = zeros(size(ddphix));
      for i = 1:size(ddphix,3)
        ddphiX(:,:,i) = Ijaco * Ijaco * ddphix(:,:,i);
      end
      B = LocalB_Dkir(dphiX,nnip,nddp,phi,idim,ddphiX);
%     Virtual node dofs have to be eliminated with discrete Kirchhoff
%     condition (see Batoz tome 2, with different notations: herein,
%     rotational dofs are true rotations)
%     tangent-to-the-edge rotation should be linear per edge,
%     normal-to-the-edge rotation is linked to deflexion U3
%     (integral version)
      if (nbno == 4)
%       DKQ: formerly 28 dofs, in fine 20 dofs
        CC = sparse(8,20);
      elseif (nbno == 3)
%       DKT: formerly 21 dofs, in fine 15 dofs
        CC = sparse(6,15);
      else
        error('only for DKQ and DKT up to now')
      end
      for iedge = 1:nbno
% On devrait sortir cette boucle de la boucle sur les ptg
        i = iedge; j = mod(iedge,nbno)+1;
        ddlU3i = (i-1)*5+3; ddlR1i = ddlU3i+1; ddlR2i = ddlR1i+1;
        ddlU3j = (j-1)*5+3; ddlR1j = ddlU3j+1; ddlR2j = ddlR1j+1;
        ddlR1  = (iedge-1)*2+1; ddlR2  = ddlR1+1;
        xcoorji = [N1 N2]' * (xcoorel(j,:)'-xcoorel(i,:)');
        xji = xcoorji(1,1); yji = xcoorji(2,1);
        Ledge = norm(xcoorji);
        Cedge = xji / Ledge; Sedge = yji / Ledge;
%       Q = [Cedge Sedge ; Sedge -Cedge];
%       Rt        R1
%       Rn   = Q  R2    , Q^-1 = Q
%       Conditions: Rt = 0,
%                   Rn = 3/(2L)(U3j-U3i) - (3/4)(S(R1i+R1j) - C(R2i+R2j))
        CC(ddlR1:ddlR1+1,ddlU3i:ddlU3i+2) = ...
            [Sedge ; -Cedge] * (3./4.) * [-2./Ledge -Sedge +Cedge];
        CC(ddlR1:ddlR1+1,ddlU3j:ddlU3j+2) = ...
            [Sedge ; -Cedge] * (3./4.) * [2./Ledge -Sedge +Cedge];
      end
%     Virtual dof elimination
      B = B * [ sparse(eye(size(CC,2))); CC];

    elseif strcmp(mode,'DSHE')
      phi    = intge1.PHI(nnip,ptg);
      dphix  = intge1.DPHI(:,nnip,ptg);
      dphiX  = Ijaco * dphix;
      ddphix = intge1.DDPHI(:,:,nnip,ptg);
      ddphiX = zeros(size(ddphix));
      for i = 1:size(ddphix,3)
        ddphiX(:,:,i) = Ijaco * Ijaco * ddphix(:,:,i);
      end
      B = LocalB_Dshe(dphiX,nnip,nddp,phi,idim,ddphiX);
%     Virtual node dofs have to be eliminated with discrete shear
%     condition (see Batoz tome 2, with different notations: herein,
%     rotational dofs are true rotations)
%     tangent-to-the-edge rotation should be linear per edge,
%     normal-to-the-edge rotation is linked to deflexion U3
%     (integral version)
	error('TO BE DONE FOR DSHE')
      if (nbno == 4)
%       DKQ: formerly 28 dofs, in fine 20 dofs
        CC = sparse(8,20);
      elseif (nbno == 3)
%       DKT: formerly 21 dofs, in fine 15 dofs
        CC = sparse(6,15);
      else
        error('only for DKQ and DKT up to now')
      end
      for iedge = 1:nbno
% On devrait sortir cette boucle de la boucle sur les ptg
        i = iedge; j = mod(iedge,nbno)+1;
        ddlU3i = (i-1)*5+3; ddlR1i = ddlU3i+1; ddlR2i = ddlR1i+1;
        ddlU3j = (j-1)*5+3; ddlR1j = ddlU3j+1; ddlR2j = ddlR1j+1;
        ddlR1  = (iedge-1)*2+1; ddlR2  = ddlR1+1;
        xcoorji = [N1 N2]' * (xcoorel(j,:)'-xcoorel(i,:)');
        xji = xcoorji(1,1); yji = xcoorji(2,1);
        Ledge = norm(xcoorji);
        Cedge = xji / Ledge; Sedge = yji / Ledge;
%       Q = [Cedge Sedge ; Sedge -Cedge];
%       Rt        R1
%       Rn   = Q  R2    , Q^-1 = Q
%       Conditions: Rt = 0,
%                   Rn = 3/(2L)(U3j-U3i) - (3/4)(S(R1i+R1j) - C(R2i+R2j))
        CC(ddlR1:ddlR1+1,ddlU3i:ddlU3i+2) = ...
            [Sedge ; -Cedge] * (3./4.) * [-2./Ledge -Sedge +Cedge];
        CC(ddlR1:ddlR1+1,ddlU3j:ddlU3j+2) = ...
            [Sedge ; -Cedge] * (3./4.) * [2./Ledge -Sedge +Cedge];
      end
%     Virtual dof elimination
      B = B * [ sparse(eye(size(CC,2))); CC];

    else
      mode
      error('2D mode not known')
    end
  else
    idim
    idimr
    error('plongement pas prevu')
  end
  if strcmp(mode,'AXIS')
%   For axisymmetry, the 'R' coordinate is present as well as 2pi
    wptg = wptg*xcoorptg(ptg,1)*2.*pi;
  end

% Operateurs elementaires
% """""""""""""""""""""""
% Rigidite elementaire
  KE = KE + ((wptg * abs(Jaco)) * (B' * (D1(:,:,ptg) * B)));
% Matrice BE associee
  BE = [BE ; ((wptg * abs(Jaco)) * B)];
% Et pour la deformation
  bbE = [bbE ; B];
end 

% Final rotations for structural elements
% """""""""""""""""""""""""""""""""""""""
if ((strcmp(mode,'BARR') || ...
     strcmp(mode,'TIMO') || strcmp(mode,'POUT') || ...
     strcmp(mode,'DKIR') || ...
     strcmp(mode,'JOIN')))
% Back into global coordinate basis
  KE = Q' * KE * Q;
  BE = BE * Q;
  bbE = bbE * Q;
  if nou1
    CCE = CC * Q;
    varargout{1} = CCE;
  end
end


if strcmp(mode,'DKIR')
% Additional spin stiffness for plates
  if (idim == 3)
%   6 global dofs (UX,UY,UZ,RX,RY,RZ)
    for ino = 1:nbno
      ddlrot = [4 5 6] + ((ino-1)*6);
      Kadd1 = spin1 * sum(diag(KE(ddlrot,ddlrot)))/length(ddlrot);
%%DD 02/11/06 FAUTE... MERCI AHMAD      KE(ddlrot,ddlrot) = KE(ddlrot,ddlrot) + (Kadd1 * (N3' * N3));
      KE(ddlrot,ddlrot) = KE(ddlrot,ddlrot) + (Kadd1 * (N3 * N3'));
    end
  else
    error('less than 3D not implemented for DKIR')
  end
end
