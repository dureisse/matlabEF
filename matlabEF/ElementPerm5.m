function [KE,BE,bbE] = ElementPerm5(modle1,D1,intge1,xcoorel,mode)
% Elemental permeability matrix and related operators of one element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 07 / 2008
%   Ajout du mode BARR (unidimensionnel)
%
% Integrale du produit des gradients des fonctions de forme
% Et matrice B associee BE
% Et pour la champ de deformation associe, matrice bbE
% (voir ElementStiffness5, ElementMass4, ElementCompress5)

% nbno : Number of primal or dual nodes (not used)
% idim : Physical dimension (real space)
[nbno idim] = size(xcoorel);

% idimr : Reference dimension (reference space)
% nbnni : Number of primal and dual basis functions (not used)
% nbptg : Number of integration points
[idimr nbnni nbptg] = size(intge1.DPHI);

if (length(modle1.DDLP) ~= length(modle1.DDLD))
  modle1.DDLP
  modle1.DDLD
  error('Different primal and dual physical components not implemented')
end

if (length(modle1.NDDP) ~= length(modle1.NDDD))
  modle1.NDDP
  modle1.NDDD
  error('Different primal and dual ddl not implemented')
end
nbddl = length(modle1.NDDP); % Number of primal and dual dof
                             % for physical field

%%% A REVOIR
nbcopC = idim; % Number of components


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

BE = zeros(nbddl,0);
KE = zeros(nbddl,nbddl);
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

if strcmp(mode,'AXIS')
  if (idim ~= 2)
    idim
    error('Dimension should be 2 for AXIS')
  end
% Real coordinates of integration points
  xcoorptg = intge1.PHI' * xcoorel;
end

% Permeability (local basis for gradient and its dual counterpart)
% """"""""""""
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

% C-matrix (operator to be integrated)
% """"""""
% (here: same for dual and primal)
  if (idim == idimr)
    Ijaco = inv(Mjaco);
    dphix = intge1.DPHI(:,nnip,ptg); % basis functions (correct dof order)
    dphiX = Ijaco * dphix;
    C = zeros(nbcopC,nbddl);
    for ddl = 1:nbddl
      C(:,ddl) = dphiX(:,ddl);
    end
    if strcmp(mode,'AXIS')
%     For axisymmetry no modification in pressure gradient, since its
%     last component is 0 in basis (R,Z,T)
%     For axisymmetry, the 'R' coordinate is present as well as 2pi
      wptg = wptg*xcoorptg(ptg,1)*2.*pi;
    end

  elseif (idim == 2 & idimr == 1)
%   Element 1D dans un espace 2D
    Ijaco = 1./Jaco;
    if strcmp(mode,'BARR')
      dphix = intge1.DPHI(:,nnip,ptg); % basis functions (correct dof order)
      dphiX = Ijaco * dphix;
      C = dphiX;
    else
      mode
      error('1D mode not known...')
    end

  elseif (idim == 3 & idimr == 1)
%   Element 1D dans un espace 3D
    Ijaco = 1./Jaco;
    if strcmp(mode,'BARR')
      dphix = intge1.DPHI(:,nnip,ptg); % basis functions (correct dof order)
      dphiX = Ijaco * dphix;
      C = dphiX;
    else
      mode
      error('1D mode not known...')
    end

  else
    idim
    idimr
    error('plongement pas prevu')
  end

% Operateurs elementaires
% """""""""""""""""""""""
% Elemental matrix KE
  KE = KE + (wptg * abs(Jaco)) * (C' * (D1(:,:,ptg) * C));
% Associated BE matrix
  BE = [BE (wptg * abs(Jaco)) * C'];
% And for associated dual components
  bbE = [bbE ; C];
end 
BE = BE';
