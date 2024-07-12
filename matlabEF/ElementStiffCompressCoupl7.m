function [NE,BE,bbE] = ElementStiffCompressCoupl7(modlce1,D1,TR1, ...
                                       intge1,xcoorel,mode,varargin)
% Elemental stiffness-compressibility coupling matrix and associated
% operators
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 11 / 2002
%   passage d'un modele couple, et rigidites elem hors diagonale
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 09 / 2007
%   Ajout mode AXIS

% (voir ElementStiffness5, ElementCompress5)

nin1 = nargin-6;
switch nin1
  case 0,
    clear xcoorptg;
  case 1,
    xcoorptg = varargin{1};
  otherwise,
    nin1
    error('Bad umber of optional input arguments')
end
if strcmp(mode,'AXIS')
  if nin1 ~= 1
    error('For AXIS mode, coordinates of integration points are required')
  end
else
  if nin1 ~= 0
    error('No optional arguments are needed, but one is provided')
  end
end

% nbno : Number of primal or dual nodes (not used)
% idim : Physical dimension (real space)
[nbno idim] = size(xcoorel);

% idimr : Reference dimension (reference space)
% nbnni : Number of primal and dual basis functions (not used)
% nbptg : Number of integration points
[idimr nbnni nbptg] = size(intge1.DPHI);

nbddls = length(modlce1.NDDP);
nbddlf = length(modlce1.NDDD);
% nbddl  = nbddls + nbddlf;

nnot1 = modlce1.NNOT;
nnit1 = modlce1.NNIT;

NE  = zeros(nbddlf,nbddls);
BE  = zeros(0,nbddls);
bbE = zeros(0,nbddls);

% Coordinate of nodes that participate to the transformation
% of the element
xcoort1 = xcoorel(nnot1,:);

% list of node numbers
nnop = modlce1.NNOP;
% list of dof names
nddp  = modlce1.NDDP;
% list of basis functions
nnip  = modlce1.NNIP;
nnid  = modlce1.NNID;

% list of basis functions for element transformation
nnit = modlce1.NNIT;

% Loop on integration points
for ptg = 1:nbptg
  wptg  = intge1.WEIGHT(ptg);

% Transformation for elements
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
  Ijaco = inv(Mjaco);

% C-matrix (operator to be integrated for the solid)
% (here: not the same for dual and primal; solid->primal)
  dphix = intge1.DPHI(:,nnip,ptg);
  dphiX = Ijaco * dphix;
  B = LocalB4(dphiX,nnip,nddp);
  if strcmp(mode,'AXIS')
    if idim ~= 2
      idim
      error('Only dimension 2 for AXIS mode')
    end
    junk = intge1.PHI(nnip(1:2:end),ptg) / xcoorptg(ptg,1);
    B(4,1:2:end) = junk';
  end
  C = TR1(:,:,ptg) * B;

% N-matrix (operator to be integrated for the fluid)
% (here: not the same for dual and primal; fluide->dual)
  N = intge1.PHI(nnid,ptg)';

  if strcmp(mode,'AXIS')
%   For axisymmetry, the 'R' coordinate is present as well as 2pi
    wptg = wptg*xcoorptg(ptg,1)*2.*pi;
  end

% Rigidite elementaire
  NE  = NE + (wptg * abs(Jaco)) * (N' * (D1(:,:,ptg) * C));
  BE  = [BE ; (wptg * abs(Jaco)) * C]; % c'est BE et non BE'
  bbE = [bbE ; C];
end
