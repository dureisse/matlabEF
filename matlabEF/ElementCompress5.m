function [KE,BE,bbE] = ElementCompress5(modle1,D1,intge1,xcoorel,mode, ...
                                        varargin)
% Elemental compressibility matrices and related operators
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 05 / 2008
%   Ajout du mode DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 07 / 2008
%   Ajout du mode BARR
%
% Elemental compressibility matrix KE
% Integrale du produit des fonctions de forme
% Et matrice B associee BE
% Et pour la champ de deformation associe, matrice bbE
% (voir ElementStiffness5, ElementMass4)
%
% Inputs
%   modle1	Part of the model for the considered zone
%   D1(:,:,nbpgt)	Behavior operator
%   intge1	Part of the integration segment for the considered zone
%   xcoorel(nbno,idim)	Real coordinates of nodes
%   mode	Analysis mode ('TRID','COPL','DEPL','AXIS','DKIR')
% Optional inputs
%   xcoorptg(nbptg,idim) For AXIS mode, real coordinates
%   			 of integration points
%   T(ino*nvec,idim)     For DKIR mode, elemental local geometry
%                        information at nodes
% Outputs
%   KE		Compressibility matrix of the element
%   BE		Generalized associated dual operator
%   bbE		Associated primal operator
%

nin1 = nargin-5;
if strcmp(mode,'AXIS')
  if nin1 ~= 1
    error('For AXIS mode, coordinates of integration points are required')
  else
    xcoorptg = varargin{1};
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

if (length(modle1.DDLP) ~= length(modle1.DDLD))
  modle1.DDLP
  modle1.DDLD
  error('Different primal and dual physical components not implemented')
end

% nbcop : Number of primal (and dual) components for physical field
nbcop = length(modle1.NCOP);
if length(modle1.NCOP) ~= length(modle1.NCOD)
  modle1.NCOP
  modle1.NCOD
  error('Not the same number of primal and dual components')
end
% Pre-treatment : who are the dofs with the same names as the comp
% (primal only: dual is supposed to be identical)
clear list_idx;
switch mode
  case {'TRID','COPL','DEPL','AXIS','DKIR','BARR'}
    for cop = 1:nbcop
%     Name of the component
      NCOP = modle1.COMP(modle1.NCOP(cop));
%     Where is the dof with the same name
      dop = findoccur(NCOP,modle1.DDLP);
%     dof contributing to this component
      idx = find(modle1.NDDP == dop);
      list_idx{cop} = idx;
    end
  otherwise,
    mode
    error('Mode not yet available')
end

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
%%DD nbddl = length(modle1.NDDP);

% Size within local (physical) space
KE = zeros(nbddl,nbddl);
BE = zeros(nbddl,0);
bbE = zeros(0,nbddl);

% Coordinate of nodes that participate to the transformation
% of the element
xcoort1 = xcoorel(modle1.NNOT,:);

% list of node numbers
nnop = modle1.NNOP;
% list of dof names
nddp = modle1.NDDP;
% list of basis functions
nnip = modle1.NNIP;
% list of basis functions for element transformation
nnit = modle1.NNIT;

% Loop on integration points
for ptg = 1:nbptg
  wptg  = intge1.WEIGHT(ptg);

% Transformation for isoparametric elements
% """""""""""""""""""""""""""""""""""""""""
  if ~isempty(nnit)
    dphixt1 = intge1.DPHI(:,nnit,ptg);
    [Mjaco,Jaco] = LocalJaco2(dphixt1,xcoort1);
    if (ptg == 1)
      Jaco0 = Jaco;
    else
      if (Jaco0 * Jaco < 0.)
        Jaco0
        ptg
        Jaco
        xcoorel
        error('Change of sign in Jacobian')
      end
    end
  end

  if ~isempty(nnip)

%   B-matrix (operator to be integrated)
%   """"""""
%   (here: same for dual and primal)
%    B = intge1.PHI(modle1.NNIP,ptg)';
    B = zeros(nbcop,nbddl);

    switch mode
      case {'TRID','COPL','DEPL','AXIS','DKIR','BARR'}
%       Loop on physical component
        for cop = 1:nbcop
%DD 07/09/11          idx = find(modle1.NDDP == cop);
          idx = list_idx{cop};
          B(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
        end
      otherwise,
        mode
        error('Mode not yet available bis')
    end
    if strcmp(mode,'AXIS')
%     For axisymmetry, the 'R' coordinate is present as well as 2pi
      wptg = wptg*xcoorptg(ptg,1)*2.*pi;
    end

%   Elemental operators
%   """""""""""""""""""
%   Elemental compressibility
    KE = KE + (wptg * abs(Jaco)) * (B' * (D1(:,:,ptg) * B));
%   And associated operators
    BE = [BE (wptg * abs(Jaco)) * B'];
    bbE = [bbE ; B];

  else
    error('Dirac integration not implemented foor compressibility')
  end

end 
BE = BE';
