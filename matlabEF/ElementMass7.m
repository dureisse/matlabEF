function [ME,BE,bbE] = ElementMass7(modle1,D1,intge1, ...
                                    xcoorel,mode,varargin);
% Mass matrix and related operators of one element
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 01 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 08 / 2002
%   ajout plongement 1D dans 2D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
%   ajout matrices BE,bbE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 01 / 2003
%   cas particulier ou il n'y a pas de fonction de forme
%   (integration discrete sur element POI1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 08 / 2005
%   Ajout du mode d'analyse POUT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 07 / 07 / 2006
%   Passage au 3D pour BARR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  LaMCoS                           le 28 / 06 / 2016
%   Passage au 3D pour POUT
%
% Inputs
%   D1			Volumic density matrix of the element
%                       D1(3,3,nbptg) for 3D
%                       D1(2,2,nbptg) for 2D
%   intge1		Integration information on the element
%   xcoorel		Element node coordinates
% Optional inputs
%   T(ino*nvec,idim)    Elemental local geometry information at nodes,
%                       mandatory for structural elements
% Outputs
%   ME(nbddl,nbddl)	Mass matrix of the element
%   BE(:,nbddl)         Dual generalized operator
%   bbE(:,nbddl)        Primal operator 
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
%   This is strange... to be checked for SEG3 TIMO or POUT elements!
%   The Hooke operator is expressed in the geometric local basis
%   i.e. the global basis for massive elements, and the local
%   basis T for structural elements

% Transformation from reference space to real space is done
% with inte1.NNOT and inte1.NNIT

nin1 = nargin-5;
if nin1 == 1
  T = varargin{1};
end

% nbno : Number of primal or dual nodes
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

if (length(modle1.NDDP) ~= length(modle1.NDDD))
  modle1.NDDP
  modle1.NDDD
  error('Different primal and dual dof not implemented')
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
% nbddl : Number of primal and dual dof for physical field
nbddl = length(modle1.NDDP);


% Size within local (physical) space
ME  = zeros(nbddl,nbddl);
BE  = zeros(nbddl,0);
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
  nbddl_loc = nbddl; % idim local dof (U1,U2) or (U1,U2,U3)
  Q = sparse(nbddl_loc,nbddl);
  switch idim
    case 2,
      iddl = 0;
      for ino = 1:nbno
        e1 = T(ino,:)';
        e2 = ProdVect(e1); % normal vector (2D: directly rotated)
        Q(iddl+1,iddl+1:iddl+idim) = e1';
        Q(iddl+2,iddl+1:iddl+idim) = e2';
        iddl = iddl + idim;
      end
    case 3,
%      disp('Warning: Mass for 3D BARR, isotropic case is assumed')
      iddl = 0;
      for ino = 1:nbno
        e1 = T(ino,:)';
%       Remainder for the local basis: no orientation is assumed
        MM = eye(idim) - e1*e1';
	[U,S,V] = svd(MM);
	e2 = U(:,1);
	e3 = U(:,2);
%       Direct basis
        if (det([e1 e2 e3]) < 0)
	  e0 = e2; e2 = e3; e3 = e0;
	end
        Q(iddl+1,iddl+1:iddl+idim) = e1';
        Q(iddl+2,iddl+1:iddl+idim) = e2';
        Q(iddl+3,iddl+1:iddl+idim) = e3';
        iddl = iddl + idim;
      end
  end
 
elseif strcmp(mode,'POUT')
  switch idim
      case 2,
%       Local basis rotation at nodes of the element
%       Assumption of correct ordering of dof
        nbddl_loc = nbddl; % local dofs
        Q = sparse(nbddl_loc,nbddl);
        for ino = 1:nbno
          e1 = T(ino,:)';
          e2 = ProdVect(e1); % normal vector (2D: directly rotated)
          lddl = [(ino-1)*3+1:(ino-1)*3+3];
          Q(lddl,lddl) = [e1'   0.
                          e2'   0.
                          0. 0. 1.];
        end
      case 3,
        nbddl_loc = nbddl; % local dofs
        Q = sparse(nbddl_loc,nbddl);
        % isotropic : only the tangent vector is given in T, complete
        for ino = 1:nbno
          e1 = T(ino,:)';
%%          tt = [e1 eye(idim)]; m = tt' * tt; bb = chol(m); tt = (B' \ tt')';
          tt = [e1 eye(idim)]; [qq,rr] = qr(tt);
          if norm(rr(:,1)/rr(1,1) - [1 0 0]') > 1.e-8
              error('pb in qr')
          end
          e2 = qq(:,2)*sign(rr(1,1));
          e3 = qq(:,3)*sign(rr(1,1));
          lddl = [(ino-1)*6+1:(ino-1)*6+6];
          Q(lddl,lddl) = [e1' 0 0 0 
                          e2' 0 0 0 
                          e3' 0 0 0
                          0 0 0 e1'
                          0 0 0 e2'
                          0 0 0 e3'];
        end
%         full(Q)
%         error('more than 2D not implemented for TIMO nor POUT')
      otherwise,
        error('Bad dimension')
  end

elseif strcmp(mode,'AXIS')
  if (idim ~= 2)
    idim
    error('Dimension should be 2 for AXIS')
  end
% Real coordinates of integration points
  xcoorptg = intge1.PHI' * xcoorel;
  disp('PASSER CA EN ARGUMENT OPTIONNEL')

end

 % Mass in local basis
% """""""""""""""""""
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
        error('Sign change of Jacobian')
      end
    end
  end

  if ~isempty(nnip)

%   N-matrix (operator to be integrated)
%   """"""""
%   (here: same for dual and primal)
    N = zeros(nbcop,nbddl);

    switch mode
      case {'TRID','COPL','DEPL','AXIS'}
%       Loop on physical component
        for cop = 1:nbcop
%         dof contributing to this component
          idx = find(modle1.NDDP == cop);
          N(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
        end
      case 'BARR',
%       Loop on physical component
        for cop = 1:nbcop
%         local dof contributing to this component
          idx = find(modle1.NDDPLOC == cop);
          N(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
	end
      case 'POUT',
	switch idim
	  case 2,
	    cop = 1; % U1 as usual
	      idx = find(modle1.NDDPLOC == cop);
	      N(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
        cop = 2; % U2 use both U2 and RZ local dof
	      idx2 = [find(modle1.NDDPLOC == cop) find(modle1.NDDPLOC == 3)];
	      N(cop,idx2) = intge1.PHI(modle1.NNIP(idx2),ptg)';
	    cop = 3; % RZ is the derivate of U2
	      N(cop,idx2) = intge1.DPHI(1,modle1.NNIP(idx2),ptg);
%       Problem of rotational dofs:
%       they have to be multiplied by transformation gradient to be
%       interpreted as derivatives in real space
	    for ino = 1:nbno
	      lddl = [(ino-1)*3+1:(ino-1)*3+3];
	      N(:,lddl(3)) = Jaco * N(:,lddl(3));
	    end
	    N(3,:) = (1 / Jaco) * N(3,:);
	  case 3,
	    cop = 1; % U1 as usual
	      idx = find(modle1.NDDPLOC == cop);
	      N(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
        cop = 2; % U2 use both U2 and R3 local dof
	      idx2 = [find(modle1.NDDPLOC == cop) find(modle1.NDDPLOC == 6)];
	      N(cop,idx2) = intge1.PHI(modle1.NNIP(idx2),ptg)';
        cop = 3; % U3 use both U3 and R2 local dof
	      idx2 = [find(modle1.NDDPLOC == cop) find(modle1.NDDPLOC == 2)];
	      N(cop,idx2) = intge1.PHI(modle1.NNIP(idx2),ptg)';
	    cop = 4; % R1 as usual
	      idx = find(modle1.NDDPLOC == cop);
	      N(cop,idx) = intge1.PHI(modle1.NNIP(idx),ptg)';
	    cop = 5; % R2 is the derivate of U3
	      idx = find(modle1.NDDPLOC == 3);
	      N(cop,idx) = intge1.DPHI(1,modle1.NNIP(idx),ptg);
	    cop = 6; % R3 is the derivate of U2
	      idx = find(modle1.NDDPLOC == 2);
	      N(cop,idx) = intge1.DPHI(1,modle1.NNIP(idx),ptg);
%       Problem of rotational dofs:
%       they have to be multiplied by transformation gradient to be
%       interpreted as derivatives in real space
	    for ino = 1:nbno
	      lddl = [(ino-1)*6+1:(ino-1)*6+6];
	      N(:,lddl([4 5 6])) = Jaco * N(:,lddl([4 5 6]));
	    end
	    N([4 5 6],:) = (1 / Jaco) * N([4 5 6],:);
%         full(N)          
% 	    error('To be checked in 3D... sorry')
	  otherwise,
	    idim
	    error('bad idim')
	end
      otherwise,
        mode
	error('Mode not yet available')
    end
    if strcmp(mode,'AXIS')
%     For axisymmetry, the 'R' coordinate is present as well as 2pi
      wptg = wptg*xcoorptg(ptg,1)*2.*pi;
    end

%   Elemental operators
%   """""""""""""""""""
%   Elemental mass
    ME  = ME + ((wptg*abs(Jaco)) * (N' * D1(:,:,ptg) * N));
%   And associated B matrices
    BE  = [BE (wptg*abs(Jaco)) * N'];
    bbE = [bbE ; N];

  else
%   Dirac integration: sum of discrete value
%   Elemental mass
    ME  = ME + (wptg * D1(:,:,ptg));
%   And associated B matrices
    BE  = [BE wptg*eye(nbddl,nbcop)];
    bbE = [bbE ; eye(nbcop,nbddl)];
  end

end 
BE = BE';

% Final rotations for structural elements
% """""""""""""""""""""""""""""""""""""""
if (strcmp(mode,'BARR') || strcmp(mode,'POUT'))
% Back to global coordinate basis
  ME = Q' * ME * Q;
  BE = BE * Q;
  bbE = bbE * Q;
end
