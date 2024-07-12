function [KE,BE,bbE] = ElementStiffnessDST( ...
                                          modle1,D1,intge1, ...
                                          xcoorel,mode,varargin)
% Stiffness matrix and related operators of DST element
% inspired by routines of Cast3M
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 13 / 10 / 2010
%
% Integrale du produit des gradients symetriques des fonctions de forme
% Et matrice B associee BE
% Et pour la champ de deformation, matrice bbE
%
% Inputs
%   modle1              Elemental model
%   D1(:,:,nptg)        Hooke's operator (material behavior)
%   intge1              Elemental integration description
%   xcoorel(nbno,idim)  Coordinates of nodes of current element
%   mode                Mode of analysis
% Optional inputs
%   T(ino*nvec,idim)    Elemental local geometry information at nodes,
%                       mandatory for structural elements
%   spin1               Percentage of artificial spin stiffness for plates
% Outputs
%   KE(nbddl,nbddl)     Stiffness matrix of the element
%   BE(:,nbddl)         Dual generalized operator
%   bbE(:,nbddl)        Primal operator
%
% Comments
%   For plate elements (DSHE), T(nvec,idim) contains the local basis,
%   constant on all the element. T(3,:) is the normal.
%
%   For structural elements,
%   the rotation matrix is not local (i.e. at material point)
%   due to the interpolation of local dofs.
%   This is strange.. to be checked for SEG3 TIMO or POUT elements!
%   The Hooke operator is expressed in the geometric local basis
%   i.e. the global basis for massive elements, and the local
%   basis T for structural elements

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
    if ~strcmp(mode,'DSHE')
      mode
      error('Artifical spin stiffness only for DSHE')
    end
  end
end

if ~strcmp(mode,'DSHE')
  error('Only DHSE in ElementStiffnessDST.m')
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

% Node coordinates (idim,nbno)
XE = xcoorel';
POIGAU = intge1.WEIGHT;

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

%   Element 2D dans un espace 3D
    N1 = T(1,:)';
    N2 = T(2,:)';
    N3 = T(3,:)';
    F = Mjaco';
    Fnew = [N1 N2]' * F;
    Ijaco = inv(Fnew');

    XEL = [N1 N2]' * XE;

    DDHOMU = D1(:,:,ptg);
% Voigt and not Cowin here
    DDHOMU(3,:) = 0.5*DDHOMU(3,:);
    DDHOMU(6,:) = 0.5*DDHOMU(6,:);
    DDHOMU(7,:) = 0.5*DDHOMU(7,:);
    DDHOMU(8,:) = 0.5*DDHOMU(8,:);
    [HS4,HS5,HS6,BGENE] = rcdst(XEL,DDHOMU);
    QSIGAU = intge1.COOR(:,1)';
    ETAGAU = intge1.COOR(:,2)';
    SHPTOT = zeros(6,3,nbptg);
    SHPTOT(1,:,:) = intge1.PHI(1:3,:);
    SHPTOT(2,:,:) = intge1.DPHI(1,1:3,:);
    SHPTOT(3,:,:) = intge1.DPHI(2,1:3,:);
    IGAU = ptg;
    [BGENE,DJAC] = bmfdst(IGAU,XEL,QSIGAU,ETAGAU,SHPTOT,HS4,HS5,HS6,BGENE);
    DJAC=DJAC*POIGAU(IGAU);
    EXCEN = 0.0;
%   ON MODIFIE LA MATRICE B EN CAS D'EXCENTREMENT
      for IJL = 1:3
        for IJC = 1:size(BGENE,2)
          BGENE(IJL,IJC)=BGENE(IJL,IJC)+EXCEN*BGENE(IJL+3,IJC);
        end
      end
%      CALL BDBS1(BGENE,DJAC,DDHOMU,LRE,NSTRS,REL,MFR,IFOUR,MATE,
%     1     IGAU,IMAT,EXCEN)
%
%
%      REL(6,6)=REL(5,5)*1.D-7
%      REL(12,12)=REL(6,6)
%      REL(18,18)=REL(6,6)
%      ICOM=0
%      on remplit les terme symetriques...

      B = BGENE;
% Back to Cowin
      B(3,:) = (sqrt(2.)/2.)*B(3,:);
      B(6,:) = (sqrt(2.)/2.)*B(6,:);
      B(7,:) = (sqrt(2.)/2.)*B(7,:);
      B(8,:) = (sqrt(2.)/2.)*B(8,:);
%truc = B' * (D1(:,:,ptg) * B);
%troc = BGENE' * DDHOMU * BGENE;
%max(max(abs(truc-troc)))
% Get rid of spin local dofs
      B = B(:,[1:5 7:11 13:17]);

% Operateurs elementaires
% """""""""""""""""""""""
% Rigidite elementaire
  KE = KE + ((wptg * abs(Jaco)) * (B' * (D1(:,:,ptg) * B)));
% Matrice BE associee
  BE = [BE ; ((wptg * abs(Jaco)) * B)];
% Et pour la deformation
  bbE = [bbE ; B];


%  troc*wptg*abs(Jaco)

end 

% Final rotations for structural elements
% """""""""""""""""""""""""""""""""""""""
if (strcmp(mode,'BARR') || ...
    strcmp(mode,'TIMO') || strcmp(mode,'POUT') || ...
    strcmp(mode,'DKIR') || strcmp(mode,'DSHE') || ...
    strcmp(mode,'JOIN'))
% Back into global coordinate basis
  KE = Q' * KE * Q;
  BE = BE * Q;
  bbE = bbE * Q;
end

%KE
%error('a suivre')

if (strcmp(mode,'DKIR') || strcmp(mode,'DSHE'))
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
    error('less than 3D not implemented for DKIR or DSHE')
  end
end
