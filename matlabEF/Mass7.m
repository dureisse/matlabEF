function [mass1,varargout] = Mass7(modl1,matr1,mail1,intg1, ...
                                   xcoor1,mode1,varargin)
% Elementary mass matrices, and associated operators
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 01 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
%   Ajout des operateurs B
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%   Ajout du mode BARR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 08 / 2005
%   Ajout du mode POUT et des arguments optionnels
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  LaMCoS                           le 28 / 06 / 2016
%   Ajout de la 3D pour la MASSE en mode POUT
%
% Inputs
%  modl1		Model
%  matr1		Material coefficients
%  mail1		Mesh
%  intg1		Integration information
%  xcoor1(nbno,idim)	Node coordinates
%  mode1		Analysis mode
% Optional inputs
%   'PrimOp'            Keyword for primal operator bu1
%   'GenDualOp'         Keyword for dual generalized operator bgamma1
%   'ConstiOp'          Keyword for constitutive operator d1
%   chamno1             Element field at nodes for information on local
%                       geometry
% Outputs
%  mass1		Elementary mass matrices
%  and a list of operators as requested by keywords
%
% Examples
%   mass1 = Mass7(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [mass1,bgamma1] = Mass7(modl1,matr1,mail1,intg1,xcoor1,mode1,'GenDualOp')
%   [mass1,bu1,bgamma1] = Mass7(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                              'PrimOp','GenDualOp')
%   [mass1,d1] = Mass7(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                      chamno1,'ConstiOp')
%
% Comments
%   Une fois assembles,
%   bgamma1' permet de passer d'un champ par element d'acceleration
%   gamma aux forces generalisees associees G : G = bgamma1' * gamma
%   bu1 permet de passer d'un champ par point de deplacement U a un
%   champ par elements de deplacement u : u = bu1 * U
%
% (see also Rigi9)

% Analyse optional input arguments,
% check order of operators and find if chamno1 is provided
nin1 = nargin-6;
nout = nargout-1;
 
clear chamno1;
listOp = zeros(1,3);
ipos1 = 0;
for in1 = 1:nin1
  argin1 = varargin{in1};
  if iscell(argin1)
    chamno1 = argin1;
  else
    ipos1 = ipos1 + 1;
    switch argin1
      case 'GenDualOp',
%        listOp(ipos1) = 1; ipos1 = ipos1 + 1;
        listOp(1) = ipos1;
      case 'PrimOp',
%        listOp(ipos1) = 2; ipos1 = ipos1 + 1;
        listOp(2) = ipos1;
      case 'ConstiOp',
%        listOp(ipos1) = 3; ipos1 = ipos1 + 1;
        listOp(3) = ipos1;
      otherwise,
        argin1
        error('unknown option')
    end
  end
end
 
if (ipos1 ~= nout)
  ipos1
  nout
  error('bad number of arguments and/or outputs')
end


GlobalVar;
idim = size(xcoor1,2); % Dimension of physical space

nbzone1 = length(modl1);
if ((nbzone1 ~= length(matr1)) || (nbzone1 ~= length(mail1)) || ...
    (nbzone1 ~= length(intg1)))
  nbzone1
  length(matr1)
  length(mail1)
  length(intg1)
  error('Bad number of sub-zones')
end

% Number of physical components in B matrix
% """""""""""""""""""""""""""""""""""""""""
if strcmp(mode1,liste_mode{1})
% TRID
  if (idim ~= 3)
    mode1
    idim
    error('Bad dimension for this mode')
  end
  nbcopB = 3;    % (UX,UY,UZ)

elseif (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || ...
        strcmp(mode1,'AXIS'))
% COPL Plane stress | DEPL Plane strain | AXIS Axisymmmetric
  if (idim ~= 2)
    mode1
    idim
    error('Bad dimension for this mode-2')
  end
  nbcopB = 2;    % (UX,UY) or (UR,UZ)

elseif strcmp(mode1,liste_mode{4})
% BARR Bar element (in local basis)
  if (idim == 2)
    nbcopB = 2;    % (U1,U2)
  elseif (idim == 3)
    nbcopB = 3;    % (U1,U2,U3)
  else
    mode1
    idim
    error('Bad dimension for this mode-1')
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

elseif strcmp(mode1,'POUT')
% POUT Euler-Bernoulli beam without curvature (in local basis)
  if (idim == 2)
    nbcopB = 3;    % (U1,U2,RZ)
  elseif (idim == 3)
    nbcopB = 6;    % (U1,U2,U3,R1,R2,R3)
  else
    mode1
    idim
    error('Bad dimension for this mode-3')
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

else
  mode1
  error('Bad mode')
end

clear mass1 bgamma1 bu1 d1;

% Loop on zones
% """""""""""""
for zo1 = 1:nbzone1
  intge1 = intg1{zo1};
  matre1 = matr1{zo1};
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};
  clear chamnoe1;

  nbed = length(modle1.NDDD);
  nbep = length(modle1.NDDP);
  nbel = size(maile1.MAIL,1);

  topo1  = maile1.MAIL;
  nbptg = length(intge1.WEIGHT);

  xval1 = zeros(nbed,nbep,nbel);
  if listOp(1); xval2 = zeros(nbcopB*nbptg,nbep,nbel);         end
  if listOp(2); xval3 = zeros(nbcopB*nbptg,nbep,nbel);         end
  if listOp(3); xval4 = zeros(nbcopB*nbptg,nbcopB*nbptg,nbel); end

% Elastic material coefficients
% """""""""""""""""""""""""""""
  [listComp1,listUnit1] = ListCompCham2(matr1,zo1);
  isotropic1 = 0;
  if strcmp(mode1,liste_mode{4})
%   BARR: adding characteristic field
    list_material_names = [{'RHO'} {'SECT'}];
%   and local geometry
    chamnoe1 = chamno1{zo1};
  elseif strcmp(mode1,'POUT')
%   POUT: adding characteristic field 
    switch idim
      case 2,
        list_material_names = [{'RHO'} {'SECT'} {'INRZ'}];
      case 3,
        list_material_names = [{'RHO'} {'SECT'} {'INR1'} {'INR2'} {'INR3'}];
      otherwise,
        idim
	error('Bad dimension')
    end
%   and local geometry
    chamnoe1 = chamno1{zo1};
  else
    list_material_names = [{'RHO'}];
  end

% Find right order of material components
  [junk,ia,ib] = intersect(list_material_names,listComp1);
% junk=list_material_names(ia)=listComp1(ib)
  [ia,junk] = sort(ia); ib=ib(junk);
  if length(ia) == length(list_material_names)
    isotropic1 = 1;
  else
    listComp1
    list_material_names
    error('Material non implemented')
  end

% Loop on elements
% """"""""""""""""
  for el1 = 1:nbel
    node1   = topo1(el1,:);
    xcoorel = xcoor1(node1,:);
    
    dd1 = sparse(0,0); 
    D1 = zeros(nbcopB,nbcopB,nbptg);

%   Loop on integration points
    for ptg = 1:nbptg
      clear D1e;
      if isotropic1
	    switch mode1
	      case 'BARR',
%           rho.S
            RHOS = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(2)}.XVAL(el1,ptg);
            D1e = RHOS * eye(idim,idim);
	      case 'POUT',
%           diag(rho.S,rho.I)
	        switch idim
	          case 2,
		        RHOS = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(2)}.XVAL(el1,ptg);
		        RHOI = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(3)}.XVAL(el1,ptg);
                D1e = [RHOS * eye(idim,idim), zeros(idim,1)
		               zeros(1,idim)        , RHOI];
	          case 3,
		        RHOS = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(2)}.XVAL(el1,ptg);
		        RHOI1 = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(3)}.XVAL(el1,ptg);
		        RHOI2 = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(4)}.XVAL(el1,ptg);
		        RHOI3 = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(5)}.XVAL(el1,ptg);
		    	D1e = [RHOS * eye(idim,idim), zeros(idim,idim)
		               zeros(1,idim)        , RHOI1 , 0     , 0
		               zeros(1,idim)        , 0     , RHOI2 , 0
		               zeros(1,idim)        , 0     , 0     , RHOI3];
%	            error('3D case to be done...')
              otherwise,
                error('Bad dimension')
            end
	      otherwise,
            D1e = matre1{ib(1)}.XVAL(el1,ptg) * eye(idim,idim);
	    end
      else
        error('Anisotropic material not yet done...')
      end
      D1(:,:,ptg) = D1e;
      aa = size(dd1,1);
      bb = size(D1e,1);
      dd1(aa+1:aa+bb,aa+1:aa+bb) = D1e;
      clear D1e;
    end

    if (strcmp(mode1,'BARR') || strcmp(mode1,'POUT'))
      switch idim
        case 2,
          T = [chamnoe1{1}.XVAL(el1,:)' chamnoe1{2}.XVAL(el1,:)'];
        case 3,
          T = [chamnoe1{1}.XVAL(el1,:)' chamnoe1{2}.XVAL(el1,:)' ...
               chamnoe1{3}.XVAL(el1,:)'];
        otherwise,
          error('bad idim')
      end
      xT = sum(T.^2,2).^-0.5; % normalize
      T = T .* repmat(xT,1,idim); clear xT;
      [ME,BE,bbE] = ElementMass7(modle1,D1,intge1,xcoorel,mode1,T);
    else
      [ME,BE,bbE] = ElementMass7(modle1,D1,intge1,xcoorel,mode1);
    end

    xval1(:,:,el1) = ME;
    if listOp(1); xval2(:,:,el1) = BE;  end
    if listOp(2); xval3(:,:,el1) = bbE; end
    if listOp(3); xval4(:,:,el1) = dd1; end
    clear D1 ME BE bbE;
  end
  mass1{zo1} = struct('XVAL',xval1);
  if listOp(1); bgamma1{zo1} = struct('XVAL',xval2); end
  if listOp(2); bu1{zo1}     = struct('XVAL',xval3); end
  if listOp(3); d1{zo1}      = struct('XVAL',xval4); end
  clear xval1 xval2 xval3 xval4;
end

if listOp(1); varargout(listOp(1)) = {bgamma1}; end
if listOp(2); varargout(listOp(2)) = {bu1};     end
if listOp(3); varargout(listOp(3)) = {d1};      end
