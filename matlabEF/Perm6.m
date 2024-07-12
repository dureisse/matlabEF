function [perm1,varargout] = Perm6(modl1,matr1,mail1,intg1, ...
                                   xcoor1,mode1,varargin)
% Permeability elementary matrices and associated operators
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2003
%   arguments optionnels
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 08 / 2006
%   ajout mode TRID
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 09 / 2007
%   ajout mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 10 / 2007
%   ajout anisotropie
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 07 / 2008
%   Ajout mode BARR (unidimensionnel)
%
% Matrices de permeabilite elementaires perm1
% Et operateur dual generalise bw1 associe
% Et operateur primal associe bz1
% Et operateur de comportement h1 associe
%
% Entrees
%   modl1               modele
%   matr1               materiau
%   mail1               maillage
%   intg1               segment d'integration
%   xcoor1(nbno,idim)   coordonnees des noeuds
%   mode1               mode d'analyse
%  et une liste d'options parmi
%   PrimOp      operateur primal b1
%   GenDualOp   operateur dual generalise bsig1
%   ConstiOp    operateur de comportement d1
%
% Sorties
%   perm1       rigidites elementaires
%  et une liste d'operateurs ranges comme demandes dans la liste des
%  options
%
% Exemples
%   perm1 = Perm6(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [perm1,bw1] = Perm6(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                       'GenDualOp')
%   [perm1,bz1,bw1] = Perm6(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                           'PrimOp','GenDualOp')
%   [perm1,bw1,bz1,h1] = Perm6(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                              'GenDualOp','PrimOp','ConstiOp')

% (voir Compress5,Rigi5)

nin1 = nargin-6;
nout = nargout-1;
if (nin1 ~= nout)
  nin1
  nout
  error('bad number of arguments and/or outputs')
else
  listRegularOp = [{'GenDualOp'} {'PrimOp'} {'ConstiOp'}];
  listOp = findoccur(listRegularOp,varargin);
  test = findoccur(varargin,listRegularOp);
  if ~isempty(find(test==0))
    varargin
    listRegularOp
    error('bad option')
  end
end

GlobalVar;

nbzone1 = length(modl1);
if (nbzone1 ~= length(matr1)) | (nbzone1 ~= length(mail1)) | ...
   (nbzone1 ~= length(intg1))
  nbzone1
  length(matr1)
  length(mail1)
  length(intg1)
  error('Bad number of sub-zones')
end

switch mode1
  case {'COPL','DEPL','AXIS'},
%   COPL Plane stress | DEPL Plane strain
    idim   = 2;   % Dimension de l'espace physique
    nbcopB = 2;   % Nombre de composantes physiques dans les matrices B
                  % (PX,PY) pour le gradient de pression, (WX,WY) pour le flux
                  % ou (PR,PZ), (WR,WZ)
  case 'TRID',
%   TRID tridimensionnel
    idim   = 3;   % Dimension de l'espace physique
    nbcopB = 3;   % Nombre de composantes physiques dans les matrices B
                  % (PX,PY,PZ) pour le gradient de pression,
                  % (WX,WY,WZ) pour le flux
  case 'BARR',
%   BARR unidimensionnel
    idim   = 1;   % Dimension de l'espace physique
    nbcopB = 1;   % Nombre de composantes physiques dans les matrices B
                  % PX pour le gradient de pression,
                  % WX pour le flux
  otherwise,
    mode1
    error('Bad mode')
end

clear perm1 bw1 bz1;

% Loop on zones
% """""""""""""
for zo1 = 1:nbzone1
  intge1 = intg1{zo1};
  matre1 = matr1{zo1};
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};

  nbed = length(modle1.NDDD);     % Number of dual dof
  nbep = length(modle1.NDDP);     % Number of primal dof
  nbel = size(maile1.MAIL,1);     % Number of elements

  topo1  = maile1.MAIL;
  nbptg = length(intge1.WEIGHT);  % Number of integration points

  xval1 = zeros(nbed,nbep,nbel);
  if listOp(1); xval2 = zeros(nbcopB*nbptg,nbep,nbel); end
  if listOp(2); xval3 = zeros(nbcopB*nbptg,nbep,nbel); end
  if listOp(3); xval4 = zeros(nbcopB*nbptg,nbcopB*nbptg,nbel); end

% Material coefficients
% """""""""""""""""""""
  list_material_names = [{'PERM'} {'VISC'}];
  switch mode1
    case {'COPL','DEPL'},
      list_material_names_aniso = [{'PEXX'} {'PEYY'} ...
                                   {'PEXY'} {'PEYX'} {'VISC'}];
    case {'TRID'},
      list_material_names_aniso = [{'PEXX'} {'PEYY'} {'PEZZ'} ...
                                   {'PEXY'} {'PEYX'} ...
                                   {'PEYZ'} {'PEZY'} ...
                                   {'PEZX'} {'PEXZ'} {'VISC'}];
    case {'AXIS'},
      list_material_names_aniso = [{'PERR'} {'PEZZ'} ...
                                   {'PERZ'} {'PEZR'} {'VISC'}];
  end
  [listComp1,listUnit1] = ListCompCham2(matr1,zo1);

% On cherche deja l'isotropie
  isotropic1 = 0;
  [junk,ia,ib] = intersect(list_material_names,listComp1);
% junk=list_material_names(ia)=listComp1(ib)
  [ia,junk] = sort(ia); ib=ib(junk);
  if length(ia) == length(list_material_names)
    isotropic1 = 1;
    disp('  Perm6: Isotropic material detected')
  elseif exist('list_material_names_aniso')
%   On cherche l'anisotropie
    [junk,ia_aniso,ib_aniso] = intersect(list_material_names_aniso,listComp1);
%   junk=list_material_names_aniso(ia_aniso)=listComp1(ib_aniso)
    [ia_aniso,junk] = sort(ia_aniso); ib_aniso=ib_aniso(junk);
    if length(ia_aniso) == length(list_material_names_aniso)
      disp('  Perm6: Anisotropic material detected')
    else
      listComp1
      list_material_names
      list_material_names_aniso
      error('Mandatory material coefficients not found')
    end
  else
    listComp1
    list_material_names
    error('Mandatory material coefficients not found')
  end

% Loop on elements
% """"""""""""""""
  for el1 = 1:nbel
    node1 = topo1(el1,:);
    xcoorel = xcoor1(node1,:);

    dd1 = sparse(0,0);
    D1 = zeros(nbcopB,nbcopB,nbptg);
    for ptg = 1:nbptg
      if isotropic1
        PERM = matre1{ib(1)}.XVAL(el1,ptg);
        VISC = matre1{ib(2)}.XVAL(el1,ptg);
        D1e = (PERM / VISC) * eye(nbcopB);
      else
        switch mode1
          case {'COPL','DEPL','AXIS'},
            PERM11 = matre1{ib_aniso(1)}.XVAL(el1,ptg);
            PERM22 = matre1{ib_aniso(2)}.XVAL(el1,ptg);
            PERM12 = matre1{ib_aniso(3)}.XVAL(el1,ptg);
            PERM21 = matre1{ib_aniso(4)}.XVAL(el1,ptg);
            VISC   = matre1{ib_aniso(5)}.XVAL(el1,ptg);
            D1e = (1. / VISC) * [PERM11 PERM12 ; PERM21 PERM22];
          case {'TRID'},
            PERMXX = matre1{ib_aniso(1)}.XVAL(el1,ptg);
            PERMYY = matre1{ib_aniso(2)}.XVAL(el1,ptg);
            PERMZZ = matre1{ib_aniso(3)}.XVAL(el1,ptg);
            PERMXY = matre1{ib_aniso(4)}.XVAL(el1,ptg);
            PERMYX = matre1{ib_aniso(5)}.XVAL(el1,ptg);
            PERMYZ = matre1{ib_aniso(6)}.XVAL(el1,ptg);
            PERMZY = matre1{ib_aniso(7)}.XVAL(el1,ptg);
            PERMZX = matre1{ib_aniso(8)}.XVAL(el1,ptg);
            PERMXZ = matre1{ib_aniso(9)}.XVAL(el1,ptg);
            VISC   = matre1{ib_aniso(10)}.XVAL(el1,ptg);
            D1e = (1. / VISC) * [PERMXX PERMXY PERMXZ
                                 PERMYX PERMYY PERMYZ
                                 PERMZX PERMZY PERMZZ];
        end
      end
      D1(:,:,ptg) = D1e;
      aa = size(dd1,1);
      bb = size(D1e,1);
      dd1(aa+1:aa+bb,aa+1:aa+bb) = D1e;
      clear D1e;
    end
    [KE,BE,bbE] = ElementPerm5(modle1,D1,intge1,xcoorel,mode1);
    xval1(:,:,el1) = KE;
    if listOp(1); xval2(:,:,el1) = BE; end
    if listOp(2); xval3(:,:,el1) = bbE; end
    if listOp(3); xval4(:,:,el1) = dd1; clear dd1;
    end
    clear D1 KE BE bbE;
  end
  perm1{zo1} = struct('XVAL',xval1);
  if listOp(1); bw1{zo1}   = struct('XVAL',xval2); end
  if listOp(2); bz1{zo1}   = struct('XVAL',xval3); end
  if listOp(3); d1{zo1}    = struct('XVAL',xval4); end
  clear xval1 xval2 xval3 xval4;
end

if listOp(1); varargout(listOp(1)) = {bw1}; end
if listOp(2); varargout(listOp(2)) = {bz1}; end
if listOp(3); varargout(listOp(3)) = {d1}; end
