function [rigi1,varargout] = Compress7(modl1,matr1,mail1,intg1, ...
                                       xcoor1,mode1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 05 / 2003
%   Modification de syntaxe pour arguments et sorties sur demande
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 12 / 2003
%   Possibilite de passer l'inverse du module de Biot
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 10 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 05 / 2008
%   Ajout du mode DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 07 / 2008
%   Ajout du mode BARR
%
% Matrices de compressibilites elementaires rigi1,
% et autres operateurs associes.
% (voir Rigi7, Mass7)
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
%   rigi1       compressibilites elementaires
%  et une liste d'operateurs ranges comme demandes dans la liste des
%  options

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

switch mode1
  case {'COPL','DEPL','AXIS'},
%   COPL Plane stress | DEPL Plane strain
    if (idim ~= 2)
      mode1
      idim
      error('Bad dimension for this mode')
    end
    nbcopB = 1;   % Nombre de composantes physiques dans les matrices B
  case 'TRID',
%   TRID Tridimensional
    if (idim ~= 3)
      mode1
      idim
      error('Bad dimension for this mode (bis)')
    end
    nbcopB = 1;   % Nombre de composantes physiques dans les matrices B
  case 'BARR',
%   BARR Unimensional
    nbcopB = 1;   % Nombre de composantes physiques dans les matrices B
  case 'DKIR',
%   DKIR
    if (idim == 3)
      nbcopB = 1;   % Nombre de composantes physiques dans les matrices B
    else
      idim
      error('pas prevu ate que 3D pour DKIR')
    end
  otherwise,
    mode1
    error('Bad mode (not yet available)')
end

clear rigi1 bsig1 b1 d1;

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
% IMOB: inverse of Biot modulus, or
% MOB: Biot modulus
  [listComp1,listUnit1] = ListCompCham2(matr1,zo1);
  isotropic1 = 0;
  [junk,ia,ib] = intersect([{'IMOB'}],listComp1);
  [ia,junk] = sort(ia); ib_IMOB=ib(junk);
  [junk,ia,ib] = intersect([{'MOB'}],listComp1);
  [ia,junk] = sort(ia); ib_MOB=ib(junk);
  if xor((length(ib_MOB) == 1),(length(ib_IMOB) == 1))
    isotropic1 = 1;
  else
    listComp1
    error('The correct material coefficient has not been found')
  end

% Loop on elements
% """"""""""""""""
  for el1 = 1:nbel
    node1 = topo1(el1,:);
    xcoorel = xcoor1(node1,:);

    dd1 = sparse(0,0);
    D1 = zeros(1,1,nbptg);

%   Loop on integration points
    for ptg = 1:nbptg
      if isotropic1
%       1./MOB
        if ib_IMOB
          D1e = matre1{ib_IMOB(1)}.XVAL(el1,ptg);
        else
          D1e = inv(matre1{ib_MOB(1)}.XVAL(el1,ptg));
        end
      else
        isotropic1
        error('Anisotropy not yet available')
      end
      D1(:,:,ptg) = D1e;
      aa = size(dd1,1);
      bb = size(D1e,1);
      dd1(aa+1:aa+bb,aa+1:aa+bb) = D1e;
      clear D1e;
    end
    if strcmp (mode1,'AXIS')
%     Real coordinates of integration points
%% METTRE CA DANS LA ROUTINE ELEMENTAIRE (CF ElementStiffness11)
      xcoorptg = intge1.PHI' * xcoorel;
      [KE,BE,bbE] = ElementCompress5(modle1,D1,intge1,xcoorel,mode1,xcoorptg);
    else
      [KE,BE,bbE] = ElementCompress5(modle1,D1,intge1,xcoorel,mode1);
    end
    xval1(:,:,el1) = KE;
    if listOp(1); xval2(:,:,el1) = BE; end
    if listOp(2); xval3(:,:,el1) = bbE; end
    if listOp(3); xval4(:,:,el1) = dd1; clear dd1; end
    clear D1 KE BE bbE;
  end
  rigi1{zo1} = struct('XVAL',xval1);
  if listOp(1); bsig1{zo1} = struct('XVAL',xval2); end
  if listOp(2); b1{zo1}    = struct('XVAL',xval3); end
  if listOp(3); d1{zo1}    = struct('XVAL',xval4); end
  clear xval1 xval2 xval3 xval4;
end

if listOp(1); varargout(listOp(1)) = {bsig1}; end
if listOp(2); varargout(listOp(2)) = {b1}; end
if listOp(3); varargout(listOp(3)) = {d1}; end
