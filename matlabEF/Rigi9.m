function [rigi1,varargout] = Rigi9(modl1,matr1,mail1,intg1, ...
                                   xcoor1,mode1,varargin)
% Elemental stiffness rigidities and associated operators
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 01 / 08 / 2002
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 05 / 08 / 2002
%   Ajout des B
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%   Ajout du mode BARR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 05 / 2003
%   Modification de syntaxe pour arguments et sorties sur demande
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   Modification composantes 2D ZZ
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 12 / 2003
%   Ajout du mode TIMO (poutre de Timoshenko a courbure negligee)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2003
%   Ajout du mode POUT (poutre d'Euler-Bernoulli a courbure negligee)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 04 / 2004
%   Ajout du champ de geometrie locale pour les elements de structure
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 04 / 2004
%   Passage au 3D pour les barres
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 07 / 2004
%   Ajout du mode DKIR (plaque Kirchhoff discret supposee plane)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 05 / 2005
%   Ajout du mode JOIN (joint d'epaisseur nulle)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 12 / 2005
%   Ajout de la rigidite de spin pour les plaques DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 17 / 01 / 2006
%   Passage au 3D pour les POUT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2006
%   Ajout de la reconstruction d'une rotation quadratique pour DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 20 / 07 / 2008
%   Ajout du materiau HookeMatrix en anisotrope complet
% DUREISSEIX David  LaMCoS                           le 08 / 10 / 2010
%   Ajout du mode DSHE (plaque avec cisaillement discret supposee plane)
%
% Entrees
%   modl1		modele
%   matr1		materiau
%   mail1		maillage
%   intg1		segment d'integration
%   xcoor1(nbno,idim)	coordonnees des noeuds
%   mode1		mode d'analyse
% Entrees optionnelles
%   PrimOp		mot cle pour l'operateur primal b1
%   GenDualOp		mot cle pour l'operateur dual generalise bsig1
%   ConstiOp		mot cle pour l'operateur de comportement d1
%   QuadRotOp		mot cle pour reconstruire les rotations quadratiques
%                       pour un element DKIR (RX,RY,RZ aux noeuds aretes)
%   chamno1		champ par elements aux noeuds contenant des
%                       informations sur la geometrie
%   'Spin'		mot cle pour ajouter une raideur artificielle
%                       dans le cas des plaques. Il faut fournir aussi
%                       un coefficient xrig1
%   xrig1		% de rigidite affecte au spin pour les plaques
%
% Sorties
%   rigi1		rigidites elementaires
%   et une liste d'operateurs ranges comme demandes dans la liste des
%   options
%
% Exemples
%   rigi1 = Rigi9(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [rigi1,bsig1] = Rigi9(modl1,matr1,mail1,intg1,xcoor1,mode1,'GenDualOp')
%   [rigi1,b1,bsig1] = Rigi9(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                            'PrimOp','GenDualOp')
%   [rigi1,d1] = Rigi9(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                      chamno1,'ConstiOp')
%   rigi1 = Rigi9(modl1,matr1,mail1,intg1,xcoor1,'DKIR','Spin',0.01)
%   [rigi1,cc1] = Rigi9(modl1,matr1,mail1,intg1,xcoor1,'DKIR', ...
%                       'Spin',0.01,'QuadRotOp')
%
% Commentaires
%   On peut ainsi avoir
%   l'operateur BSIGMA bsig1 dont le transpose donne un champ
%   d'efforts generalise a partir d'un champ de contraintes,
%   et l'operateur BEPSILON b1 qui donne un champ de deformation a partir
%   d'un champ de deplacement.
%   Le champ chamno1 peut etre utile pour les raccords d'elements de
%   structure : par defaut, les tangentes aux poutres (normales aux
%   plaques) sont construites independemment element par element, ce
%   qui pose probleme en cas de brisure due a la discretisation
%   (rigidification excessive) ; chamno1 les surcharge.

% Analyse optional input arguments,
% check order of operators and find if chamno1 is provided
nin1 = nargin-6;
nout = nargout-1;

clear chamno1;
listOp = zeros(1,4);
ipos1 = 0;
clear spin1;
in1 = 0;
while (in1 < nin1)
  in1 = in1 + 1;
  argin1 = varargin{in1};
  if iscell(argin1)
    chamno1 = argin1;
  else
    switch argin1
      case 'GenDualOp',
%       Operateur dual generalise de type BSIGMA
        ipos1 = ipos1 + 1;
        listOp(1) = ipos1;
      case 'PrimOp',
%       Operateur primal de type BEPSILON
        ipos1 = ipos1 + 1;
        listOp(2) = ipos1;
      case 'ConstiOp',
%       Operateur de comportement
        ipos1 = ipos1 + 1;
        listOp(3) = ipos1;
      case 'Spin',
%       Ajout d'une rigidite articielle de spin pour les plaques
        in1 = in1 + 1;
        spin1 = varargin{in1};
        disp(['  Warning in Rigi9: adding artificial spin '...
	      'stiffness, coefficient ' num2str(spin1)])
      case 'QuadRotOp',
        ipos1 = ipos1 + 1;
        listOp(4) = ipos1;
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

if strcmp(mode1,'DSHE')
  disp('  Rigi9: Warning DSHE for DST only')
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
if strcmp(mode1,'TRID')
% TRID Tridimensional
  if (idim ~= 3)
    mode1
    idim
    error('Bad dimension for this mode')
  end
  nbcopB = 6;    % EPXX,EPYY,EPZZ,GAYZ,GAZX,GAXY

elseif (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
% COPL Plane stress | DEPL Plane strain | AXIS Axisymmetric
  if (idim ~= 2)
    mode1
    idim
    error('Bad dimension for this mode-2')
  end
  nbcopB = 4;    % (EPXX,EPYY,GAXY,EPZZ) ou (EPRR,EPZZ,GARZ,EPTT)

elseif strcmp(mode1,liste_mode{4})
% BARR Bar element (in local basis)
  nbcopB = 1;    % EPS1
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

elseif strcmp(mode1,liste_mode{6})
% POUT Euler-Bernoulli beam without curvature (in local basis)
  if (idim == 2)
    nbcopB = 2;   % EPS1,C3
  elseif (idim == 3)
    nbcopB = 4;   % EPS1,C1,C2,C3
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

elseif strcmp(mode1,liste_mode{7})
% DKIR Discrete Kirchhoff plane plate (in local basis)
  if (idim == 3)
    nbcopB = 6;  % EP11,EP22,GA12,CP11,CP22,CG12
                 % EF11,EF22,EG12,MF11,MF22,MG12
  else
    idim
    error('pas prevu autre que 3D pour DKIR')
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

elseif strcmp(mode1,'DSHE')
% DSHE Discrete sshear plane plate (in local basis)
  if (idim == 3)
    nbcopB = 8;  % EP11,EP22,GA12,CP11,CP22,CG12,G1,G2
                 % EF11,EF22,EG12,MF11,MF22,MG12,T1,T2
  else
    idim
    error('pas prevu autre que 3D pour DSHE')
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

elseif strcmp(mode1,liste_mode{8})
% JOIN Joint d'epaisseur nulle (in local basis)
  if (idim == 2)
    nbcopB = 2;   % DRSN,DRN
  elseif (idim == 3)
    error('JOIN not yet in 3D')
  end
  if ~exist('chamno1')
    [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
  end

else
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

  if strcmp(mode1,liste_mode{4})
%   BARR: adding characteristic field
    list_material_names = [{'YOUN'} {'NU'} {'SECT'}];
%   and local geometry
    chamnoe1 = chamno1{zo1};

  elseif (strcmp(mode1,liste_mode{5}) || strcmp(mode1,liste_mode{6}))
%   TIMO or POUT: adding characteristic fields
    if (idim == 2)
      list_material_names = [{'YOUN'} {'NU'} {'SECT'} {'INRZ'}];
    elseif (idim == 3)
      disp('Warning Rigi9 TIMO POUT 3D: section circulaire uniquement')
%      list_material_names = [{'YOUN'} {'NU'} {'SECT'} ...
%                             {'INR1'} {'INR2'} ...
%			     {'TORS'}];
      list_material_names = [{'YOUN'} {'NU'} {'SECT'} ...
                             {'INR'} {'TORS'}];
    end
%   and local geometry
    chamnoe1 = chamno1{zo1};

  elseif strcmp(mode1,liste_mode{7})
%   DKIR
    list_material_names = [{'YOUN'} {'NU'} {'EPAI'}];
    list_material_names_aniso = [{'HookeMatrixPlate'}];
%   and local geometry
    chamnoe1 = chamno1{zo1};

  elseif strcmp(mode1,'DSHE')
%   DSHE
    list_material_names = [{'YOUN'} {'NU'} {'EPAI'}];
    list_material_names_aniso = [{'HookeMatrixPlate'}];
%   and local geometry
    chamnoe1 = chamno1{zo1};

  elseif strcmp(mode1,liste_mode{8})
%   JOIN
    list_material_names = [{'KS'} {'KN'}];
%   and local geometry
    chamnoe1 = chamno1{zo1};

  elseif strcmp(mode1,'TRID') || strcmp(mode1,'COPL') || ...
         strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS')
    list_material_names = [{'YOUN'} {'NU'}];
    list_material_names_aniso = [{'YG1'} {'YG2'} {'YG3'} ...
                                 {'NU12'} {'NU23'} {'NU13'} ...
                                 {'G12'} {'G23'} {'G13'} ...
                                 {'V1X'} {'V1Y'} {'V1Z'} ...
                                 {'V2X'} {'V2Y'} {'V2Z'}];
    list_material_names_aniso2 = [{'HookeMatrix'}];

  else
    mode1
    error('Mode not recognized')
  end

% Find right order of material components
  ia_aniso2 = [];
  ia_aniso = [];
  ia = [];
  if exist('list_material_names_aniso2')
    [junk,ia_aniso2,ib_aniso2] = intersect(list_material_names_aniso2, ...
                                           listComp1);
%   junk=list_material_names_aniso2(ia_aniso2)=listComp1(ib_aniso2)
    [ia_aniso2,junk] = sort(ia_aniso2); ib_aniso2=ib_aniso2(junk);
  end
  if exist('list_material_names_aniso')
    [junk,ia_aniso,ib_aniso] = intersect(list_material_names_aniso,listComp1);
%   junk=list_material_names_aniso(ia_aniso)=listComp1(ib_aniso)
    [ia_aniso,junk] = sort(ia_aniso); ib_aniso=ib_aniso(junk);
  end
  [junk,ia,ib] = intersect(list_material_names,listComp1);
% junk=list_material_names(ia)=listComp1(ib)
  [ia,junk] = sort(ia); ib=ib(junk);

  clear isotropic1;
  if length(ia) == length(list_material_names)
    isotropic1 = 1;
    disp('  Rigi9: Isotropic material detected')
  elseif exist('list_material_names_aniso')
    if length(ia_aniso) == length(list_material_names_aniso)
      isotropic1 = 0;
      disp('  Rigi9: Anisotropic material detected')
    elseif exist('list_material_names_aniso2')
      if length(ia_aniso2) == length(list_material_names_aniso2)
        isotropic1 = -1;
        disp('  Rigi9: Full anisotropic material detected')
      else
        listComp1
        list_material_names
        list_material_names_aniso
        error('Full anisotropic material coefficients not found')
      end
    else
      listComp1
      list_material_names
      list_material_names_aniso
      error('Anisotropic material coefficients not found')
    end
  else
    listComp1
    list_material_names
    list_material_names_aniso
    error('Isotropic material coefficients not found')
  end

% Loop on elements
% """"""""""""""""
  for el1 = 1:nbel
    node1 = topo1(el1,:);
    xcoorel = xcoor1(node1,:);

    dd1 = sparse(0,0);
    D1 = zeros(nbcopB,nbcopB,nbptg);
    for ptg = 1:nbptg

      if isotropic1==1

%       Isotropic case
%       """"""""""""""
        if strcmp(mode1,'BARR')
%         BARR: E.S
          D1e = matre1{ib(1)}.XVAL(el1,ptg) * matre1{ib(3)}.XVAL(el1,ptg);
        elseif (strcmp(mode1,'TIMO') || strcmp(mode1,'POUT'))
%         TIMO or POUT
          if idim==2
%           YOUN NU SECT INR
            D1e = HookeBeam(matre1{ib(1)}.XVAL(el1,ptg), ...
                            matre1{ib(2)}.XVAL(el1,ptg), ...
                            matre1{ib(3)}.XVAL(el1,ptg), ...
                            matre1{ib(4)}.XVAL(el1,ptg), ...
                            mode1);
          else
%           YOUN NU SECT INR TORS
            D1e = HookeBeam3D(matre1{ib(1)}.XVAL(el1,ptg), ...
                              matre1{ib(2)}.XVAL(el1,ptg), ...
                              matre1{ib(3)}.XVAL(el1,ptg), ...
                              matre1{ib(4)}.XVAL(el1,ptg), ...
                              matre1{ib(5)}.XVAL(el1,ptg), ...
                              mode1);
          end
        elseif strcmp(mode1,'DKIR')
%         DKIR: YOUN NU EPAI
          D1e = HookePlate(matre1{ib(1)}.XVAL(el1,ptg), ...
                           matre1{ib(2)}.XVAL(el1,ptg), ...
                           matre1{ib(3)}.XVAL(el1,ptg), ...
                           mode1);
        elseif strcmp(mode1,'DSHE')
%         DSHE: YOUN NU EPAI
          D1e = HookePlateS(matre1{ib(1)}.XVAL(el1,ptg), ...
                            matre1{ib(2)}.XVAL(el1,ptg), ...
                            matre1{ib(3)}.XVAL(el1,ptg), ...
                            mode1);
        elseif strcmp(mode1,'JOIN')
%         JOIN: KS KN
          D1e = [matre1{ib(1)}.XVAL(el1,ptg) 0.
	         0. matre1{ib(2)}.XVAL(el1,ptg)];
        elseif strcmp(mode1,'TRID') || strcmp(mode1,'COPL') || ...
               strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS')
%         TRID: YOUN NU
          D1e = HookeIsotropic3(matre1{ib(1)}.XVAL(el1,ptg), ...
                                matre1{ib(2)}.XVAL(el1,ptg), ...
                                mode1);
        else
          mode1
          error('Mode not implemented for isotropic case')
        end

      elseif isotropic1==0

%       Anisotropic case
%       """"""""""""""""
%       matre1 contient la matrice de Hooke dans la base geometrique
%       locale pour les elements de structure, dans la base globale
%       pour les elements massifs
        if strcmp(mode1,'DKIR')
          if ~isempty(ib_aniso)
%           HookeMatrixPlate
            junk = size(matre1{ib_aniso(1)}.XVAL);
            D1e = zeros(junk(3),junk(4));
            D1e(:,:) = matre1{ib_aniso(1)}.XVAL(el1,ptg,:,:);
          elseif ~isempty(ia_aniso)
            YG1 = matre1{ia_aniso(1)}.XVAL(el1,ptg);
            YG2 = matre1{ia_aniso(2)}.XVAL(el1,ptg);
            YG3 = matre1{ia_aniso(3)}.XVAL(el1,ptg);
            NU12 = matre1{ia_aniso(4)}.XVAL(el1,ptg);
            NU23 = matre1{ia_aniso(5)}.XVAL(el1,ptg);
            NU13 = matre1{ia_aniso(6)}.XVAL(el1,ptg);
            G12 = matre1{ia_aniso(7)}.XVAL(el1,ptg);
            G23 = matre1{ia_aniso(8)}.XVAL(el1,ptg);
            G13 = matre1{ia_aniso(9)}.XVAL(el1,ptg);
            V1X = matre1{ia_aniso(10)}.XVAL(el1,ptg);
            V1Y = matre1{ia_aniso(11)}.XVAL(el1,ptg);
            V1Z = matre1{ia_aniso(12)}.XVAL(el1,ptg);
            V2X = matre1{ia_aniso(13)}.XVAL(el1,ptg);
            V2Y = matre1{ia_aniso(14)}.XVAL(el1,ptg);
            V2Z = matre1{ia_aniso(15)}.XVAL(el1,ptg);
            V1 = [V1X V1Y V1Z ; V2X V2Y V2Z]';
            NU31 = NU13 * YG3/YG1;
            G31 = G13;
            D1e = HookOrthotropic3D(YG1,YG2,YG3,NU23,NU31,NU12, ...
                                    G23,G31,G12,mode1,V1);
            error('Appeler HookePlateAniso')
            D1e = epai1 * D1e;
          else
            error('PB avec coef mat')
          end
        else
          mode1
          error('Mode not implemented for anisotropic case')
        end

      elseif isotropic1==-1

%       Full anisotropic case
%       """""""""""""""""""""
        if strcmp(mode1,'TRID')
%         TRID: HookeMatrix
          junk = size(matre1{ib_aniso2(1)}.XVAL);
          D1e = zeros(junk(3),junk(4));
          D1e(:,:) = matre1{ib_aniso2(1)}.XVAL(el1,ptg,:,:);
        else
          mode1
          error('Mode not implemented for full anisotropic case')
        end
      end
      D1(:,:,ptg) = D1e;
      aa = size(dd1,1);
      bb = size(D1e,1);
      dd1(aa+1:aa+bb,aa+1:aa+bb) = D1e;
      clear D1e;
    end

    if (strcmp(mode1,'BARR') || strcmp(mode1,'TIMO') || ...
        strcmp(mode1,'POUT') || strcmp(mode1,'JOIN'))
%     Local basis for BARR, TIMO, POUT, JOIN
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
      [KE,BE,bbE] = ElementStiffness11(modle1,D1,intge1,xcoorel,mode1,T);

    elseif (strcmp(mode1,'DKIR') || strcmp(mode1,'DSHE'))
%     Local basis for DKIR, DSHE
      switch idim
        case 3,
          junk = [chamnoe1{1}.XVAL(el1,:)
                  chamnoe1{2}.XVAL(el1,:)
                  chamnoe1{3}.XVAL(el1,:)];
          N1 = junk(:,1:3:end);
          N2 = junk(:,2:3:end);
          N3 = junk(:,3:3:end);
%         Average to get a local basis constant per element
          N1 = sum(N1,2) / size(N1,2);
          N2 = sum(N2,2) / size(N2,2);
          N3 = sum(N3,2) / size(N3,2);
          N3 = N3 / norm(N3);
          N1 = ProdVect([N2 N3]); N1 = N1 / norm(N1);
          N2 = ProdVect([N3 N1]); N2 = N2 / norm(N2);
          T = [N1 N2 N3]';
          clear junk N1 N2 N3;
        otherwise,
          error('bad idim for DKIR or DSHE')
      end
      if strcmp(mode1,'DKIR')
        if exist('spin1')
          [KE,BE,bbE,CCE] = ElementStiffness11(modle1,D1,intge1,xcoorel, ...
                                               mode1,T,spin1);
        else
          [KE,BE,bbE,CCE] = ElementStiffness11(modle1,D1,intge1,xcoorel, ...
                                               mode1,T);
        end
      elseif strcmp(mode1,'DSHE')
        if exist('spin1')
          [KE,BE,bbE] = ElementStiffnessDST(modle1,D1,intge1,xcoorel, ...
                                           mode1,T,spin1);
        else
          [KE,BE,bbE] = ElementStiffnessDST(modle1,D1,intge1,xcoorel, ...
                                           mode1,T);
        end
      end

    else
%     Massive case: global basis is used
      [KE,BE,bbE] = ElementStiffness11(modle1,D1,intge1,xcoorel,mode1);
    end

    if el1 == 1
      if listOp(4); xval5 = zeros(size(CCE,1),nbep,nbel); end
    end
    xval1(:,:,el1) = KE;
    if listOp(1); xval2(:,:,el1) = BE;  end
    if listOp(2); xval3(:,:,el1) = bbE; end
    if listOp(3); xval4(:,:,el1) = dd1; end
    if listOp(4); xval5(:,:,el1) = CCE; end
    clear D1 KE BE bbE dd1 CCE;
  end
  rigi1{zo1} = struct('XVAL',xval1);
  if listOp(1); bsig1{zo1} = struct('XVAL',xval2); end
  if listOp(2); b1{zo1}    = struct('XVAL',xval3); end
  if listOp(3); d1{zo1}    = struct('XVAL',xval4); end
  if listOp(4); cc1{zo1}   = struct('XVAL',xval5); end
  clear xval1 xval2 xval3 xval4 xval5;
end

if listOp(1); varargout(listOp(1)) = {bsig1}; end
if listOp(2); varargout(listOp(2)) = {b1};    end
if listOp(3); varargout(listOp(3)) = {d1};    end
if listOp(4); varargout(listOp(4)) = {cc1};   end
