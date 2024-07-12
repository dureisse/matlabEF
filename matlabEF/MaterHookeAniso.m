function [hooke1] = MaterHookeAniso(coef1,local1,modl1,intg1, ...
                                    idim,mode2,varargin)
% Build element-field of Hooke matrix for anisotropic models
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2004
%
% For massive elements
%   hooke1 = MaterHookeAniso(coef1,local1,modl1,intg1,idim,mode2)
% For structural elements (beam, plate...)
%   hooke1 = MaterHookeAniso(coef1,local1,modl1,intg1,idim,mode2, ...
%                            chamno1,mail1,intgno1,mode1)
%
% Cette routine construit la matride de Hooke relative au modele.
% Elle est identique a celle relative qu materiau pour les elements
% finis massifs ; pour les elements de structure, elle est modifiee
% (les composantes "deformation"-"contraintes" sont celles du modele)
%
% Inputs
%   coef1	Champ par element des coefficients materiau
%   local1	Champ par element des caracteristiques geometriques
%               (repere local matiere, epaisseur des plaques...)
%   modl1	Modele sous-tendant coef1 ET local1
%   intg1	Segment d'integration de coef1 ET local1
%   idim	Dimension de l'espace (2 ou 3)
%   mode2	Type de symetries
%		QUEL : quelconque
%               ORTH : orthotrope
%               ORTE : orthotrope equilibre
%               ISTR : orthotrope de revolution ou isotrope transverse
%               ISOT : isotrope
% Optional inputs
%   chamno1	Champ par element aux noeuds du repere local geometrique
%   mail1	Son maillage
%   intgno1	Segment d'integration correspondant
%   mode1	Mode d'analyse (DKIR...)
% Output
%   hooke1	Champ par element de matrices de Hooke du modele
%               (pas du materiau : elle est adaptee pour les elements
%                de structure, ou elle est definie dans la base locale
%                geometrique)
%
% Pour le repere local materiau
%   local1{zo1}{1}.COMP = 'Repere'
%   local1{zo1}{1}.XVAL(nbel,nbpgau,nbvec,idim)
%   en 2D, il y a 1 vecteur
%   en 3D, il y a 2 vecteurs
% Eventuellement pour les plaques : epaisseur (materiau) et 
% repere local element
%   local1{zo1}{1}.COMP = 'EPAI'
%   local1{zo1}{1}.XVAL(nbel,nbpgau)
%   chamno1{zo1}{1}.COMP = 'Repere'
%   chamno1{zo1}{1}.XVAL(nbel,nbno,nbvec,idim)
%   en poutre, il y a 1 vecteur (tangent)
%   en plaque, il y a 3 vecteurs (le 3e est la normale)
%
% Pour la matrice de Hooke
%   hooke1{zo1}{1}.COMP = 'HookeMatrix'
%   hooke1{zo1}{1}.XVAL(nbel,nbpgau,nbc,nbc)

% Arguments optionnels
narg = nargin - 6;
switch narg
  case 0,
    clear chamno1 mail1 intgno1 mode1;
  case 4,
    chamno1 = varargin{1};
    mail1   = varargin{2};
    intgno1 = varargin{3};
    mode1   = varargin{4};
%   On passe le repere geometrique local au points d'integration
%   du materiau
    geom1 = ChamnoToCham(chamno1,mail1,intgno1,modl1,intg1);
  otherwise,
    narg
    error('Bad bumber of optional arguments')
end

% Nombre de zones
nbzo1 = length(local1);

% Boucle sur les zones
% """"""""""""""""""""
clear hooke1;
for zo1 = 1 : nbzo1
  locale1 = local1{zo1}{1};
  intge1  = intg1{zo1};
  modle1  = modl1{zo1};
  if exist('mode1')
    geome1 = geom1{zo1}{1};
  end

% Nombre d'elements dans la zone
  nbel = size(locale1.XVAL,1);
% Nombre de points d'integration dans la zone
  nbpgau = length(intge1.WEIGHT);
% Nombre de composantes primales/duales dans le modele
  nbc = length(modle1.COMP);
% Liste des noms de coefficients materiau dans la zone
  [listComp1,listUnit1] = ListCompCham2(coef1,zo1); 
% Liste des noms des parametres geometriques dans la zone
  [listComp2,listUnit2] = ListCompCham2(local1,zo1); 

  xval = zeros(nbel,nbpgau,nbc,nbc);

% On recherche les coefficients materiau necessaires dans la zone
% '''''''''''''''''''''''''''''''''''''' 
  switch idim
    case 2,
%
%     Cas 2D
%     """"""
      error('Cas 2D pas encore implante, desole...')

    case 3,
%
%     Cas 3D
%     """"""
      switch mode2
        case 'QUEL',
%         quelconque
          Lcoef = [{'D11'} {'D12'} {'D13'} {'D14'} {'D15'} {'D16'} ...
                   {'D21'} {'D22'} {'D23'} {'D24'} {'D25'} {'D26'} ...
                   {'D31'} {'D32'} {'D33'} {'D34'} {'D35'} {'D36'} ...
                   {'D41'} {'D42'} {'D43'} {'D44'} {'D45'} {'D46'} ...
                   {'D51'} {'D52'} {'D53'} {'D54'} {'D55'} {'D56'} ...
                   {'D61'} {'D62'} {'D63'} {'D64'} {'D65'} {'D66'}];
          Icoef = findoccur(Lcoef,listComp1);
          D11 = coef1{zo1}{Icoef(1)}.XVAL;
          D12 = coef1{zo1}{Icoef(2)}.XVAL;
          D13 = coef1{zo1}{Icoef(3)}.XVAL;
          D14 = coef1{zo1}{Icoef(4)}.XVAL;
          D15 = coef1{zo1}{Icoef(5)}.XVAL;
          D16 = coef1{zo1}{Icoef(6)}.XVAL;
          D21 = coef1{zo1}{Icoef(7)}.XVAL;
          D22 = coef1{zo1}{Icoef(8)}.XVAL;
          D23 = coef1{zo1}{Icoef(9)}.XVAL;
          D24 = coef1{zo1}{Icoef(10)}.XVAL;
          D25 = coef1{zo1}{Icoef(11)}.XVAL;
          D26 = coef1{zo1}{Icoef(12)}.XVAL;
          D31 = coef1{zo1}{Icoef(13)}.XVAL;
          D32 = coef1{zo1}{Icoef(14)}.XVAL;
          D33 = coef1{zo1}{Icoef(15)}.XVAL;
          D34 = coef1{zo1}{Icoef(16)}.XVAL;
          D35 = coef1{zo1}{Icoef(17)}.XVAL;
          D36 = coef1{zo1}{Icoef(18)}.XVAL;
          D41 = coef1{zo1}{Icoef(19)}.XVAL;
          D42 = coef1{zo1}{Icoef(20)}.XVAL;
          D43 = coef1{zo1}{Icoef(21)}.XVAL;
          D44 = coef1{zo1}{Icoef(22)}.XVAL;
          D45 = coef1{zo1}{Icoef(23)}.XVAL;
          D46 = coef1{zo1}{Icoef(24)}.XVAL;
          D51 = coef1{zo1}{Icoef(25)}.XVAL;
          D52 = coef1{zo1}{Icoef(26)}.XVAL;
          D53 = coef1{zo1}{Icoef(27)}.XVAL;
          D54 = coef1{zo1}{Icoef(28)}.XVAL;
          D55 = coef1{zo1}{Icoef(29)}.XVAL;
          D56 = coef1{zo1}{Icoef(30)}.XVAL;
          D61 = coef1{zo1}{Icoef(31)}.XVAL;
          D62 = coef1{zo1}{Icoef(32)}.XVAL;
          D63 = coef1{zo1}{Icoef(33)}.XVAL;
          D64 = coef1{zo1}{Icoef(34)}.XVAL;
          D65 = coef1{zo1}{Icoef(35)}.XVAL;
          D66 = coef1{zo1}{Icoef(36)}.XVAL;
        case 'ORTH',
%         orthotrope
          Lcoef = [{'YG1'} {'YG2'} {'YG3'} ...
                   {'NU23'} {'NU31'} {'NU12'} ...
	           {'G23'} {'G31'} {'G12'}];
          Icoef = findoccur(Lcoef,listComp1);
          YG1 = coef1{zo1}{Icoef(1)}.XVAL;
          YG2 = coef1{zo1}{Icoef(2)}.XVAL;
          YG3 = coef1{zo1}{Icoef(3)}.XVAL;
          NU23 = coef1{zo1}{Icoef(4)}.XVAL;
          NU31 = coef1{zo1}{Icoef(5)}.XVAL;
          NU12 = coef1{zo1}{Icoef(6)}.XVAL;
          G23 = coef1{zo1}{Icoef(7)}.XVAL;
          G31 = coef1{zo1}{Icoef(8)}.XVAL;
          G12 = coef1{zo1}{Icoef(9)}.XVAL;
        case 'ISTR',
%         orthotrope de revolution ou isotrope transverse,
%         comme cas particulier de l'orthotrope
          Lcoef = [{'YG1'} {'YG2'} {'NU23'} {'NU12'} {'G12'}];
	  Icoef = findoccur(Lcoef,listComp1);
	  mode2 = 'ORTH';
	  YG1 = coef1{zo1}{Icoef(1)}.XVAL;
	  YG2 = coef1{zo1}{Icoef(2)}.XVAL;
	  YG3 = coef1{zo1}{Icoef(2)}.XVAL;
	  NU23 = coef1{zo1}{Icoef(3)}.XVAL;
	  NU31 = coef1{zo1}{Icoef(4)}.XVAL;
	  NU12 = coef1{zo1}{Icoef(4)}.XVAL;
	  G23 = 0.5 * (YG2 ./ (1 + NU23));
	  G31 = coef1{zo1}{Icoef(5)}.XVAL;
	  G12 = coef1{zo1}{Icoef(5)}.XVAL;
        case 'ORTE',
%         orthotrope equilibre, comme cas particulier de l'orthotrope
          Lcoef = [{'YOUN'} {'NU'} {'G'}];
          Icoef = findoccur(Lcoef,listComp1);
          mode2 = 'ORTH';
	  YG1 = coef1{zo1}{Icoef(1)}.XVAL;
	  YG2 = coef1{zo1}{Icoef(1)}.XVAL;
	  YG3 = coef1{zo1}{Icoef(1)}.XVAL;
	  NU23 = coef1{zo1}{Icoef(2)}.XVAL;
	  NU31 = coef1{zo1}{Icoef(2)}.XVAL;
	  NU12 = coef1{zo1}{Icoef(2)}.XVAL;
	  G23 = coef1{zo1}{Icoef(3)}.XVAL;
	  G31 = coef1{zo1}{Icoef(3)}.XVAL;
	  G12 = coef1{zo1}{Icoef(3)}.XVAL;
        case 'ISOT',
%         isotrope
	  Lcoef = [{'YOUN'} {'NU'}];
          Icoef = findoccur(Lcoef,listComp1);
	  YOUN = coef1{zo1}{Icoef(1)}.XVAL;
	  NU   = coef1{zo1}{Icoef(2)}.XVAL;
        otherwise,
	  mode2
	  error('symetrie materielle non reconnue')
      end
      switch mode1
        case 'DKIR',
	  Lcoef = [{'EPAI'}];
	  Icoef = findoccur(Lcoef,listComp2);
	  EPAI = local1{zo1}{Icoef(1)}.XVAL;
	case 'TRID',
	otherwise,
	  mode1
	  error('mode d analyse non reconnu')
      end

    otherwise,
      idim
      error('Bad idim')
  end

%   Boucle sur les points d'integration
%   ''''''''''''''''''''''''''''''''''' 
    for el1 = 1:nbel

      for ptg1 = 1:nbpgau

%       Repere local materiau Q = [v1 v2 v3] eventuel
        if ~strcomp(mode2,'ISOT')
	  Q      = zeros(idim,idim); 
	  Q(:,1) = locale1.XVAL(el1,ptg1,1,:)';
	  switch idim
	    case 2,
	      error('Pas encore fait en 2D')
	    case 3,
	      Q(:,2) = locale1.XVAL(el1,ptg1,2,:)';
%             reconstruit
	      vec3 = ProdVect(Q(:,1:2));
	      Q(:,3) = vec3 / norm(vec3);
              vec2 = ProdVect(Q(:,[3 1]));
              Q(:,2) = vec2 / norm(vec2);
	      vec1 = Q(:,1);
	      Q(:,1) = vec1 / norm(vec1);
          end
%         Modification pour retourner au repere local geometrique
%         R = [v1 v2 v3] eventuel
	  if exist('mode1')
	    R      = zeros(idim,idim); 
	    R(:,1) = geome1.XVAL(el1,ptg1,1,:)';
keyboard
%Separer si POUT ou DKIR
	    switch idim
	      case 2,
	        error('Pas encore fait en 2D bis')
	      case 3,
	        R(:,2) = geom1.XVAL(el1,ptg1,2,:)';
%             reconstruit
	        vec3 = ProdVect(R(:,1:2));
	        R(:,3) = vec3 / norm(vec3);
                vec2 = ProdVect(R(:,[3 1]));
                R(:,2) = vec2 / norm(vec2);
	        vec1 = R(:,1);
	        R(:,1) = vec1 / norm(vec1);
            end
	    Q = R' * Q;
	  end
        end

	switch mode2
          case 'QUEL',
%           Matrice de Hooke en repere local materiau
	    HookeLocale1 = zeros(6,6);
	    HookeLocale1(1,1) = D11(el1,ptg1);
	    HookeLocale1(1,2) = D12(el1,ptg1);
	    HookeLocale1(1,3) = D13(el1,ptg1);
	    HookeLocale1(1,4) = D14(el1,ptg1);
	    HookeLocale1(1,5) = D15(el1,ptg1);
	    HookeLocale1(1,6) = D16(el1,ptg1);
	    HookeLocale1(2,1) = D21(el1,ptg1);
	    HookeLocale1(2,2) = D22(el1,ptg1);
	    HookeLocale1(2,3) = D23(el1,ptg1);
	    HookeLocale1(2,4) = D24(el1,ptg1);
	    HookeLocale1(2,5) = D25(el1,ptg1);
	    HookeLocale1(2,6) = D26(el1,ptg1);
	    HookeLocale1(3,1) = D31(el1,ptg1);
	    HookeLocale1(3,2) = D32(el1,ptg1);
	    HookeLocale1(3,3) = D33(el1,ptg1);
	    HookeLocale1(3,4) = D34(el1,ptg1);
	    HookeLocale1(3,5) = D35(el1,ptg1);
	    HookeLocale1(3,6) = D36(el1,ptg1);
	    HookeLocale1(4,1) = D41(el1,ptg1);
	    HookeLocale1(4,2) = D42(el1,ptg1);
	    HookeLocale1(4,3) = D43(el1,ptg1);
	    HookeLocale1(4,4) = D44(el1,ptg1);
	    HookeLocale1(4,5) = D45(el1,ptg1);
	    HookeLocale1(4,6) = D46(el1,ptg1);
	    HookeLocale1(5,1) = D51(el1,ptg1);
	    HookeLocale1(5,2) = D52(el1,ptg1);
	    HookeLocale1(5,3) = D53(el1,ptg1);
	    HookeLocale1(5,4) = D54(el1,ptg1);
	    HookeLocale1(5,5) = D55(el1,ptg1);
	    HookeLocale1(5,6) = D56(el1,ptg1);
	    HookeLocale1(6,1) = D61(el1,ptg1);
	    HookeLocale1(6,2) = D62(el1,ptg1);
	    HookeLocale1(6,3) = D63(el1,ptg1);
	    HookeLocale1(6,4) = D64(el1,ptg1);
	    HookeLocale1(6,5) = D65(el1,ptg1);
	    HookeLocale1(6,6) = D66(el1,ptg1);
%           On change de repere : repere global pour les elements
%           massifs, repere local geometrique pour les elements
%           de structure
	    HookeGlobale1 = HookeLocalToGlobal(HookeLocale1,Q);
	    clear HookeLocale1;
          case 'ORTH',
%           Matrice de Hooke en repere local materiau
	    HookeLocale1 = HookeOrthotropic3D(...
	       YG1(el1,ptg1) ,YG2(el1,ptg1) ,YG3(el1,ptg1), ...
	       NU23(el1,ptg1),NU31(el1,ptg1),NU12(el1,ptg1), ...
	       G23(el1,ptg1) ,G31(el1,ptg1) ,G12(el1,ptg1));
%           On change de repere : repere global pour les elements
%           massifs, repere local geometrique pour les elements
%           de structure
	    HookeGlobale1 = HookeLocalToGlobal(HookeLocale1,Q);
	    clear HookeLocale1;
	  case 'ISOT',
	    HookeGlobale1 = HookeIsotropic3D(YOUN(el1,ptg1), ...
	                                     NU(el1,ptg1));
	end

%       Modification de la matrice de Hooke materielle pour les
%       elements de structure, afin d'avoir la matrice de Hooke
%       du modele
        switch mode1
          case 'DKIR',
	    HookeGlobale1 = HookePlateAniso(HookeGlobale1, ...
	                                    EPAI(el1,ptg1),mode1);
        end

        xval(el1,ptg1,:,:) = HookeGlobale1;
	clear HookeGlobale1 Q;
      end
    end

  hookel1{1} = struct('COMP','HookeMatrix','UNIT','','XVAL',xval);
  hooke1{zo1} = hookel1;
  clear hookel1 locale1 xval;
  clear intge1;
end
