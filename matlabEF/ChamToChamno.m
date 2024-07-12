function [chamno2,intg2] = ChamToChamno(cham1,mail1,intg1,xcoor1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 30 / 03 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 05 / 2008
%   On traite le cas DKIR
%
% MAUVAIS NOM : ChamToChamno
%
% Transforme un champ par elements (cham1,mail1,intg1) en un champ par
% elements defini aux noeuds (chamno2,mail1,intg2)
% par extrapolation.
%
% Entrees
%  cham1	Champ par elements a extrapoler
%  mail1	Son maillage sous tendant
%  intg1	Les points de Gauss sous tendant
%  xcoor1(nbnn,idim)	Pile de noeuds
%
% Sorties
%  chamno2	Champ par elements aux noeuds
%  intg2	Son segment d'integration

idim = size(xcoor1,2);
clear chamno2 intg2;

% On a besoin des matrices de compressibilite elementaires sur mail1
% (integrales des produit des fct de base)
% Le segment d'integration devant etre intg1, on ecrase le segment
% usuel... en supposant que les fcts de base sont rangees dans le
% meme ordre...
if (idim == 2)
  mode1 = 'COPL';
else
  zo1 = 1;
  idimr = size(intg1{zo1}.DPHI,1);
  if (idimr == 3)
    mode1 = 'TRID';
  elseif (idimr == 2)
    mode1 = 'DKIR';
  elseif (idimr == 1)
    mode1 = 'BARR';
  else
    error('Bad idimr')
  end
end
%% DD 08/05/2003 [modl2,intg2] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim); 
%% DD 08/05/2003 intg2 = intg1; % On ecrase...
[modl2,intg2] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim,intg1); 
mater2 = ManuChml(mail1,'MOB','',1.);
mater2 = ChmlToCham(mater2,mail1,intg2);
[m2,bg2] = Compress7(modl2,mater2,mail1,intg2,xcoor1,mode1,'GenDualOp');
clear mater2 modl2 intg2 mode1;

% Boucle sur les zones
nbzone1 = length(mail1);
for zo1 = 1:nbzone1

  chame1 = cham1{zo1};
  intge1 = intg1{zo1};
  maile1 = mail1{zo1};
  me2  = m2{zo1};
  bge2 = bg2{zo1};
  [nbel,nbno] = size(maile1.MAIL);
  nptg = size(intge1.PHI,2);
  nbcomp = length(chame1);

  M = me2.XVAL;
  B = bge2.XVAL;

  xval1 = zeros(nptg,nbcomp,nbel);
  xval2 = zeros(nbno,nbcomp,nbel);

% Boucle sur les composantes
  for i = 1:nbcomp
    xval = chame1{i}.XVAL;
    [nbel1,nbvalt] = size(xval);
    if (nbvalt ~= nptg)
      nbvalt
      nptg
      error('nbvalt different de nptg pas encore prevu')
    end
    if (nbel1 ~= nbel)
      zo1
      i
      nbel1
      nbel
      error('bad number of elements')
    end
    xval1(:,i,:) = xval';

    clear xval;
  end

  toto = eye(nbno,nbno) - (1./nbno)*ones(nbno,nbno) ;
% Boucle sur les elements
  for el1 = 1:nbel
%   Valeurs aux ptg
    f = xval1(:,:,el1);
%   Valeurs aux noeuds (cas singulier : on limite les oscillations)
%    G = M(:,:,el1) \ (B(:,:,el1)' * f);
%    G = MM * (B(:,:,el1)' * f) + RR * alpha;
    MM = pinv(M(:,:,el1));
    RR = null(M(:,:,el1));
    G0 = MM * (B(:,:,el1)' * f);
    alpha = - (RR' * toto * RR) \ (RR' * toto * G0);
    xval2(:,:,el1) = G0 + RR * alpha;

    clear f G;
  end
  clear toto;

% Boucle sur les composantes
  clear chamnoe2;
  for i = 1:nbcomp
    toto = zeros(nbno,nbel);
    toto(:,:) = xval2(:,i,:);
    chamnoe2{i} = struct('COMP',chame1{i}.COMP,'UNIT',chame1{i}.UNIT, ...
                         'XVAL',toto');
    clear toto;
  end
  chamno2{zo1} = chamnoe2;

  clear chamnoe2 xval2 xval1 M B;
  clear chame1 intge1 maile1 me2 bge2;
end

clear m2 bg2;

% On cree maintenant intg2 aux noeuds
% Pseudo integration segment associated to mail1
intg2 = SegmentIntgNo(mail1);

%clear intg2;
%nbzone1 = length(mail1);
%for zo1 = 1:nbzone1
%  maile1 = mail1{zo1};
%  [nbel,nbno] = size(maile1.MAIL);
%  type1 = maile1.TYPE;
%  Xcorel1 = EF_CoorRefNod(type1);
%  if strcmp(type1,'TRI3')
%    w = 0.5 / nbno * ones(1,nbno);
%  elseif strcmp(type1,'QUA4')
%    w = 4. / nbno * ones(1,nbno);
%  else
%    type1
%    error('Type of element not yet implemented')
%  end
%  intge1 = struct('PHI',eye(nbno),'COOR',Xcorel1,'WEIGHT',w);
%  intg2{zo1} = intge1;
%  clear intge1 Xcorel1 maile1;
%end
