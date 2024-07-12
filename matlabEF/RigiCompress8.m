function [rigi1,varargout] = RigiCompress8(modl1,matr1,mail1,intg1,...
                                           xcoor1,mode1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 11 / 2002
%   modification des modeles version 1.3 passe a version 1.4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 12 / 2003
%   Possibilite de passer l'inverse du module de Biot
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 10 / 11 / 2007
%   Ajout mode TRID
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 18 / 08 / 2009
%   Modification de syntaxe pour arguments et sorties sur demande
%
% Matrices elementaires du systeme couple rigidite-compressibilite
% Et matrices B elementaire associees
% (voir Rigi9, Compress5, RigiCompress8)
% Ici l'operateur primal correspond a la deformation, puis la pression
% le comportement a Hooke et a l'oppose de la compressibilite
% Attention : l'oppose de la compressibilite est utilise pour
% compatibilite avec rigi1

% Analyse optional input arguments,
% check order of operators 
nin1 = nargin-6;
nout = nargout-1;

listOp = zeros(1,4);
ipos1 = 0;
in1 = 0;
while (in1 < nin1)
  in1 = in1 + 1;
  argin1 = varargin{in1};
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
      otherwise,
        argin1
        error('unknown option')
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

if strcmp(mode1,'TRID')
% TRID tridimensional
  if (idim ~= 3)
    mode1
    idim
    error('Bad dimension for this mode')
  end
  nbcopBs = 6;   % Number of physical components for solid
%                  'EPXX'  'EPYY'  'EPZZ'  'GAYZ'  'GAZX'  'GAXY'   or
%                  'SMXX'  'SMYY'  'SMZZ'  'TAYZ'  'TAZX'  'TAXY'
  nbcopBf = 1;   % Number of physical components for fluid (P)
  nbcopBc = 1;   % Number of physical components for trace of solid
elseif (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
% COPL Plane stress | DEPL Plane strain | AXIS Axisymmetric
  if (idim ~= 2)
    mode1
    idim
    error('Bad dimension for this mode-2')
  end
  nbcopBs = 4;   % Number of physical components for solid
%                  (EPXX,EPYY,GAXY,EPZZ) or (EPRR,EPZZ,GARZ,EPTT)
  nbcopBf = 1;   % Number of physical components for fluid (P)
  nbcopBc = 1;   % Number of physical components for trace of solid
else
  mode1
  error('Bad mode')
end
nbcopB = nbcopBs + nbcopBf;

clear rigi1 Bt1 bt1 d1;

% Separation of models for solid and fluid
modls1 = TireModlIntg8(modl1,intg1,'ELASTIQUE','ELASTIQUE',mode1);
modlf1 = TireModlIntg8(modl1,intg1,'FLUIDE','FLUIDE',mode1);
modlc1 = TireModlIntg8(modl1,intg1,'FLUIDE','ELASTIQUE',mode1);

% Loop on zones
for zo1 = 1:nbzone1
  intge1 = intg1{zo1};
  matre1 = matr1{zo1};
  modle1 = modl1{zo1};
  modlse1 = modls1{zo1};
  modlfe1 = modlf1{zo1};
  modlce1 = modlc1{zo1};
  maile1 = mail1{zo1};

  nbeds = length(modlse1.NDDD);
  nbeps = length(modlse1.NDDP);
  nbedf = length(modlfe1.NDDD);
  nbepf = length(modlfe1.NDDP);
  nbed = nbeds + nbedf;
  nbep = nbeps + nbepf;
  nbel = size(maile1.MAIL,1);

  nbcomds = length(modlse1.NCOD);
  nbcomps = length(modlse1.NCOP);
  nbcomdf = length(modlfe1.NCOD);
  nbcompf = length(modlfe1.NCOP);
  nbcomd  = length(modle1.NCOD);
  nbcomp  = length(modle1.NCOP);

  topo1 = maile1.MAIL;
  nbptg = length(intge1.WEIGHT);

  xval1 = zeros(nbed,nbep,nbel);
  if listOp(1); xval2 = zeros(nbcopB*nbptg,nbed,nbel);         end
  if listOp(2); xval3 = zeros(nbcopB*nbptg,nbep,nbel);         end
  if listOp(3); xval4 = zeros(nbcopB*nbptg,nbcopB*nbptg,nbel); end

%  KtE = zeros(nbed,nbep);
%  BtE = zeros(nbcopB*nbptg,nbed);
%  btE = zeros(nbcopB*nbptg,nbep);

% Material coefficient
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

% Material coefficient
% COB: Biot coefficient
  [listComp1,listUnit1] = ListCompCham2(matr1,zo1);
  isotropic1 = 0;
  [junk,ia,ib] = intersect(liste_rigicompr_isotropic,listComp1);
% junk=liste_rigicompr_isotropic(ia)=listComp1(ib)
  [ia,junk] = sort(ia); ib_COB=ib(junk);
  if length(ia) == length(liste_rigicompr_isotropic)
    isotropic1 = 1;
  else
    listComp1
    error('Material non implemented')
  end

% Material coefficient (YOUN,NU)
  [junk,ia,ib] = intersect(liste_rigi_isotropic,listComp1);
% junk=liste_rigi_isotropic(ia)=listComp1(ib)
  [ia,junk] = sort(ia); ib_YN=ib(junk);
  if length(ia) == length(liste_rigi_isotropic)
  else
    listComp1
    error('Material non implemented2')
  end

% For correct assembly with the order in modle1
% (a little bit hard to follow)
% One must have both the same node and the same name of dof
% for dof names
%  [listDdlDual1,listDdlPrimal1] = ListDdlModl2(modl1,zo1);
%  [listDdlDuals1,listDdlPrimals1] = ListDdlModl2(modls1,zo1);
%  [listDdlDualf1,listDdlPrimalf1] = ListDdlModl2(modlf1,zo1);
%  ddlps1 = findoccur(listDdlPrimals1,listDdlPrimal1);
%  ddlpf1 = findoccur(listDdlPrimalf1,listDdlPrimal1);
%  ddlds1 = findoccur(listDdlDuals1,listDdlDual1);
%  ddldf1 = findoccur(listDdlDualf1,listDdlDual1);
  ddlps1 = findoccur(modlse1.DDLP,modle1.DDLP);
  ddlpf1 = findoccur(modlfe1.DDLP,modle1.DDLP);
  ddlds1 = findoccur(modlse1.DDLD,modle1.DDLD);
  ddldf1 = findoccur(modlfe1.DDLD,modle1.DDLD);
% for both names and nodes
  full_ddlp = [modle1.NNOP' modle1.NDDP']; % pairs primal (name,node)
  full_ddld = [modle1.NNOD' modle1.NDDD']; % pairs dual (name,node)
  [junk,nddps1,nddps0] = intersect([modlse1.NNOP' ddlps1(modlse1.NDDP)'], ...
                                   full_ddlp,'rows');
  [junk,nddpf1,nddpf0] = intersect([modlfe1.NNOP' ddlpf1(modlfe1.NDDP)'], ...
                                   full_ddlp,'rows');  
  [junk,nddds1,nddds0] = intersect([modlse1.NNOD' ddlps1(modlse1.NDDD)'], ...
                                   full_ddld,'rows');
  [junk,ndddf1,ndddf0] = intersect([modlfe1.NNOD' ddlpf1(modlfe1.NDDD)'], ...
                                   full_ddld,'rows');
  clear full_ddlp full_ddld ddlps1 ddlpf1 ddlds1 ddldf1;

% For components (much simpler, cause all are used; but repeated nbptg times)
%  [listCompDual1,listCompPrimal1] = ListCompModl2(modl1,zo1);
%  [listCompDuals1,listCompPrimals1] = ListCompModl2(modls1,zo1);
%  [listCompDualf1,listCompPrimalf1] = ListCompModl2(modlf1,zo1);
%  comps1 = findoccur(listCompPrimals1,listCompPrimal1);
%  compf1 = findoccur(listCompPrimalf1,listCompPrimal1);
%  comds1 = findoccur(listCompDuals1,listCompDual1);
%  comdf1 = findoccur(listCompDualf1,listCompDual1);
  comps1 = findoccur(modlse1.COMP(modlse1.NCOP),modle1.COMP(modle1.NCOP));
  compf1 = findoccur(modlfe1.COMP(modlfe1.NCOP),modle1.COMP(modle1.NCOP));
  comds1 = findoccur(modlse1.COMD(modlse1.NCOD),modle1.COMD(modle1.NCOD));
  comdf1 = findoccur(modlfe1.COMD(modlfe1.NCOD),modle1.COMD(modle1.NCOD));
  comp1  = [1:length(modle1.NCOP)];
  comd1  = [1:length(modle1.NCOD)];
  comps2 = zeros(1,0); compf2 = zeros(1,0);
  comds2 = zeros(1,0); comdf2 = zeros(1,0);
  for ptg = 1:nbptg
    comps2 = [comps2 comps1+(ptg-1)*nbcomp];
    compf2 = [compf2 compf1+(ptg-1)*nbcomp];
    comds2 = [comds2 comds1+(ptg-1)*nbcomd];
    comdf2 = [comdf2 comdf1+(ptg-1)*nbcomd];
  end
  clear comps1 compf1 comds1 comdf1;

% Loop on elements
% """"""""""""""""
  for el1 = 1:nbel
    node1 = topo1(el1,:);
    xcoorel = xcoor1(node1,:);
    dd1 = sparse(0,0);
    Ds1 = zeros(nbcopBs,nbcopBs,nbptg);
    Df1 = zeros(nbcopBf,nbcopBf,nbptg);
    Dc1 = zeros(nbcopBf,nbcopBc,nbptg);
    TR1 = zeros(nbcopBf,nbcopBs,nbptg);
    for ptg = 1:nbptg
      if isotropic1
%	Hooke matrix
        Ds1e = HookeIsotropic3(matre1{ib_YN(1)}.XVAL(el1,ptg), ...
                               matre1{ib_YN(2)}.XVAL(el1,ptg), ...
                               mode1);
      end
      Ds1(:,:,ptg) = Ds1e;
      aa = size(dd1,1);
      bb = size(Ds1e,1);
      dd1(aa+1:aa+bb,aa+1:aa+bb) = Ds1e;
      clear Ds1e;

      if isotropic1
%       1./MOB
        Df1e = inv(matre1{ib_MOB(1)}.XVAL(el1,ptg));
      end
      Df1(:,:,ptg) = Df1e;
      aa = size(dd1,1);
      bb = size(Df1e,1);
%     Attention : oppose de 1./MOB
      dd1(aa+1:aa+bb,aa+1:aa+bb) = -1. * Df1e;
      clear Df1e;

      if isotropic1
%	COB
        Dc1e = matre1{ib_COB(1)}.XVAL(el1,ptg);
      end
      Dc1(:,:,ptg) = Dc1e;
      clear Dc1e;

      if strcmp(mode1,liste_mode{2})
%% DD 25/11/03        TR1e = TraceIsotropic(matre1{ib_YN(1)}.XVAL(el1,ptg), ...
%% DD 25/11/03                              matre1{ib_YN(2)}.XVAL(el1,ptg),mode1);
        TR1e = TraceIsotropic2(matre1{ib_YN(1)}.XVAL(el1,ptg), ...
                               matre1{ib_YN(2)}.XVAL(el1,ptg),mode1);
      else
%% DD 25/11/03        TR1e = TraceIsotropic([],[],mode1);
        TR1e = TraceIsotropic2([],[],mode1);
      end
      TR1(:,:,ptg) = TR1e;
      clear TR1e;
    end

%   Create with the order of dof in modls1 and modlf1
    if strcmp(mode1,'AXIS')
      intges1 = intge1;
      intges1.PHI = intge1.PHI(1:size(xcoorel,1),:);
% Real coordinates of integration points
disp('Warning: RigiCompress7 verifier transformation')
      xcoorptg = intge1.PHI(1:size(xcoorel,1),:)' * xcoorel;
      [KE,BsigmE,bepsiE] = ElementStiffness11(modlse1,Ds1,intges1,xcoorel,mode1);
      [SE,BqE,bpE] = ElementCompress5(modlfe1,Df1,intge1,xcoorel,mode1,xcoorptg);
      [NE,BPE,beE] = ElementStiffCompressCoupl7(modlce1,Dc1,TR1,intge1,...
                                              xcoorel,mode1,xcoorptg);
    else
      [KE,BsigmE,bepsiE] = ElementStiffness11(modlse1,Ds1,intge1,xcoorel,mode1);
      [SE,BqE,bpE] = ElementCompress5(modlfe1,Df1,intge1,xcoorel,mode1);
      [NE,BPE,beE] = ElementStiffCompressCoupl7(modlce1,Dc1,TR1,intge1,...
                                                xcoorel,mode1);
    end

%   Assemble with the order in modle1
    KtE([nddds0;ndddf0],[nddps0;nddpf0]) = ...
      [KE(nddds1,nddps1) -NE(nddpf1,nddds1)'
      -NE(ndddf1,nddps1) -SE(ndddf1,nddpf1)];
    BtE(comds2,nddds0) = BsigmE(:,nddds1);
    BtE(comdf2,ndddf0) = BqE(:,ndddf1);
    btE(comps2,nddps0) = bepsiE(:,nddps1);
    btE(compf2,nddpf0) = bpE(:,nddpf1);

    xval1(:,:,el1) = KtE;
    if listOp(1); xval2(:,:,el1) = BtE;  end
    if listOp(2); xval3(:,:,el1) = btE; end
    if listOp(3); xval4(:,:,el1) = dd1; end
    clear D1 TR1 KE BsigmE bepsiE SE BqE bpE NE BPE beE;
  end
  rigi1{zo1} = struct('XVAL',xval1);
  if listOp(1); Bt1{zo1} = struct('XVAL',xval2); end
  if listOp(2); bt1{zo1} = struct('XVAL',xval3); end
  if listOp(3); d1{zo1}  = struct('XVAL',xval4); end
  clear xval1 xval2 xval3 xval4;
end

if listOp(1); varargout(listOp(1)) = {Bt1}; end
if listOp(2); varargout(listOp(2)) = {bt1};    end
if listOp(3); varargout(listOp(3)) = {d1};    end
