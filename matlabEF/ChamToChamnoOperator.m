function [T1] = ChamToChamnoOperator(mail1,intg1,numerptg,numerptg2,xcoor1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 04 / 2003
%
% Construit l'operateur discretise qui permet de faire passer d'un
% vecteur representant un champ par elements a une composante,
% au vecteur representant un champ par elements defini aux noeuds
% par extrapolation.
%
% Entrees
%  mail1	Maillage sous tendant le premier champ par elements
%  intg1	Les points de Gauss associes
%  numerptg(Nptg1)	Numerotation locale des points support
%  numerptg2(Nptg2)	Numerotation locale des points support cibles
%                       (i.e. les noeuds des elements)
%  xcoor1(nbnn,idim)	Pile de noeuds
%
% Sorties
%  T1(Nptg2,Nptg1)	Matrice de passage
%
% (voir ChamToChamno)

% Number of target points
Nptg2 = length(numerptg2);
% Number of source points
Nptg1 = length(numerptg);
% Inverse renumbering
numerptg_inv = InverseList(numerptg,max(numerptg));
numerptg2_inv = InverseList(numerptg2,max(numerptg2));

T1 = sparse(Nptg2,Nptg1);

% On a besoin des matrices de compressibilite elementaires sur mail1
% (integrales des produit des fct de base)
% Le segment d'integration devant etre intg1, on ecrase le segment
% usuel... en supposant que les fcts de base sont rangees dans le
% meme ordre...
idim = size(xcoor1,2);
if (idim == 2)
  mode1 = 'COPL';
else
  mode1 = 'TRID';
end
%% DD 08/05/2003 [modl2,intg2] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim); 
%% DD 08/05/2003 intg2 = intg1; % On ecrase...
[modl2,intg2] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim,intg1); 
mater2 = ManuChml(mail1,'MOB','',1.);
mater2 = ChmlToCham(mater2,mail1,intg2);
[m2,bg2] = Compress7(modl2,mater2,mail1,intg2,xcoor1,mode1,'GenDualOp');
clear bu2 mater2 modl2 intg2 mode1;

% Boucle sur les zones
ind1 = 0;
ind2 = 0;
nbzone1 = length(mail1);
for zo1 = 1:nbzone1

  intge1 = intg1{zo1};
  maile1 = mail1{zo1};
  me2  = m2{zo1};
  bge2 = bg2{zo1};
  [nbel,nbno] = size(maile1.MAIL);

  nptg = size(intge1.PHI,2);
  nbcomp = 1;

  M = me2.XVAL;
  B = bge2.XVAL;

  toto = eye(nbno,nbno) - (1./nbno)*ones(nbno,nbno) ;
% Boucle sur les elements
  for el1 = 1:nbel
%   Valeurs aux noeuds (cas singulier : on limite les oscillations)
%    G = M(:,:,el1) \ (B(:,:,el1)' * f);
%    G = MM * (B(:,:,el1)' * f) + RR * alpha;
    MM = pinv(M(:,:,el1));
    RR = null(M(:,:,el1));
    G0 = MM * B(:,:,el1)';
    alpha = - (RR' * toto * RR) \ (RR' * toto * G0);

    ptg1 = [ind1+1:ind1+nptg];
    ptg1 = numerptg_inv(ptg1);
    ptg2 = [ind2+1:ind2+nbno];
    ptg2 = numerptg2_inv(ptg2);

    T1(ptg2,ptg1) = G0 + RR * alpha;

    ind1 = ind1 + nptg;
    ind2 = ind2 + nbno;
    
    clear MM RR G0 alpha ptg1 ptg2;
  end
  clear toto;

  clear M B me2 bge2 maile1 intge1;
end

clear m2 bg2;
