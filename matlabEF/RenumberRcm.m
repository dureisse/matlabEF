function [numer1] = RenumberRcm(mail1,xcoor1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 18 / 10 / 2003
% NERON David       L.M.T.   STRUCTURES ET SYSTEMES  le 18 / 10 / 2003
%
% Donne la renumerotation de Cuthill-McKee inverse d'un maillage
%
% Inputs
%   mail1		Maillage
%   xcoor1(nbnot1,idim)	Pile de noeuds
%                       vecteur resultat
% Outputs
%   numer1(nbno1)	renumerotation des noeuds
%

nmail1 = ChangeMesh2(mail1,'POI1');
numer1bof = nmail1{1}.MAIL';

idim = size(xcoor1,2);
switch idim
  case 2,
    mode1 = 'DEPL';
  case 3,
    mode1 = 'TRID';
end
[modl1,intg1] = ModlIntg13(mail1,'THERMIQUE','CAPACITE',mode1,idim);
%% DD 03/11/25 matr1 = ManuChml(mail1,[{'K'}],[{''}],1.);
matr1 = ManuChml(mail1,[{'C'}],[{''}],1.,[{'RHO'}],[{''}],1.);
mater1 = ChmlToCham(matr1,mail1,intg1);
%%DD 06/08/06 [rigi1,bsig1,b1] = Capa2(modl1,mater1,mail1,intg1,xcoor1,mode1, ...
%%DD 06/08/06                          'GenDualOp','ConstiOp');
rigi1 = Capa2(modl1,mater1,mail1,intg1,xcoor1,mode1);
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1bof,listDdlDual1, ...
                                        numer1bof,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1bof,mapddlDual1,listDdlDual1, ...
                  numer1bof,mapddlPrim1,listDdlPrim1);
%%DD 04/01/27 FAUTE  numer1 = symrcm(K1);
p = symrcm(K1);
numer1 = numer1bof(p);

clear K1 mapddlPrim1 mapddlDual1 numer1bof rigi1 modl1 bsig1 b1 nmail1;
