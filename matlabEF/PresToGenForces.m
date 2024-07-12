function [FD1] = PresToGenForces(cfext1,nfext1,mail1, ...
                                 numerd1,mapddlDual1,listDdlDual1, ...
                                 numerp1,mapddlPrim1,listDdlPrim1, ...
                                 xcoort1,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 08 / 2003
%
% Passe d'un champ par point de pressions (de composantes listDdlPrim1)
% a un vecteur de forces generalisees vect1
% (entre les deux, il y a la matrice des fonctions de forme croisees)
%
% Entrees
%   cfext1		Champ par point de pression (interpolee)
%   nfext1		Maillage POI1 qui le sous-tend
%   mail1		Maillage du bord sur lequel calculer les forces
%   numerd1		Numerotation des noeuds duaux
%   mapddlDual1		Matrice d'assemblage des ddl duaux
%   listDdlDual1	liste des noms de ddl duaux
%   numerp1		Numerotation des noeuds primaux
%   mapddlPrim1		Matrice d'assemblage des ddl primaux
%   listDdlPrim1	liste des noms de ddl primaux
%   xcoort1(nbno,idim)	Coordonnees des noeuds
%   mode1		Mode d'analyse
% Sorties
%   FD1(nddl,1)		Vecteur de forces generalisees
%

  idim = size(xcoort1,2);

% Masse unitaire sur le bord
  [modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode1,idim);
  matr1 = ManuChml(mail1,'RHO','',1.);
  matr1 = ChmlToCham(matr1,mail1,intg1);
%  [mass1,bgamma1,bu1] = Mass6(modl1,matr1,mail1,intg1,xcoort1,mode1);
%%DD 05/08/13
  [mass1,bgamma1,bu1] = Mass7(modl1,matr1,mail1,intg1,xcoort1,mode1, ...
                              'GenDualOp','PrimOp');

% On passe en matriciel
  M1 = RigiToMatrix(mass1,modl1,mail1, ...
                    numerd1,mapddlDual1,listDdlDual1, ...
                    numerp1,mapddlPrim1,listDdlPrim1);
  F1 = ChpoToVect3(cfext1,nfext1, ...
                   numerd1,mapddlDual1,listDdlDual1);
%%DD 04/01/27 FAUTE                   numerp1,mapddlPrim1,listDdlPrim1);
  FD1 = M1 * F1;

  clear modl1 intg1 matr1 mass1 bgamma1 bu1 M1 F1;
