function [rigi1,BP1,be1] = RigiCompressCoupl7(modlc1,matr1,mail1,...
                                              intg1,xcoor1,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 11 / 2002
%  ajout des matrices BP1,be1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 11 / 2002
%  passage d'un modele couple, et rigidites elem hors diagonale


% Matrices elementaires de couplage rigidite-compressibilite
% (voir Rigi5, Compress5)
% Et matrices B elementaire associees
% BP1' permet de passer d'un champ par element de pression p,
% aux forces generalisees P qui l'equilibrent en volume :
% U' * P = \int p \Trace \epsilon(u)
% avec P = BP1' * p (une fois BP1 mise sous forme matricielle)
% be1 permet de passer d'un champ par point de deplacement U au
% champ par element de variation de volume, e = \Trace \epsilon(u) :
% e = be1 * U (une fois be1 mise sous forme matricielle)

% Attention, ces dernieres sont peut etre les transposees suivant
% le modele couple modlc1

disp('RigiCompressCoupl7: Warning OBSOLETE, utiliser RigiCompressCoupl8')
[rigi1,BP1,be1] = RigiCompressCoupl8(modlc1,matr1,mail1,...
                                     intg1,xcoor1,mode1,'GenDualOp','PrimOp');
