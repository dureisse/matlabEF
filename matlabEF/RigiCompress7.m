function [rigi1,Bt1,bt1] = RigiCompress7(modl1,matr1,mail1,...
                                         intg1,xcoor1,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 11 / 2002
%   modification des modeles version 1.3 passe a version 1.4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 12 / 2003
%   Possibilite de passer l'inverse du module de Biot
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 10 / 11 / 2007
%   Ajout mode TRID
%
% Matrices elementaires du systeme couple rigidite-compressibilite
% Et matrices B elementaire associees
% (voir Rigi5, Compress5, RigiCompress5)

disp('Warning: RigiCompress7 OBSOLETE, utiliser RigiCompress8')
 [rigi1,Bt1,bt1] = RigiCompress8(modl1,matr1,mail1,intg1, ...
                                 xcoor1,mode1,'GenDualOp','PrimOp');

