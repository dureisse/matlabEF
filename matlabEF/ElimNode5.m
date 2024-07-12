function [pile_mail2,ListUsedNodes] = ElimNode5(pile_mail,xcoor,lxcrit);
% DUREISSEIX David   L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 22 / 01 / 2003
%   autant de valeurs de critere que de maillages dans la pile
%   extension au 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 07 / 04 / 2003
%   correction de bug
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 26 / 12 / 2003
%   correction de bug
%
% Elimination des noeuds a meme emplacement dans tous les maillages
% d'une pile de maillages, entre eux (on n'elimine pas a priori les
% noeuds en interne d'un meme maillage).
%
% Input
%   pile_mail{nbmail1}	pile de maillages a traiter
%   xcoor(nbno,idim)	coordonnees des noeuds
%   lxcrit(nbmail1)		criteres de proximite de noeuds
%
% Output
%  pile_mail2{nbmail1}	pile de maillages traites
%  ListUsedNodes(nbnou)	liste des noeuds effectivement utilises dans
%                       pile_mail2


% On boucle sur les maillages croissant en augmentant le maillage
% de reference (nmail3) au fur et a mesure
clear pile_mail2;
ListUsedNodes = [];
imail1 = 1;
  mail3 = pile_mail{imail1};
  pile_mail2{imail1} = mail3;
  nmail3 = ChangeMesh2(mail3,'POI1');
  numer3 = unique(nmail3{1}.MAIL');   %  numeros tries

nbmail1 = length(pile_mail);
for imail1 = 2:nbmail1
  mail1 = pile_mail{imail1};
  xcrit = lxcrit(imail1);
  [mail2,ListUsedNod] = ElimNodeMails5(xcoor,xcrit,mail1,nmail3);

  pile_mail2{imail1} = mail2;
  nmail2 = ChangeMesh2(mail2,'POI1');
  numer3 = unique([numer3 nmail2{1}.MAIL']);
  nmail3{1}.MAIL = numer3';

  clear nmail2 ListUsedNod mail1;
end
ListUsedNodes = numer3';
