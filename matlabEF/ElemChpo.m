function [nmail2,listelem1] = ElemChpo(chpo1,nmail1,val1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 08 / 2003

% Extrait le maillage du champ par point (chpo1,nmail1)
% sur lequel il vaut val1
% Retourne le maillage nmail2 et la correspondance d'elements dans
% nmail1: l'element i de nmail2 est l'element listelem1(i) de nmail1
% (numerotations globales des elements)

% List of global numbers of elements
listelem1 = zeros(1,0);
decal1 = 0; % for global numbering of elements over subzones
nbzo1 = length(chpo1);
for zo1 = 1:nbzo1
  chpoe1 = chpo1{zo1};
  nbcomp1 = length(chpoe1);
  if (nbcomp1 ~= 1)
    zo1
    nbcomp1
    error('More than 1 component in the field')
  end
  xval1 = chpoe1{nbcomp1}.XVAL;
  if (size(xval1,2) ~= 1)
    size(xval1)
    error('Complex field component not usable')
  end
  listelem1 = [listelem1 (find(xval1(:,1)' == val1) + decal1)];
  decal1 = decal1 + size(xval1,1);
end

% Element extraction
nmail2 = ElemMesh(nmail1,listelem1);
