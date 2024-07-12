function [mail2,listelem1] = ElemChml(chml1,mail1,val1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2002

% Extrait le maillage du champ par element constant (chml1,mail1)
% sur lequel il vaut val1
% Retourne le maillage mail2 et la correspondance d'elements dans
% mail1: l'element i de mail2 est l'element listelem1(i) de mail1
% (numerotations globales des elements)

nbcomp1 = length(chml1);
if (nbcomp1 ~= 1)
  nbcomp1
  error('More than 1 component in the field')
end

xval1 = chml1{1}.XVAL;
[nbelt,nbval] = size(xval1);
if (nbval ~= 1)
  nbval
  error('More than 1 value for the component')
end

% List of global numbers of elements
listelem1 = find(xval1'==val1);

% Element extraction
mail2 = ElemMesh(mail1,listelem1);
