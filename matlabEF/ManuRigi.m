function [rigi1] = ManuRigi(mail1,Kel1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2002

% Construit la rigidite rigi1 manuellement,
% a partir d'un maillage mail1 a une seule zone
% (pour avoir le nombre d'elements)
% et d'une matrice elementaire Kel1
% (identique pour tous les elements)

nbel1 = Nelements2(mail1);
nbzone1 = length(mail1);
if (nbzone1 ~= 1)
  nbzone1
  error('Only 1 zone expected')
end

clear rigi1;
for el1 = 1:nbel1
  xval1(:,:,el1) = Kel1;
end
rigi1{1} = struct('XVAL',xval1);

clear xval1;
