function [ListMeshExt] = DDExt(ListContSdm,mail1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2002

% Avec la liste des maillages des contours des sous domaines
% ListContSdm(nsdm1), et le maillage d'un bord exterieur mail1,
% on decoupe mail1
% en la liste des maillages des exterieurs locaux aux sous domaines
% ListMeshExt(nsdm1).

nsdm1 = length(ListContSdm);

% Avec les maillages des contours, on cherche l'interface exterieure
for sdm1 = 1:nsdm1
  mext1 = IntersectMesh(ListContSdm{sdm1},mail1);
  ListMeshExt{sdm1} = mext1;
  clear mext1;
end
