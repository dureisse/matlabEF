function [ListMeshSdm,ListContSdm,varargout] = DD2(chml1,mail1,idim,varargin)
% Construit les sous-domaines a partir d'un champ les reperants
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 02 / 2003
%   Cas des SEG2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 12 / 2004
%   Ajout du 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 05 / 2005
%   Limitation possible du nombre de sous domaines
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 07 / 2006
%   Possibilite d'avoir la correspondance d'elements
%
% Entrees
%   chml1		Champ constant par element reperant le sous-domaine
%   mail1		Maillage associe
%   idim		Dimension de l'espace
% Entree optionnelle
%   nsdm1		Restriction sur le nb de sous domaines
% Sorties
%   ListMeshSdm{sdm2}	Liste des maillages des sous domaines
%   ListContSdm{sdm2}	Liste des maillages des contours des sous domaines
% Sortie optionnelle
%   ListElemSdm{sdm2}   Liste des correspondances d'elements

% Avec le champ par element (chml1,mail1) reperant l'appartenance
% des elements a un sous domaine (composante de nom SDM),
% on decoupe le maillage mail1, pour obtenir
% la liste des maillages des sous domaines ListMeshSdm(nsdm1)
% et en prime la liste des maillages de leur contour ListContSdm(nsdm1)
%
% L'element i de ListMeshSdm{sdm2} est l'element ListElemSdm{sdm2}(i)
% de mail1.

nout = nargout - 2;

chmls1 = ExtrChml2(chml1,[{'SDM'}]);
nsdm1 = MaxiChml(chmls1);

nin1 = nargin-3;
if nin1 == 0
  disp(['Number of subdomains ' int2str(nsdm1)])
elseif nin1 == 1
  disp(['Total number of subdomains ' int2str(nsdm1)])
  nsdm1 = varargin{1};
  disp(['Modified number of subdomains ' int2str(nsdm1)])
else
  nin1
  error('Wrong number of optional arguments')
end

clear ListElemSdm;
for sdm1 = 1:nsdm1
  [msdm1,listelem1] = ElemChml(chmls1,mail1,sdm1);
  csdm1 = ContourMesh3(msdm1,idim);
  ListMeshSdm{sdm1} = msdm1;
  ListContSdm{sdm1} = csdm1;
  ListElemSdm{sdm1} = listelem1;
  clear msdm1 csdm1;
end

clear chmls1;

switch nout
  case 0,
  case 1,
    varargout(1) = {ListElemSdm};
  otherwise,
    nout
    error('Bad number of output variables')
end
