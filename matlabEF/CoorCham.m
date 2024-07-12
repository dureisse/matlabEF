function [cham2] = CoorCham(mail1,intg1,xcoor1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 07 / 04 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 05 / 2008
%   Using transformation functions only (this requires the model)
%
% Retourne le champ par elements (cham2,mail1,intg1) des coordonnees
% des points d'integration d'un champ par elements defini sur
% (mail1,intg1)
%
% Entrees
%   mail1		Son maillage support
%   intg1		Son support d'integration
%   xcoor1(nbno,idim)	Pile des noeuds (coordonnees)
% Optional entries
%   modl1		Son modele
%     (si non fourni : on suppose que toutes les fonctions de forme
%      sont utilisees pour la transformation geometrique)
% Sortie
%   cham2		Champ par elements des coordonnees

GlobalVar;
clear cham2;

narg = nargin - 3;
switch narg
  case 0,
    disp('  CoorCham: Warning: all basis functions are assumed to transform')
    clear modl1;
  case 1,
    modl1 = varargin{1};
  otherwise,
    error(['Bad number of optional inputs ' int2str(narg)])
end

nbzone1 = length(mail1);
if (nbzone1 ~= length(intg1))
  nbzone1
  length(intg1)
  disp('Bad number of zones')
  error('stopped in CoorCham')
end

idim = size(xcoor1,2);

% Loop on zones
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  intge1 = intg1{zo1};
  topo1  = maile1.MAIL;
  [nbel,nbno] = size(topo1);
  nbptg  = length(intge1.WEIGHT);
  if exist('modl1');
    modle1 = modl1{zo1};
    nnot1 = modle1.NNOT; % numeros des noeuds de l'element pour la transf
    nnit1 = modle1.NNIT; % numeros des fct de forme de l'element pour la transf
  else
    nnot1 = [1:nbno];
    nnit1 = [1:size(intge1.PHI,1)];
  end

  xval2 = zeros(nbel,nbptg,idim);
% Loop on elements
  for el1 = 1:nbel
    lno1 = topo1(el1,:);
    lno1 = lno1(nnot1);
    lphi1 = intge1.PHI(nnit1,:);
%%DD    xcorel1 = xcoor1(topo1(el1,:)',:);
%%DD    xval2(el1,:,:) = intge1.PHI' * xcorel1;
    xcorel1 = xcoor1(lno1',:);
    xval2(el1,:,:) = lphi1' * xcorel1;
    clear xcorel1;
  end

  clear chame2;
% Loop on components
  for i = 1:idim
    chame2{i} = struct('COMP',liste_ddlp(i),'UNIT','', ...
                       'XVAL',xval2(:,:,i));
  end

  cham2{zo1} = chame2;
  clear xval2 maile1 intge1;
end
