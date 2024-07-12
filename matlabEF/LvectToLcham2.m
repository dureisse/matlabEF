function [Lcham1] = LvectToLcham2(LS1,mail1,intg1, ...
                                  numer1,mapComp1,ListComp1,varargin)
% Disassemble a set of vectors into a list of element-based fields
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 18 / 07 / 2008
%
% Entrees
%   LS1(n1,n)			Liste de n vecteurs
%   mail1			Mesh of the field
%   intg1			Integration points of the field
%   numer1(nbptg1)		Ordering of integration points used
%   mapComp1(nbptg1,nbcomp1)	Mapping matrix
%   ListComp1{nbcomp1}		List of component names
% Entrees optionnelles
%   ListUnit1{nbcomp1}		Liste des unites
% Sorties
%   Lcham1{n}			Array of element-based fields
%
%   Tous les champs par element d'une liste ont la meme structure.

% Treatment of optional arguments
nin1 = nargin-6;
if nin1 == 0
  ListUnit1 = repmat([{''}],length(ListComp1));
elseif nin1 == 1
  ListUnit1 = varargin{1};
else
  nin1
  error('Wrong number of optional arguments')
end

nbval = 1; % On fait l'hypothese que les composantes sont scalaires !

% Sizes and consistency
[n1,n] = size(LS1);
[nbptg1,nbcomp1] = size(mapComp1);
if (nbcomp1 ~= length(ListComp1)) | (nbcomp1 ~= length(ListUnit1))
  nbcomp1
  length(ListDdl1)
  length(ListUnit1)
  error('no consistent lengths')
end
if nbptg1 ~= length(numer1)
  nbptg1
  length(numer1)
  error('no consistent lengths bis')
end

numerptg1 = [1:nbptg1];

% Boucle sur les champs
% """""""""""""""""""""
clear Lcham1;
for ic = 1:n
  clear cham1;
  S1 = LS1(:,ic);

% Boucle sur les sous-zones du champ par element
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  intge1  = intg1{zo1};
  maile1  = mail1{zo1};
  nbpt1   = length(intge1.WEIGHT); clear intge1;
  [nbel1,nbno1] = size(maile1.MAIL); clear maile1;

  clear chamel1;

% Boucle sur les composantes du champ par element
  for i = 1:nbcomp1
    name = ListComp1(i);
    xval1 = zeros(nbel1,nbval*nbpt1);

    j = findoccur(numerptg1,numer1);
    k = find(j); j = j(k);
    xval1(:,:) = reshape(S1(mapComp1(j,i),:),nbpt1,nbel1)';

    chamel1{i} = struct('COMP',name,'UNIT','','XVAL',xval1); clear xval1;
  end

  cham1{zo1} = chamel1; clear chamel1;
  clear S1;
end

  Lcham1{ic} = cham1;
end
