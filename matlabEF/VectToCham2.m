function [cham1] = VectToCham2(S1,mail1,intg1,numer1,mapComp1,listComp1)
% Disassemble a vector into an element-based field
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 01 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 04 / 2004
%  Correction de bug dans le rangement des valeurs de cham1
%
% cham1 = VectToCham2(S1,mail1,intg1,numer1,mapComp1,listComp1);
%
% Built the element-based field (not constant per element a priori)
% (cham1,mail1,intg1)
% from a given vector of values S1, according to the mapping
% (numer1,mapComp1,listComp1).
%
% Inputs
%   S1(n1,1)			Vector of values
%   mail1			Mesh of the field
%   intg1			Integration points of the field
%   numer1(nbptg1)		Ordering of integration points used
%   mapComp1(nbptg1,nbcomp1)	Mapping matrix
%   listComp1{nbcomp1}		List of component names
% Outputs
%   cham1			Element-based field

clear cham1;
nbcomp1 = length(listComp1);
nbptg1  = Nptg(mail1,intg1);
numerptg1 = [1:nbptg1];

nbval = 1; % Cas particulier ou 1 composante n'a que 1 valeur


% Boucle sur les sous-zones du champ par element
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  intge1  = intg1{zo1};
  maile1  = mail1{zo1};
  nbpt1   = length(intge1.WEIGHT);
  [nbel1,nbno1] = size(maile1.MAIL);

  clear chamel1;

% Boucle sur les composantes du champ par element
  for i = 1:nbcomp1
    name = listComp1(i);
%%DD13/04/04    xval1 = zeros(nbptg1,nbval);
    xval1 = zeros(nbel1,nbval*nbpt1);

    j = findoccur(numerptg1,numer1);
    k = find(j); j = j(k);
%%DD13/04/04    xval1(k,1) = S1(mapComp1(j,i),:);
    xval1(:,:) = reshape(S1(mapComp1(j,i),:),nbpt1,nbel1)';

    chamel1{i} = struct('COMP',name,'UNIT','','XVAL',xval1);
  end

  cham1{zo1} = chamel1;
end
