function [r1,modlr1] = IntegrOperator2(modl1,mail1,intg1,xcoor1,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 03 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 03 / 2003
%  Ajout du modele modlr1
%
% Matrices elementaires de l'operateur d'integration
% (s'applique sur un champ scalaire)
%
% Entrees
%   modl1	Modele
%   mail1	Maillage
%   intg1	Segment d'integration
%   xcoor1	Pile des noeuds
%   mode1	Modele de l'analyse
%
% Sorties
%   r1		Operateur d'integration
%   modlr1	Modele associe

nbzone1 = length(mail1);
if (nbzone1 ~= length(intg1)) | (nbzone1 ~= length(modl1))
  nbzone1
  length(intg1)
  length(modl1)
  error('Bad number of sub-zones')
end


nbcopB = 1;    % Number of physical components (SCAL)

clear r1 modlr1;

% Loop on zones
for zo1 = 1:nbzone1
  intge1 = intg1{zo1};
  maile1 = mail1{zo1};
  modle1 = modl1{zo1};

  topo1 = maile1.MAIL;
  nbel  = size(topo1,1);
  nbptg = length(intge1.WEIGHT);

  xval1 = zeros(nbcopB*nbptg,nbcopB*nbptg,nbel);

% Loop on elements
  for el1 = 1:nbel
    node1 = topo1(el1,:);
    xcoorel = xcoor1(node1,:);
    [KE] = ElementIntegrOperator(modle1,intge1,xcoorel,mode1);
    xval1(:,:,el1) = KE;
    clear KE;
  end
  r1{zo1} = struct('XVAL',xval1);
  modlre1 = struct('NCOP',[1],'NCOD',[1]);
  modlre1.COMP = [{'SCAL'}];
  modlre1.COMD = [{'SCAL'}];
  modlr1{zo1} = modlre1;
  clear xval1;
end
