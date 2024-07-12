function [chpo1] = ManuChpo(nmail1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 06 / 12 / 2002          

% Construction d'un champ par point constant manuellement
% On donne le maillage nmail1 (de type POI1)
% puis les triplet (nomi,uniti,valuei) i=1,...
% de noms de composante, unite et valeur
% (voir ManuChml)

narg = nargin; % Nombre d'arguments
if (narg >= 1)
  nvalue1 = fix((narg-1)/3);
  if ~(nvalue1*3 == narg-1)
    narg
    nvalue1
    error('Bar number of arguments')
  end
else
  error('Bar number of arguments2')
end

% Number of elements
nbzone0 = TestMeshType(nmail1,'POI1');
if (nbzone0 == 0)
  error('Mesh should be POI1')
end

% Loop on zones
clear chpo1;
for zo1 = 1:nbzone0
  nbpt1 = size(nmail1{zo1}.MAIL,1);
% Fill in
  clear chpoe1;
  for c1=1:nvalue1
    nomi   = varargin{(c1-1)*3+1};
    uniti  = varargin{(c1-1)*3+2};
    valuei = varargin{(c1-1)*3+3};
    chpoe1{c1} = struct('COMP',nomi,'UNIT',uniti, ...
                        'XVAL',valuei*ones(nbpt1,1));
  end
  chpo1{zo1} = chpoe1;
end
