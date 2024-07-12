function [mater1] = ManuChml(mail1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002          

% Construction d'un champ constant par element manuellement
% (ce champ est en fait constant sur tout le maillage)
% On donne le maillage mail1
% puis les triplet (nomi,uniti,valuei) i=1,...
% de noms de composante, unite et valeur

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
nbelt1 = Nelements2(mail1);

% Fill in material field
clear mater1;
for c1=1:nvalue1
  nomi = varargin{(c1-1)*3+1};
  uniti = varargin{(c1-1)*3+2};
  valuei = varargin{(c1-1)*3+3};
  mater1{c1} = struct('COMP',nomi,'UNIT',uniti, ...
                      'XVAL',valuei*ones(nbelt1,1));
end
