function [mail1] = ManuMesh(varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 12 / 2002

% Construction manuelle d'un maillage
% On donne les couples (type1,topo1) ou
%   type1 est le nom du type d'element,
%   topo1(nbel1,nbnn1) est la connectivite

narg = nargin; % Nombre d'arguments
nbzone1 = fix(narg/2);
if ~(nbzone1*2 == narg)
  narg
  nvalue1
  error('Bar number of arguments')
end

clear mail1;
for zo1 = 1:nbzone1
  type1 = varargin{(zo1-1)*2+1};
  topo1 = varargin{(zo1-1)*2+2};
  maile1 = struct('TYPE',type1,'MAIL',topo1);
  mail1{zo1} = maile1;
  clear maile1 type1 topo1;
end
