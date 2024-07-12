function [lnod2,lind2,varargout] = InverseGraph2(lnod1,lind1,varargin)
% Graph inversion (with optional order preservation)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 04 / 2004
%  Ajout sortie optionnelle pour conserver l'ordre
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 17 / 06 / 2005
%  traitement d'un graphe vide
%
% [lnod2,lind2] = InverseGraph2(lnod1,lind1);
% Inversion of a graph (lno1,lind1) into the dual graph (lnod2,lind2).
% [lnod2,lind2,lord2] = InverseGraph2(lnod1,lind1);
% Same operation, but get the information about initial ordering in
% such a way that
% [lnod1,lind1] = InverseGraph2(lnod2,lind2,lord2);
% gives again the same (ordered) (lnod1,lind1)
%
% Inputs
%   lnod1(n1)		Graph to be inverted
%   lind1(nbel+1)	Indirection vector if this graph
% Optional input
%   lord1(n1)		Information on edge number order
% Outputs
%   lnod2(n2)		Inverted graph
%   lind2(nbno+1)	Indirection vector if this graph
% Optional outputs
%   lord2(n2)		Information on initial edge number,
%                       usefull to rebuilt initial graph with the
%                       same ordering
% Comments
%   Beware that the numbering of vertices of the first graph must be
%   local, i.e. continuous from 1 to nbno.
% About graphs
%   The edge iel of the first graph is linked to the vertices
%   lnod1(lind1(iel):lind1(iel+1)-1)
%   The edge ino of the second graph is linked to the vertices
%   lnod2(lind2(ino):lind2(ino+1)-1)
% About optional numberings
%   Optional arguments lord1 and lord2 are mutually exclusive.
%   With the numbering lord2,
%   lord2(lind2(ino):lind2(ino+1)-1) is the occurence number of
%   the edge ino of (lnod2,lind2) in the vertices list
%   {lnod1(lind1(iel):lind1(iel+1)-1)}
%   Therefore, one gets the property
%     if    ielk = lnod2(lind2(ino)+k-1)  (1 <= k <= lind2(ino+1)-lind2(ino))
%     and   iloc = lord2(lind2(ino)+k-1)
%     then  lnod1(lind1(iel)+iloc-1) = ino
%   The following sequence gives the same ordered graphs
%   (lnod3,lind3) = (lnod1,lind1)
%     [lnod2,lind2,lord2] = InverseGraph2(lnod1,lind1);
%     [lnod3,lind3] = InverseGraph2(lnod2,lind2,lord2);
%   The following sequence gives the same unordered graphs but
%   maybe not with the same order
%     [lnod2,lind2] = InverseGraph2(lnod1,lind1);
%     [lnod3,lind3] = InverseGraph2(lnod2,lind2);

nin1 = nargin-2;
switch nin1
  case 0,
    clear lord1;
  case 1,
    lord1 = varargin{1};
  otherwise,
    nin1
    error('Bar number of input arguments')
end
nout = nargout-2;

if ((nin1 == 1) && (nout == 1))
  error('Input and output orderings are mutually exclusive')
end

% Matrix of connectivity topo2
nbno = max(lnod1);
nbel = length(lind1)-1;
if nbel ~= 0
  topo2 = sparse(nbno,nbel);
  for iel = 1:nbel
    nodes1 = lnod1(lind1(iel):lind1(iel+1)-1);
    if exist('lord1')
      ordre1 = lord1(lind1(iel):lind1(iel+1)-1);
    else
      ordre1 = [1:length(nodes1)];
    end
%  topo2(nodes1,iel) = 1;
    topo2(nodes1,iel) = ordre1';
    clear nodes1 ordre1;
  end
else
  topo2 = [];
end

% Building the dual graph
lnod2 = []; ind2 = 0;
lind2 = [1];
lord2 = [];

for ino = 1:nbno
  listel = find(topo2(ino,:));
  nbel2  = length(listel);
  if exist('lord1')
    listel = listel(topo2(ino,listel));
  end
  lnod2(ind2+1:ind2+nbel2) = listel;
  lord2(ind2+1:ind2+nbel2) = topo2(ino,listel);
  ind2 = ind2 + nbel2;
  lind2(ino+1) = ind2 + 1;
  clear listel;
end

if nbel == 0
  lind2 = [1 1];
end

switch nout
  case 0,
  case 1,
    varargout(1) = {lord2};
  otherwise
    nout
    error('Bad number of output arguments')
end

clear topo2;
