function [numer1] = LocateSuperposedNodes(xcoor1,xcrit1);
% Find nodes that share the same location inspace
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 24 / 02 / 2005
%
% Use the coordinates of a liste of n nodes xcoor1(n,idim)
% to find superposed nodes.
%
% Inputs
%   xcoor1(n,idim)	Coordinates of the nodes
%   xcrit1		Proximity criteria to decide of superposition
% Outputs
%   numer1(n)		Renumbering
%
% if numer1(i) = i, the node i is unique
% if numer1(i) = j < i, the node i is superposed with the node i
% superposed means that max(abs(xcoor1(i,:) - xcoor1(j,:))) < xcrit1

[n,idim] = size(xcoor1);

% Initial value of numer1
numer1 = [1:n];

% Find a local renumbering: lexicographic order with xcrit1 as
% threshold

% sort along 1st coordinate
[junk,i1] = sort(xcoor1(:,1));
xcoor2 = xcoor1(i1,:);

lblocks = [];
inod    = 1;
iblock  = 1;

% Loop on blocks
while (inod < n)

% Find block of potential neighbours of the current node
  toto = xcoor2(inod:n,1) - xcoor2(inod,1);
  ii = find (toto > xcrit1);
  if isempty(ii); ii = n+1; end
  inod:ii(1)-1




end
