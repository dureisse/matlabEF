function [mapddl2] = FindDof(numer1,mapddl1,listddl1,numer2,listddl2)
% Find degrees of freedom in an assembling mapping
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 10 / 04 / 2004
%
% Inputs
%  numer1(nbno1)		List of nodes in the mapping
%  mapddl1(nbno1,nbname1)	Mapping matrix
%  listddl1{nbname1}		Array of dof names
%  numer2(nbno2)		List of searched nodes
%  listddl2{nbname2}		Array of searched dof names
%
% Outputs
%  mapddl2(nbno2,nbname2)	Mapping matrix of searched dofs
%
% Comments
%   mapddl1(no1,name1) gives the number of the dof associated to the
%   node no1 and with the dof name name1 (0 if not found).

[nbno1,nbname1] = size(mapddl1);
if (nbno1 ~= length(numer1) | nbname1 ~= length(listddl1))
  nbno1
  length(numer1)
  nbname1
  length(listddl1)
  error('Unconsistent datas')
end

nbno2 = length(numer2);
nbname2 = length(listddl2);
mapddl2 = zeros(nbno2,nbname2);

% Find occurences of node numbers 
in = findoccur(numer2,numer1);
% Find occurences of dof names
jn = findoccur(listddl2,listddl1);

% Avoid null values (no occurences)
il = find(in); in = in(il);
jl = find(jn); jn = jn(jl);

mapddl2(il,jl) = mapddl1(in,jn);

clear in il jn jl;
