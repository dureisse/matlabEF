function [lnod1,lind1] = MeshToGraph(mail1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2002

% Transform a mesh into a connectivity graph.
% Some information is lost during the conversion (type of elements...)
% The edges of the graph are the elements, the vertices are the nodes:
% the element iel (local numbering of elements) is linked to the nodes
% lnod1(lind1(iel):lind1(iel+1)-1)

lnod1 = []; ind1 = 0;
lind1 = [1];

% Loop on elements
iel = 0;
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  topo1 = mail1{zo1}.MAIL;
  [nbel1,nbno1] = size(topo1);
  for el1 = 1:nbel1
    iel = iel + 1;
    lnod1(ind1+1:ind1+nbno1) = topo1(el1,:);
    ind1 = ind1 + nbno1;
    lind1(iel+1) = ind1 + 1; 
  end
end
