function [mail3] = ContourMesh3(mail1,idim)
% Contour (en 2D) ou peau (en 3D) du maillage mail1
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 02 / 2003
%   Cas du SEG2 : tous les noeuds sont sur le bord si idim > 1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 12 / 2004
%   Ajout du CUB8 en 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 05 / 2008
% ICETA Damien      L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 05 / 2008
%   Ajout du POI1
%
% Construit le maillage mail3, contour du maillage mail1 en 2D,
% peau du maillage mail1 en 3D.


% In order to place element types correctly,
% build an empty mesh mail3 with all element types possible
% (2D or POI1 only up to now).
clear mail3;
mail3{1} = struct('TYPE','SEG2','MAIL',zeros(0,2));
mail3{2} = struct('TYPE','SEG3','MAIL',zeros(0,3));
mail3{3} = struct('TYPE','POI1','MAIL',zeros(0,1));
mail3{4} = struct('TYPE','QUA4','MAIL',zeros(0,4));

% From here, we use only local node numbering
nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = nmail1{1}.MAIL';
numer_inv1 = InverseList(numer1,max(numer1));
clear nmail1;
mail2 = RenumMesh(mail1,numer_inv1);

% With local node numbering, it is easier to work
% on graphs rather than on meshes.
% The element i possesses the nodes lnod1(lind1(i):lind1(i+1)-1)
[lnod1,lind1] = MeshToGraph(mail2);

% Inverse the graph.
% The node j is linked to the elements lnod2(lind2(j):lind2(j+1)-1)
[lnod2,lind2] = InverseGraph2(lnod1,lind1);


% Find the contour mail3
% Loop on elements
nbzone1 = length(mail2);
for zo1 = 1:nbzone1
  maile1 = mail2{zo1};
  typ1   = maile1.TYPE;
  topo1  = maile1.MAIL;
  [nbel1,nbnn1] = size(topo1);

% For each element, look for the neighbors of each edge (face in 3D)
  if strcmp(typ1,'TRI3')
%   3 SEG2 for edges
    facZone = [1 1 1];
    facTemplate{1} = [1 2];
    facTemplate{2} = [2 3];
    facTemplate{3} = [3 1];
  elseif strcmp(typ1,'QUA4')
%   4 SEG2 for edges
    facZone = [1 1 1 1];
    facTemplate{1} = [1 2];
    facTemplate{2} = [2 3];
    facTemplate{3} = [3 4];
    facTemplate{4} = [4 1];
  elseif strcmp(typ1,'TRI6')
%   3 SEG3 for edges
    facZone = [2 2 2];
    facTemplate{1} = [1 2 4];
    facTemplate{2} = [2 3 5];
    facTemplate{3} = [3 1 6];
  elseif strcmp(typ1,'QUA8')
%   4 SEG3 for edges
    facZone = [2 2 2 2];
    facTemplate{1} = [1 2 5];
    facTemplate{2} = [2 3 6];
    facTemplate{3} = [3 4 7];
    facTemplate{4} = [4 1 8];
  elseif strcmp(typ1,'SEG2')
%   2 POI1 for edges
    facZone = [3 3];
    facTemplate{1} = [1];
    facTemplate{2} = [2];
  elseif strcmp(typ1,'POI1')
%   1 POI1 for edges
    facZone = [3];
    facTemplate{1} = [1];
  elseif strcmp(typ1,'CUB8')
%   6 QUA4 as faces
    facZone = [4 4 4 4 4 4];
    facTemplate{1} = [1 2 3 4];
    facTemplate{2} = [5 6 7 8];
    facTemplate{3} = [1 2 6 5];
    facTemplate{4} = [2 3 7 6];
    facTemplate{5} = [3 4 8 7];
    facTemplate{6} = [4 1 5 8];
  else
    typ1
    error('element not implemented yet')
  end

  if (strcmp(typ1,'SEG2') && (idim > 1))
%   Special case where all nodes are in the contour
    fac1 = 1; 
    zo3 = facZone(fac1);
    topo3 = mail3{zo3}.MAIL;
    mail3{zo3}.MAIL = [topo3 ; unique([topo1(:,1) ; topo1(:,2)])];
    clear topo3;
  else
%
  for el1 = 1:nbel1 
%   Loop on faces
    for fac1 = 1:length(facZone)
      nodes1 = topo1(el1,facTemplate{fac1});
      elvois1 = FindNeighborElement(nodes1,lnod2,lind2);
      if ~elvois1
        el1
        nodes1
        error('PB: a face has even not this element as neighbor...')
      end
      if length(elvois1) == 1
%       If no neighbors, the edge is in the contour 
        zo3 = facZone(fac1);
        topo3 = mail3{zo3}.MAIL;
        mail3{zo3}.MAIL = [topo3 ; nodes1];
        clear topo3;
      end
      clear nodes1;
    end
  end
%
  end

end
clear mail2 lnod1 lind1 lnod2 lind2;

% Clean up mail3 (get rid of void zones)
clear mail4; zo4 = 0;
nbzone3 = length(mail3);
for zo3 = 1:nbzone3
  if mail3{zo3}.MAIL
    zo4 = zo4 + 1;
    mail4{zo4} = mail3{zo3};
  end
end
clear mail3;

% Coming back to global node numbering
mail3 = RenumMesh(mail4,numer1);





%-----------------------------------------------------------------------
function [elvois1] = FindNeighborElement(nodes1,lnod2,lind2)

nbnode1 = length(nodes1);
i1 = 1;
  node1 = nodes1(i1);
  elvois1 = lnod2(lind2(node1):lind2(node1+1)-1);
for i1 = 2:nbnode1
  node1 = nodes1(i1);
  elvois1 = intersect(elvois1,lnod2(lind2(node1):lind2(node1+1)-1));
end
