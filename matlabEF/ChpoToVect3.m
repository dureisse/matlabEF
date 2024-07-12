function [U1] = ChpoToVect3(chpo1,nmail1,numer1,mapddl1,listComp1)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 02 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 04 / 12 / 2002
%   Cas ou on assemble pas tout (j=0)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2003
%   Modification pour assemblages speciaux
%
% Assemble un champ par point (chpo1,nmail1) dans un vecteur U1
% selon le mapping mapddl1 (voir MapddlChpo)

nbddl1 = max(max(mapddl1));
U1 = zeros(nbddl1,1);

% Boucle sur les sous-zones du champ par point
nbzone1 = length(chpo1);
for zo1 = 1:nbzone1
  chpoel1 = chpo1{zo1};
  nmaile1 = nmail1{zo1};
% Boucle sur les composantes du champ par point
  nbcomp1 = size(chpoel1,2);
  for i = 1:nbcomp1
    xval1 = chpoel1{i}.XVAL;
    if (size(xval1,2) ~= 1)
      error('multiple value component not implemented')
    end

%   Component in listComp1
    j = findoccur({chpoel1{i}.COMP},listComp1);

    if j~=0
%     Nodes in nmail1, in numer1
      list_node1 = find(mapddl1(:,j));
      list_node3 = numer1(list_node1); % nodes in numer1
      list_node2 = nmaile1.MAIL';      % nodes in nmail1
      [junk,ia,ib] = intersect(list_node3,list_node2);
%     junk=list_node3(ia)=list_node2(ib)
%     list_node3(ia)=numer1(list_node1(ia))

      list_ddl1  = mapddl1(list_node1(ia),j);
%%DD 03/12/25      U1(list_ddl1,1) = U1(list_ddl1,1) + xval1(ib,1);
%     For special cases where several values go at the same position:
%     the result should sum the values and not override them
      for k = 1:length(ib)
        U1(list_ddl1(k),1) = U1(list_ddl1(k),1) + xval1(ib(k),1);
      end
    end

  end
end
