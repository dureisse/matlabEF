function [bE] = AssembleOperator(num1,map1,nam1,num2,map2,nam2);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2002

% Use two mappings of dof, to build the mapping matrix from
% one dof ordering to the other

% Input
%   num1(nbnode1)         node numbers of the first mapping
%   nam1(nbname1)         dof names of the first mapping
%   map1(nbnode1,nbname1) dof number in the first mapping
%   num2(nbnode2)         node numbers of the second mapping
%   nam2(nbname2)         dof names of the second mapping
%   map2(nbnode2,nbname2) dof number in the second mapping

% Output
%   bE(nbddl1,nbddl2)     boolean mapping operator

% The dof number of node occurence ino1, i.e. of node number
% num1(ino1), and name occurence ina1, i.e. of name nam1(ina1),
% is dof1 = map1(ino1,ina1)
% if dof1 = 0, the dof is not present in the assembled vector.
% Same thing for the second mapping.
% If V1 is the vector of the first mapping, and V2 of the second one,
% V1 = bE * V2

nbddl1 = max(max(map1));
nbddl2 = max(max(map2));
bE = sparse(nbddl1,nbddl2);

% Find num1 in num2, and nam1 in nam2
  in = findoccur(num1,num2);
  id = findoccur(nam1,nam2);
  for ino1 = 1:length(num1)
    for ina1 = 1:length(nam1)
      ddl1 = map1(ino1,ina1);
      ino2 = in(ino1);
      ina2 = id(ina1);
%if (num1(ino1) ~= num2(ino2) | ~all(nam1{ina1} == nam2{ina2}))
%  error('PBPB')
%end
      if (ino2 && ina2 && ddl1)
        ddl2 = map2(ino2,ina2);
        if ddl2
          bE(ddl1,ddl2) = 1.;
        end
      end 
    end
  end
