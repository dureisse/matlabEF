function [ddl1] = FindDdl(mapddl1,listddl1,listddl2)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 11 / 2002

% Trouve ou sont les ddl de noms listddl1 dans listddl2

ddl1 = findoccur(listddl2,listddl1);
ddl1 = mapddl1(:,ddl1);
ddl1 = reshape(ddl1',1,size(ddl1,1)*size(ddl1,2));
ddl1 = ddl1(find(ddl1));

