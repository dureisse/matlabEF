function [nbzone1] = TestMeshType(mail1,typ1)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002

% Test if the mesh mail1 has only elements of the type typ1,
% If yes, nbzone1 is the number of sub-zones
% If no, nbzone1=0

nbzone0 = length(mail1);

nbzone1 = nbzone0;
for zo0 = 1:nbzone0
  if ~strcmp(mail1{zo0}.TYPE,typ1)
    nbzone1 = 0;
  end
end
