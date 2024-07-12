function [list1] = ListTypeMesh(mail1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2002

% Retourne la liste des types d'elements de mail1
% (autant que de zones, meme s'il y a des memes types dans
% des zones differentes)

list1 = [];

% Loop on zones
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  typ1 = mail1{zo1}.TYPE;
  if isempty(list1)
    list1 = [{typ1}];
  else
    list1 = [list1 {typ1}];
  end
end
