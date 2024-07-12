function [mail2] = ElemMesh(mail1,listelem1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2002

% Extrait le maillage mail2 du mailllage mail1 qui contient les
% elements listelem1 de mail1
% (numerotations globales des elements)

listelem1 = sort(listelem1);
nbelt1  = length(listelem1);
nbzone1 = length(mail1);

clear mail2; zo2 = 0;
i1 = 1;
nbelc1 = 0;
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  topo1  = maile1.MAIL;
  [nbel1,nbnn1] = size(topo1);
  nbelc1 = nbelc1 + nbel1;

  ielem2 = 0; topo2 = zeros(0,nbnn1);
  while ((i1 <= nbelt1) && (listelem1(i1) <= nbelc1))
%   The element is in this zone
    ielem2 = ielem2 + 1;
    ielem1 = nbel1 - nbelc1 + listelem1(i1);
    topo2(ielem2,:) = topo1(ielem1,:);
    i1 = i1 + 1;
  end

  if ielem2
    zo2 = zo2 + 1;
    mail2{zo2} = struct('TYPE',maile1.TYPE,'MAIL',topo2);
  end
  clear topo2;
end

% disp(['Nb elements ' int2str(ielem2) ' Nb zones ' int2str(zo2)])
