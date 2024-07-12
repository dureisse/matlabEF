function [mail2] = ElemMeshLeadingStrictly(mail1,nmail1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 18 / 01 / 2003

% Get the elements of mail1 that leads strictly on the nodes of nmail1
% (i.e. for which all the nodes are in nmail1)

nmail1 = ChangeMesh2(nmail1,'POI1'); % with only 1 subzone
numer1 = nmail1{1}.MAIL';

clear mail2; zo2 = 0;

% Loop on subzones of mail1
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  topo1  = maile1.MAIL;
  [nbel1,nbnn1] = size(topo1);
  topo2  = zeros(0,nbnn1);
  for el1 = 1:nbel1
    lnode1 = findoccur(topo1(el1,:),numer1);
    if (length(find(lnode1)) == nbnn1)
%     An element has been found
      topo2(end+1,:) = topo1(el1,:);
    end
  end
  if topo2
%   A subzone has been found
    zo2 = zo2 + 1;
    mail2{zo2} = struct('TYPE',maile1.TYPE,'MAIL',topo2);
  end
  clear topo2 topo1 maile1;
end

clear numer1 nmail1;
