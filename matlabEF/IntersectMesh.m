function [mail3] = IntersectMesh(mail1,mail2)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2002

% Construit le maillage mail3 intersection des maillages
% mail1 et mail2
% Un element est commun quand tous les noeuds sont communs
% (dans l'ordre ou dans l'ordre inverse)
% et que les elements sont de meme type


nbzone1 = length(mail1);
nbzone2 = length(mail2);

% List of elements type of mail2
ListType2 = ListTypeMesh(mail2);

clear mail3; zo3 = 0;

% Loop on elements of mail1
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  typ1   = maile1.TYPE;
  topo1  = maile1.MAIL;
  [nbel1,nbnn1] = size(topo1);
  listzo2 = findoccur([{typ1}],ListType2);

  topo3 = zeros(0,nbnn1);
  for i2 = 1:length(listzo2)
    zo2 = listzo2(i2);
    if zo2
      maile2 = mail2{zo2};
      topo2  = maile2.MAIL;
      nbel2 = size(topo2,1);

%     Inefficient but robust...
      for el1 = 1:nbel1
        nel1 = topo1(el1,:);
        for el2 = 1:nbel2
          nel2 = topo2(el2,:);
          if all(ismember(nel1,nel2))
            topo3 = [topo3 ; nel1];
            break
          end 
        end
      end

%     More efficient, but with a pb...
%      nel2 = sort(topo2,2);
%      for el1 = 1:nbel1
%        nel0 = topo1(el1,:);
%        nel1 = sort(nel0,2);
%        if ismember(nel1,nel2,'rows')
%          topo3 = [topo3 ; nel0];
%          break
%        end
%      end

      clear nel1 nel2 topo2 maile2;

%     Wrong
%      topo2i = topo2(:,[nbnn1:-1:1]); % reverse order is considered
%%     Keep as possible the same order of elements of mail1
%      [t1,ia,ib]  = intersect(topo1,topo2,'rows'); t1 = topo1(sort(ia),:);
%      [t1i,ia,ib] = intersect(topo1,topo2i,'rows'); t1i = topo1(sort(ia),:);
%      t1 = [t1;t1i];
%      [t3,ia,ib] = unique(t1,'rows'); t3 = t1(sort(ia),:);
%      topo3 = [topo3 ; t3];
%      clear maile2 topo2 topo2i t1 t1i t3;
    end

% PLUS SIMPLE (A REMPLACER, DONC):
%    if zo2
%      maile2 = mail2{zo2};
%      topo2  = maile2.MAIL;
%      topo2i = topo2(:,[nbnn1:-1:1]); % reverse order is considered
%%     Keep as possible the same order of elements of mail1
%      [t1,ia1,ib]  = intersect(topo1,topo2,'rows');
%      [t1i,ia1i,ib] = intersect(topo1,topo2i,'rows');
%      listel1 = unique([ia1 ia1i]); % result is sorted
%      topo3 = [topo3 ; topo1(listel1,:)];
%      clear maile2 topo2 topo2i t1 t1i listel1;
%    end

  end
  
  if topo3
%   There is an intersection as a new zone
    zo3 = zo3 + 1;
    mail3{zo3} = struct('TYPE',typ1,'MAIL',topo3);
  end

  clear maile1 topo1 topo3;
end

if ~exist('mail3')
  mail3 = [];
end
