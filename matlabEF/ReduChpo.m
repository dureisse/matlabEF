function [chpo1,nmail1] = ReduChpo(chpo2,nmail2,nmail3)
% Reduction of a nodal field onto a mesh
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 12 / 2005
%  Treatment of void chpo2

% Reduit un champ par point (chpo2,nmail2)
% sur un maillage nmail3,
% en un champ par point (chpo1,nmail1).
% (voir IntersectMesh.m ReduChml.m)


nbzone3 = TestMeshType(nmail3,'POI1');
nbzone2 = TestMeshType(nmail2,'POI1');

if (nbzone3 ~= 1)
  nbzone3
  error('More than 1 zone not yet implemented')
end

numer3 = nmail3{1}.MAIL';

chpo1 = [];
nmail1 = [];
if length(chpo2);
  zo1 = 0;
  for zo2 = 1:nbzone2
    numer2 = nmail2{zo2}.MAIL';
    chpoe2 = chpo2{zo2};
    [numer1,ia,ib] = intersect(numer2,numer3);
%   numer1=numer2(ia)=numer3(ib)
    if numer1
%     If there is an intersection, all components are selected
      clear chpoe1;
      for i = 1:length(chpoe2)
        chpoe1{i} = struct('COMP',chpoe2{i}.COMP,'UNIT',chpoe2{i}.UNIT);
        xval2 = chpoe2{i}.XVAL;
        xval1 = xval2(ia,:);
        chpoe1{i}.XVAL = xval1;
        clear xval1 xval2,
      end
      zo1 = zo1 + 1;
      nmail1{zo1} = struct('TYPE','POI1','MAIL',numer1');
      chpo1{zo1} = chpoe1;
      clear chpoe1;
    end
    clear chpoe2 numer2 numer1;
  end
else
  zo1 = 0;
  for zo2 = 1:nbzone2
    numer2 = nmail2{zo2}.MAIL';
    [numer1,ia,ib] = intersect(numer2,numer3);
%   numer1=numer2(ia)=numer3(ib)
    if numer1
      zo1 = zo1 + 1;
      nmail1{zo1} = struct('TYPE','POI1','MAIL',numer1');
    end
    clear numer1 numer2;
  end
end
