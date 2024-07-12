function [error1] = WriteMergeAVS(xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
                                  ListChml1,ListCara1,file1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 04 / 2004

% Ecriture d'une base de donnees en plusieurs passes dans un fichier AVS - UCD

fid = fopen(file1,'w');

n1 = length(ListMesh1);

for i = 1:n1
  mail1 = ListMesh1{i};
  nmail1 = ListnMesh1{i};
  chpo1 = ListChpo1{i};
  chml1 = ListChml1{i};
  cara1 = ListCara1{i};
%DD 08/01/03  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS(xcoort1,mail1,nmail1,chpo1);
%DD 26/04/04  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS2(xcoort1,mail1,nmail1,chpo1);
  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS3(xcoort1,mail1,nmail1,chpo1);
  error1 = Write1AVS5(xcoor1a,mail1a,chpo1a,chml1,cara1,fid);
end

fclose(fid);
