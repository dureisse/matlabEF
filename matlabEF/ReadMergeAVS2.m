function [xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
          ListChml1,ListCara1,error1] = ReadMergeAVS2(xthresh1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 19 / 12 / 2002
%   plusieurs fichiers
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 01 / 2003
%   plusieurs seuils xthresh1
%
% Lecture de n1 passes dans un fichier AVS - UCD, et fusion
% des bases de donnees correspondantes
% (les noeuds sont fusionnes si distants de moins de xthresh1)
% On donne le critere xthresh1, puis les couples (fichier,nombre de passe)

clear ListXcoor1 ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

narg = nargin; % Nombre d'arguments
if (narg >= 3)
  nfiles1 = fix((narg-1)/2);
  if ~(nfiles1*2 == narg-1)
    narg
    nfiles1
    error('Bar number of arguments')
  end
else
  error('Bar number of arguments 2')
end

% Total number of readings
nread1 = 0;
for ifile1 = 1:nfiles1
  n1 = varargin{ifile1*2};
  nread1 = nread1 + n1;
end
if length(xthresh1) == 1
  xthresh1 = xthresh1 * ones(1,nread1);
end
if length(xthresh1) ~= nread1
  nread1
  xthresh1
  error('Bad number of thresholds')
end

% Loop on files
iread1 = 0;
for ifile1 = 1:nfiles1
  file1 = varargin{ifile1*2-1};
  n1    = varargin{ifile1*2};
  disp(['  File number ' int2str(ifile1) ': ' file1])

  fid = fopen(file1,'rt');

%   Reading
  for i = 1:n1
  disp(['  Reading number ' int2str(i)])
    [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
if error1
  error(error1)
end
    iread1 = iread1 + 1;
    ListXcoor1{iread1} = xcoor1;
    ListMesh1{iread1} = mail1;
    ListChpo1{iread1} = chpo1;
    clear nmail1;
    nbno1 = size(xcoor1,1);
    nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');
    ListnMesh1{iread1} = nmail1;
    ListChml1{iread1} = chml1;
    ListCara1{iread1} = cara1;
    clear xcoor1 mail1 chpo1 chml1 cara1;
  end

  fclose(fid);
end

% Merging
clear PileMail1;
clear ListXthresh1;
xcoort1 = zeros(0,3);
for i = 1:nread1
  decal1 = size(xcoort1,1);
  xcoort1 = [xcoort1 ; ListXcoor1{i}];
  mail1 = ListMesh1{i};
  mail1 = ShiftNodeNumber(mail1,decal1);
  nmail1 = ListnMesh1{i};
  nmail1 = ShiftNodeNumber(nmail1,decal1);
  PileMail1{2*(i-1)+1} = mail1;
  PileMail1{2*(i-1)+2} = nmail1;
  ListXthresh1(2*(i-1)+1) = xthresh1(i);
  ListXthresh1(2*(i-1)+2) = xthresh1(i);
  clear mail1 nmail1;
end
clear ListMesh1 ListnMesh1 ListXcoor1;

% Elimination of useless nodes
[PileMail1,ListUsedNodes] = ElimNode5(PileMail1,xcoort1,ListXthresh1);
numer_inv1 = InverseList(ListUsedNodes',max(ListUsedNodes));
for i = 1:nread1
  mail1 = PileMail1{2*(i-1)+1};
  mail1 = RenumMesh(mail1,numer_inv1); 
  nmail1 = PileMail1{2*(i-1)+2};
  nmail1 = RenumMesh(nmail1,numer_inv1); 
  ListMesh1{i} = mail1;
  ListnMesh1{i} = nmail1;
  clear mail1 nmail1;
end
clear PileMail1;

xcoort1 = xcoort1(ListUsedNodes,:);
