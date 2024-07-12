function [error1] = WriteMeshGMSH(xcoor1,ListMesh1,fid,varargin)
% Write a list of meshes in a gmsh ascii meshing file (format 1.0)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 14 / 04 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 04 / 08 / 2006
%  Ajout CUB8, TET4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 08 / 2006
%  Possibilite de ne pas renumeroter
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2006
%  Ajout PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 02 / 2007
%  Le TRI6 est fait en tant que tel (on ne sous-decoupe plus en TRI3)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 03 / 2007
%  Ajout QUA8
% DUREISSEIX David  LaMCoS                           le 27 / 11 / 2015
%  Ajout TE10
%
% Inputs
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   ListMesh1{nmesh1}	List of meshes
%   fid			file identifier (file must be opened, see fopen)
% Optional input
%   option1		Option: 'None' is default,
%                               'PreserveNodeNumber'
% Output
%   error1		

error1 = 0;

nin1 = nargin - 3;
switch nin1
  case 0,
    option1 = 'None';
  case 1,
    option1 = varargin{1};
  otherwise,
    nin1
    error('Bad number of optional input arguments')
end

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('3D assumed')
end

nmesh1 = length(ListMesh1);
disp([' Number of physical entities: ' int2str(nmesh1)])

% Local node numbering for the meshes
% """""""""""""""""""""""""""""""""""

numer1 = []; % usefull nodes
nbelt1 = 0; % Total number of elements if all the meshes
for imesh1 = 1:nmesh1
  mail1  = ListMesh1{imesh1};
  nmail1 = ChangeMesh2(mail1,'POI1');
  numer1 = [numer1 nmail1{1}.MAIL'];
  nbelt1 = nbelt1 + Nelements2(mail1);
  clear mail1 nmail1;
end
numer1 = unique(numer1);

if ~strcmp(option1,'PreserveNodeNumber')

  numer_inv1 = InverseList(numer1,max(numer1));
  clear ListMesh2;
  for imesh1 = 1:nmesh1
    mail1  = ListMesh1{imesh1};
    lmail1 = RenumMesh(mail1,numer_inv1); % Locally numbered mesh
    ListMesh2{imesh1} = lmail1;
    clear lmail1 mail1;
  end
  clear numer_inv1;

else

  ListMesh2 = ListMesh1;
  numer1 = [1:size(xcoor1,1)];

end

% Nodes
% """""
  count = fprintf(fid,'$NOD\n');
  nbnot1 = length(numer1);
  count = fprintf(fid,'%d\n',nbnot1);
  count = fprintf(fid,'%d %e %e %e\n',[[1:nbnot1]; xcoor1(numer1,:)']);
  count = fprintf(fid,'$ENDNOD\n');

% Elements
% """"""""
  ielt1 = 0;
  count = fprintf(fid,'$ELM\n');
  count = fprintf(fid,'%d\n',nbelt1);

% Loop on meshes (physical entities)
  for imesh1 = 1:nmesh1
    lmail1 = ListMesh2{imesh1};

  nbzone1 = length(lmail1);
  for izo1 = 1:nbzone1
    maile1 = lmail1{izo1};
    type1 = maile1.TYPE;
    topo1 = maile1.MAIL;
    [nel1,nno1] = size(topo1);
%   GMSH type of element
    switch type1
      case 'SEG2',
        permut1 = [1 2];
        ityp1 = 1;
      case 'SEG3',
        permut1 = [1 2 3];
        ityp1 = 8;
      case 'TRI3',
        permut1 = [1 2 3];
        ityp1 = 2;
      case 'TRI6',
        permut1 = [1 2 3 4 5 6];
        ityp1 = 9;
      case 'QUA4',
        permut1 = [1 2 3 4];
        ityp1 = 3;
      case 'QUA8',
        permut1 = [1 2 3 4 5 6 7 8];
        ityp1 = 16;
      case 'TET4',
        permut1 = [1 2 3 4];
        ityp1 = 4;
      case 'CUB8',
        permut1 = [1 2 3 4 5 6 7 8];
        ityp1 = 5;
      case 'CU20',
        permut1 = [1:8 9 12 17 10 18 11 19 20 13 16 14 15];
        ityp1 = 17;
      case 'PRI6',
        permut1 = [1 2 3 4 5 6];
        ityp1 = 6;
      case 'PR15',
        permut1 = [1:6 7 9 13 8 14 15 10 12 11];
        ityp1 = 18;
      case 'POI1',
        permut1 = [1];
        ityp1 = 15;
      case 'TE10',
 % DD 08/12/2015       permut1 = [1 3 5 10 2 4 6 7 9 8];
 % DD 08/12/2015       permut1 = [1:4 5 6 7 8 9 10];
        permut1 = [1:10];
        ityp1 = 11;
      otherwise,
        type1
        error('Element not implemented yet')
    end
    RegPhys = imesh1;
    RegElem = 1;
    nno1 = size(permut1,2);
    format1 = strcat('%d %d %d %d %d',repmat(' %d',1,nno1),'\n');
    for iel1 = 1:nel1
      for el2 = 1:size(permut1,1)
        ielt1 = ielt1 + 1;
        line1 = [ielt1,ityp1,RegPhys,RegElem,nno1,topo1(iel1,permut1(el2,:))];
        count = fprintf(fid,format1,line1);
        clear line1;
      end
    end
    clear topo1 maile1;
  end
  clear lmail1;

% End loop on meshes
  end
  count = fprintf(fid,'$ENDELM\n');

clear ListMesh2 numer1;
