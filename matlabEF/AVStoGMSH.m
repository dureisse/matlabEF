function [error1] = AVStoGMSH(avsfile1,nrec1,gmshfile1,varargin)
%  Traduce one record in an AVS UCD file to GMSH files (msh, pos)
%  Only scalar components are managed for the various fields
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 07 / 2006
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 03 / 2007
%  Allows to read a particular record irec1
%
% Inputs
%   avsfile1		Name of the avs file to be opened
%   nrec1		Number of records to be read in the avs file
%   gmshfile1		Root name of the gmsh file to be written
% Optional input
%   irec1		Record to manage (default: 1)
% Outputs
%   error1		Integer (0 if no problem)
%
% error1 = AVStoGMSH(avsfile1,nrec1,gmshfile1,varargin)
%
% For fields, up to now, each component is written as an independant
% scalar field (not yet vector nor tensor components)

narg = nargin - 3;
switch narg
  case 0,
    irec1 = 1;
  case 1,
    irec1 = varargin{1};
  otherwise,
    error(['Bad number of optional inputs ' int2str(narg)])
end

% Read the avs file
xcrit1 = -1.; % No node merging
[xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
          ListChml1,ListCara1,error1] = ReadMergeAVS2(xcrit1,avsfile1,nrec1);
% xcoort1(nbnot1,idim)  Node coordinates
% ListMesh1{nrec1}      List of meshes
% ListChpo1{nrec1}	List of nodal fields
% ListnMesh1{nrec1}	List of their clouds of nodes
% ListChml1{nrec1}	List of elemental fields
% ListCara1{nrec1}	list of characteristics

if error1
  error(['Error in avs file reading ' int2str(error1)])
end

% Mesh is saved in gmshfile1.msh file
% """""""""""""""""""""""""""""""""""
mail2 = ListMesh1{irec1};

if ~isempty(mail2)
  clear ListMesh2;
  ListMesh2{1} = mail2;

  gmshfilem1 = strcat(gmshfile1,'.msh');
  disp(['Writing Mesh in ' gmshfilem1])
  fid = fopen(gmshfilem1,'w');
    if (fid == -1)
      error(['Error in opening mesh file ' gmshfilem1])
    end
    error1 = WriteMeshGMSH(xcoort1,ListMesh2,fid);
    if error1
      error(['Error in writing mesh ' int2str(error1)])
    end
  fclose(fid);

  clear ListMesh2;
end


% Node-based fields are saved in gmshfile1_chpo.pos file
% """"""""""""""""""""""""""""""""""""""""""""""""""""""

chpo1 = ListChpo1{irec1};
nbzo0 = length(chpo1);
if nbzo0
  nmail1 = ListnMesh1{irec1};
  if ~isempty(mail2)
    disp('Converting node-based field')
    disp('into element-based field defined at nodes')
    [chamno1,intg1] = ChpoToChamno3(chpo1,nmail1,mail2);
  else

    chamno1 = [];

    gmshfilep0 = strcat(gmshfile1,'_chpo.pos');
    disp(['Writing node-based field in ' gmshfilep0])
    fid = fopen(gmshfilep0,'w');
      if (fid == -1)
        error(['Error in opening node-based post file ' gmshfilep0])
      end
      [listComp1,listUnit1] = ListCompChpo2(chpo1,[1:nbzo0]);
      error1 = WriteChpoGMSH(xcoort1,nmail1,chpo1,fid, ...
                             'AVStoGMSH_chpo',listComp1);
      if error1
        error(['Error in writing node-based field ' int2str(error1)])
      end
      clear listComp1 listUnit1;
    fclose(fid)

  end
else
  chamno1 = [];
end
clear chpo1;

% Element-based fields are saved in gmshfile1.pos file
% """"""""""""""""""""""""""""""""""""""""""""""""""""

%DD 08/06/17 chml1 = ListChml1{1};
chml1 = ListChml1{irec1};
if (length(chml1) ~= 0)
  disp('Converting element-based uniform field')
  disp('into element-based field defined at nodes')
%DD 08/06/17   mail2 = ListMesh1{1};
  mail2 = ListMesh1{irec1};
%  intg1 = SegmentIntgNo(mail2);
  chamno2 = ChmlToChamno(chml1,mail2);
else
  chamno2 = [];
end
clear chml1;

chamno3 = FusChamno(chamno1,chamno2);

if ~isempty(chamno3)
  gmshfilef1 = strcat(gmshfile1,'.pos');
  disp(['Writing element-oriented fields in ' gmshfilef1])
  fid = fopen(gmshfilef1,'w');
    if (fid == -1)
      error(['Error in opening field file ' gmshfilef1])
    end
    error1 = WriteChamnoGMSH(xcoort1,mail2,chamno3,fid);
    if error1
      error(['Error in writing element-based field ' int2str(error1)])
    end
  fclose(fid)
end

clear chamno3;
clear mail2;
