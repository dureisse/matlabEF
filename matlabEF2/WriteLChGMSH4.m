function [error1] = WriteLChGMSH4(fid,xcoor1,lt1, ...
                                  mail1,Lchamno1,listComp1,name1, ...
                                  nmail2,Lchpo2,listComp2,name2)
% Write a time-evolution field in a gmsh ascii msh file (version 4.1)
% 2 possible fields: 
% either nodal ones, or fields defined by elements at nodes
%
% DUREISSEIX David  LaMCoS                     le 01 / 07 / 2022
% DUREISSEIX David  LaMCoS                     le 13 / 10 / 2022
%   Correction de bug pour la renumérotation, maintenant tous les 
%   maillages sont renumérotés au début
%
% To write a scalar/vector/tensor time nodal or element-wise field 
%   fid = fopen(file1,'w'); 
%     error1 = WriteLChGMSH4(fid,xcoor1,lt1, ...
%                            mail1,Lchamno1,listComp1,name1, ...
%                            nmail2,Lchpo2,listComp2,name2)
%   fclose(fid)
%
% Inputs
%   fid			    file identifier (file must be opened, see fopen)
%   xcoor1(nbno,3)	3D coordinates
%   lt1(nbf1)		liste of time steps
%   mail1		    mesh
%   Lchamno1{nbf1}  liste des champs par element définis aux noeuds
%   listComp1		liste de composantes du champ par point
%   name1		    name of the view for the nodal data
%   nmail2		    maillage nuage de points (POI1) sous-tendant
%   Lchpo2{nbf1}	liste de champs par point
%   listComp2		liste de composantes du champ par element
%                   (1 : champ scalaire, 3 : champ vectoriel en 3D, 
%                    9 : champ de tenseur en 3D)
%                        sigx sigxy sigxz sigxy sigy sigyz sigxz sigyz sigz 
%   name2		    name of the view for the element data
% Output
%   error1		    en cas d'erreur si <> 0
%
% Tous les champs de la liste doivent etre formattes identiquement,
% seules les valeurs changent.
% The file has to be opened before, and closed after.
% For gmsh updates, see http://www.geuz.org/gmsh

error1 = 0;
impi1 = 0;

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('Dimension should be 3')
end

% Number of time steps
nbf1 = length(lt1);
if impi1; nbf1
end

% DD 2022/10/13
% Pre-processing : renumerotation des maillages
% """""""""""""""""""""""""""""""""""""""""""""
mail2 = ChangeMesh2(mail1,'POI1');
numer1 = mail2{1}.MAIL'; 
numer2 = nmail2{1}.MAIL';
numer2 = unique([numer1 numer2]); % renumerotation
clear numer1 mail2;
numer2_inv = InverseList(numer2,max(numer2)); % liste inverse
mail1 = RenumMesh(mail1,numer2_inv);
nmail2 = RenumMesh(nmail2,numer2_inv);
localxcoor2 = xcoor1(numer2,:); % corresponding coordinates

% Mesh mail1
% """"
  nbzone1 = length(mail1);
  if (nbzone1 ~= 1)
    error('More than 1 zone not yet implemented for the mesh')
  end
  nelt1 = Nelements2(mail1); % total number of elements
  maile1 = mail1{1}; % only 1 zone assumed
  topo1 = maile1.MAIL;
  type1 = maile1.TYPE; 
  % for elementType in gmsh see https://gitlab.onelab.info/gmsh/gmsh/blob/master/src/common/GmshDefines.h
  % for node ordering in gmsh see https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
  % ityp1 gmsh type of element
  % permut1 node permutation to match matlabEF/AVS and gmsh node ordering
  % topo_matlabEF(:,permut1) = topo_gmsh(:,1:nbno)
  switch type1
      case 'SEG2',
        permut1 = [1 2];
        ityp1 = 1;
        entityDim = 1;
      case 'SEG3',
        permut1 = [1 2 3];
        ityp1 = 8;
        entityDim = 1;
      case 'TRI3',
        permut1 = [1 2 3];
        ityp1 = 2;
        entityDim = 2;
      case 'TRI6',
        permut1 = [1 2 3 4 5 6];
        ityp1 = 9;
        entityDim = 2;
      case 'QUA4',
        permut1 = [1 2 3 4];
        ityp1 = 3;
        entityDim = 2;
      case 'QUA8',
        permut1 = [1 2 3 4 5 6 7 8];
        ityp1 = 16;
        entityDim = 2;
      case 'TET4',
        permut1 = [1 2 3 4];
        ityp1 = 4;
        entityDim = 3;
      case 'CUB8',
        permut1 = [1 2 3 4 5 6 7 8];
        ityp1 = 5;
        entityDim = 3;
      case 'CU20',
        permut1 = [1:8 9 12 17 10 18 11 19 20 13 16 14 15];
        ityp1 = 17;
        entityDim = 3;
      case 'PRI6',
        permut1 = [1 2 3 4 5 6];
        ityp1 = 6;
        entityDim = 3;
      case 'PR15',
        permut1 = [1:6 7 9 13 8 14 15 10 12 11];
        ityp1 = 18;
        entityDim = 3;
      case 'POI1',
        permut1 = [1];
        ityp1 = 15;
        entityDim = 0; % to be tested
      case 'TE10',
        permut1 = [1:10];
        ityp1 = 11;
        entityDim = 3;
      otherwise,
        type1
        error('Element not implemented yet')
    end

% Node data size
% """"""""""""""
nScalarPoints = 0; nVectorPoints = 0; nTensorPoints = 0;
if ~isempty(nmail2)
  disp('Nodal values')

% Local node numbering (1 zone assumed) for nodal field chpo1
  nbzone2 = length(nmail2);
  if (nbzone2 ~= 1)
    error('More than 1 zone not implemented for the point fields')
  end
  zo2 = 1;
  numer2 = nmail2{zo2}.MAIL'; % node ordering
%  localxcoor2 = xcoor1(numer2,:); % corresponding coordinates
  nbpts = length(numer2); % total number of nodes

% Lists of components
  [LListComp2,LListUnit2] = ListCompChpo2(Lchpo2{1});
  ind2 = findoccur(listComp2,LListComp2); % listComp2(i) = LListComp2(ind2(i))
  switch length(ind2)
      case 1,
          nScalarPoints = nbpts;
      case 3,
          nVectorPoints = nbpts;
      case 9,
          nTensorPoints = nbpts;
      otherwise
          listComp2
          LListComp2
          error('Not a correct list of nodal components')
  end
  disp(['Number of points for scalar field ' int2str(nScalarPoints)])
  disp(['Number of points for vector field ' int2str(nVectorPoints)])
  disp(['Number of points for tensor field ' int2str(nTensorPoints)])
end

% Element data size
% """""""""""""""""
nScalarElements = 0; nVectorElements = 0; nTensorElements = 0;
if ~isempty(Lchamno1)
  disp('Element values')

% Lists of components for the element data
  nbzone1 = length(mail1);
  list1   = ListTypeMesh(mail1);
  lnelt1  = Nelements2(mail1,[1:nbzone1]);
  nelt1   = Nelements2(mail1);
  [LListComp1,LListUnit1] = ListCompCham2(Lchamno1{1});
  ind1 = findoccur(listComp1,LListComp1); % listComp1(i) = LListComp1(ind1(i))
  switch length(ind1)
      case 1,
          nScalarElements = nelt1;
      case 3,
          nVectorElements = nelt1;
      case 9,
          nTensorElements = nelt1;
      otherwise
          listComp1
          LListComp1
          error('Not a correct list of element components')
  end
  disp(['Number of elements for scalar field ' int2str(nScalarElements)])
  disp(['Number of elements for vector field ' int2str(nVectorElements)])
  disp(['Number of elements for tensor field ' int2str(nTensorElements)])
end


% Header
% """"""
  count = fprintf(fid,'$MeshFormat\n');
  disp('Current msh format version is 4.1')
  count = fprintf(fid,'%g %d %d\n',[4.1 0 8]);
  count = fprintf(fid,'$EndMeshFormat\n');

% Nodes
% """""
  count = fprintf(fid,'$Nodes\n');
%   numEntityBlocks = 1;
%   numNodes = nbpts;
%   minNodeTag = 1;
%   maxNodeTag = nbpts;
  count = fprintf(fid,'%d %d %d %d\n',[1 nbpts 1 nbpts]);
%   entityDim
%   entityTag = 1 ??
%   parametric = 0;
%   numNodesInBlock = nbpts;
  count = fprintf(fid,'%d %d %d %d\n',[entityDim 1 0 nbpts]);
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',ino);
  end
  for ino = 1:nbpts
    count = fprintf(fid,'%d %d %d\n',[localxcoor2(ino,:)]);
  end
  count = fprintf(fid,'$EndNodes\n');

% Elements
% """"""""
  count = fprintf(fid,'$Elements\n');
%   numEntityBlocks = 1;
%   numElements = nelt1;
%   minElementTag = 1;
%   maxElementTag = nelt1;
  count = fprintf(fid,'%d %d %d %d\n',[1 nelt1 1 nelt1]);
%   entityDim 
%   entityTag = 1 ??
%   elementType 
%   numElementsInBlock = nelt1;
  count = fprintf(fid,'%d %d %d %d\n',[entityDim 1 ityp1 nelt1]);
  nno1 = length(permut1);
  format1 = strcat('%d',repmat(' %d',1,nno1),'\n');
  for iel1 = 1:nelt1
      count = fprintf(fid,format1,[iel1 topo1(iel1,permut1)]);
  end
  count = fprintf(fid,'$EndElements\n');

% Loop on fields and Node Data
% """"""""""""""""""""""""""""
  if ~isempty(Lchpo2)
    for it1 = 1:nbf1
      time1 = lt1(it1);
      chpo1 = Lchpo2{it1};
      chpoe1 = chpo1{1}; % 1 zone
%      xval1 = chpoe1{1}.XVAL; % valeurs de la premiere composante (scalaire)
      XXVAL1 = zeros(nbpts,length(ind2));
      for ii = 1:length(ind2)
          XXVAL1(:,ii) = chpoe1{ind2(ii)}.XVAL;
      end
      format2 = strcat('%d',repmat(' %d',1,length(ind2)),'\n');

      count = fprintf(fid,'$NodeData\n');
      count = fprintf(fid,'%d\n',[1]);
      count = fprintf(fid,'"%s"\n',name2); % name of the view
      count = fprintf(fid,'%d\n',[1]);
      count = fprintf(fid,'%d\n',[time1]); % time step value
      count = fprintf(fid,'%d\n',[3]);
      count = fprintf(fid,'%d\n',[it1-1]); % time steps number
      count = fprintf(fid,'%d\n',[length(ind2)]); % nb of components (1,3,9)
      count = fprintf(fid,'%d\n',[nbpts]); % nb of nodes
      for ino = 1:nbpts
%        count = fprintf(fid,'%d %d\n',[ino xval1(ino,:)]);
        count = fprintf(fid,format2,[ino XXVAL1(ino,:)]);
      end
      count = fprintf(fid,'$EndNodeData\n');

      clear XXVAL1 chpoe1 chpo1;
    end
  end

% Loop on fields and Element Data
% """""""""""""""""""""""""""""""
if ~isempty(Lchamno1)
%   Boucle sur les pas de temps
    for it1 = 1:nbf1
      time1 = lt1(it1);
      chamno1 = Lchamno1{it1};
      chamnoe1 = chamno1{1}; % 1 zone
%      xval1 = chamnoe1{1}.XVAL; % (nb element,nb nodes per element)
%      valeurs de la premiere composante (scalaire), 
%     XXVAL1(nbel,nbnoel*nbcomp)
      XXVAL1 = zeros(nelt1,nno1*length(ind1)); 
% %       ii = 1;
% %           PPERMUT1 = [permut1(ii):nno1:nno1*length(ind1)];
% %       for ii = 2:nno1
% %           PPERMUT1 = [PPERMUT1 permut1(ii):nno1:nno1*length(ind1)];
% %           % 1:nno1:end 2:nno1:end ... length(ind1):nno1:end
% %           % mais avec la permutation des noeuds permut1
% %       end
%     Boucle sur les composantes
      for ii = 1:length(ind1)
%         XX1(nbel,nbnoel) nod1valii nod2valii ...
          XX1 = chamnoe1{ind1(ii)}.XVAL;
%         on permute les noeuds matlabEF/AVS vers gmsh
%         topo_matlabEF(:,permut1) = topo_gmsh(:,1:nbno)
          XX2 = XX1(:,permut1);
%         on range les noeuds/composantes
%         ligne = nod1val1 nod1val2 nod3val1 nod2val1 nod2val2 ...
          PPERMUT2 = [ii:length(ind1):nno1*length(ind1)];
          XXVAL1(:,PPERMUT2) = XX2;
          clear XX1 XX2;
      end
%      format1 = strcat('%d',repmat(' %d %d',1,1*nnoel),'\n'); % scalaire 1 composante
      format1 = strcat('%d %d',repmat(' %d',1,nno1*length(ind1)),'\n');

      count = fprintf(fid,'$ElementNodeData\n');
      count = fprintf(fid,'%d\n',[1]);
      count = fprintf(fid,'"%s"\n',name1); % name of the view
      count = fprintf(fid,'%d\n',[1]);
      count = fprintf(fid,'%d\n',[time1]); % time step value
      count = fprintf(fid,'%d\n',[4]);
      count = fprintf(fid,'%d\n',[it1-1]); % time step index
      count = fprintf(fid,'%d\n',[length(ind1)]); % number of components (scalar : 1)
      count = fprintf(fid,'%d\n',[nelt1]); % nb of elements
      count = fprintf(fid,'%d\n',[0]); % partition index (no partition)
      nnoel = length(permut1);
      for iel = 1:nelt1
        count = fprintf(fid,format1,[iel nno1 XXVAL1(iel,:)]);
      end
      count = fprintf(fid,'$EndElementNodeData\n');

      clear chamno1 chamnoe1 XXVAL1;    
    end
end

% End of file
% """""""""""
  count = fprintf(fid,'$EndView\n');
