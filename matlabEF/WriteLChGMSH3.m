function [error1] = WriteLChGMSH3(fid,name1,xcoor1,lt1, ...
                                   mail1,Lchamno1,listComp1, ...
                                   nmail2,Lchpo2,listComp2)
% Write a time-evolution scalar field in a gmsh ascii post-pro file
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 12 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 03 / 2005
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 10 / 2006
%   Vecteurs et tenseurs pour les champs par elements aux noeuds
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 02 / 2007
%   Elements quadratiques et format 1.4 de fichier
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 08 / 2007
%   Ajout PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2007
%   Ajout TE10
% DUREISSEIX David  L.M.G.C. THERMOMECANIQUE DES MAT le 18 / 03 / 2010
%   Change in format interpretation
%
% To only write a scalar-element field
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                           mail1,Lchamno1,LlistComp1,[],[],[])
%
% To only write a scalar/vector/tensor-point field
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                           [],[],[],nmail2,Lchpo2,,LlistComp2)
%
% To write both
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                 mail1,Lchamno1,LlistComp1,nmail2,Lchpo2,LlistComp2)
%
% Inputs
%   fid			file identifier (file must be opened, see fopen)
%   name1		name of the view
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   lt1(nbf1)		liste des pas de temps
%   mail1		maillage
%   Lchamno1{nbf1}	liste de champs par element definis aux noeuds
%   LlistComp1		Liste de listes de composantes du champ
%                       1 component for scalar point field
%                       3 components for vector point field
%                       9 components for tensor point field
%   if different lists are provided, they should be all length different
%   nmail2		maillage nuage de points (POI1) sous-tendant
%   Lchpo2{nbf1}	liste de champs par point
%   LlistComp2		Liste de listes de composantes du champ
% Optional inputs
%   ListC1{nbcomp1}	list of components to be written
% Output
%   error1		en cas d'erreur si <> 0
%
% Tous les champs de la liste doivent etre formattes identiquement,
% seules les valeurs changent.
% Only one view is saved in the file.
% The file has to be opened before, and closed after.
% For gmsh updates, see http://www.geuz.org/gmsh

error1 = 0;
impi1 = 0;
disp('DD:2010/03/18 - Due to a format change, all should be changed according to nVectorTriangles - To be done')

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('Dimension should be 3')
end

% Number of time steps
nbf1 = length(lt1);
if impi1; nbf1
end

% Node data size
% """"""""""""""
nScalarPoints = 0; nVectorPoints = 0; nTensorPoints = 0;
if ~isempty(nmail2)
  disp('Nodal values')

% Local node numbering (1 zone assumed) for node field chpo1
  nbzone2 = length(nmail2);
  if (nbzone2 ~= 1)
    error('More than 1 zone not implemented for the point fields')
  end
  zo2 = 1;
  numer2 = nmail2{zo2}.MAIL';
  localxcoor2 = xcoor1(numer2,:);
  nbpts = length(numer2);

% Lists of components
  [ListComp2,ListUnit2] = ListCompChpo2(Lchpo2{1})
  list1 = listComp2;
  if length(list1) ~= 1
    list1 
    error('long list not implemented')
  end
       nScalarPoints = nbpts;
       ListCompScalar = list1;

  disp(['Number of points for scalar field ' int2str(nScalarPoints)])
  disp(['Number of points for vector field ' int2str(nVectorPoints)])
  disp(['Number of points for tensor field ' int2str(nTensorPoints)])
end
if impi1
  nScalarPoints
  nVectorPoints
  nTensorPoints
end

% Element data size
% """""""""""""""""
nScalarElements = 0; nVectorElements = 0; nTensorElements = 0;

% For v1.2
  nScalarLines        = 0; nVectorLines        = 0; nTensorLines        = 0;
  nScalarTriangles    = 0; nVectorTriangles    = 0; nTensorTriangles    = 0;
  nScalarQuadrangles  = 0; nVectorQuadrangles  = 0; nTensorQuadrangles  = 0;
  nScalarTetrahedra   = 0; nVectorTetrahedra   = 0; nTensorTetrahedra   = 0;
  nScalarHexahedra    = 0; nVectorHexahedra    = 0; nTensorHexahedra    = 0;
  nScalarPrisms       = 0; nVectorPrisms       = 0; nTensorPrisms       = 0;
  nScalarPyramids     = 0; nVectorPyramids     = 0; nTensorPyramids     = 0;
% For v1.4
  nScalarLines2       = 0; nVectorLines2       = 0; nTensorLines2       = 0;
  nScalarTriangles2   = 0; nVectorTriangles2   = 0; nTensorTriangles2   = 0;
  nScalarQuadrangles2 = 0; nVectorQuadrangles2 = 0; nTensorQuadrangles2 = 0;
  nScalarTetrahedra2  = 0; nVectorTetrahedra2  = 0; nTensorTetrahedra2  = 0;
  nScalarHexahedra2   = 0; nVectorHexahedra2   = 0; nTensorHexahedra2   = 0;
  nScalarPrisms2      = 0; nVectorPrisms2      = 0; nTensorPrisms2      = 0;
  nScalarPyramids2    = 0; nVectorPyramids2    = 0; nTensorPyramids2    = 0;

if ~isempty(mail1)
  disp('Element values')

% Lists of components
  nbzone1 = length(mail1);
  list1   = ListTypeMesh(mail1);
  lnelt1  = Nelements2(mail1,[1:nbzone1]);
  nelt1   = Nelements2(mail1);
  [ListComp1,ListUnit1] = ListCompCham2(Lchamno1{1});
  list2 = listComp1;
  switch length(list2)
    case 1,
      disp('  Scalar values')
      nScalarElements = nelt1;
      ListCompScalar2 = list2;
    case 3,
      disp('  Vector values')
      nVectorElements = nelt1;
      ListCompVector2 = list2;
      if impi1
        nVectorElements
        ListCompVector2
      end
    case 9,
      disp('  Tensor values')
      nTensorElements = nelt1;
      ListCompTensor2 = list2;
      if impi1
        nTensorElements
        ListCompTensor2
      end
    otherwise,
      listComp1
      error('Bad number of components')
  end
  disp(['Number of elements for scalar field ' int2str(nScalarElements)])
  disp(['Number of elements for vector field ' int2str(nVectorElements)])
  disp(['Number of elements for tensor field ' int2str(nTensorElements)])

% Number of scalar values
  if nScalarElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [ListCompe1,ListUnite1] = ListCompCham2(Lchamno1{1},[izo1]);
    icomp1 = findoccur(ListCompScalar2,listComp1);
    if icomp1
    switch type1
      case 'SEG2',
        nScalarLines = nScalarLines + nbel1;
      case 'TRI3',
        nScalarTriangles = nScalarTriangles + nbel1;
      case 'TRI6',
        nScalarTriangles2 = nScalarTriangles2 + nbel1;
      case 'QUA4',
        nScalarQuadrangles = nScalarQuadrangles + nbel1;
      case 'PRI6',
        nScalarPrisms = nScalarPrisms + nbel1;
      case 'TET4',
        nScalarTetrahedra  = nScalarTetrahedra  + nbel1;
      case 'TE10',
        nScalarTetrahedra2  = nScalarTetrahedra2  + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented 1')
    end
    end
  end
  end

% Number of vector values
  if nVectorElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(Lchamno1{1},[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if impi1
      izo1
      type1
    end
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
        nVectorLines = nVectorLines + nbel1;
      case 'TRI3',
        nVectorTriangles = nVectorTriangles + nbel1;
      case 'TRI6',
        nVectorTriangles2 = nVectorTriangles2 + nbel1;
      case 'QUA4',
        nVectorQuadrangles = nVectorQuadrangles + nbel1;
      case 'PRI6',
        nVectorPrisms = nVectorPrisms + nbel1;
      case 'TET4',
        nVectorTetrahedra = nVectorTetrahedra + nbel1;
      case 'TE10',
        nVectorTetrahedra2 = nVectorTetrahedra2 + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented 2')
    end
    end
  end
  end

% Number of tensor values
  if nTensorElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(Lchamno1{1},[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if impi1
      izo1
      type1
      icomp1
    end
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
        nTensorLines = nTensorLines + nbel1;
      case 'TRI3',
        nTensorTriangles = nTensorTriangles + nbel1;
      case 'TRI6',
        nTensorTriangles2 = nTensorTriangles2 + nbel1;
      case 'QUA4',
        nTensorQuadrangles = nTensorQuadrangles + nbel1;
      case 'PRI6',
        nTensorPrisms = nTensorPrisms + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented 3')
    end
    end
  end
  end

end

% Up to now, no text provided
% """""""""""""""""""""""""""
  nText2d = 0; nText2dChars = 0; nTest3d = 0; nText3dChars = 0;

% Header
% """"""
  count = fprintf(fid,'$PostFormat\n');
  disp('Current format version is 1.4')
  count = fprintf(fid,'%g %d %d\n',[1.4 0 8]);
  count = fprintf(fid,'$EndPostFormat\n');

% Only one view up to now
% """""""""""""""""""""""
iview = 1;
  disp(['  writing view ' int2str(iview) ': ' name1])
  count = fprintf(fid,'$View\n');

% Beginning of field
  count = fprintf(fid,'%s\n',name1);
  nTimeSteps = nbf1;   % Number of time steps
  TimeStepValues = lt1; % List of time step values
  count = fprintf(fid,'%d\n',nTimeSteps);

  count = fprintf(fid,'%d %d %d\n',[nScalarPoints nVectorPoints nTensorPoints]);
  count = fprintf(fid,'%d %d %d\n',[nScalarLines        nVectorLines        nTensorLines]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTriangles    nVectorTriangles    nTensorTriangles]);
  count = fprintf(fid,'%d %d %d\n',[nScalarQuadrangles  nVectorQuadrangles  nTensorQuadrangles]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTetrahedra   nVectorTetrahedra   nTensorTetrahedra]);
  count = fprintf(fid,'%d %d %d\n',[nScalarHexahedra    nVectorHexahedra    nTensorHexahedra]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPrisms       nVectorPrisms       nTensorPrisms]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPyramids     nVectorPyramids     nTensorPyramids]);
  count = fprintf(fid,'%d %d %d\n',[nScalarLines2       nVectorLines2       nTensorLines2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTriangles2   nVectorTriangles2   nTensorTriangles2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarQuadrangles2 nVectorQuadrangles2 nTensorQuadrangles2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTetrahedra2  nVectorTetrahedra2  nTensorTetrahedra2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarHexahedra2   nVectorHexahedra2   nTensorHexahedra2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPrisms2      nVectorPrisms2      nTensorPrisms2]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPyramids2    nVectorPyramids2    nTensorPyramids2]);
  count = fprintf(fid,'%d %d %d %d\n',[nText2d nText2dChars nTest3d nText3dChars]);

  count = fprintf(fid,'%.3e ',TimeStepValues); count = fprintf(fid,'\n');

% ScalarPointValues
% """""""""""""""""
if nScalarPoints
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  localvalues = zeros(nbpts,1);
  ind1 = findoccur(ListCompScalar,ListComp2);
  for t1 = 1:nbf1
    chpo2 = Lchpo2{t1};
    localvalues(:,1) = chpo2{zo2}{ind1}.XVAL;
    for ino = 1:nbpts
      count = fprintf(fid,'%d\n',localxcoor2(ino,:));
      count = fprintf(fid,'%d\n',localvalues(ino,:));
    end
    clear chpo2 localvalues;
  end
end

% VectorPointValues
% """""""""""""""""
if nVectorPoints
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  localvalues = zeros(nbpts,3);
  ind1 = findoccur(ListCompVector,ListComp2);
  for t1 = 1:nbf1
    chpo2 = Lchpo2{t1};
    if ind1(1); localvalues(:,1) = chpo2{zo2}{ind1(1)}.XVAL; end
    if ind1(2); localvalues(:,2) = chpo2{zo2}{ind1(2)}.XVAL; end
    if ind1(3); localvalues(:,3) = chpo2{zo2}{ind1(3)}.XVAL; end
    for ino = 1:nbpts
      count = fprintf(fid,'%d\n',localxcoor2(ino,:));
      count = fprintf(fid,'%d %d %d\n',localvalues(ino,:));
    end
    clear chpo2 localvalues;
  end
end

% TensorPointValues
% """""""""""""""""
if nTensorPoints
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  localvalues = zeros(nbpts,9);
  ind1 = findoccur(ListCompTensor,ListComp2);
  for t1 = 1:nbf1
    chpo2 = Lchpo2{t1};
    for i = 1:9
      if ind1(i); localvalues(:,i) = chpo2{zo2}{ind1(i)}.XVAL; end
    end
    for ino = 1:nbpts
      count = fprintf(fid,'%d\n',localxcoor2(ino,:));
      count = fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',localvalues(ino,:));
    end
    clear chpo2 localvalues;
  end
end

% ScalarLineValue
% """""""""""""""
if nScalarLines
t1 = 1;
  chamno1 = Lchamno1{t1};

  for izo1 = 1:nbzone1
    maile1 = mail1{izo1};
    type1  = list1{izo1};
    nbel1  = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'SEG2') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chmlnoe1 = Lchamno1{t1}{izo1};
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e\n',xvalue1);
          clear xvalue1 chmlnoe1;
        end
        clear xcoore1 topo1;
      end
    end
    clear maile1;
  end
  clear chamno1;
end

% VectorLineValue 
% """""""""""""""
if nVectorLines
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if (strcmp(type1,'SEG2') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = zeros(2,3);
        for i = 1:3
          if icomp1(i)
            xvalue(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
          end
        end
        count = fprintf(fid,'%.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% TensorLineValue 
% """""""""""""""
if nTensorLines
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if (strcmp(type1,'SEG2') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = zeros(2,9);
        for i = 1:9
          if icomp1(i)
            xvalue(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
          end
        end
        count = fprintf(fid,'%.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3 %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% ScalarTriangleValue 
% """""""""""""""""""
if nScalarTriangles
t1 = 1;
  chamno1 = Lchamno1{t1};

  for izo1 = 1:nbzone1
    maile1 = mail1{izo1};
    type1  = list1{izo1};
    nbel1  = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'TRI3') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chmlnoe1 = Lchamno1{t1}{izo1};
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e\n',xvalue1);
          clear xvalue1 chmlnoe1;
        end
        clear xcoore1 topo1;
      end
    end
    clear maile1;
  end
  clear chamno1;

end

% VectorTriangleValue
% """""""""""""""""""
if nVectorTriangles

% DD 2010/03/18 New format interpretation
% Loop on elements
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    if strcmp(type1,'TRI3')
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
%       Node coordinates
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
%       Loop on time steps
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1};
          chmlnoe1 = chamno1{izo1};
          [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
          icomp1 = findoccur(ListCompVector2,listCompe1);
          if ~all(icomp1==0)
            xvalue1 = zeros(3,3); % (nbno,ncomp)
            for i = 1:3
              if icomp1(i)
                xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
              end
            end
%           Field values
            count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
            clear xvalue1 
          end
          clear chmlnoe1 chamno1;
        end
      end
    end
    clear maile1;
  end

end

% TensorTriangleValue
% """""""""""""""""""
if nTensorTriangles

  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    if strcmp(type1,'TRI3')
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1};
          chmlnoe1 = chamno1{izo1};
          [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
          icomp1 = findoccur(ListCompTensor2,listCompe1);
          if ~all(icomp1==0)
            xvalue1 = zeros(3,9); % (nbno,ncomp)
            for i = 1:9
              if icomp1(i)
                xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
              end
            end
            count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
            clear xvalue1;
          end
          clear chmlnoe1 chamno1;
        end
      end
    end
    clear maile1;
  end

end

% ScalarQuadrangleValue 
% """""""""""""""""""""
if nScalarQuadrangles
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'QUA4') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% VectorQuadrangleValue 
% """""""""""""""""""""
if nVectorQuadrangles
% disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
% DD 2016/05/19 New format interpretation
% Loop on elements
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    if strcmp(type1,'QUA4')
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
%       Node coordinates
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
%       Loop on time steps
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1};
          chmlnoe1 = chamno1{izo1};
          [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
          icomp1 = findoccur(ListCompVector2,listCompe1);
          if ~all(icomp1==0)
            xvalue1 = zeros(4,3); % (nbno,ncomp)
            for i = 1:3
              if icomp1(i)
                xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
              end
            end
%           Field values
            count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
            clear xvalue1 
          end
          clear chmlnoe1 chamno1;
        end
          end
    end
    clear maile1;
  end
  
end

% TensorQuadrangleValue 
% """""""""""""""""""""
if nTensorQuadrangles
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if (strcmp(type1,'QUA4') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = zeros(4,9); % (nbno,ncomp)
        for i = 1:9
          if icomp1(i)
            xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
          end
        end
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% ScalarTetrahedronValue 
% """"""""""""""""""""""
if nScalarTetrahedra
disp('Il faudra decider une fois pour toutes si on met')
disp('  1. les elements, 2. les pas de temps')
disp('  ou le contraire')
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    t1 = 1; chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'TET4') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e\n',xvalue1);
          clear xvalue1 chmlnoe1 chamno1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

% VectorTetrahedronValue 
% """"""""""""""""""""""
if nVectorTetrahedra
disp('Il faudra decider une fois pour toutes si on met')
disp('  1. les elements, 2. les pas de temps')
disp('  ou le contraire')
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    t1 = 1; chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if (strcmp(type1,'TET4') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
          xvalue1 = zeros(4,3); % (nbno,ncomp)
          for i = 1:3
            if icomp1(i)
              xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
          clear xvalue1 chmlnoe1 chamno1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

% TensorTetrahedronValue 
% """"""""""""""""""""""
if nTensorTetrahedra
  error('Not yet implemented 3')
end

% ScalarHexahedronValue 
% """""""""""""""""""""
if nScalarHexahedra
  error('Not yet implemented 4')
end

% VectorHexahedronValue 
% """""""""""""""""""""
if nVectorHexahedra
  error('Not yet implemented 5')
end

% TensorHexahedronValue 
% """""""""""""""""""""
if nTensorHexahedra
  error('Not yet implemented 6')
end

% ScalarPrismValue 
% """"""""""""""""
if nScalarPrisms
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'PRI6') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% VectorPrismValue 
% """"""""""""""""
if nVectorPrisms
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if (strcmp(type1,'PRI6') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = zeros(6,3); % (nbno,ncomp)
        for i = 1:3
          if icomp1(i)
            xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
          end
        end
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% TensorPrismValue 
% """"""""""""""""
if nTensorPrisms
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
for t1 = 1:nbf1
  chamno1 = Lchamno1{t1};
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if (strcmp(type1,'PRI6') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = zeros(6,9); % (nbno,ncomp)
        for i = 1:9
          if icomp1(i)
            xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
          end
        end
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
        clear xvalue1 xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end
end

% ScalarPyramidValue 
% """"""""""""""""""
if nScalarPyramids
  error('Not yet implemented 7')
end

% VectorPyramidValue 
% """"""""""""""""""
if nVectorPyramids
  error('Not yet implemented 8')
end

% TensorPyramidValue 
% """"""""""""""""""
if nTensorPyramids
  error('Not yet implemented 9')
end

% Lines2Value 
% """""""""""
if nScalarLines2
  error('SEG3 not yet implemented... 1')
end
if nVectorLines2
  error('SEG3 not yet implemented... 2')
end
if nTensorLines2
  error('SEG3 not yet implemented... 3')
end

% ScalarTriangles2Value 
% """""""""""""""""""""
if nScalarTriangles2
t1 = 1;
  chamno1 = Lchamno1{t1};

  for izo1 = 1:nbzone1
    maile1 = mail1{izo1};
    type1  = list1{izo1};
    nbel1  = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'TRI6') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chmlnoe1 = Lchamno1{t1}{izo1};
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e\n',xvalue1);
          clear xvalue1 chmlnoe1;
        end
        clear xcoore1 topo1;
      end
    end
    clear maile1;
  end
  clear chamno1;
end

% VectorTriangles2Value 
% """""""""""""""""""""
if nVectorTriangles2
t1 = 1;
  chamno1 = Lchamno1{t1};

  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if (strcmp(type1,'TRI6') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chmlnoe1 = Lchamno1{t1}{izo1};
          xvalue1 = zeros(6,3); % (nbno,ncomp)
          for i = 1:3
            if icomp1(i)
              xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
          clear xvalue1 chmlnoe1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

% TensorTriangles2Value 
% """""""""""""""""""""
if nTensorTriangles2
t1 = 1;
  chamno1 = Lchamno1{t1};

  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if (strcmp(type1,'TRI6') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chmlnoe1 = Lchamno1{t1}{izo1};
          xvalue1 = zeros(6,9); % (nbno,ncomp)
          for i = 1:9
            if icomp1(i)
              xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
          clear xvalue1 chmlnoe1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

if nScalarQuadrangles2
  error('Not yet done...')
end
if nVectorQuadrangles2
  error('Not yet done...')
end
if nTensorQuadrangles2
  error('Not yet done...')
end

if nScalarTetrahedra2 
disp('nScalarTetrahedra2: Il faudra decider une fois pour toutes si on met')
disp('  1. les elements, 2. les pas de temps')
disp('  ou le contraire')
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    t1 = 1; chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if (strcmp(type1,'TE10') && icomp1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e\n',xvalue1);
          clear xvalue1 chmlnoe1 chamno1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

if nVectorTetrahedra2 
disp('nScalarTetrahedra2: Il faudra decider une fois pour toutes si on met')
disp('  1. les elements, 2. les pas de temps')
disp('  ou le contraire')
disp('DD:2010/03/18 - Due to a format change, this should be rewritten')
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    t1 = 1; chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if (strcmp(type1,'TE10') && ~all(icomp1==0))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        for t1 = 1:nbf1
          chamno1 = Lchamno1{t1}; chmlnoe1 = chamno1{izo1};
          xvalue1 = zeros(10,3); % (nbno,ncomp)
          for i = 1:3
            if icomp1(i)
              xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
          clear xvalue1 chmlnoe1 chamno1;
        end
        clear xcoore1 topo1;
      end
    end
    clear chmlnoe1 maile1;
  end
  clear chamno1;
end

if nTensorTetrahedra2
  error('Not yet done...')
end
if nScalarHexahedra2  
  error('Not yet done...')
end
if nVectorHexahedra2  
  error('Not yet done...')
end
if nTensorHexahedra2
  error('Not yet done...')
end
if nScalarPrisms2     
  error('Not yet done...')
end
if nVectorPrisms2     
  error('Not yet done...')
end
if nTensorPrisms2
  error('Not yet done...')
end
if nScalarPyramids2   
  error('Not yet done...')
end
if nVectorPyramids2   
  error('Not yet done...')
end
if nTensorPyramids2
  error('Not yet done...')
end

% text2d  text2d-chars 
% text3d  text3d-chars 

% End of field
  count = fprintf(fid,'$EndView\n');
