function [error1] = WriteChGMSH2(fid,name1,xcoor1, ...
                                mail1,chamno1,LlistComp1, ...
                                nmail2,chpo2,LlistComp2)
% Write a field in a gmsh ascii post-pro file
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 12 / 2004
%
% To only write a scalar-element field
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                           mail1,chamno1,LlistComp1,[],[],[])
% To only write a scalar/vector/tensor-point field
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                           [],[],[],nmail2,chpo2,,LlistComp2)
% To write both
%   error1 = WriteChpoGMSH2(fid,name1,xcoor1, ...
%                 mail1,chamno1,LlistComp1,nmail2,chpo2,LlistComp2)
%
% Inputs
%   fid			file identifier (file must be opened, see fopen)
%   name1		name of the view
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   mail1		maillage
%   chamno1		champ par element defini aux noeuds
%   LlistComp1		Liste de listes de composantes pour chamno1
%                       1 component for scalar point field
%                       3 components for vector point field
%                       9 components for tensor point field
%   if different lists are provided, they should be all length different
%   nmail2		maillage nuage de points (POI1) sous-tendant
%   chpo2		champ par point
%   LlistComp2		Liste de listes de composantes pour chpo2
% Optional inputs
%   ListC1{nbcomp1}	list of components to be written
% Output
%   error1		en cas d'erreur si <> 0
%
% Only one view is saved in the file.
% Only one time step is considered (if not, there should be a list of
% fields as arguments, with a list of time step values, and the routine
% should be WriteLchGMSH... most of the content should be shared with
% this routine... to be done)
% The file has to be opened before, and closed after.
% For gmsh updates, see http://www.geuz.org/gmsh

error1 = 0;

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('Dimension should be 3')
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
  [ListComp2,ListUnit2] = ListCompChpo2(chpo2); 
  nlist1 = length(LlistComp2);
  for ilist1 = 1:nlist1
    list1 = LlistComp2{ilist1};
    switch length(list1)
      case 1,
       if nScalarPoints
         error('Scalar list of component already provided')
       end
       nScalarPoints = nbpts;
       ListCompScalar = list1;
      case 3,
       if nVectorPoints
         error('Vector list of component already provided')
       end
       nVectorPoints = nbpts;
       ListCompVector = list1;
      case 9,
       if nTensorPoints
         error('Tensor list of component already provided')
       end
       nTensorPoints = nbpts;
       ListCompTensor = list1;
      otherwise,
        list1
        error('Not a recognized list of comp. (bad number of comp.)')
    end
  end
  disp(['Number of points for scalar field ' int2str(nScalarPoints)])
  disp(['Number of points for vector field ' int2str(nVectorPoints)])
  disp(['Number of points for tensor field ' int2str(nTensorPoints)])
end

% Element data size
% """""""""""""""""
nScalarElements = 0; nVectorElements = 0; nTensorElements = 0;
nScalarLines = 0; nVectorLines = 0; nTensorLines = 0;
nScalarTriangles = 0; nVectorTriangles = 0; nTensorTriangles = 0;
nScalarQuadrangles = 0; nVectorQuadrangles = 0; nTensorQuadrangles = 0;
nScalarTetrahedra = 0; nVectorTetrahedra = 0; nTensorTetrahedra = 0;
nScalarHexahedra = 0; nVectorHexahedra = 0; nTensorHexahedra = 0;
nScalarPrisms = 0; nVectorPrisms = 0; nTensorPrisms = 0;
nScalarPyramids = 0; nVectorPyramids = 0; nTensorPyramids = 0;
if ~isempty(mail1)
  disp('Element values')

% Lists of components
  nbzone1 = length(mail1);
  list1   = ListTypeMesh(mail1);
  lnelt1  = Nelements2(mail1,[1:nbzone1]);
  nelt1   = Nelements2(mail1);
  [listComp1,listUnit1] = ListCompCham2(chamno1);
  nlist2 = length(LlistComp1);
  for ilist2 = 1:nlist2
    list2 = LlistComp1{ilist2};
    switch length(list2)
      case 1,
       if nScalarElements
         error('Scalar list of component already provided')
       end
       nScalarElements = nelt1;
       ListCompScalar2 = list2;
      case 3,
       if nVectorElements
         error('Vector list of component already provided')
       end
       nVectorElements = nelt1;
       ListCompVector2 = list2;
      case 9,
       if nTensorElements
         error('Tensor list of component already provided')
       end
       nTensorElements = nelt1;
       ListCompTensor2= list2;
      otherwise,
        list2
        error('Not a recognized list of comp.-2 (bad number of comp.)')
    end
  end
  disp(['Number of elements for scalar field ' int2str(nScalarElements)])
  disp(['Number of elements for vector field ' int2str(nVectorElements)])
  disp(['Number of elements for tensor field ' int2str(nTensorElements)])

% Number of scalar values
  if nScalarElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if icomp1
    switch type1
      case 'SEG2',
        nScalarLines = nScalarLines + nbel1;
      case 'TRI3',
        nScalarTriangles = nScalarTriangles + nbel1;
      case 'TRI6',
error('TRI6 do not have to be subdivised into TRI3... to be changed')
%       subtriangles
        TRI6toTRI3 = [1 4 6
                      4 2 5
                      6 5 3
                      4 5 6];
        nScalarTriangles = nScalarTriangles + size(TRI6toTRI3,1)*nbel1;
      case 'QUA4',
        nScalarQuadrangles = nScalarQuadrangles + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented')
    end
    end
  end
  disp(['Number of lines for scalar field ' int2str(nScalarLines)])
  disp(['Number of triangles for scalar field ' int2str(nScalarTriangles)])
  disp(['Number of quadrangles for scalar field ' int2str(nScalarQuadrangles)])
  end

% Number of vector values
  if nVectorElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
        nVectorLines = nVectorLines + nbel1;
      case 'TRI3',
        nVectorTriangles = nVectorTriangles + nbel1;
      case 'TRI6',
%       subtriangles
        TRI6toTRI3 = [1 4 6
                      4 2 5
                      6 5 3
                      4 5 6];
        nVectorTriangles = nVectorTriangles + size(TRI6toTRI3,1)*nbel1;
      case 'QUA4',
        nVectorQuadrangles = nVectorQuadrangles + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented')
    end
    end
  end
  disp(['Number of lines for vector field ' int2str(nVectorLines)])
  disp(['Number of triangles for vector field ' int2str(nVectorTriangles)])
  disp(['Number of quadrangles for vector field ' int2str(nVectorQuadrangles)])
  end

% Number of tensor values
  if nTensorElements
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
        nTensorLines = nTensorLines + nbel1;
      case 'TRI3',
        nTensorTriangles = nTensorTriangles + nbel1;
      case 'TRI6',
%       subtriangles
        TRI6toTRI3 = [1 4 6
                      4 2 5
                      6 5 3
                      4 5 6];
        nTensorTriangles = nTensorTriangles + size(TRI6toTRI3,1)*nbel1;
      case 'QUA4',
        nTensorQuadrangles = nTensorQuadrangles + nbel1;
      otherwise,
        type1
        error('Element type not yet implemented')
    end
    end
  end
  disp(['Number of lines for tensor field ' int2str(nTensorLines)])
  disp(['Number of triangles for tensor field ' int2str(nTensorTriangles)])
  disp(['Number of quadrangles for tensor field ' int2str(nTensorQuadrangles)])
  end

end

% Up to now, no text provided
% """""""""""""""""""""""""""
  nText2d = 0; nText2dChars = 0; nTest3d = 0; nText3dChars = 0;

% Header
% """"""
  count = fprintf(fid,'$PostFormat\n');
  count = fprintf(fid,'%g %d %d\n',[1.2 0 8]);
  count = fprintf(fid,'$EndPostFormat\n');

% Only one view up to now, and one time step
% """"""""""""""""""""""""""""""""""""""""""
iview = 1;
  disp(['  writing view ' int2str(iview) ': ' name1])
  count = fprintf(fid,'$View\n');

% Beginning of field
  count = fprintf(fid,'%s\n',name1);
  nTimeSteps = 1;      % Number of time steps
  TimeStepValues = 1.; % List of time step values
  count = fprintf(fid,'%d\n',nTimeSteps);


  count = fprintf(fid,'%d %d %d\n',[nScalarPoints nVectorPoints nTensorPoints]);
  count = fprintf(fid,'%d %d %d\n',[nScalarLines nVectorLines nTensorLines]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTriangles nVectorTriangles nTensorTriangles]);
  count = fprintf(fid,'%d %d %d\n',[nScalarQuadrangles nVectorQuadrangles nTensorQuadrangles]);
  count = fprintf(fid,'%d %d %d\n',[nScalarTetrahedra nVectorTetrahedra nTensorTetrahedra]);
  count = fprintf(fid,'%d %d %d\n',[nScalarHexahedra nVectorHexahedra nTensorHexahedra]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPrisms nVectorPrisms nTensorPrisms]);
  count = fprintf(fid,'%d %d %d\n',[nScalarPyramids nVectorPyramids nTensorPyramids]);
  count = fprintf(fid,'%d %d %d %d\n',[nText2d nText2dChars nTest3d nText3dChars]);

  count = fprintf(fid,'%.3e\n',TimeStepValues);

% ScalarPointValues
% """""""""""""""""
if nScalarPoints
  localvalues = zeros(nbpts,1);
  ind1 = findoccur(ListCompScalar,ListComp2);
  localvalues(:,1) = chpo2{zo2}{ind1}.XVAL;
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor2(ino,:));
    count = fprintf(fid,'%d\n',localvalues(ino,:));
  end
end

% VectorPointValues
% """""""""""""""""
if nVectorPoints
  localvalues = zeros(nbpts,3);
  ind1 = findoccur(ListCompVector,ListComp2);
  if ind1(1); localvalues(:,1) = chpo2{zo2}{ind1(1)}.XVAL; end
  if ind1(2); localvalues(:,2) = chpo2{zo2}{ind1(2)}.XVAL; end
  if ind1(3); localvalues(:,3) = chpo2{zo2}{ind1(3)}.XVAL; end
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor2(ino,:));
    count = fprintf(fid,'%d %d %d\n',localvalues(ino,:));
  end
end

% TensorPointValues
% """""""""""""""""
if nTensorPoints
  localvalues = zeros(nbpts,9);
  ind1 = findoccur(ListCompTensor,ListComp2);
  for i = 1:9
    if ind1(i); localvalues(:,i) = chpo2{zo2}{ind1(i)}.XVAL; end
  end
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor2(ino,:));
    count = fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',localvalues(ino,:));
  end
end

% ScalarLineValue
% """""""""""""""
if nScalarLines
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar2,listCompe1);
    if icomp1
    switch type1
      case 'SEG2',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e\n',xvalue1);
        end
      case 'SEG3',
        error('SEG3 not yet implemented')
      otherwise
        error('not yet implemented')
    end
    end
  end
end

% VectorLineValue 
% """""""""""""""
if nVectorLines
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector2,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
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
        end
      case 'SEG3',
        error('SEG3 not yet implemented')
      otherwise
        error('not yet implemented')
    end
    end
  end
end

% TensorLineValue 
% """""""""""""""
if nTensorLines
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'SEG2',
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
        end
      case 'SEG3',
        error('SEG3 not yet implemented')
      otherwise
        error('not yet implemented')
    end
    end
  end
end

% ScalarTriangleValue 
% """""""""""""""""""
if nScalarTriangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompScalar,listCompe1);
    if icomp1
    switch type1
      case 'TRI3',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e\n',xvalue1);
        end
      case 'TRI6',
        for el1 = 1:nbel1
          topo1   = maile1.MAIL(el1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          for el2 = 1:size(TRI6toTRI3,1)
            topo2   = topo1(1,TRI6toTRI3(el2,:));
            xcoore2 = xcoor1(topo2,:);
            xvalue2 = xvalue1(1,TRI6toTRI3(el2,:));
            count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore2);
            count = fprintf(fid,'%.3e\n',xvalue2);
          end
        end
    end
    end
  end
end

% VectorTriangleValue
% """""""""""""""""""
if nVectorTriangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'TRI3',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = zeros(3,3); % (nbno,ncomp)
          for i = 1:3
            if icomp1(i)
              xvalue(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
        end
      case 'TRI6',
        for el1 = 1:nbel1
          topo1   = maile1.MAIL(el1,:);
          xvalue1 = zeros(3,6); % (nbcomp,nbno)
          for i = 1:3
            if icomp1(i)
              xvalue(i,:) = chmlnoe1{icomp1(i)}.XVAL(el1,:);
            end
          end
          for el2 = 1:size(TRI6toTRI3,1)
            topo2   = topo1(1,TRI6toTRI3(el2,:));
            xcoore2 = xcoor1(topo2,:);
            xvalue2 = xvalue1(:,TRI6toTRI3(el2,:));
            count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore2);
            count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue2);
          end
        end
    end
    end
  end
end

% TensorTriangleValue
% """""""""""""""""""
if nTensorTriangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor2,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'TRI3',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = zeros(3,9); % (nbno,ncomp)
          for i = 1:9
            if icomp1(i)
              xvalue1(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
        end
      case 'TRI6',
        for el1 = 1:nbel1
          topo1   = maile1.MAIL(el1,:);
          xvalue1 = zeros(9,6); % (nbcomp,nbno)
          for i = 1:9
            if icomp1(i)
              xvalue(i,:) = chmlnoe1{icomp1(i)}.XVAL(el1,:);
            end
          end
          for el2 = 1:size(TRI6toTRI3,1)
            topo2   = topo1(1,TRI6toTRI3(el2,:));
            xcoore2 = xcoor1(topo2,:);
            xvalue2 = xvalue1(:,TRI6toTRI3(el2,:));
            count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore2);
            count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue2);
          end
        end
    end
    end
  end
end

% ScalarQuadrangleValue 
% """""""""""""""""""""
if nScalarQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur({name1},listCompe1);
    if icomp1
    switch type1
      case 'QUA4',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e\n',xvalue1);
        end
      case 'QUA8',
        error('QUA8 not yet implemented')
    end
    end
  end
end

% VectorQuadrangleValue 
% """""""""""""""""""""
if nVectorQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompVector,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'QUA4',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = zeros(4,3); % (nbno,ncomp)
          for i = 1:3
            if icomp1(i)
              xvalue(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1');
        end
      case 'QUA8',
        error('QUA8 not yet implemented')
    end
    end
  end
end

% TensorQuadrangleValue 
% """""""""""""""""""""
if nTensorQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
    icomp1 = findoccur(ListCompTensor,listCompe1);
    if ~all(icomp1==0)
    switch type1
      case 'QUA4',
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = zeros(4,9); % (nbno,ncomp)
          for i = 1:9
            if icomp1(i)
              xvalue(:,i) = chmlnoe1{icomp1(i)}.XVAL(el1,:)';
            end
          end
          count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xvalue1');
        end
      case 'QUA8',
        error('QUA8 not yet implemented')
    end
    end
  end
end

% ScalarTetrahedronValue 
% """"""""""""""""""""""
if nScalarTetrahedra
  error('Not yet implemented')
end

% VectorTetrahedronValue 
% """"""""""""""""""""""
if nVectorTetrahedra
  error('Not yet implemented')
end

% TensorTetrahedronValue 
% """"""""""""""""""""""
if nTensorTetrahedra
  error('Not yet implemented')
end

% ScalarHexahedronValue 
% """""""""""""""""""""
if nScalarHexahedra
  error('Not yet implemented')
end

% VectorHexahedronValue 
% """""""""""""""""""""
if nVectorHexahedra
  error('Not yet implemented')
end

% TensorHexahedronValue 
% """""""""""""""""""""
if nTensorHexahedra
  error('Not yet implemented')
end

% ScalarPrismValue 
% """"""""""""""""
if nScalarPrisms
  error('Not yet implemented')
end

% VectorPrismValue 
% """"""""""""""""
if nVectorPrisms
  error('Not yet implemented')
end

% TensorPrismValue 
% """"""""""""""""
if nTensorPrisms
  error('Not yet implemented')
end

% ScalarPyramidValue 
% """"""""""""""""""
if nScalarPyramids
  error('Not yet implemented')
end

% VectorPyramidValue 
% """"""""""""""""""
if nVectorPyramids
  error('Not yet implemented')
end

% TensorPyramidValue 
% """"""""""""""""""
if nTensorPyramids
  error('Not yet implemented')
end

% text2d  text2d-chars 
% text3d  text3d-chars 

% End of field
  count = fprintf(fid,'$EndView\n');
