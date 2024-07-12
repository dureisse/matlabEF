function [error1] = WriteChamnoGMSH2(xcoor1,mail1,chamno1, ...
                                     LlistComp1,listnames,fid)
% Write a field in a gmsh ascii post-pro file v1.4
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 19 / 07 / 2008
%
% Inputs
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   mail1		maillage
%   chamno1		champ par element defini aux noeuds
%   LlistComp1{n}	list of list of components to write
%   listnames{n}	list of view names
%   fid			file identifier (file must be opened, see fopen)
% Output
%   error1		error indicator
%
% On ecrit n vues. Dans chaque vue i, un champ est ecrit.
% Il s'agit d'un champ scalaire si length(LlistComp1{i}) = 1, 
%                      vecteur si length(LlistComp1{i}) = 3 (en 3D), 
%                      tenseur si length(LlistComp1{i}) = 9 (en 3D).
% Si la composante de la liste n'est pas dans le champ, elle est mise
% a zero, ou elle n'est pas ecrite... !!!

error1 = 0;

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('Dimension should be 3')
end

% Header
  count = fprintf(fid,'$PostFormat\n');
  disp('Current format version is 1.4')
  count = fprintf(fid,'%g %d %d\n',[1.4 0 8]);
  count = fprintf(fid,'$EndPostFormat\n');

% Data size
  nbzone1 = length(mail1);
  list1   = ListTypeMesh(mail1);
  lnelt1  = Nelements2(mail1,[1:nbzone1]);  
  [listComp1,listUnit1] = ListCompCham2(chamno1);
  nview1 = length(LlistComp1);

% Loop on views
% """""""""""""
  for iview = 1:nview1
    listComp1 = LlistComp1{iview};
    name1     = listnames{iview};
    switch length(listComp1)
      case 1,
        disp(['  writing scalar view ' int2str(iview) ': ' name1])
      case 3,
        disp(['  writing vector view ' int2str(iview) ': ' name1])
      case 9,
        disp(['  writing tensor view ' int2str(iview) ': ' name1])
      otherwise,
        iview
        name1
        listComp1
        error('Number of components not understood')
    end

% Beginning of the view
  count = fprintf(fid,'$View\n');
  nTimeSteps = 1;
  TimeStepValues = 1.;
  count = fprintf(fid,'%s\n',name1);
  count = fprintf(fid,'%d\n',nTimeSteps);

  nScalarPoints       = 0; nVectorPoints       = 0; nTensorPoints       = 0;
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

  nText2d = 0; nText2dChars = 0; nTest3d = 0; nText3dChars = 0;

% Loop on zones
% """""""""""""
  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);

    if all(icomp1~=0)
%   Toutes les composantes ont ete trouvees
    switch type1
      case 'POI1',
        switch length(listComp1)
          case 1,
            nScalarPoints = nScalarPoints + nbel1;
          case 3,
            nVectorPoints = nVectorPoints + nbel1;
          case 9,
            nTensorPoints = nTensorPoints + nbel1;
        end
      case 'SEG2',
        switch length(listComp1)
          case 1,
            nScalarLines = nScalarLines + nbel1;
          case 3,
            nVectorLines = nVectorLines + nbel1;
          case 9,
            nTensorLines = nTensorLines + nbel1;
        end
      case 'SEG3',
        switch length(listComp1)
          case 1,
            nScalarLines2 = nScalarLines2 + nbel1;
          case 3,
            nVectorLines2 = nVectorLines2 + nbel1;
          case 9,
            nTensorLines2 = nTensorLines2 + nbel1;
        end
      case 'TRI3',
        switch length(listComp1)
          case 1,
            nScalarTriangles = nScalarTriangles + nbel1;
          case 3,
            nVectorTriangles = nVectorTriangles + nbel1;
          case 9,
            nTensorTriangles = nTensorTriangles + nbel1;
        end
      case 'TRI6',
        switch length(listComp1)
          case 1,
            nScalarTriangles2 = nScalarTriangles2 + nbel1;
          case 3,
            nVectorTriangles2 = nVectorTriangles2 + nbel1;
          case 9,
            nTensorTriangles2 = nTensorTriangles2 + nbel1;
        end
      case 'QUA4',
        switch length(listComp1)
          case 1,
            nScalarQuadrangles = nScalarQuadrangles + nbel1;
          case 3,
            nVectorQuadrangles = nVectorQuadrangles + nbel1;
          case 9,
            nTensorQuadrangles = nTensorQuadrangles + nbel1;
        end
      case 'QUA8',
        switch length(listComp1)
          case 1,
            nScalarQuadrangles2 = nScalarQuadrangles2 + nbel1;
          case 3,
            nVectorQuadrangles2 = nVectorQuadrangles2 + nbel1;
          case 9,
            nTensorQuadrangles2 = nTensorQuadrangles2 + nbel1;
        end
      case 'CUB8',
        switch length(listComp1)
          case 1,
            nScalarHexahedra = nScalarHexahedra + nbel1;
          case 3,
            nVectorHexahedra = nVectorHexahedra + nbel1;
          case 9,
            nTensorHexahedra = nTensorHexahedra + nbel1;
        end
      case 'CU20',
        switch length(listComp1)
          case 1,
            nScalarHexahedra2 = nScalarHexahedra2 + nbel1;
          case 3,
            nVectorHexahedra2 = nVectorHexahedra2 + nbel1;
          case 9,
            nTensorHexahedra2 = nTensorHexahedra2 + nbel1;
        end
      case 'TET4',
        switch length(listComp1)
          case 1,
            nScalarTetrahedra = nScalarTetrahedra + nbel1;
          case 3,
            nVectorTetrahedra = nVectorTetrahedra + nbel1;
          case 9,
            nTensorTetrahedra = nTensorTetrahedra + nbel1;
        end
      case 'TE10',
        switch length(listComp1)
          case 1,
            nScalarTetrahedra2 = nScalarTetrahedra2 + nbel1;
          case 3,
            nVectorTetrahedra2 = nVectorTetrahedra2 + nbel1;
          case 9,
            nTensorTetrahedra2 = nTensorTetrahedra2 + nbel1;
        end
      case 'PRI6',
        switch length(listComp1)
          case 1,
            nScalarPrisms = nScalarPrisms + nbel1;
          case 3,
            nVectorPrisms = nVectorPrisms + nbel1;
          case 9,
            nTensorPrisms = nTensorPrisms + nbel1;
        end
      case 'PR15',
        switch length(listComp1)
          case 1,
            nScalarPrisms2 = nScalarPrisms2 + nbel1;
          case 3,
            nVectorPrisms2 = nVectorPrisms2 + nbel1;
          case 9,
            nTensorPrisms2 = nTensorPrisms2 + nbel1;
        end
      otherwise,
        type1
        error('Element type not implemented')
    end
    end
  end
  count = fprintf(fid,'%d %d %d\n',[nScalarPoints       nVectorPoints       nTensorPoints]);
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

  count = fprintf(fid,'%.3e\n',TimeStepValues);

% ======================================================================
% Support = point
% ======================================================================

% ScalarPointValues
% """""""""""""""""
  if nScalarPoints
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'POI1') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorPointValues
% """""""""""""""""
  if nVectorPoints
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'POI1') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...1')
    end
  end
  end

% TensorPointValues
% """""""""""""""""
  if nTensorPoints
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'POI1') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...2')
    end
  end
  end

% ======================================================================
% Support = line (2 nodes)
% ======================================================================

% ScalarLineValue
% """""""""""""""
  if nScalarLines
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG2') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e\n',xvalue1);
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
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG2') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...3')
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
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG2') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...4')
    end
  end
  end

% ======================================================================
% Support = triangle (3 nodes)
% ======================================================================

% ScalarTriangleValue
% """"""""""""""""""" 
  if nScalarTriangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI3') && all(icomp1~=0) && length(listComp1)==1)
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1);
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
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI3') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...5')
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
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI3') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...6')
    end
  end
  end

% ======================================================================
% Support = quagrangle (4 nodes)
% ======================================================================

% ScalarQuadrangleValue 
% """"""""""""""""""" 
  if nScalarQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA4') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end
  
% VectorQuadrangleValue 
% """"""""""""""""""" 
  if nVectorQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA4') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...7')
    end
  end
  end
  
% TensorQuadrangleValue 
% """"""""""""""""""" 
  if nTensorQuadrangles
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA4') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...8')
    end
  end
  end

% ======================================================================
% Support = tetrahedron (4 nodes)
% ======================================================================

% ScalarTetrahedronValue 
% """""""""""""""""""""" 
  if nScalarTetrahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TET4') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorTetrahedronValue 
% """""""""""""""""""""" 
  if nVectorTetrahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TET4') && all(icomp1~=0) && length(listComp1)==3)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
%       4 points, 3 valeurs par point
        xvalue1 = xvalue1([1 5 9 2 6 10 3 7 11 4 8 12]); % Convention de numeros
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xvalue1);
      end
    end
  end
  end

% TensorTetrahedronValue 
% """""""""""""""""""""" 
  if nTensorTetrahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TET4') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...10')
    end
  end
  end

% ======================================================================
% Support = hexahedron (6 nodes)
% ======================================================================

% ScalarHexahedronValue 
% """""""""""""""""""""
  if nScalarHexahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CUB8') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorHexahedronValue 
% """""""""""""""""""""
  if nVectorHexahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CUB8') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...11')
    end
  end
  end

% TensorHexahedronValue 
% """""""""""""""""""""
  if nTensorHexahedra
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CUB8') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...12')
    end
  end
  end

% ======================================================================
% Support = Prism (6 nodes)
% ======================================================================

% ScalarPrismValue 
% """"""""""""""""
  if nScalarPrisms
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PRI6') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorPrismValue 
% """"""""""""""""
  if nVectorPrisms
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PRI6') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...13')
    end
  end
  end

% TensorPrismValue 
% """"""""""""""""
  if nTensorPrisms
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PRI6') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...14')
    end
  end
  end

% ======================================================================
% Support = pyramid (5 nodes)
% ======================================================================

% ScalarPyramidValue
% """""""""""""""""" 
  if nScalarPyramids
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PYR5') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorPyramidValue 
% """""""""""""""""" 
  if nVectorPyramids
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PYR5') && all(icomp1~=0) && length(listComp1)==3)
      error('Pas encore fait...15')
    end
  end
  end

% TensorPyramidValue 
% """""""""""""""""" 
  if nTensorPyramids
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PYR5') && all(icomp1~=0) && length(listComp1)==9)
      error('Pas encore fait...16')
    end
  end


% ======================================================================
% Support = line (3 nodes)
% ======================================================================

% ScalarLineValue2
% """""""""""""""
  if nScalarLines2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG3') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorLineValue2
% """""""""""""""
  if nVectorLines2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG3') && all(icomp1~=0) && length(listComp1)==3)
      error('SEG3 not yet implemented 2')
    end
  end
  end

% TensorLineValue2
% """""""""""""""
  if nTensorLines2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'SEG3') && all(icomp1~=0) && length(listComp1)==9)
      error('SEG3 not yet implemented 3')
    end
  end
  end

% ======================================================================
% Support = triangle (6 nodes)
% ======================================================================

 

% ScalarTriangleValue2
% """"""""""""""""""""
  if nScalarTriangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI6') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorTriangleValue2
% """"""""""""""""""""
  if nScalarTriangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI6') && all(icomp1~=0) && length(listComp1)==3)
      error('To be done... 1')
    end
  end
  end

% TensorTriangleValue2
% """"""""""""""""""""
  if nTensorTriangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TRI6') && all(icomp1~=0) && length(listComp1)==9)
      error('To be done... 2')
    end
  end
  end

% ======================================================================
% Support = quagrangle (8 nodes)
% ======================================================================

% ScalarQuadrangles2
% """"""""""""""""""
  if nScalarQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA8') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorQuadrangles2
% """"""""""""""""""
  if nVectorQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA8') && all(icomp1~=0) && length(listComp1)==3)
      error('to be done... 3')
    end
  end
  end

% TensorQuadrangles2
% """"""""""""""""""
  if nTensorQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'QUA8') && all(icomp1~=0) && length(listComp1)==9)
      error('to be done... 4')
    end
  end
  end

% ======================================================================
% Support = tetrahedron (10 nodes)
% ======================================================================

% ScalarTetrahedra2
% """"""""""""""""" 
  if nScalarTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TE10') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorTetrahedra2
% """"""""""""""""" 
  if nVectorTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TE10') && all(icomp1~=0) && length(listComp1)==3)
      error('to be done... 6')
    end
  end
  end

% TensorTetrahedra2
% """"""""""""""""" 
  if nTensorTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'TE10') && all(icomp1~=0) && length(listComp1)==9)
      error('to be done... 6')
    end
  end
  end

% ======================================================================
% Support = hexahedron (20 nodes)
% ======================================================================

% ScalarHexahedra2
% """"""""""""""""
  if nScalarHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CU20') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorHexahedra2
% """"""""""""""""
  if nVectorHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CU20') && all(icomp1~=0) && length(listComp1)==3)
      error('to be done... 10')
    end
  end
  end

% TensorHexahedra2
% """"""""""""""""
  if nTensorHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'CU20') && all(icomp1~=0) && length(listComp1)==9)
      error('to be done... 11')
    end
  end
  end

% ======================================================================
% Support = Prism (15 nodes)
% ======================================================================

% ScalarPrisms2
% """""""""""""
  if nScalarPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PR15') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorPrisms2
% """""""""""""
  if nVectorPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PR15') && all(icomp1~=0) && length(listComp1)==3)
      error('to be done... 12')
    end
  end
  end

% TensorPrisms2
% """""""""""""
  if nTensorPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PR15') && all(icomp1~=0) && length(listComp1)==9)
      error('to be done... 13')
    end
  end
  end

% ======================================================================
% Support = pyramid (13 nodes)
% ======================================================================

% ScalarPyramids2
% """"""""""""""" 
  if nScalarPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PY13') && all(icomp1~=0) && length(listComp1)==1)
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  end

% VectorPyramids2
% """"""""""""""" 
  if nVectorPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PY13') && all(icomp1~=0) && length(listComp1)==3)
      error('to be done... 14')
    end
  end
  end

% TensorPyramids2
% """"""""""""""" 
  if nTensorPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
%   Liste des composantes du champ pour la zone
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%   On cherche ou sont les composantes demandees dans la zone
    icomp1 = findoccur(listComp1,listCompe1);
    if (strcmp(type1,'PY13') && all(icomp1~=0) && length(listComp1)==9)
      error('to be done... 15')
    end
  end
  end

% End loop on zones
% """""""""""""""""
end

% text2d  text2d-chars 
% text3d  text3d-chars 

% End of the view
  count = fprintf(fid,'$EndView\n');

% End loop on views
% """""""""""""""""
end

