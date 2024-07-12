function [error1] = WriteChamnoGMSH(xcoor1,mail1,chamno1,fid)
disp('Warning: WriteChamnoGMSH obsolete... should be replaced with WriteChamnoGMSH2')
% Write a field in gmsh ascii post-pro file (scalar element) v1.4
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 04 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 20 / 05 / 2004
%   Implantation du TRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 04 / 08 / 2006
%   Implantation du CUB8, TET4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2006
%   Implantation du PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 01 / 2007
%   Champs a plusieurs composantes... ce n'est pas un succes !
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 02 / 2007
%   Passage au format 1.4 (supporte les elements quadratiques)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 03 / 2007
%   Implantation du QUA8
%
% Inputs
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   mail1		maillage
%   chamno1		champ par element defini aux noeuds
%   fid			file identifier (file must be opened, see fopen)
% Output
%   error1		

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
  ncomp1 = length(listComp1);

% Loop on views (i.e. components)
% """""""""""""
  for iview = 1:ncomp1
    name1 = {listComp1{iview}};
    switch length(name1)
      case 1,
        disp(['  writing scalar view ' int2str(iview) ': '])
        disp([name1])
      case 3,
        disp(['  writing vector view ' int2str(iview) ': '])
        disp([name1])
      case 9,
        disp(['  writing tensor view ' int2str(iview) ': '])
        disp([name1])
      otherwise,
        name1
        length(name1)
        error('Number of components not understood')
    end

% Beginning of field
  count = fprintf(fid,'$View\n');
  nTimeSteps = 1;
  TimeStepValues = 1.;
  switch length(name1)
    case 1,
      count = fprintf(fid,'%s\n',name1{1});
    otherwise,
      count = fprintf(fid,'%s\n',cell2mat(name1));
  end
  count = fprintf(fid,'%d\n',nTimeSteps);

  nScalarPoints = 0; nVectorPoints = 0; nTensorPoints = 0;
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

  for izo1 = 1:nbzone1
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);

    if icomp1
    switch type1
      case 'POI1',
        switch length(name1)
          case 1,
            nScalarPoints = nScalarPoints + nbel1;
          case 3,
            nVectorPoints = nVectorPoints + nbel1;
          case 9,
            nTensorPoints = nTensorPoints + nbel1;
        end
      case 'SEG2',
        switch length(name1)
          case 1,
            nScalarLines = nScalarLines + nbel1;
          case 3,
            nVectorLines = nVectorLines + nbel1;
          case 9,
            nTensorLines = nTensorLines + nbel1;
        end
      case 'SEG3',
        switch length(name1)
          case 1,
            nScalarLines2 = nScalarLines2 + nbel1;
          case 3,
            nVectorLines2 = nVectorLines2 + nbel1;
          case 9,
            nTensorLines2 = nTensorLines2 + nbel1;
        end
      case 'TRI3',
        switch length(name1)
          case 1,
            nScalarTriangles = nScalarTriangles + nbel1;
          case 3,
            nVectorTriangles = nVectorTriangles + nbel1;
          case 9,
            nTensorTriangles = nTensorTriangles + nbel1;
        end
      case 'TRI6',
%       subtriangles
%        TRI6toTRI3 = [1 4 6
%                      4 2 5
%                      6 5 3
%                      4 5 6];
%        switch length(name1)
%          case 1,
%            nScalarTriangles = nScalarTriangles + size(TRI6toTRI3,1)*nbel1;
%          case 3,
%            nVectorTriangles = nVectorTriangles + size(TRI6toTRI3,1)*nbel1;
%          case 9,
%            nTensorTriangles = nTensorTriangles + size(TRI6toTRI3,1)*nbel1;
%        end
        switch length(name1)
          case 1,
            nScalarTriangles2 = nScalarTriangles2 + nbel1;
          case 3,
            nVectorTriangles2 = nVectorTriangles2 + nbel1;
          case 9,
            nTensorTriangles2 = nTensorTriangles2 + nbel1;
        end
      case 'QUA4',
        switch length(name1)
          case 1,
            nScalarQuadrangles = nScalarQuadrangles + nbel1;
          case 3,
            nVectorQuadrangles = nVectorQuadrangles + nbel1;
          case 9,
            nTensorQuadrangles = nTensorQuadrangles + nbel1;
        end
      case 'QUA8',
        switch length(name1)
          case 1,
            nScalarQuadrangles2 = nScalarQuadrangles2 + nbel1;
          case 3,
            nVectorQuadrangles2 = nVectorQuadrangles2 + nbel1;
          case 9,
            nTensorQuadrangles2 = nTensorQuadrangles2 + nbel1;
        end
      case 'CUB8',
        switch length(name1)
          case 1,
            nScalarHexahedra = nScalarHexahedra + nbel1;
          case 3,
            nVectorHexahedra = nVectorHexahedra + nbel1;
          case 9,
            nTensorHexahedra = nTensorHexahedra + nbel1;
        end
      case 'CU20',
        switch length(name1)
          case 1,
            nScalarHexahedra2 = nScalarHexahedra2 + nbel1;
          case 3,
            nVectorHexahedra2 = nVectorHexahedra2 + nbel1;
          case 9,
            nTensorHexahedra2 = nTensorHexahedra2 + nbel1;
        end
      case 'TET4',
        switch length(name1)
          case 1,
            nScalarTetrahedra = nScalarTetrahedra + nbel1;
          case 3,
            nVectorTetrahedra = nVectorTetrahedra + nbel1;
          case 9,
            nTensorTetrahedra = nTensorTetrahedra + nbel1;
        end
      case 'PRI6',
        switch length(name1)
          case 1,
            nScalarPrisms = nScalarPrisms + nbel1;
          case 3,
            nVectorPrisms = nVectorPrisms + nbel1;
          case 9,
            nTensorPrisms = nTensorPrisms + nbel1;
        end
      case 'PR15',
        switch length(name1)
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

  count = fprintf(fid,'%.3e\n',TimeStepValues);

% ScalarPointValues
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'POI1') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorPointValues
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'POI1') && icomp1 && (length(name1)==3))
      error('Pas encore fait...1')
    end
  end

% TensorPointValues
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'POI1') && icomp1 && (length(name1)==9))
      error('Pas encore fait...2')
    end
  end

% ScalarLineValue
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG2') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e\n',xvalue1);
      end
    end
  end

% VectorLineValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG2') && icomp1 && (length(name1)==3))
      error('Pas encore fait...3')
    end
  end

% TensorLineValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG2') && icomp1 && (length(name1)==9))
      error('Pas encore fait...4')
    end
  end

% ScalarTriangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI3') && icomp1 && (length(name1)==1))
        for el1 = 1:nbel1
          topo1 = maile1.MAIL(el1,:);
          xcoore1 = xcoor1(topo1,:);
          xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
          count = fprintf(fid,'%.3e %.3e %.3e\n',xvalue1);
        end
     end
  end

% VectorTriangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI3') && icomp1 && (length(name1)==3))
      error('Pas encore fait...5')
    end
  end

% TensorTriangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI3') && icomp1 && (length(name1)==9))
      error('Pas encore fait...6')
    end
  end

% ScalarQuadrangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA4') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end
  
% VectorQuadrangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA4') && icomp1 && (length(name1)==3))
      error('Pas encore fait...7')
    end
  end
  
% TensorQuadrangleValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA4') && icomp1 && (length(name1)==9))
      error('Pas encore fait...8')
    end
  end

% ScalarTetrahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TET4') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorTetrahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TET4') && icomp1 && (length(name1)==3))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
%       4 points, 3 valeurs par point
        xvalue1 = xvalue1([1 5 9 2 6 10 3 7 11 4 8 12]);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e\n',xvalue1);
      end
    end
  end

% TensorTetrahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TET4') && icomp1 && (length(name1)==9))
      error('Pas encore fait...10')
    end
  end

% ScalarHexahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CUB8') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorHexahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CUB8') && icomp1 && (length(name1)==3))
      error('Pas encore fait...11')
    end
  end

% TensorHexahedronValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CUB8') && icomp1 && (length(name1)==9))
      error('Pas encore fait...12')
    end
  end

% ScalarPrismValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PRI6') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorPrismValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    if (strcmp(type1,'PRI6') && icomp1 && (length(name1)==3))
      error('Pas encore fait...13')
    end
  end

% TensorPrismValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    if (strcmp(type1,'PRI6') && icomp1 && (length(name1)==9))
      error('Pas encore fait...14')
    end
  end

% ScalarPyramidValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    if (strcmp(type1,'PYR5') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorPyramidValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    if (strcmp(type1,'PYR5') && icomp1 && (length(name1)==3))
      error('Pas encore fait...15')
    end
  end

% TensorPyramidValue 
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    if (strcmp(type1,'PYR5') && icomp1 && (length(name1)==9))
      error('Pas encore fait...16')
    end
  end

% ScalarLineValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG3') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorLineValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG3') && icomp1 && (length(name1)==3))
      error('SEG3 not yet implemented 2')
    end
  end

% TensorLineValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'SEG3') && icomp1 && (length(name1)==9))
      error('SEG3 not yet implemented 3')
    end
  end

% ScalarTriangleValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI6') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorTriangleValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI6') && icomp1 && (length(name1)==3))
      error('To be done... 1')
    end
  end

% TensorTriangleValue2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TRI6') && icomp1 && (length(name1)==9))
      error('To be done... 2')
    end
  end

% ScalarQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA8') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA8') && icomp1 && (length(name1)==3))
      error('to be done... 3')
    end
  end

% TensorQuadrangles2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'QUA8') && icomp1 && (length(name1)==9))
      error('to be done... 4')
    end
  end

% ScalarTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TE10') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TE10') && icomp1 && (length(name1)==3))
      error('to be done... 6')
    end
  end

% TensorTetrahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'TE10') && icomp1 && (length(name1)==9))
      error('to be done... 6')
    end
  end

% ScalarHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CU20') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CU20') && icomp1 && (length(name1)==3))
      error('to be done... 10')
    end
  end

% TensorHexahedra2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'CU20') && icomp1 && (length(name1)==9))
      error('to be done... 11')
    end
  end

% ScalarPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PR15') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PR15') && icomp1 && (length(name1)==3))
      error('to be done... 12')
    end
  end

% TensorPrisms2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PR15') && icomp1 && (length(name1)==9))
      error('to be done... 13')
    end
  end

% ScalarPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PY13') && icomp1 && (length(name1)==1))
      for el1 = 1:nbel1
        topo1 = maile1.MAIL(el1,:);
        xcoore1 = xcoor1(topo1,:);
        xvalue1 = chmlnoe1{icomp1}.XVAL(el1,:);
        count = fprintf(fid,'%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n',xcoore1);
        count = fprintf(fid,'%.3e\n',xvalue1);
      end
    end
  end

% VectorPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PY13') && icomp1 && (length(name1)==3))
      error('to be done... 14')
    end
  end

% TensorPyramids2
  for izo1 = 1:nbzone1
    maile1   = mail1{izo1};
    chmlnoe1 = chamno1{izo1};
    type1 = list1{izo1};
    nbel1 = lnelt1(izo1);
    [listCompe1,listUnite1] = ListCompCham2(chamno1,[izo1]);
%%    icomp1 = findoccur({name1},listCompe1);
%%    icomp1 = findoccur1cell(name1,listCompe1);
    icomp1 = findoccur(name1,listCompe1);
    if (strcmp(type1,'PY13') && icomp1 && (length(name1)==9))
      error('to be done... 15')
    end
  end

% text2d  text2d-chars 
% text3d  text3d-chars 

% End of field
  count = fprintf(fid,'$EndView\n');

% End loop on views
% """""""""""""""""
end







function [icomp1] = findoccur1cell(name1,listCompe1)
    if iscell(name1)
      icomp1 = 0;
      for j = 1:length(listCompe1)
        listCompe2 = listCompe1(j);
        if iscell(listCompe2)
           listCompe2 = listCompe2{1};
           if length(name1) == length(listCompe2)
             for k = 1:length(name1)
               icomp1 = j;
               name2 = name1(k);
               Compe2 = listCompe2(k);
               if ~strcmp(name2,Compe2)
                 icomp1 = 0;
                 break
               end
             end
           end
        end
      end
    else
      icomp1 = findoccur({name1},listCompe1);
    end
