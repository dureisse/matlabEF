function [error1] = WriteChpoGMSH(xcoor1,nmail1,chpo1,fid,name1,varargin)
% Write a field in gmsh ascii post-pro file (scalar/vector/tensor point)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 12 / 2004
%
% error1 = WriteChpoGMSH(xcoor1,nmail1,chpo1,fid,ListC1)
% error1 = WriteChpoGMSH(xcoor1,nmail1,chpo1,fid,ListC1,ListC2)
% error1 = WriteChpoGMSH(xcoor1,nmail1,chpo1,fid,ListC1,ListC2,ListC3)
%
% Inputs
%   xcoor1(nbno,3)	coordonnees (si 2D la troisieme est nulle)
%   nmail1		maillage nuage de points (POI1) sous-tendant
%   chpo1		champ par point
%   name1		view name
%   fid			file identifier (file must be opened, see fopen)
% Optional inputs
%   ListC1{nbcomp1}	list of components to be written
%                       1 component for scalar point field
%                       3 components for vector point field
%                       9 components for tensor point field
%   if different lists are provided, they should be all length different
% Output
%   error1		en cas d'erreur si <> 0
%
% Only one view is saved in the file.
% The file has to be opened before calling the routine, and closed after

error1 = 0;

% Checking arguments
% """"""""""""""""""
narg = nargin - 5; % Nombre d'arguments optionnels
nScalarPoints = 0; nVectorPoints = 0; nTensorPoints = 0;
nbpts = Nelements2(nmail1);
[ListComp1,ListUnit1] = ListCompChpo2(chpo1); 
for arg1 = 1:narg
  List1 = varargin{arg1};
  switch length(List1)
    case 1,
     if nScalarPoints
       error('Scalar list of component already provided')
     end
     nScalarPoints = nbpts;
     ListCompScalar = List1;
    case 3,
     if nVectorPoints
       error('Vector list of component already provided')
     end
     nVectorPoints = nbpts;
     ListCompVector = List1;
    case 9,
     if nTensorPoints
       error('Tensor list of component already provided')
     end
     nTensorPoints = nbpts;
     ListCompTensor = List1;
    otherwise,
      List1
      error('Not a recognized list of components (bad number of comp)')
  end
end
disp(['Number of points for scalar field ' int2str(nScalarPoints)])
disp(['Number of points for vector field ' int2str(nVectorPoints)])
disp(['Number of points for tensor field ' int2str(nTensorPoints)])

idim = size(xcoor1,2);
if (idim ~= 3)
  idim
  error('Dimension should be 3')
end

% Local node numbering (1 zone assumed)
% """"""""""""""""""""
nbzone1 = length(nmail1);
if (nbzone1 ~= 1)
  error('More than 1 zone not implemented')
end
zo1 = 1;
numer1 = nmail1{zo1}.MAIL';
localxcoor1 = xcoor1(numer1,:);

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

%  nScalarPoints = 0; nVectorPoints = 0; nTensorPoints = 0;
  nScalarLines = 0; nVectorLines = 0; nTensorLines = 0;
  nScalarTriangles = 0; nVectorTriangles = 0; nTensorTriangles = 0;
  nScalarQuadrangles = 0; nVectorQuadrangles = 0; nTensorQuadrangles = 0;
  nScalarTetrahedra = 0; nVectorTetrahedra = 0; nTensorTetrahedra = 0;
  nScalarHexahedra = 0; nVectorHexahedra = 0; nTensorHexahedra = 0;
  nScalarPrisms = 0; nVectorPrisms = 0; nTensorPrisms = 0;
  nScalarPyramids = 0; nVectorPyramids = 0; nTensorPyramids = 0;
  nText2d = 0; nText2dChars = 0; nTest3d = 0; nText3dChars = 0;

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
  ind1 = findoccur(ListCompScalar,ListComp1);
  localvalues(:,1) = chpo1{zo1}{ind1}.XVAL;
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor1(ino,:));
    count = fprintf(fid,'%d\n',localvalues(ino,:));
  end
end

% VectorPointValues
% """""""""""""""""
if nVectorPoints
  localvalues = zeros(nbpts,3);
  ind1 = findoccur(ListCompVector,ListComp1);
  if ind1(1); localvalues(:,1) = chpo1{zo1}{ind1(1)}.XVAL; end
  if ind1(2); localvalues(:,2) = chpo1{zo1}{ind1(2)}.XVAL; end
  if ind1(3); localvalues(:,3) = chpo1{zo1}{ind1(3)}.XVAL; end
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor1(ino,:));
    count = fprintf(fid,'%d %d %d\n',localvalues(ino,:));
  end
end

% TensorPointValues
% """""""""""""""""
if nTensorPoints
  localvalues = zeros(nbpts,9);
  ind1 = findoccur(ListCompTensor,ListComp1);
  for i = 1:9
    if ind1(i); localvalues(:,i) = chpo1{zo1}{ind1(i)}.XVAL; end
  end
  for ino = 1:nbpts
    count = fprintf(fid,'%d\n',localxcoor1(ino,:));
    count = fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',localvalues(ino,:));
  end
end

% ScalarLineValue
% VectorLineValue 
% TensorLineValue 
% ScalarTriangleValue 
% VectorTriangleValue 
% TensorTriangleValue 
% ScalarQuadrangleValue 
% VectorQuadrangleValue 
% TensorQuadrangleValue 
% ScalarTetrahedronValue 
% VectorTetrahedronValue 
% TensorTetrahedronValue 
% ScalarHexahedronValue 
% VectorHexahedronValue 
% TensorHexahedronValue 
% ScalarPrismValue 
% VectorPrismValue 
% TensorPrismValue 
% ScalarPyramidValue 
% VectorPyramidValue 
% TensorPyramidValue 
% text2d  text2d-chars 
% text3d  text3d-chars 

% End of field
  count = fprintf(fid,'$EndView\n');
