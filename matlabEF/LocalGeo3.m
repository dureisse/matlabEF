function [chamno1,intgno1,varargout] = LocalGeo3(modl1,mail1,intg1,...
                                                 xcoor1,mode1,varargin)
% Element-based field, defined at nodes, of local geometry (oriented)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 04 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 04 / 2004
%   Possibilite d'orienter les elements, passage au 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 07 / 2004
%   Modeles de plaque DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 01 / 2005
%   Changement de structure dans le champ chamno1
%
% [chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
% For structural elements (TIMO,POUT) find local geometry vector
% (tangent) independantly for each element.
% For structural elements (DKIR) find local geometry vector
% (all vectors) independantly for each element, constant per element
% (i.e. same values at each node of the same element)
%
% [chamno1,intgno1,mail2] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1,xcoor2);
% If coordinates of a point are provided, elements are reoriented and
% a new topology mail2 is provided.
%
% Local elemental basis :
%   chamno1{zo1}{1}.COMP = 'Repere'
%   chamno1{zo1}{1}.XVAL(nbel,nbno,nbvec,idim)

%
% Inputs
%   modl1		Model informations
%   mail1		Mesh
%   intg1		Integration informations
%   xcoor1(nbnot1,idim)	Coordinates of nodes
%   mode1		Mode of analyse (BARR,TIMO,POUT,DKIR)
% Optional inputs
%   xcoor2(1,idim)	Coordinates of a point to orient elements
% Outputs
%   chamno1		Element-based field
%   intgno1		Its integration informations
% Optional outputs
%   mail2		Oriented mesh if xcoor2 has been provided

nin1 = nargin - 5;
nout = nargout - 2;
switch nin1
  case 0,
    clear xcoor2;
  case 1,
    xcoor2 = varargin{1};
    nperm2 = 0;
  otherwise,
    nin1
    error('Wrong number of optional inputs')
end
if (nin1 ~= nout)
  nin1
  nout
  error('Incompatible numbers of input/output arguments')
end
clear mail2;

idim = size(xcoor1,2); % Dimension of the physical space
nbzo1 = length(intg1); % Number of zones

if nbzo1 ~= length(modl1) | nbzo1 ~= length(mail1)
  nbzo1
  length(modl1)
  length(mail1)
  error('Unconsistent data ')
end

% Loop on zones
% """""""""""""
for zo1 = 1:nbzo1
  intge1 = intg1{zo1};
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};
  clear maile2;

% list of basis functions for element transformation
  nnit = modle1.NNIT;
% list of nodes for element transformation
  nnot = modle1.NNOT;

  topo1 = maile1.MAIL;
  [nbel1,nbno1] = size(topo1);
  if strcmp(mode1,'BARR') | strcmp(mode1,'TIMO') | strcmp(mode1,'POUT')
%   Straight beam elements: 1 vector
    nbvec = 1;
  elseif strcmp(mode1,'DKIR')
%   Plane plate elements: 2 vectors
    nbvec = 2;
  end
  xvalt1 = zeros(nbel1,nbno1,nbvec,idim);

% Loop on elements
% """"""""""""""""
  for iel1 = 1:nbel1

%   Coordinate of nodes that participate to the transformation
    xcoortr1 = xcoor1(topo1(iel1,nnot),:);

    if strcmp(mode1,'BARR') | strcmp(mode1,'TIMO') | strcmp(mode1,'POUT')
%     Beam elements: look for tangent local vector at nodes

      T = zeros(idim,nbno1);

%     Loop on nodes
      for ino = 1:nbno1
%       Derivates of basis function used for element transformation
        dphix = intge1.DPHIN(:,nnit,ino);
%       Local transformation gradient and Jacobian
        [Mjaco,Jaco] = LocalJaco2(dphix,xcoortr1);
%       Tangent at node ino
        T(:,ino) = (1. / Jaco) * Mjaco';
      end

      xvalt1(nbel1,:,1,:) = T';
      clear T;

    elseif strcmp(mode1,'DKIR')
%     Plate elements: look for a unique local basis (N1,N2,N3)
%     N3 is the normal
      if (idim ~= 3)
        error('Only 3D for DKIR')
      end

      N1 = zeros(idim,nbno1);
      N2 = zeros(idim,nbno1);
      N3 = zeros(idim,nbno1);

%     Loop on nodes
      for ino = 1:nbno1
%       Derivates of basis function used for element transformation
        dphix = intge1.DPHIN(:,nnit,ino);
%       Local transformation gradient and Jacobian
        [Mjaco,Jaco] = LocalJaco2(dphix,xcoortr1);
%       Local basis by projection of gradient transformation
        F = Mjaco';
%%        Nnew = orth(F);
        Nnew = GramS(F);
        N3(:,ino) = ProdVect(Nnew); % normal
        Fnew = Nnew' * F; % projected gradient
        [Vnew,Rnew] = PolarDec(Fnew); % polar decomp. Fnew = Vnew*Rnew
        NN = Nnew * Rnew; N1(:,ino) = NN(:,1); N2(:,ino) = NN(:,2); clear NN;
      end
%     Average to get a local basis constant per element
      N1 = sum(N1,2) / nbno1;
      N2 = sum(N2,2) / nbno1;
      N3 = sum(N3,2) / nbno1;
      N3 = N3 / norm(N3);
      N1 = ProdVect([N2 N3]); N1 = N1 / norm(N1);
      N2 = ProdVect([N3 N1]); N2 = N2 / norm(N2);
%      xval1(iel1,:) = repmat([N1(1,1) N2(1,1) N3(1,1)],1,nbno1);
%      xval2(iel1,:) = repmat([N1(2,1) N2(2,1) N3(2,1)],1,nbno1);
%      xval3(iel1,:) = repmat([N1(3,1) N2(3,1) N3(3,1)],1,nbno1);

      for ino = 1:nbno1
        xvalt1(nbel1,ino,1,:) = N1';
        xvalt1(nbel1,ino,2,:) = N2';
        xvalt1(nbel1,ino,3,:) = N3';
      end
      clear N1 N2 N3;
    else
      mode1
      error('Not yet available for other modes than BARR, TIMO, POUT, DKIR')
    end

%   Permutation of nodes for correct orientation
    if exist('xcoor2')
      if (idim ~= 2)
        error('pb of orientation to be check for 3D')
      end
      ipos = 1; ineg = 1;
      for ino1 = 1:nbno1
%        v1 = [(xcoor1(topo1(iel1,ino),:) - xcoor2(1,:))
%              [xval1(iel1,ino) xval2(iel1,ino)]]';
        v1 = [(xcoor1(topo1(iel1,ino),:) - xcoor2(1,:))
	      xvalt1(iel1,ino,1,:)]';
        test = ProdVect(v1) / norm(v1(:,1));
        if abs(test) < 1.e-2
          test
          error('Orientation origin not suited')
        elseif test < 0
          ipos = 0;
        elseif test > 0
          ineg = 0;
        end
      end
%     ipos = ineg = 0: the sign is changing
%     ipos = 1, ineg = 0: all >0
%     ipos = 0, ineg = 1: all <0
      if ~xor(ipos,ineg)
        ipos
        ineg
        error('Orientation origin not suited 2')
      end
      if ineg
%       Permutation required
        nperm2 = nperm2 + 1;
        topo1(iel1,:) = topo1(iel1,[nbno1:-1:1]);
        xvalt1(iel1,:,:,:) = xvalt1(iel1,[nbno1:-1:1],:,:);
%        xval1(iel1,:) = xval1(iel1,[nbno1:-1:1]);
%        xval2(iel1,:) = xval2(iel1,[nbno1:-1:1]);
%        if (idim == 3); xval3(iel1,:) = xval3(iel1,[nbno1:-1:1]); end
      end
    end

    clear xcoortr1;
% End loop on elements
  end

% Store the result in the element-wise field
  clear chamnoe1;
%  chamnoe1{1} = struct('COMP','X1','UNIT','','XVAL',xval1);
%  chamnoe1{2} = struct('COMP','Y1','UNIT','','XVAL',xval2);
%  if (idim == 3)
%    chamnoe1{3} = struct('COMP','Z1','UNIT','','XVAL',xval3);
%  end
  chamnoe1{1} = struct('COMP','Repere','UNIT','','XVAL',xvalt1);

  chamno1{zo1} = chamnoe1;
  %clear chamnoe1;
  if exist('xcoor2')
    maile2 = maile1;
    maile2.MAIL = topo1;
    mail2{zo1} = maile2;
    clear maile2;
    disp(['  Number of reoriented elements in the zone: ' ...
          int2str(nperm2) ' / ' int2str(nbel1)])
  end

  clear xvalt1;
  clear topo1 maile1 intge1 modle1;
% End loop on zones
end

intgno1 = SegmentIntgNo(mail1);


if (nin1 == 1)
  varargout(1) = {mail2};
end