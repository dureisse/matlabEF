function cham1 = MeshPartitioning0(mail1,lnsdm1,Q,xcoor1);
% Coarse partitioning of a mesh in subdomains

idim = size(xcoor1,2);
if length(lnsdm1)~=idim
  idim
  error('bad list of sdm numbers')
end
nx1 = lnsdm1(1);
ny1 = lnsdm1(2);
if idim == 3
  nz1 = lnsdm1(3);
end

% Coordinates of element centroids
intg1 = SegmentIntgNo(mail1,'centroid');
cham2 = CoorCham(mail1,intg1,xcoor1);

% Rotated box
el1 = 0;
nbzo1 = length(mail1);
for izo1 = 1 : nbzo1
  nbel = size(mail1{izo1}.MAIL,1);
  chame2 = cham2{izo1};
  X2 = chame2{1}.XVAL;
  Y2 = chame2{2}.XVAL;
  if idim == 3
    Z2 = chame2{3}.XVAL;
  end
  for iel = 1 : nbel
    el1 = el1 + 1;
    if idim == 2
      xel1 = [X2(iel,1) Y2(iel,1)] * Q';
    elseif idim == 3
      xel1 = [X2(iel,1) Y2(iel,1) Z2(iel,1)] * Q';
    end
    if el1 == 1
      boxmin = xel1;
      boxmax = xel1;
    else
      boxmin = min([boxmin;xel1]);
      boxmax = max([boxmax;xel1]);
    end
  end
end
nelt1 = el1;

% subdomain attribution
lsdm1 = zeros(nelt1,1);
el1 = 0;
nbzo1 = length(mail1);
for izo1 = 1 : nbzo1
  nbel = size(mail1{izo1}.MAIL,1);
  chame2 = cham2{izo1};
  X2 = chame2{1}.XVAL;
  Y2 = chame2{2}.XVAL;
  if idim == 3
    Z2 = chame2{3}.XVAL;
  end
  for iel = 1 : nbel
    el1 = el1 + 1;
    if idim == 2
      xel1 = [X2(iel,1) Y2(iel,1)] * Q';
      X = 1. + floor((xel1(1)-boxmin(1))/(boxmax(1)-boxmin(1))*nx1);
      Y = 1. + floor((xel1(2)-boxmin(2))/(boxmax(2)-boxmin(2))*ny1);
      if (X > nx1); X = nx1; end
      if (Y > ny1); Y = ny1; end
      sdm1 = (Y-1)*nx1 + X;
    elseif idim == 3
      xel1 = [X2(iel,1) Y2(iel,1) Z2(iel,1)] * Q';
      X = 1. + floor((xel1(1)-boxmin(1))/(boxmax(1)-boxmin(1))*nx1);
      Y = 1. + floor((xel1(2)-boxmin(2))/(boxmax(2)-boxmin(2))*ny1);
      Z = 1. + floor((xel1(3)-boxmin(3))/(boxmax(3)-boxmin(3))*nz1);
      if (X > nx1); X = nx1; end
      if (Y > ny1); Y = ny1; end
      if (Z > nz1); Z = nz1; end
      sdm1 = (Z-1)*(nx1*ny1) + (Y-1)*nx1 + X;
    end
%    [el1 X Y sdm1]
    lsdm1(el1) = sdm1;
  end
end


clear cham1;
cham1{1} = struct('COMP','SDM','UNIT','','XVAL',lsdm1);
