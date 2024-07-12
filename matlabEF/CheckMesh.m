function [val1] = CheckMesh(mail1,xcoor1,option1)
% Check integrity of a mesh
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 07 / 10 / 2007
%
% Inputs
%   mail1		Mesh
%   xcoor1(nbno,idim)	Coordinates
%   option1		Type of verifiation: 'Surface' 'EdgeLength'
% Output
%   val1		

lmax1 = []; % maximal edge length
lmin1 = []; % minimal edge length
lsur1 = []; % area
lmoy1 = []; % average edge length (equivalence with resp. to surface)

nbzo1 = length(mail1);
for zo1 = 1 : nbzo1
  maile1 = mail1{zo1};
  type1  = maile1.TYPE;
  topo1  = maile1.MAIL;
  nbel1  = size(topo1,1);
  nbno1  = size(topo1,2);
  switch type1
    case {'TRI3','TRI6'},
      for el1 = 1 : nbel1
        lnod1 = topo1(el1,1:3); % vertex only for TRI6
        lnod2 = circshift(lnod1,[0 1]);
        xvec1 = xcoor1(lnod1',:) - xcoor1(lnod2',:);
%       maximal edge length
        lmax1(end+1) = max(sum(xvec1.^2,2)');
%       minimal edge length
        lmin1(end+1) = min(sum(xvec1.^2,2)');
%       area
        lsur1(end+1) = 0.5*abs(ProdVect(xvec1(1:2,:)'));
%       surf = 0.5*lon1*lon1*sqrt(3.)/2.; % for equilateral triangle
        lmoy1(end+1) = sqrt((4./sqrt(3.)) * lsur1(end));
      end
    otherwise,
      type1
      error('Mesh type not implemented')
   end
end

switch option1
  case 'EdgeLength',
    val1 = min(lmin1 ./ lmax1); % ratio between edge lengths
  case 'Surface',
    val1 = min(lmin1 ./ lmoy1);
  otherwise,
    option1
    error('option not known')
end
