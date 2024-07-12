function [n] = MediumPlane(xcoor)
% Find the medium plane of a set of points
%
% Inputs
%   xcoor(npts,idim)		Coordinates of the points
% Outputs
%   n(idim,1)			Vector defining the plane
%
% The plane has n as a normal (not unitary)
% and pass by the point A such that A = n


%b = sum(xcoor); % b(1,idim)
%idim = size(xcoor,2);
%M = zeros(idim,idim);
%for i = 1:size(xcoor,1)
%  M = M + xcoor(i,:)' * xcoor(i,:);
%end
%s = svd(M);
%if (s(1)-s(end))/s(1) < 1.e-6
%  error('No plane can be defined')
%end
%n = M \ b';


[npts,idim] = size(xcoor);

% 'Mass' matrix (or weighting matrix)
M1 = eye(idim*npts);
xddl = [1:idim:idim*npts];
yddl = [2:idim:idim*npts];
zddl = [3:idim:idim*npts];
[XG1,S1,I1] = mM_CaractMass3(M1,xcoor',xddl,yddl,zddl);

% Maximal eigenvalue: normal to the plane
error('to be done')

