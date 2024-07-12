function [lnod1,lind1,ListMeshInt] = DDInt2(ListContSdm,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 06 / 11 / 2005
%  Treatment of weak interfaces and optional arguments

% Avec la liste des maillages des contours des sous domaines
% ListContSdm(nsdm1), on trouve :
% - la liste des maillages des interfaces    ListMeshInt(nint1)
% - le graphe de connectivite de la decomposition (lnod1,lind1)
%   (les aretes sont les interfaces, les sommets les sous domaines)

% Inputs
%  ListContSdm{nsdm1}	list of meshes of boundary of subdomains
% Optional inputs
%  WeakIntOption	keyword for weak interfaces detection
%      'None':    no detection of weak interfaces
%      'Skip':    do not consider weak interfaces (beware that
%                 this may change the reference problem)
%      'Micro':   let weak interface be micro only
%      'Macro':   let weak interface be macro only
%      'Rbms':    macro part on weak interfaces is rbm per subdomain
%      'Corners': consider weak interfaces as an augmentation of
%                 corner nodes (see FETI-DP)
%  nbming1		minimum number of nodes to consider an interface
%                       as strong
% Outputs
%  (lnod1,lind1)	decomposition graph connectivity
%                       (edges are interfaces, vertices are subdomains)
%  ListMeshInt{nint1}	list of meshes of interfaces

nsdm1 = length(ListContSdm);
ListMeshInt = {};

nin1 = nargin-1;
if nin1 == 0
  WeakIntOption = 'None';
elseif nin1 == 2
  WeakIntOption = varargin{1};
  nbming1       = varargin{2};
  disp(['Treatment of weak interfaces: ' WeakIntOption])
  disp(['Minimum number of nodes:      ' int2str(nbming1)])
else
  nin1
  error('Wrong number of optional arguments')
end


% Avec les maillages des contours, on cherche les interfaces
% et leurs maillages dans la liste ListMeshInt(nint1),
% on place de plus la topologie de la decomposition dans
% TopoDec(nint1,2)
nint1 = 0;
nint0 = 0;
TopoDec = zeros(0,2);
for sdm1 = 1:nsdm1-1
  for sdm2 = sdm1+1:nsdm1
    mint1 = IntersectMesh(ListContSdm{sdm1},ListContSdm{sdm2});
    if length(mint1)
%     A new interface

      nmail1 = ChangeMesh2(mint1,'POI1');
      nnod1 = size(nmail1{1}.MAIL,1);
      if strcmp(WeakIntOption,'Skip') & (nnod1 <= nbming1)
%       Get rid of this interface
        nint0 = nint0 + 1;
      else
%       Classical treatment
        nint1 = nint1 + 1;
        TopoDec(nint1,:) = [sdm1 sdm2];
        ListMeshInt{nint1} = mint1;
      end
      
    end
  end
end

% Graphe de connectivite entre sous domaines
clear MeshDec;
MeshDec{1} = struct('TYPE','SEG2','MAIL',TopoDec);
[lnod1,lind1] = MeshToGraph(MeshDec);
clear TopoDec MeshDec;

disp(['Number of detected interfaces ' int2str(nint1)])
if strcmp(WeakIntOption,'Skip')
  disp(['Number of skipped weak interfaces ' int2str(nint0)])
end
