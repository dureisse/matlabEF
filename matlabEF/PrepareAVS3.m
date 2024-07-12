function [xcoor2,mail2,nmail3,varargout] = PrepareAVS3...
	  (xcoor1,mail1,nmail0,varargin)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 01 / 2003
%   Bug numerotation nmail3 corrige
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2003
%   Utilisation des listes de chpo en option
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 04 / 2003
%   On permute les valeurs du chpo aussi avec la renumerotation
%
% Verifie et rectifie la coherence des donnees pour ecrire en AVS UCD
%
% Entrees
%   xcoor1		Pile des noeuds
%   mail1		Maillage numerote en accord avec xcoor1
%   nmail0		Maillage POI1 numerote en accord avec xcoor1
% Entree optionnelle
%   chpo1		Champ par point
%   Lchpo1		Liste de champs par point
% Sorties
%   xcoor2		Pile des noeuds utiles (etendue au 3D)
%   mail2		Maillage renumerote en accord avec xcoor2
%   nmail3		Maillage POI1 numerote en accord avec xcoor2
% Sortie optionnelle
%   chpo2		Champ par point etendu a tous les noeuds utiles
%   Lchpo2		Liste de champs par point
%
% Les entrees optionnelles chpo1 et Lchpo1 sont mutuellement exclusives.
% Si chpo1 est donne, chpo2 est retourne,
% si Lchpo1 est donne, Lchpo2 est retourne.

nin1 = nargin - 3;
switch nin1
  case 0,
    Lchpo1 = []; chpo1 = [];
    test_list = -1;
  case 1,
    chpo1 = varargin{1};
    if isempty(chpo1)
      Lchpo1 = []; chpo1 = [];
%%DD ajout 28/02/2006
      test_list = -2;
    else
      if iscell(chpo1{1}{1})
        Lchpo1 = varargin{1}; clear chpo1;
        test_list = 1;
      elseif isstruct(chpo1{1}{1})
        Lchpo1{1} = chpo1; clear chpo1;
        test_list = 0;
      else
        error('Unknown type for optional argument')
      end
    end
  otherwise,
    nin1
    error('Bad number of optional argument')
end
if nargout ~= nargin
  nargin
  nargout
  error('Bad number of ouputs')
end


[nbnot1,idim] = size(xcoor1);

% Test and update space dimension
% """""""""""""""""""""""""""""""
if (idim == 2)
  disp('Extension from 2D to 3D with null 3rd coordinate')
  xcoor1 = [xcoor1 zeros(nbnot1,1)];
elseif (idim == 3)
else
  idim
  error('Bad space dimension')
end

% Local node numbering
% """"""""""""""""""""
% Select nodes
list_node1 = zeros(1,0);
if ~isempty(mail1)
  nmail1 = ChangeMesh2(mail1,'POI1');
  list_node1 = [list_node1 nmail1{1}.MAIL'];
  clear nmail1;
end
if ~isempty(nmail0)
  nbzone0 = TestMeshType(nmail0,'POI1');
  if nbzone0 == 0
    error('Mesh is not POI1')
  end
  for zo0 = 1:nbzone0
    list_node1 = [list_node1 nmail0{zo0}.MAIL'];
  end
end
list_node1 = unique(list_node1);
nbnot2 = length(list_node1);

disp(['Used ' int2str(nbnot2) ' nodes among ' int2str(nbnot1)])

% Update xcoor2 with new numbering
xcoor2 = xcoor1(list_node1,:);

% Update meshes with new numbering
inv_list_node1 = InverseList(list_node1,max(list_node1));
mail2 = [];
if ~isempty(mail1)
  disp('Renumbering Mesh')
  mail2 = RenumMesh(mail1,inv_list_node1);
end
if ~isempty(nmail0)
  disp('Renumbering Mesh')
  nmail2 = RenumMesh(nmail0,inv_list_node1);
end
clear inv_list_node1;
list_node1new = [1:length(list_node1)];

% Extend chpo to nodes and only one subzone
Lchpo2 = [];
chpo2 = [];
if ~isempty(Lchpo1)
  nchpo1 = length(Lchpo1);
  for ichpo1 = 1:nchpo1
    chpo1i = Lchpo1{ichpo1};
    comp2 = 0;
    clear chpo2i chpoe2;
    nbzone1 = length(chpo1i);
    for zo1 = 1:nbzone1
      list_node2 = nmail2{zo1}.MAIL';
      chpoe1 = chpo1i{zo1};
      nbcomp1 = length(chpoe1);
      for comp1 = 1:nbcomp1
        comp2 = comp2 + 1;
        nbval = size(chpoe1{comp1}.XVAL,2);
        xval2 = zeros(nbnot2,nbval);
        [junk,ia,ib] = intersect(list_node2,list_node1new);
%       junk=list_node2(ia)=list_node1new(ib)
        xval2(ib,:) = chpoe1{comp1}.XVAL(ia,:);
        chpoe2{comp2} = struct('COMP',chpoe1{comp1}.COMP, ...
                               'UNIT',chpoe1{comp1}.UNIT, ...
                               'XVAL',xval2);
        clear xval2;
      end
      clear chpoe1;
    end
    chpo2i{1} = chpoe2;
    Lchpo2{ichpo1} = chpo2i; clear chpo2i;
    clear chpo1i;
  end
end
if ~test_list
  chpo2 = Lchpo2{1}; Lchpo2 = [];
end

% extend the nmail accordingly
nmail3 = [];
if ~isempty(nmail0)
  clear nmail3;
  zo1 = 1;
    nmail3{zo1} = struct('TYPE','POI1','MAIL',list_node1new');
end

if ~isempty(Lchpo2)
  varargout(1) = {Lchpo2};
end

switch test_list
  case -2,
%%DD Ajout 28/02/2006
    varargout(1) = {[]};
  case -1,
  case 0,
    varargout(1) = {chpo2};
  case 1,
    varargout(1) = {Lchpo2};
end
