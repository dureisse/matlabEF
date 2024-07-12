function [xcoor,mail,chpo1,chml,cara,error1] = Read1AVS6(fid)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 28 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2002
%  pb de fin de fichier et echo


% Lecture d'une passe d'un fichier AVS - UCD
% Attention : pour AVS, meme en 2D il y a 3 coordonnees par noeud,
%   et la troisieme est nulle. Il faut donc faire ensuite 
%   xcoor = xcoor(:,1:2)
% Lecture d'elements quadratiques

% Input
%  fid : file identifier (file must be opened, see fopen)

% Output
%   xcoor(nbno,3) : coordonnees (si 2D la troisieme est nulle)
%   mail : maillage
%     mail{i} : maillage elementaire de la sous-zone i
%     mail{i}.TYPE : type d'elements
%     mail{i}.MAIL(nbeli,nbnni) : connectivite
%   chpo1 : champ par point (a 1 seule sous-zone par limitation d'AVS)
%     chpo1{1}{i}.COMP : nom de la composante numero i
%     chpo1{1}{i}.UNIT : nom de l'unite de la composante numero i
%     chpo1{1}{i}.XVAL(nbno,nbval) : valeurs
%   chml : champ constant par element
%     chml{i}.COMP : nom de la composante numero i
%     chml{i}.UNIT : nom de l'unite de la composante numero i
%     chml{i}.XVAL(nbel,nbval) : valeurs
%   cara : champ constant de caracteristiques (modele selon AVS)
%     cara{i}.COMP : nom de la composante numero i
%     cara{i}.UNIT : nom de l'unite de la composante numero i
%     cara{i}.XVAL(nbel,nbval) : valeurs
%
% Remarques sur la numerotation : 
%   les noeuds sont numerotes dans leur ordre d'apparition dans xcoor
%   l'ensemble de ces noeuds forme un nuage de points qui sous-tend les
%     champs par point
%   les elements sont numerotes dans leur ordre d'apparition dans les
%     sous-zones successives
%   le maillage mail sous-tend les champs par elements

echo = 0;

error1 = 0;
GlobalVar

% Skip comments
  line1 = '#'; while (line1(1) == '#'),
    line1 = fgetl(fid);
    if size(line1,2) == 0
      line1 = fgetl(fid);
      if size(line1,2) == 0
%       End of file
        disp('EOF')
        error1 = 1;
        return;
      end
    end
  end

% Header
  format = '%i %i %i %i %i\n'; size1 = 5;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line')
end
  if (size(errmsg,1) ~= 0)
    error1 = errmsg;
    return;
  end
  num_nodes = line(1,1); % number of nodes
  num_cells = line(2,1); % number of cells
  num_ndata = line(3,1); % length of the data associated with nodes
  num_cdata = line(4,1); % length of the data associated with cells
  num_mdata = line(5,1); % length of the data associated with the model
if echo
disp(line')
end

% Loop on nodes
  xcoor = zeros(0,3);
%%BUG DD 2023/10/05  format = '%5i  %e  %e  %e\n'; size1 = 4;
  if (num_nodes > 999999)
      error('pb format fichier')
  else
      format = '%6i  %e  %e  %e\n'; size1 = 4;
  end
  for ino=1:num_nodes
    [line,count] = fscanf(fid,format,size1);
if echo
disp(line')
end
    node_id          = line(1,1); % node-id
    xcoor(node_id,:) = line(2:4,1)'; % coordinates
  end

% Loop on cells
  nb_zone1 = size(list_type_AVS,2);
%%  nb_zones = 2 * nb_zone1; % For eventual quadratic elements
  nb_zones = size(list_type_C3M1,2) + size(list_type_C3M2,2);
  if (nb_zones ~= 2 * nb_zone1)
    error('problem in common data for element names')
  end
%  Initialize zones to void
  for izo=1:nb_zones
    mail_elem{izo} = [];
  end
  list_ZoneNum = zeros(num_cells,2);
  for icel=1:num_cells
    format = '%i %i'; size1 = 2;
    [line,count] = fscanf(fid,format,size1);
if echo
disp(line')
end
    cell_id = line(1,1); % cell-id
    mat_id  = line(2,1); % material
    format = '%s'; size1 = 1;
    [line2,count] = fscanf(fid,format,size1);
    cell_type = line2; % type
    [junk,junk1,itype] = intersect(cell_type,list_type_AVS);
    if isempty(itype)
      error1 = ['type of cell not implemented: ' cell_type];
      disp(error1)
      return;
    end
%
    linetext = fgets(fid);
if echo
disp([int2str(line') ' ' line2 linetext])
end
    format = list_format1{itype};
    size1  = list_size1{itype}; % For "linear" elements
    size2  = list_size2{itype}; % For "quadratic" elements, size2 > size1
%DD    [line,count] = fscanf(fid,format,size2);
    [line,count] = sscanf(linetext,'%i');
if echo
disp(line')
end
    if (count == size1)
%     Nothing to do      
    elseif (count == size2)
      itype = itype + nb_zone1;
      size1 = size2;
    else
      error1 = ['bad number of nodes for cell: ' cell_type];
      error(error1)
      return;
    end
    list_ZoneNum(icel,1) = itype;
    inum = size(mail_elem{itype},1) + 1;
    list_ZoneNum(icel,2) = inum;
    nbnn = size1;
    if ((itype == 13) && (size1 == 10))
%     Pour TE10, les noeuds ne sont pas dans le bon ordre
%     AVS, Cast3M et gmsh codent differemment. On choisi ici gmsh.
%% DD modif 08/12/2015      recode = [1 6 10 2 5 9 3 8 4 7]; % AVS to gmsh
      recode = [1 2 3 4 5 8 6 7 9 10]; % AVS to gmsh
    else
      recode = [1:nbnn];
    end
    for ino = 1:nbnn
      mail_elem{itype}(inum,ino) = line(recode(ino),1);
    end
  end
% New global numbering of cells (according to ordered zones)
  clear new_zonenum;
  new_sum = 0;
  for izo=1:nb_zones
    new_zonenum(izo) = new_sum;
    new_sum = new_sum + size(mail_elem{izo},1);
  end
  clear new_num;
  for icel=1:num_cells
    izo = list_ZoneNum(icel,1);
    inu = list_ZoneNum(icel,2);
    new_num(icel) = new_zonenum(izo) + inu;
  end
  clear new_zonenum list_ZoneNum;
% Mesh
  mail = [];
  ielem = 0;
  for izo=1:nb_zones
    if size(mail_elem{izo},1)
      ielem = ielem + 1;
      if (izo <= size(list_type_C3M1,2))
	 el_type = list_type_C3M1{izo}; 
      else
	 el_type = list_type_C3M2{izo-size(list_type_C3M1,2)};
      end
      mail{ielem} = struct('MAIL',mail_elem{izo},'TYPE',el_type); 
    end
  end

% Data vector associated with nodes
  if (num_ndata)
    clear chpo1 chpo;
    [node_size_comps,node_comp_label,node_units_label,node_xval] = ...
      Read2AVS2(fid,num_ndata,num_nodes,echo);
    node_num_comp = size(node_size_comps,2);
    itemp0 = 0;
    for icomp=1:node_num_comp
      itemp = itemp0 + node_size_comps(icomp);
      temp = node_xval(:,itemp0+1:itemp);
      chpoi = struct('XVAL',temp);
      chpoi.COMP = node_comp_label{icomp};
      chpoi.UNIT = node_units_label{icomp};
      itemp0 = itemp;
      chpo{icomp} = chpoi;
    end
    chpo1{1} = chpo;
  else
    chpo1 = [];
  end

% Data vector associated with cells
  if (num_cdata)
    [cell_size_comps,cell_comp_label,cell_units_label,cell_xval1] = ...
      Read2AVS2(fid,num_cdata,num_cells,echo);
%   new cell numbering
    cell_xval = zeros(size(cell_xval1));
    cell_xval(new_num,:) = cell_xval1;
    clear cell_xval1;
    cell_num_comp = size(cell_size_comps,2);
    itemp0 = 0;
    for icomp=1:cell_num_comp
      itemp = itemp0 + cell_size_comps(icomp);
      temp = cell_xval(:,itemp0+1:itemp);
      chmli = struct('XVAL',temp);
%DD 28/07/2002      chmli.COMP = node_comp_label{icomp};
%DD 28/07/2002      chmli.UNIT = node_units_label{icomp};
      chmli.COMP = cell_comp_label{icomp};
      chmli.UNIT = cell_units_label{icomp};
      itemp0 = itemp;
      chml{icomp} = chmli;
    end
  else
    chml = [];
  end

% Single-model based data
  if (num_mdata)
    [modl_size_comps,modl_comp_label,modl_units_label,modl_xval] = ...
      Read2AVS2(fid,num_mdata,1,echo);
    modl_num_comp = size(modl_size_comps,2);
    itemp0 = 0;
    for icomp=1:modl_num_comp
      itemp = itemp0 + modl_size_comps(icomp);
      temp = modl_xval(:,itemp0+1:itemp);
      carai = struct('XVAL',temp);
      carai.COMP = modl_comp_label{icomp};
      carai.UNIT = modl_units_label{icomp};
      itemp0 = itemp;
      cara{icomp} = carai;
    end
  else
    cara = [];
  end

% Should be empty string
% (il faut tester si cette ligne est necessaire... DD 28/07/2002)
%  line1 = fgetl(fid);
%if echo
%disp(line1)
%end

