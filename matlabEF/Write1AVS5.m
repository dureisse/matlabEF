function [error1] = Write1AVS5(xcoor,mail,chpo1,chml,cara,fid)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002


% Ecriture d'une passe d'un fichier AVS

% Input
%   xcoor(nbno,3) : coordonnees (si 2D la troisieme est nulle)
%   mail : maillage
%     mail{i} : maillage elementaire de la sous-zone i
%     mail{i}.TYPE : type d'elements
%     mail{i}.MAIL(nbeli,nbnni) : connectivite
%   chpo1 : champ par point
%     Il doit n'avoir qu'une zone, et etre defini sur tous les
%     noeuds de xcoor (limitation due a AVS)
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
%  fid : file identifier (file must be opened, see fopen)

% Output

% Remarques sur la numerotation : 
%   les noeuds sont numerotes dans leur ordre d'apparition dans xcoor
%   l'ensemble de ces noeuds forme un nuage de points qui sous-tend le
%     champ par point
%   les elements sont numerotes dans leur ordre d'apparition dans les
%     sous-zones successives
%   le maillage mail sous-tend les champs par elements

error1 = 0;
GlobalVar

% number of nodes
  num_nodes = size(xcoor,1);
% number of cells
  num_cells = 0; 
  icel = 0; list_ZoneNum = zeros(0,2);
  for ielem=1:size(mail,2)
    nelem = size(mail{ielem}.MAIL,1);
%%    list_ZoneNum(1:nelem,1) = ielem;
%%    list_ZoneNum(1:nelem,2) = [1:nelem]';
    list_ZoneNum(num_cells+1:num_cells+nelem,1) = ielem;
    list_ZoneNum(num_cells+1:num_cells+nelem,2) = [1:nelem]';
    num_cells = num_cells + nelem;
  end
% number of models
  num_modls = 1;
% length of the data associated with nodes
  num_ndata = 0;
  nbzone1 = length(chpo1);
  if nbzone1 == 0
    chpo = [];
  elseif nbzone1 == 1
    chpo = chpo1{1};
  else
    error('Only 1 subzone admitted for CHPO')
  end
  num_comp = size(chpo,2);
  for icomp=1:num_comp
    num_ndata = num_ndata + size(chpo{icomp}.XVAL,2);
  end
% length of the data associated with cells
  num_cdata = 0;
  num_comp = size(chml,2);
  for icomp=1:num_comp
    num_cdata = num_cdata + size(chml{icomp}.XVAL,2);
  end
% length of the data associated with the model
  num_mdata = 0;
  num_comp = size(cara,2);
  for icomp=1:num_comp
    num_mdata = num_mdata + size(cara{icomp}.XVAL,2);
  end

% Header
  line = zeros(5,1);
  line(1,1) = num_nodes; % number of nodes
  line(2,1) = num_cells; % number of cells
  line(3,1) = num_ndata; % length of the data associated with nodes
  line(4,1) = num_cdata; % length of the data associated with cells
  line(5,1) = num_mdata; % length of the data associated with the model
%  format = '%i %i %i %i %i\n'; size1 = 5;
  format = '%6i %5i %5i %5i %5i\n'; size1 = 5;
  count = fprintf(fid,format,line);

% Nodes
%  format = '%5i  %e  %e  %e\n'; size1 = 4;
%DD 29/07/2002  format = '%5i  %12.7E  %12.7E  %12.7E\n'; size1 = 4;
  format = '%5i % 12.7E % 12.7E % 12.7E\n'; size1 = 4;
  for ino=1:num_nodes
    line = zeros(4,1);
    node_id = ino;
    line(1,1) = node_id; % node-id
    line(2:4,1) = xcoor(node_id,:)'; % coordinates
    count = fprintf(fid,format,line);
  end

% Loop on cells
  for icel=1:num_cells
    line = zeros(2,1);
    cell_id = icel; % cell-id
    line(1,1) = cell_id; 
    mat_id = 1;
    line(2,1) = mat_id; % material
%    format = '%i %i'; size1 = 2;
%DD 28/07/2002    format = '%6i %3i'; size1 = 2;
    format = '%6i %2i'; size1 = 2;
    count = fprintf(fid,format,line);
%%%%    itype = list_ZoneNum(icel,1);
    izo = list_ZoneNum(icel,1);
    inum  = list_ZoneNum(icel,2);
    cell_type = mail{izo}.TYPE;
%DD 09/04/2002    line = cell_type;
    [junk,junk1,itype] = intersect(mail{izo}.TYPE,list_type_C3M1);
    if isempty(itype)
      [junk,junk1,itype] = intersect(mail{izo}.TYPE,list_type_C3M2);
    end
    line = list_type_AVS{itype};
%    format = '%s'; size1 = 1;
    format = '%6s'; size1 = 1;
    count = fprintf(fid,format,line);
%DD 09/04/2002    itype = IsWordMember(mail{izo}.TYPE,list_type_AVS);
%DD 28/07/2002    itype = IsWordMember(mail{izo}.TYPE,list_type_C3M1);
    [junk,junk1,itype] = intersect(mail{izo}.TYPE,list_type_C3M1);
    if isempty(itype)
      [junk,junk1,itype] = intersect(mail{izo}.TYPE,list_type_C3M2);
      format = list_format2{itype};
      size1  = list_size2{itype};
    else
      format = list_format1{itype};
      size1  = list_size1{itype};
    end
    nbnn = size1;
    line = zeros(nbnn,1);
    for ino = 1:nbnn
      line(ino,1) = mail{izo}.MAIL(inum,ino);
    end
    count = fprintf(fid,format,line);
  end

% Data vector associated with nodes
  if (num_ndata)
    node_xval = zeros(num_nodes,num_ndata);
    clear node_size_comps node_comp_label node_units_label;
    node_num_comp = size(chpo,2);
    itemp0 = 0;
    for icomp=1:node_num_comp
      chpoi = chpo{icomp};
      itemp = itemp0 + size(chpoi.XVAL,2);
      temp = chpoi.XVAL;
      node_xval(:,itemp0+1:itemp) = temp;
      node_size_comps(icomp)  = size(chpoi.XVAL,2);
      node_comp_label{icomp}  = chpoi.COMP;
      node_units_label{icomp} = chpoi.UNIT;
      itemp0 = itemp;
    end
    Write2AVS(fid, ...
      node_size_comps,node_comp_label,node_units_label,node_xval);
  end

% Data vector associated with cells
  if (num_cdata)
    cell_xval = zeros(num_cells,num_cdata);
    clear cell_size_comps cell_comp_label cell_units_label;
    cell_num_comp = size(chml,2);
    itemp0 = 0;
    for icomp=1:cell_num_comp
      chmli = chml{icomp};
      itemp = itemp0 + size(chmli.XVAL,2);
      temp = chmli.XVAL;
      cell_xval(:,itemp0+1:itemp) = temp;
      cell_size_comps(icomp)  = size(chmli.XVAL,2);
      cell_comp_label{icomp}  = chmli.COMP;
      cell_units_label{icomp} = chmli.UNIT;
%% DD 08/05/2003 il manquait la ligne suivante
      itemp0 = itemp;
    end
    Write2AVS(fid, ...
      cell_size_comps,cell_comp_label,cell_units_label,cell_xval);
  end

% Single-model based data
  if (num_mdata)
    modl_xval = zeros(num_modls,num_mdata);
    clear modl_size_comps modl_comp_label modl_units_label;
    modl_num_comp = size(cara,2);
    itemp0 = 0;
    for icomp=1:modl_num_comp
      carai = cara{icomp};
      itemp = itemp0 + size(carai.XVAL,2);
      temp = carai.XVAL;
      modl_xval(:,itemp0+1:itemp) = temp;
      modl_size_comps(icomp)  = size(carai.XVAL,2);
      modl_comp_label{icomp}  = carai.COMP;
      modl_units_label{icomp} = carai.UNIT;
    end
    Write2AVS(fid, ...
      modl_size_comps,modl_comp_label,modl_units_label,modl_xval);
  end
