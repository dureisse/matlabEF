function [xcoor,mail,error1] = ReadMeshGMSH(fid)
% Read a list of meshes in a gmsh ascii meshing file (format 1.0)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2006
%
% Inputs
%   fid			file identifier (file must be opened, see fopen)
% Output
%   xcoor1(nbno,3)	coordonnees
%   ListMesh1{nmesh1}	List of meshes
%   error1		

GlobalVar
error1 = 0;
echo = 1;
echo = 0;

% Nodes
% """""
% Header
  line1 = fgetl(fid);
  format = '%s\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  if (size(errmsg,1) ~= 0)
    error1 = errmsg;
    return;
  end
  if ~strcmp(line,'$NOD')
    line
    error('Unable to read node header')
  end

  line1 = fgetl(fid);
  format = '%i\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  num_nodes = line;

% Loop on nodes
  xcoor = zeros(num_nodes,3);
  format = '%i  %e  %e  %e\n'; size1 = 4;
  for ino=1:num_nodes
    [line,count] = fscanf(fid,format,size1);
if echo
disp(line')
end
    node_id          = line(1,1); % node-id
    xcoor(node_id,:) = line(2:4,1)'; % coordinates
  end

% Trailer
  line1 = fgetl(fid);
  format = '%s\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  if (size(errmsg,1) ~= 0)
    error1 = errmsg;
    return;
  end
  if ~strcmp(line,'$ENDNOD')
    line
    error('Unable to read node header')
  end

% Elements
% """"""""
% Header
  line1 = fgetl(fid);
  format = '%s\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  if (size(errmsg,1) ~= 0)
    error1 = errmsg;
    return;
  end
  if ~strcmp(line,'$ELM')
    line
    error('Unable to read element header')
  end

  line1 = fgetl(fid);
  format = '%i\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  num_cells = line;

% Loop on cells
  nb_zone1 = size(list_type_AVS,2);
  nb_zones = size(list_type_C3M1,2) + size(list_type_C3M2,2);
  if (nb_zones ~= 2 * nb_zone1)
    error('problem in common data for element names')
  end
% Initialize zones to void
  clear mail_elem;
  for izo=1:nb_zones
    mail_elem{izo} = [];
  end
  list_ZoneNum = zeros(num_cells,2);
  for icel=1:num_cells
    format = '%i %i %i %i %i'; size1 = 5;
    [line,count] = fscanf(fid,format,size1);
if echo
line0 = line;
end
    cell_id = line(1,1); % cell-id
    cell_type_GMSH = line(2,1); % cell-type
    cell_type = list_type_GMSH{cell_type_GMSH};
    reg_phys = line(3,1); % not yet used
    reg_elem = line(4,1); % not yet used
    nbnn = line(5,1); % number of nodes

    itype = findoccur({cell_type},list_type_C3M1);
    if itype
      format = list_format1{itype}; size1 = list_size1{itype};
    else
      itype = findoccur({cell_type},list_type_C3M2);
      format = list_format2{itype}; size1 = list_size2{itype};
      if itype
        itype = itype + size(list_type_C3M1,2);
      else
        cell_type
        error('Cell type not implemented')
      end
    end
    if nbnn~=size1
      cell_type
      nbnn
      size1
      error('Pas compris')
    end
    [line,count] = fscanf(fid,format,size1);
if echo
disp([line0' line'])
end
    list_ZoneNum(icel,1) = itype;
    inum = size(mail_elem{itype},1) + 1;
    list_ZoneNum(icel,2) = inum;
    nbnn = size1;
    for ino = 1:nbnn
      mail_elem{itype}(inum,ino) = line(ino,1);
    end
  end

% Trailer
  line1 = fgetl(fid);
  format = '%s\n'; size1 = 1;
  [line,count,errmsg,nextindex] = sscanf(line1,format,size1);
if echo
disp(line)
end
  if (size(errmsg,1) ~= 0)
    error1 = errmsg;
    return;
  end
  if ~strcmp(line,'$ENDELM')
    line
    error('Unable to read node header')
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

