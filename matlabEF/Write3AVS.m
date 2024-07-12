function [error1] = Write3AVS(xcoor,mail,Lchpo1,nmail, ...
                              Lchml,Lcara,ListInd,fid)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David    L.M.G.C. SYSTEMES MULTICONTACT le 20 / 10 / 2003
%   Cas des listes de champs
%
% Ecriture de plusieurs passes d'un fichier AVS,
% associees a des listes de champs.
%
% Input
%   xcoor(nbno,3) : coordonnees (si 2D la troisieme est nulle)
%   mail : maillage
%     mail{i} : maillage elementaire de la sous-zone i
%     mail{i}.TYPE : type d'elements
%     mail{i}.MAIL(nbeli,nbnni) : connectivite
%   Lchpo1(nl) : liste de champs par point
%     Il doit n'avoir qu'une zone, et etre defini sur tous les
%     noeuds de xcoor (limitation due a AVS)
%     chpo1 = Lchpo1{j} : champ par point numero j
%     chpo1{1}{i}.COMP : nom de la composante numero i
%     chpo1{1}{i}.UNIT : nom de l'unite de la composante numero i
%     chpo1{1}{i}.XVAL(nbno,nbval) : valeurs
%   nmail : maillage POI1 sous-tendant tous les champs par point
%   Lchml(nl) : liste de champs constants par element
%     chml = Lchml{j} : champ constant par element numero j
%     chml{i}.COMP : nom de la composante numero i
%     chml{i}.UNIT : nom de l'unite de la composante numero i
%     chml{i}.XVAL(nbel,nbval) : valeurs
%   Lcara(nl) : liste de champs constants de caracteristiques
%               (modele selon AVS)
%     cara = Lcara{j} : champ constants de caracteristiques j
%     cara{i}.COMP : nom de la composante numero i
%     cara{i}.UNIT : nom de l'unite de la composante numero i
%     cara{i}.XVAL(nbel,nbval) : valeurs
%  ListInd(nl) : liste d'indices (temps par exemple)
%  fid : file identifier (file must be opened, see fopen)
%
% Output
%  error1	retour non nul en cas d'erreur
%
% Remarques sur la numerotation : 
%   les noeuds sont numerotes dans leur ordre d'apparition dans xcoor
%   l'ensemble de ces noeuds forme un nuage de points qui sous-tend la
%     liste de champs par point.
%   les elements sont numerotes dans leur ordre d'apparition dans les
%     sous-zones successives
%   le maillage mail sous-tend les champs par elements
%   Tous les champs d'une liste ont la meme structure.
% Voir Write1AVS5

error1 = 0;
GlobalVar

% Controls for consistency of datas
% """""""""""""""""""""""""""""""""

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

% number of fields in the list
  nl = length(ListInd);

% length of the data associated with nodes
  num_ndata = 0;

  if length(Lchpo1) == 0
    nbzone1 = 0;  % number of subzones
    num_comp = 0; % number of components
    chpo = [];
    Lnum_ndata = [];
  elseif nl == length(Lchpo1)
    j = 1;
      chpo1 = Lchpo1{j};
      nbzone1 = length(chpo1);
      num_comp = size(chpo1,2);
      chpo = chpo1{1};
      for icomp = 1:num_comp
        Lnum_ndata(icomp) = size(chpo{icomp}.XVAL,2);
      end
    for j = 2:nl
      chpo1 = Lchpo1{j};
      if nbzone1 ~= length(chpo1);
        nbzone1
        j
        length(chpo1)
        error('All point-wise fields do not have same subzone number')
      end
      if num_comp ~= size(chpo1,2);
        num_comp
        j
        size(chpo1,2)
        error('All point-wise fields do not have same component number')
      end
      chpo = chpo1{1};
      for icomp = 1:num_comp
        if Lnum_ndata(icomp) ~= size(chpo{icomp}.XVAL,2)
          Lnum_ndata(icomp)
          j
          icomp
          size(chpo{icomp}.XVAL,2)
          error('All point-wise fields do not have the same size')
        end
      end
    end
    if nbzone1 ~= 1
      error('Only 1 subzone admitted for CHPO')
    end
  else
    nl
    length(Lchpo1)
    error('Bad list CHPO length')
  end

  for icomp=1:num_comp
    num_ndata = num_ndata + Lnum_ndata(icomp);
  end

% length of the data associated with cells
  num_cdata = 0;

  if length(Lchml) == 0
    num_comp = 0; % number of components
    Lnum_cdata = [];
  elseif nl == length(Lchml)
    j = 1;
      chml = Lchml{j};
      num_comp = size(chml,2);
      for icomp = 1:num_comp
        Lnum_cdata(icomp) = size(chml{icomp}.XVAL,2);
      end
    for j = 2:length(Lchml)
      chml = Lchml{j};
      if num_comp ~= size(chml,2);
        num_comp
        size(chml,2)
        error('All element-wise cst fields do not have the same comp')
      end
      for icomp = 1:num_comp
        if Lnum_cdata(icomp) ~= size(chml{icomp}.XVAL,2)
          Lnum_cdata(icomp)
          j
          icomp
          size(chml{icomp}.XVAL,2)
          error('All element-wise fields do not have the same size')
        end
      end
    end
  else
    nl
    length(Lchml)
    error('Bad list CHML length')
  end

  for icomp=1:num_comp
    num_cdata = num_cdata + Lnum_cdata(icomp);
  end

% length of the data associated with the model
  num_mdata = 0;

  if length(Lcara) == 0
    num_comp = 0; % number of components
    Lnum_mdata = [];
  elseif nl == length(Lcara)
    j = 1;
      cara = Lcara{j};
      num_comp = size(cara,2);
      for icomp = 1:num_comp
        Lnum_mdata(icomp) = size(cara{icomp}.XVAL,2);
      end
    for j = 2:length(Lchml)
      cara = Lcara{j};
      if num_comp ~= size(cara,2);
        num_comp
        size(cara,2)
        error('All caracteristic cst fields do not have the same comp')
      end
      for icomp = 1:num_comp
        if Lnum_mdata(icomp) ~= size(cara{icomp}.XVAL,2)
          Lnum_mdata(icomp)
          j
          icomp
          size(cara{icomp}.XVAL,2)
          error('All caracteristic fields do not have the same size')
        end
      end
    end
  else
    nl
    length(Lcara)
    error('Bad list CARA length')
  end

% Calls to the writing routine
% """"""""""""""""""""""""""""
j = 1;
  disp(['Writing for element ' int2str(j) ' / ' int2str(nl) ...
        ' in the list'])
  if length(Lchpo1) == 0
    chpo1 = [];
  else
    chpo1 = Lchpo1{j};
  end
  if length(Lchml) == 0
    chml1 = [];
  else
    chml1 = Lchml{j};
  end
  clear cara1;
  cara1{1} = struct('COMP','time','UNIT','','XVAL',ListInd(j));
  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS3(xcoor,mail,nmail,chpo1);
  error1 = Write1AVS5(xcoor1a,mail1a,chpo1a,chml1,cara1,fid);
  if error1
    disp('Problem in Write1AVS5')
    return
  end
for j = 2:nl
  disp(['Writing for element ' int2str(j) ' / ' int2str(nl) ...
        ' in the list'])
  if length(Lchpo1) == 0
    chpo1 = [];
  else
    chpo1 = Lchpo1{j};
  end
  if length(Lchml) == 0
    chml1 = [];
  else
    chml1 = Lchml{j};
  end
  clear cara1;
  cara1{1} = struct('COMP','time','UNIT','','XVAL',ListInd(j));
  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS3(xcoor,[],nmail,chpo1);
  error1 = Write1AVS5(xcoor1a,mail1a,chpo1a,chml1,cara1,fid);
  if error1
    disp('Problem in Write1AVS5')
    return
  end
end
