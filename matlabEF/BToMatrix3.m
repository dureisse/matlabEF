function [B1] = BToMatrix3(b1,modl1,mail1,intg1, ...
                           numer1,mapddl1,listDdl1, ...
                           numer2,mapcom2,listCom2, ...
                           opti1)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 31 / 08 / 2002
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 15 / 11 / 2002
%   modification des modeles version 1.3 a version 1.4
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 16 / 11 / 2002
%   correction des ddl primal/dual
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2003
%   Modification pour assemblages speciaux
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 10 / 2007
%   Amelioration des performances
%
% Assemble l'operateur sous forme elementaire b1 qui fait passer
% par element d'un vecteur de champ par point de deplacement
% a un vecteur de champ par element de deformation.
%
% La numerotation des valeurs dans le premier vecteur est donnee
% par (numer1,mapddl1,listDdl1), et dans le deuxieme par
% (numer2,mapcom2,listCom2)
%
% Si des valeurs du premier sont necessaires mais n'existent
% pas, elles sont prises a 0, si des valeurs du deuxieme sont
% demandees mais n'existent pas, elles sont mises a 0.



nddl1 = max(max(mapddl1));
ncom1 = max(max(mapcom2));
B1 = sparse(ncom1,nddl1);


% For point numbering
ind = 0;

% On inverse une fois pour toutes les listes
  numer1_inv = InverseList(numer1,max(numer1));
  numer2_inv = InverseList(numer2,max(numer2));

% Loop on zones
nbzone1 = length(modl1);
for zo1 = 1:nbzone1
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};
  intge1 = intg1{zo1};
  be1    = b1{zo1};

  nbptg = length(intge1.WEIGHT);

% On ne prend que les composantes qui existent dans listCom2
  if strcmp(opti1,'PRIMAL')
    namec1 = modle1.COMP(modle1.NCOP);
    named1 = modle1.DDLP(modle1.NDDP);
    nnox = modle1.NNOP;
  elseif strcmp(opti1,'DUAL')
    namec1 = modle1.COMD(modle1.NCOD);
    named1 = modle1.DDLD(modle1.NDDD);
    nnox = modle1.NNOD;
  else
    opti1
    erreur('option unknown')
  end

% On ne prend que les composantes de NCO? qui existent dans listCom2
  name1f    = repmat(namec1,1,nbptg);
  masq_pos4 = findoccur(name1f,listCom2);
  nbcom2 = length(masq_pos4) / nbptg;
% On ne prend que les ddl de NDD? qui existent dans listDdl1
  masq_pos2 = findoccur(named1,listDdl1);

% Loop on elements
  nbel1 = size(maile1.MAIL,1);
  for el1 = 1:nbel1
    topo1 = maile1.MAIL(el1,:);
    topo2 = [ind+1:ind+nbptg];
    ind = ind + nbptg;

%   On ne prend que les noeuds de NNO? qui existent dans numer1
%%DD 07/10/08     masq_pos1 = findoccur(topo1(nnox),numer1);
    ii = max(topo1(nnox)); if (ii > length(numer1_inv)); numer1_inv(ii) = 0; end
    masq_pos1 = numer1_inv(topo1(nnox));

%   On ne prend que l'intersection avec les composantes existantes
    masq3 = masq_pos1 & masq_pos2;
    lp1 = find(masq3);
    pos_in_numer1 = masq_pos1(lp1);
    pos_in_compp1 = masq_pos2(lp1);
    list_ddlp1    = diag(mapddl1(pos_in_numer1,pos_in_compp1));
    clear pos_in_numer1 pos_in_compp1;

%   On ne prend que les points qui existent dans numer2
    topo2f = reshape(repmat(topo2,nbcom2,1),1,nbcom2*nbptg);
%%DD 07/10/08     masq_pos1 = findoccur(topo2f,numer2);
    ii = max(topo2f); if ii > length(numer2_inv); numer2_inv(ii) = 0; end
    masq_pos1 = numer2_inv(topo2f);

%   On ne prend que l'intersection des deux
    masq3 = masq_pos1 & masq_pos4;
    l2 = find(masq3);
    pos_in_numer2 = masq_pos1(l2);
    pos_in_comp2  = masq_pos4(l2);
    list_comp2    = diag(mapcom2(pos_in_numer2,pos_in_comp2));
    clear pos_in_numer2 pos_in_comp2;

    BE = be1.XVAL(l2,lp1,el1);;
%%DD 03/12/25    B1(list_comp2,list_ddlp1) = B1(list_comp2,list_ddlp1) + BE;
%   For special cases where several values go at the same position:
%   the result should sum the values and not override them
    for i = 1:size(BE,2)
      for j = 1:size(BE,1)
        B1(list_comp2(j),list_ddlp1(i)) = ...
          B1(list_comp2(j),list_ddlp1(i)) + BE(j,i);
      end
    end

  end

end
