function [K1] = RigiToMatrix(rigi1,modl1,mail1, ...
		             numerd1,mapddld1,listDdlDual1, ...
		             numerp1,mapddlp1,listDdlPrim1);
% Assembling of elementary matrices into a global one
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 02 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2003
%   Modification pour assemblages speciaux
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 10 / 2007
%   Amelioration des performances
%
% Assemble des rigidites elementaires (rigi1,modl1,mail1) dans une
% matrice K1 en utilisant les mappings
% mapddld1 pour les lignes (composantes duales) et
% mapddlp1 pour les colonnes (composantes primales)
%
% Inputs
%   rigi1		Elementary matrices
%   modl1		Model
%   mail1		Mesh
%   numerd1(nbnod)	Numbering of dual nodes
%   mapddld1(nbnod,nbcd)Mapping if dual dofs
%   listDdlDual1(nbcd)	List of dual dof names
%   numerp1(nbnop)	Numbering of primal nodes
%   mapddlp1(nbnop,nbcp)Mapping if primal dofs
%   listDdlPrim1(nbcp)	List of primal dof names
% Outputs
%   K1(nddld1,nddlp1)	Global matrix

nddld1 = max(max(mapddld1));
nddlp1 = max(max(mapddlp1));
K1 = sparse(nddld1,nddlp1);

% Boucle sur les sous-zones de la rigidite
nbzone1 = length(rigi1);
for zo1 = 1:nbzone1
  rigie1 = rigi1{zo1};
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};

% On ne prend que les ddl de NDDD qui existent dans listDdlDual1
  masq_pos2d = findoccur(modle1.DDLD(modle1.NDDD),listDdlDual1);
% On ne prend que les ddl de NDDP qui existent dans listDdlPrim1
  masq_pos2p = findoccur(modle1.DDLP(modle1.NDDP),listDdlPrim1);
% On inverse une fois pour toute les listes
  numerd1_inv = InverseList(numerd1,max(numerd1));
  numerp1_inv = InverseList(numerp1,max(numerp1));

% Boucle sur les elements de la rigidite
  nbel1 = size(rigie1.XVAL,3);
  for el1 = 1:nbel1
    topo1 = maile1.MAIL(el1,:);

%   On ne prend que les noeuds de NNOD qui existent dans numerd1
%%DD 07/10/08    masq_pos1d = findoccur(topo1(modle1.NNOD),numerd1);
    ii = max(topo1(modle1.NNOD)); if ii > length(numerd1_inv); numerd1_inv(ii) = 0; end
    masq_pos1d = numerd1_inv(topo1(modle1.NNOD));

%%DD 07/10/08 %   On ne prend que les ddl de NDDD qui existent dans listDdlDual1
%%DD 07/10/08     masq_pos2d = findoccur(modle1.DDLD(modle1.NDDD),listDdlDual1);

%   On ne prend que l'intersection des deux pour le dual
    masq3 = masq_pos1d & masq_pos2d;
    ld1 = find(masq3);
    pos_in_numerd1 = masq_pos1d(ld1);
    pos_in_compd1  = masq_pos2d(ld1);
    list_ddld1     = diag(mapddld1(pos_in_numerd1,pos_in_compd1));
    clear pos_in_numerd1 pos_in_compd1;

%   On ne prend que les noeuds de NNOP qui existent dans numerp1
%%DD 07/10/08    masq_pos1p = findoccur(topo1(modle1.NNOP),numerp1);
    ii = max(topo1(modle1.NNOP)); if ii > length(numerp1_inv); numerp1_inv(ii) = 0; end
    masq_pos1p = numerp1_inv(topo1(modle1.NNOP));

%%DD 07/10/08 %   On ne prend que les ddl de NDDP qui existent dans listDdlPrim1
%%DD 07/10/08     masq_pos2p = findoccur(modle1.DDLP(modle1.NDDP),listDdlPrim1);

%   On ne prend que l'intersection des deux pour le primal
    masq3 = masq_pos1p & masq_pos2p;
    lp1 = find(masq3);
    pos_in_numerp1 = masq_pos1p(lp1);
    pos_in_compp1  = masq_pos2p(lp1);
    list_ddlp1     = diag(mapddlp1(pos_in_numerp1,pos_in_compp1));
    clear pos_in_numerp1 pos_in_compp1;

    KE = rigie1.XVAL(ld1,lp1,el1);
%%DD 03/12/25    K1(list_ddld1,list_ddlp1) = K1(list_ddld1,list_ddlp1) + KE;
%   For special cases where several values go at the same position:
%   the result should sum the values and not override them
    for i = 1:length(lp1)
      for j = 1:length(ld1)
        K1(list_ddld1(j),list_ddlp1(i)) = ...
          K1(list_ddld1(j),list_ddlp1(i)) + KE(j,i);
      end
    end

  end

end
