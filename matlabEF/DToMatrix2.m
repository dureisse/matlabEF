function [D1] = DToMatrix2(d1,modl1,mail1,intg1, ...
                           numerp1,mapComp1,listComp1, ...
                           numerd2,mapComp2,listComp2);
% Assemble l'operateur de Hooke sous forme elementaire d1
%
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 03 / 09 / 2002
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 15 / 11 / 2002
%   modification des modeles version 1.3 a version 1.4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 12 / 2003
%   Modification pour assemblages speciaux
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 10 / 2007
%   Amelioration des performances
%
% La numerotation des valeurs dans le vecteur des deformations est donnee
% par (numerp1,mapComp1,listComp1), et dans le vecteur des contraintes
% par (numerd2,mapComp2,listComp2)
%
% Si des valeurs du premier sont necessaires mais n'existent
% pas, elles sont prises a 0, si des valeurs du deuxieme sont
% demandees mais n'existent pas, elles sont mises a 0.

ncomp1 = max(max(mapComp1));
ncomp2 = max(max(mapComp2));
D1 = sparse(ncomp2,ncomp1);

nbcomp1 = length(listComp1);
nbcomp2 = length(listComp2);

% For point numbering
ind = 0;

% On inverse une fois pour toute les listes
  numerd2_inv = InverseList(numerd2,max(numerd2));
  numerp1_inv = InverseList(numerp1,max(numerp1));

% Loop on zones
nbzone1 = length(modl1);
for zo1 = 1:nbzone1
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};
  intge1 = intg1{zo1};
  de1    = d1{zo1};

  nbptg  = length(intge1.WEIGHT); % Number of integration points

% On ne prend que les composantes qui existent dans listComp1
  name1 = modle1.COMP(modle1.NCOP);
  name1f = repmat(name1,1,nbptg);
  masq_pos2 = findoccur(name1f,listComp1);

% On ne prend que les composantes qui existent dans listComp2
  name1 = modle1.COMD(modle1.NCOD);
  name1f = repmat(name1,1,nbptg);
  masq_pos4 = findoccur(name1f,listComp2);


% Loop on elements
  nbel1 = size(maile1.MAIL,1);
  for el1 = 1:nbel1
    topo1 = [ind+1:ind+nbptg];
    topo2 = [ind+1:ind+nbptg];
    ind = ind + nbptg;


%   On ne prend que les points qui existent dans numerp1
    topo1f = reshape(repmat(topo1,nbcomp1,1),1,nbcomp1*nbptg);
%%DD 07/10/08     masq_pos1 = findoccur(topo1f,numerp1);
    ii = max(topo1f); if ii > length(numerp1_inv); numerp1_inv(ii) = 0; end
    masq_pos1 = numerp1_inv(topo1f);

%   On ne prend que l'intersection avec les composantes dans listComp1
    masq3 = masq_pos1 & masq_pos2;
    l1 = find(masq3);
    pos_in_numerp1 = masq_pos1(l1);
    pos_in_comp1   = masq_pos2(l1);
    list_comp1     = diag(mapComp1(pos_in_numerp1,pos_in_comp1));
    clear pos_in_numerp1 pos_in_comp1;


%   On ne prend que les points qui existent dans numerd2
    topo2f = reshape(repmat(topo2,nbcomp2,1),1,nbcomp2*nbptg);
%%DD 07/10/08    masq_pos1 = findoccur(topo2f,numerd2);
    ii = max(topo2f); if ii > length(numerd2_inv); numerd2_inv(ii) = 0; end
    masq_pos1 = numerd2_inv(topo2f);

%   On ne prend que l'intersection des deux
    masq3 = masq_pos1 & masq_pos4;
    l2 = find(masq3);
    pos_in_numerd2 = masq_pos1(l2);
    pos_in_comp2   = masq_pos4(l2);
    list_comp2     = diag(mapComp2(pos_in_numerd2,pos_in_comp2));
    clear pos_in_numerd2 pos_in_comp2;


    DE = de1.XVAL(l2,l1,el1);
%%DD 03/12/25    D1(list_comp2,list_comp1) = D1(list_comp2,list_comp1) + DE;
%   For special cases where several values go at the same position:
%   the result should sum the values and not override them
    for i = 1:size(DE,2)
      for j = 1:size(DE,1)
        D1(list_comp2(j),list_comp1(i)) = ...
          D1(list_comp2(j),list_comp1(i)) + DE(j,i);
      end
    end

  end

end
