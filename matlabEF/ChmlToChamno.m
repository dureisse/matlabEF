function chamno1 = ChmlToChamno(chml1,mail1)
% Transforms an element-constant field into a field defined at nodes
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 12 / 12 / 2006
%
% Inputs
%   chml1	Constant-per-element field
%   mail1	Mesh
% Ouputs
%   chamno1	same field but defined at node of elements
%
% Transforme un champ constant par element chml1 en un
% champ par element chamno1 aux noeuds
% Evite l'interpolation necessaire si on passe d'un Chml en un Cham
% puis en un Chamno

nbcomp1 = length(chml1);
nbzone1 = length(mail1);

clear chamno1;

elt1 = 0;
for zo1 = 1:nbzone1
  nbel1 = size(mail1{zo1}.MAIL,1); % Nombre d'elements de la sous-zone
  nbno1 = size(mail1{zo1}.MAIL,2); % Nombre de noeuds par element
  for comp1 = 1:nbcomp1
    nbval  = size(chml1{comp1}.XVAL,2);
    nbvalt = nbval * nbno1;
    xval1  = zeros(nbel1,nbvalt);
    xval2  = chml1{comp1}.XVAL(elt1+1:elt1+nbel1,:);
    for ptg1 = 1:nbno1
      xval1(:,nbval*(ptg1-1)+1:nbval*(ptg1-1)+nbval) = xval2;
    end
    chame1{comp1} = struct('COMP',chml1{comp1}.COMP, ...
			   'UNIT',chml1{comp1}.UNIT, ...
			   'XVAL',xval1);
  end
  elt1 = elt1 + nbel1;
  chamno1{zo1} = chame1;
end
