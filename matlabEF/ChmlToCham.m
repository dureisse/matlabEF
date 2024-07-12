function cham1 = ChmlToCham(chml1,mail1,intg1)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002


% Transforme un champ constant par element chml1 en un
% champ par element cham1 qui s'appuie sur le segment
% d'integration intg1 et le maillage mail1

nbcomp1 = length(chml1);
nbzone1 = length(intg1);
if (nbzone1 ~= length(mail1))
  nbzone1
  length(mail1)
  error('Bad number of subzones')
end


clear cham1;

elt1 = 0;
for zo1 = 1:nbzone1
  clear chame1;
  nbptg = length(intg1{zo1}.WEIGHT);
  nbel1 = size(mail1{zo1}.MAIL,1); % Nombre d'elements de la sous-zone
  for comp1 = 1:nbcomp1
    nbval  = size(chml1{comp1}.XVAL,2);
    nbvalt = nbval * nbptg;
    xval1  = zeros(nbel1,nbvalt);
    xval2  = chml1{comp1}.XVAL(elt1+1:elt1+nbel1,:);
    for ptg1 = 1:nbptg
      xval1(:,nbval*(ptg1-1)+1:nbval*(ptg1-1)+nbval) = xval2;
    end
    chame1{comp1} = struct('COMP',chml1{comp1}.COMP, ...
			      'UNIT',chml1{comp1}.UNIT, ...
			      'XVAL',xval1);
  end
  elt1 = elt1 + nbel1;
  cham1{zo1} = chame1;
end
