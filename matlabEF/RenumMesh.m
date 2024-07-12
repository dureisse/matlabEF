function mail2 = RenumMesh(mail1,numer_inv)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002

% Renumerotation des noeuds d'un maillage

% Input
%   mail1 : maillage
%   numer_inv : liste de numerotation inverse
% Output
%   mail2 : maillage renumerote

% Si le numero du noeud ino de l'element iel de la zone izo
% est num = mail1{izo}.MAIL(iel,ino),
% alors il devient numer_inv(num)

nb_zone = size(mail1,2);
mail2 = [];
%clear mail2;
for izo = 1:nb_zone
  type1 = mail1{izo}.TYPE;
  num1 = mail1{izo}.MAIL;
  num2 = 0 * num1;
  nbel = size(num1,1);
  for iel = 1:nbel
    num2(iel,:) = numer_inv(num1(iel,:));
  end
  mail2{izo} = struct('MAIL',num2,'TYPE',type1);
end
