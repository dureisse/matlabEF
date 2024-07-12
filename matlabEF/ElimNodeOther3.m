function [pile_mail2,ListUsedNodes] = ElimNodeOther3(pile_mail1,pile_mail, ...
                                      xcoor,lxcrit);
% DUREISSEIX David   L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS le 22 / 01 / 2003
%   autant de valeurs de critere que de maillages dans la pile
%   extension au 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS le 25 / 02 / 2003
%   elimination d'une pile par rapport a une autre
%
% Elimination des noeuds a meme emplacement entre les maillages de 
% la pile de maillages pile_mail1 et ceux de pile_mail.
% Seule la premiere pile est modifiee.
% (voir aussi ElimNode3)
%
% Input
%   pile_mail1 : pile de maillages a modifier
%     pile_mail1{i} : maillage numero i
%   pile_mail : pile de maillages
%     pile_mail{i} : maillage numero i
%   xcoor(nbno,idim) : coordonnees des noeuds
%   lxcrit(i) : critere de proximite de noeuds pour le maillage i
%
% Output
%  pile_mail2 : pile de maillages a noeuds elimines
%               (en remplacement de pile_mail1)
%  ListUsedNodes : liste des noeuds effectivement utilises dans pile_mail2
%
% Remarque : si le probleme se pose (lxcrit non decroissant), la
% regle pour trouver le remplacant d'un noeud parmi plusieurs est :
% le plus recent ancien pris en compte.

A FAIRE

% On construit la liste des noeuds liste
% et l'index de remplacement replace_node, simultanement
liste = zeros(0,1);
clear replace_node;

% Boucle sur les maillages
nb_mail = size(pile_mail,2);
for imail = 1:nb_mail
  mail = pile_mail{imail};
  nmail = ChangeMesh2(mail,'POI1');
  nlist = nmail{1}.MAIL;
  liste = [liste ; nlist];

  xcrit = lxcrit(imail);
  nbno1 = length(liste);
  replace_node(nbno1-length(nlist)+1:nbno1) = 0;
  if xcrit > 0
  xcoorlocal = xcoor(liste,:);
% Boucle sur les noeuds de nlist
  for ino = nbno1:-1:nbno1-length(nlist)+1
%   Boucle sur les noeuds plus anciens
    for jno = 1:ino-1
      if max(abs(xcoorlocal(ino,:) - xcoorlocal(jno,:))) < xcrit
	replace_node(ino) = jno;
	break;
      end
    end
  end
  end

  clear xcoorlocal nlist nmail mail;
end


disp('PB avec ElimNodes A VOIR')


disp(['  Nb eliminated nodes: ' num2str(length(find(replace_node)))])
ListUsedNodes = liste(find(replace_node == 0));

global_replace_node = zeros(size(xcoor,1),1);
global_replace_node(liste) = replace_node';
for ino = 1:length(global_replace_node)
  if global_replace_node(ino) == 0
    global_replace_node(ino) = ino;
  end
end

for imail = 1:nb_mail
  mail = pile_mail{imail};
  mail2 = [];
%  clear mail2;
  for izo = 1:size(mail,2)
    email = mail{izo};
%% subtil
    toto = [email.MAIL ; ones(1,size(email.MAIL,2))];
    toto = global_replace_node(toto);
    toto = toto(1:end-1,:);
    email2 = struct('MAIL',toto,'TYPE',email.TYPE);
    mail2{izo} = email2;
  end
  pile_mail2{imail} = mail2;
end
