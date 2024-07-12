function [mail2,ListUsedNodes] = ElimNodeMails5(xcoor,xcrit,mail1,varargin);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 26 / 12 / 2003
%
% Elimination des noeuds dans un maillage mail1 qui sont au meme
% emplacement que des noeuds d'un maillage donne mail3 (qui n'est pas
% change) ou de facon interne si mail3 n'est pas fourni.
% On ne touche pas a la pile des noeuds xcoor.
%
% Entrees
%  xcoor(nbno,idim)	coordonnees des noeuds
%  xcrit		critere de proximite de noeuds
%  mail1		maillage a traiter
% Entree optionnelle
%  mail3		maillage de reference
% Sorties
%  mail2		maillage traite
%  ListUsedNodes	liste des noeuds effectivement utilises
%                       dans mail2
%
% Remarque : si xcrit < 0, on n'elimine pas de noeuds, mais on construit
% quand meme ListUsedNodes

if nargin == 3
  mail3 = mail1;
elseif nargin == 4
  mail3 = varargin{1};
else
  nargin
  error('Bad number of arguments')
end

nmail3 = ChangeMesh2(mail3,'POI1');
numer3 = unique(nmail3{1}.MAIL');   %  numeros tries


% liste des noeuds
% list_node(i) = i : le noeud i est utilise et ne sera pas remplace
% list_node(i) = j : le noeud i n'est plus utilise car il sera 
%                    remplace par le noeud j
% list_node(i) = 0 : noeud non utilise,

nbno = size(xcoor,1);
list_node = zeros(1,nbno);

% On regarde les noeuds du maillage en partant par la fin
  nmail1 = ChangeMesh2(mail1,'POI1');
  numer1 = unique(nmail1{1}.MAIL');   %  numeros tries
  for i1 = length(numer1):-1:1
    ino = numer1(i1);
    ino_new = ino;

    if xcrit > 0
%     Boucle sur les noeuds en stock en partant du debut et en sortant
%     avec un break (pour assurer qu'une seule passe suffit)
      xcoorino = xcoor(ino,:);
      for j1 = 1:length(numer3)
        jno = numer3(j1);
        if max(abs(xcoorino - xcoor(jno,:))) < xcrit
          ino_new = jno;
          break;
        end
      end
    end

    list_node(ino) = ino_new;
    clear xcoorino;
  end
  clear numer1 nmail1 numer3 nmail3;


% On repere les noeuds utiles list_node(i) == i
% """""""""""""""""""""""""""
ListUsedNodes = find(list_node == [1:nbno]);
disp(['  Nb eliminated nodes: ' num2str(nbno - length(ListUsedNodes))])


% On procede au remplacement
% """"""""""""""""""""""""""
mail2 = [];
nbzo1 = size(mail1,2);
for izo1 = 1:nbzo1
  maile1 = mail1{izo1};
  topo1 = maile1.MAIL;
% (:,:) necessaire pour traiter correctement un vecteur colonne
  topo1(:,:) = list_node(topo1);
  maile1.MAIL = topo1;
  mail2{izo1} = maile1;
  clear maile1 topo1;
end
