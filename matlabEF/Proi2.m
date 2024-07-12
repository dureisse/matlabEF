function [chpo2,nmail] = Proi2(chamno1,intg1,mail1,nmail2,xcoor1,xcrit1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 27 / 01 / 2003
%   Passage par l'assemblage matriciel
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 22 / 03 / 2003
%   Utilisation du segment d'integration intg1
%
% Projette un champ par elements defini aux noeuds (chamno1,mail1)
% en un champ par point (chpo2,nmail2)
% en utilisant les fonctions de forme de l'element geometrique
% pour interpoler.
%
% Entrees
%   chamno1 : champ par elements defini aux noeuds
%   intg1   : son segment d'integration
%   mail1   : son maillage support
%   nmail2  : maillage cible
%   xcoor1(nbno,idim)  : coordonnees des noeuds
%   xcrit1  : critere relatif de proximite (sera multiplie par la taille
%             de l'element courant pour avoir le critere absolu)
% Sorties
%   chpo2  : champ par point (de maillage support nmail2)
%
% Remarque : en les noeuds de nmail2 exterieurs a mail1,
% on extrapole par 0.

nbzone2 = TestMeshType(nmail2,'POI1');
if (nbzone2 ~= 1)
  nbzone2
  error('more than 1 zone for nmail2 not yet implemented')
end

% List of component names
[listComp1,listUnit1] = ListCompCham2(chamno1);
% Assembly
Nptg1 = Nptg(mail1,intg1);
nbcomp1 = length(listComp1);
numerptg = [1:Nptg1];
mapcomp1 = reshape([1:Nptg1*nbcomp1],nbcomp1,Nptg1)';

S1 = ChamToVect3(chamno1,mail1,intg1,numerptg,mapcomp1,listComp1);
S1 = S1(mapcomp1);

% Transformation matrix between vectors  U1 = T1 * S1
numer2 = nmail2{1}.MAIL';
nbno2  = length(numer2);
xcoor2 = xcoor1(numer2,:);
% reperage des noeuds deja trouves
dejatrouv2 = zeros(size(numer2));

T1 = sparse(nbno2,Nptg1);
%
% Loop on elements of mail1
iptg = 0;
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  type1  = maile1.TYPE;
  Xcorel1 = EF_CoorRefNod(type1);
  topo1 = maile1.MAIL;
  [nbel,nbno] = size(topo1);
  for el1 = 1:nbel
    disp(['  element ' int2str(el1) ' / ' int2str(nbel) ...
          ' in zone ' int2str(zo1)])

%   Find not-yet-found candidates (inside the box)
    xcorel1 = xcoor1(topo1(el1,:),:);
    max1 = max(xcorel1);
    min1 = min(xcorel1);
    eps1 = xcrit1 * min(max1 - min1);
    test1 = (xcoor2 > repmat(min1-eps1,nbno2,1)) & ...
            (xcoor2 < repmat(max1+eps1,nbno2,1));
    lpot1 = find(all(test1,2)' & ~dejatrouv2);

%   Find candidates that are really inside the element
    for ipot1 = 1:length(lpot1)
      locno2 = lpot1(ipot1);
      ino2 = numer2(locno2);
      xcor2 = xcoor1(ino2,:);

      [idans,xcor1,shp1,xcor2a] = EF_CoorRef2( ...
                         type1,Xcorel1,xcorel1,xcor2,eps1,xcrit1);
      if idans
        T1(locno2,iptg+1:iptg+nbno) = shp1;
        dejatrouv2(locno2) = 1;
      end
      clear xcor1 shp1 xcor2a xcor2 ino2 locno2;
    end

    iptg = iptg + nbno;
    clear lpot1 test1 xcorel1;
  end
  clear maile1 topo1;
end
n_reste = length(find(dejatrouv2 == 0));
if n_reste
  disp(['Proi: Warning ' int2str(n_reste) ' nodes not found'])
end

% The vector of node field
U1 = T1 * S1;
clear T1 S1;

% Disassembly of the vector
mapddl1 = reshape([1:nbno2*nbcomp1],nbcomp1,nbno2)';
V1 = zeros(nbno2*nbcomp1,1); V1(mapddl1) = U1;
[chpo2,nmail] = VectToChpo2(V1,numer2,mapddl1,listComp1,listUnit1);

clear U1 V1;
