function [nds1] = NansonFreeNDS(mail1s,modl1s,intg1s, ...
                       X1s,numerinv1s)
% DUREISSEIX David  LaMCoS                     le 23 / 06 / 2023
%
% Calcule le vecteur normal surface elementaire 
% d'une surface 2D dans un espace 3D
%
% Inputs
%   mail1s             Objet maillage de la surface
%   modl1s             Modele de la surface
%   intg1s             Segment d'integration de la surface
%   X1s(nbno1,3)       Vecteur position des noeuds de la surface
%                      dans la configuration consideree
%   numerinv1s(nbno1)  Numerotation inverse des noeuds
%                      le noeud local no1 a pour numero global numerinv1s(no1)
% %   X1s(nddlt1,1)      Vecteur position des noeuds de la surface
% %                      dans la configuration consideree
% %   mapddl1s(nbno1,idim)  Mapping des ddl (no noeuds locaux, idim composantes)
% Output
%   nds1(nbno1,3)      Vecteur normale*surface
%
% nds1(no1,:) contient pour le noeud local no1=1:nbno1
% un vecteur orienté sur la normale (lissée au noeud)
% une norme correspondant au poids d'integration sur la
% surface.
% Le code utilise un produit vectoriel et pas la formule de Nanson.

% %nbnot1 = max(numerinv1s);
% nbno1 = size(mapddl1s,1);
[nbno1,idim] = size(X1s);
nds1 = zeros(nbno1,idim);

% Boucle sur les zones
nbzo1 = length(mail1s);
for izo1 = 1:nbzo1

    NNIP1s = modl1s{izo1}.NNIP; % liste des numeros de fct de forme
    DPHI1s = intg1s{izo1}.DPHI; % valeurs des derivees des fonctions de forme
    [idimr,nbfct1,nbptg1] = size(DPHI1s);
    if ((idim ~= 3) || (idimr ~= 2))
        idim
        idimr
        error('mauvaises dimensions')
    end
    W1s = intg1s{izo1}.WEIGHT; % poids d'integration
    topo1 = mail1s{izo1}.MAIL; % maillage de la zone
    [nbel1,nbnoe1] = size(topo1);
    nddl = idim * nbnoe1; % nb de ddl par element
    l1 = [1:3:nddl]; % liste des ddl UX
    l2 = NNIP1s(l1); % liste des fonctions de forme pour UX
    if (~all(l2 == NNIP1s(l1+1)) || ~all(l2 == NNIP1s(l1+2)))
        error('numerotation pas standard')
    end
    if (length(l1) ~= nbnoe1) 
        error('pb coherence')
    end
    if (nbptg1 ~= nbnoe1)
        error('pas integration aux noeuds')
    end

    % Boucle sur les elements de la zone
    for iel1 = 1:nbel1
        lno1 = topo1(iel1,:); % numeros globaux des noeuds
        lnu1 = numerinv1s(lno1); % numeros locaux des noeuds

        % On desassemble X1s(nddlt1,1) en X1e(nbno1,3)
        % lddl1 = mapddl1s(lnu1,:); % (nbno1,3) numeros des ddl 
        % X1e = X1s(lddl1);
        X1e = X1s(lnu1,:);

        % Boucle sur les points d'integration de l'element
        for iptg1 = 1:nbptg1
            B1s = DPHI1s(:,l2,iptg1);
            XX = X1e' * B1s'; % (3,2) vecteurs tangents transformes
            nds1ptg = W1s(iptg1)*ProdVect(XX);
            % Lissage par sommation aux noeuds nds1(nbno1,3)
            ino1 = iptg1; % on suppose ptg = noeud (integration aux noeuds)
            inu1 = lnu1(ino1); % numero local du noeud
            nds1(inu1,:) = nds1(inu1,:) + nds1ptg';
        end

    end
    
end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% function [FD1] = Pres1ToGenForces(cfext1,nfext1,mail1, ...
%                                  numerd1,mapddlDual1,listDdlDual1, ...
%                                  numerp1,mapddlPrim1,listDdlPrim1, ...
%                                  xcoort1,mode1)
% % DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 08 / 2003
% % DUREISSEIX David  LaMCoS                           le 14 / 06 / 2023
% %
% % calcule un vecteur de forces generalisees FD1
% % issu d'un champ de pression uniforme
% % (entre les deux, il y a la matrice des fonctions de forme croisees)
% %
% % Entrees
% 
% 
% %   cfext1		Champ par point de pression (interpolee)
% %   nfext1		Maillage POI1 qui le sous-tend
% %   mail1		Maillage du bord sur lequel calculer les forces
% %   numerd1		Numerotation des noeuds duaux
% %   mapddlDual1		Matrice d'assemblage des ddl duaux
% %   listDdlDual1	liste des noms de ddl duaux
% %   numerp1		Numerotation des noeuds primaux
% %   mapddlPrim1		Matrice d'assemblage des ddl primaux
% %   listDdlPrim1	liste des noms de ddl primaux
% %   xcoort1(nbno,idim)	Coordonnees des noeuds
% %   mode1		Mode d'analyse
% % Sorties
% %   FD1(nddl,1)		Vecteur de forces generalisees
% %
% 
%   idim = size(xcoort1,2);
% 
%   % ON FAIT CA 1 SEULE FOIS
%   % SI ON FAIT INTGNO MATRICE DIAGONALE
% 
% % "Masses" elementaires sur le bord
%   [modl1s,intg1s] = ModlIntg13(mail1s,'ELASTIQUE','MASSE',mode1,idim);
%   intg1s = SegmentIntgNo(mail1s,'nodes');
%   matr1s = ManuChml(mail1s,'RHO','',1.);
%   matr1s = ChmlToCham(matr1s,mail1s,intg1s);
%   mass1s = Mass7(modl1s,matr1s,mail1s,intg1s,xcoort1,mode1);
%   clear modl1s matr1s;
% 
%   % Avec des matrices diagonales, on ne passe plus en matriciel
%   S1 = zeros(1,nbno1s); % vecteur des poids d'integration aux noeuds
%   n1 = zeros(idim,nbno1s); % normales aux noeuds
% 
%   for izo1s = 1:length(mass1s)
%       m1 = mass1s{izo1s}.XVAL;
%       for iel1s = 1:size(m1,3)
%           m1(:,:,iel1s)
%       end
%   end
