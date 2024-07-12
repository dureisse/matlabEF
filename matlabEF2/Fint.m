function [Fi1,varargout] = Fint(Sig1,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                                intg1,modl1, ...
                                Option1,U1)

% Calcule les forces generalisees interieures associees a
% un champ de contrainte 
%
% DUREISSEIX David  LaMCoS                     le 02 / 04 / 2023
% DUREISSEIX David  LaMCoS                     le 04 / 05 / 2023
%   Ajout Material-Hencky-IsotrElas
%
% Inputs
%   Sig1(ns,1)      Vecteur assemble des contraintes
%                   en notation de Voigt/Robin car symetrique
%   mail1           Objet maillage
%   numerinv1(no)   Numerotation inverse des noeuds
%                   le noeud local no1 a pour numero global numerinv1(no1)
%   mapddlPrim1(no,idim) mapping des ddl (no noeuds locaux, idim
%                   composantes)
%   Option1         Option du calcul :
%                   'HPP' en petites perturbations (U1 pas necessaire)
%                   'Green-Lagrange' en Lagrangien total, la 
%                     configuration est donnee par U1
%                   'Material-Hencky-IsotrElas' cas particulier elastique
%                     isotrope...
%   U1(nddl,1)      Vecteur des deplacements nodaux pour la configuration
% Output
%   Fi1(nddl,1)		Vecteur des forces generalisees
% Optional outputs
%   W1(nbptgt,1)	Vecteur des poids d'integration (poids*Jacobien)

switch Option1,
    case 'HPP',
    case 'Green-Lagrange',
    case 'Material-Hencky-IsotrElas',
    otherwise,
        Option1
        error('Bad option')
end

idim = size(xcoor1,2);
if (idim ~= 3)
    idim
    error('Pas prevu autre que 3D')
end
if (size(U1,2) ~= 1)
    error('Vecteur deplacement unique attendu pour la configuration')
end
if (size(Sig1,2) ~= 1)
    error('Vecteur colonne contrainte attendu')
end

switch idim
    case 2,
        % Mapping vecteur contrainte symetrique 3,1 (Voigt/Robin) dans matrice 2x2
        ri2 = 2 ^ -0.5;
        map = [1 3; 3 2]; 
        scale1 = [1 ri2; ri2 1]; 
        nbco = 3;
    case 3,
        % Mapping vecteur contrainte symetrique 6,1 (Voigt/Robin) dans matrice 3x3
        ri2 = 2 ^ -0.5;
        map = [1 6 5; 6 2 4; 5 4 3]; 
        scale1 = [1 ri2 ri2; ri2 1 ri2; ri2 ri2 1]; 
        nbco = 6;
    otherwise,
        idim
        error('Bad idim')
end
nbptgt = size(Sig1,1)/nbco; % nb total de points d'integration

nou1 = nargout-1; % number of optional outputs
if (nou1 > 1)
    error('Bad number of output arguments')
end
if (nou1 == 1)
    W1 = zeros(nbptgt,1);
end

% Boucle sur tous les elements
% """"""""""""""""""""""""""""
Fi1 = 0 * U1;
elt1 = 0; % numero d'element total courant
ptgt1 = 0; % numero de point d'integration total courant
jj = 1; % position courante dans Sig1
nzo1 = length(mail1); % nombre de zones
for zo1 = 1:nzo1
    intge1 = intg1{zo1}; % segment d'integration
    % idimr : Reference dimension (reference space)
    % nbnni : Number of primal and dual basis functions
    % nbptg : Number of integration points
    [idimr nbnni nbptg] = size(intge1.DPHI);
    % checks
    if (idim ~= idimr)
        error('elements massifs seulement')
    end

    modle1 = modl1{zo1}; % modele
    % checks
    if (length(modle1.DDLP) ~= length(modle1.DDLD))
      modle1.DDLP
      modle1.DDLD
      error('Different primal and dual physical components not implemented')
    end
    if (length(modle1.NDDP) ~= length(modle1.NDDD))
      modle1.NDDP
      modle1.NDDD
      error('Different primal and dual ddl not implemented')
    end
    if ~all(modle1.NNOP == modle1.NNOD)
      modle1.NNOP
      modle1.NNOD
      error('Different primal and dual nodes not implemented')
    end
    if ~all(modle1.NNIP == modle1.NNID)
      modle1.NNIP
      modle1.NNID
      error('Different primal and dual basis functions not implemented')
    end

    % nbddl: Number of primal and dual ddl on the element
    nbddl = length(modle1.NNIP);
    nbcomp = length(modle1.COMD); % Nombre de composantes de la contrainte (symetrique)
    % checks
    if (length(modle1.COMP) ~= length(modle1.COMD))
      modle1.COMP
      modle1.COMD
      error('Different primal and dual numbers of components not implemented')
    end
    if (nbcomp ~= nbco)
        error('pb pas bon nb de composantes')
    end

    % list of node numbers
    nnop = modle1.NNOP;
    % list of dof names
    nddp  = modle1.NDDP;
    % list of basis functions
    nnip  = modle1.NNIP;
    % list of basis functions for element transformation
    nnit = modle1.NNIT;

    maile1 = mail1{zo1}; % maillage
    topo1  = maile1.MAIL;
    nbel   = size(topo1,1); % nombre d'elements dans la zone
    nbno   = size(topo1,2); % norm de noeuds par element
    % Check
    if (nbno ~= nbddl/idim)
        error('Bad dimension')
    end
%   disp('Verifier element massif')
    ind1 = findoccur({maile1.TYPE},[{'TET4'} {'CUB8'} {'PRI6'}]);
    if (ind1 == 0)
        maile1.TYPE
        error('Bad element')
    end

    % matrices de mapping des ddl pour les 3 directions
    % on suppose idim = 3, ddl ranges dans l'ordre
    % u(j:idim:end,1) = Aj*u
    A1 = zeros(nbddl/idim,nbddl);
    A2 = zeros(nbddl/idim,nbddl);
    A3 = zeros(nbddl/idim,nbddl);
    for i = 1:nbno
        A1(i,1+(i-1)*idim) = 1;
        A2(i,2+(i-1)*idim) = 1;
        A3(i,3+(i-1)*idim) = 1;
    end

    for el1 = 1:nbel
        elt1 = elt1 + 1; % numero global de l'element

        % On desassemble par element U1 en Ue1
        % """"""""""""""
        lno1 = topo1(el1,:); % numeros globaux des noeuds
        xcoorel = xcoor1(lno1,:); % coordonnees des noeuds
        % coordonnees des noeuds qui participent a la transformation de l'element naturel
        xcoort1 = xcoorel(modle1.NNOT,:);
        % nnop = modle1.NNOP; % list of node numbers
        % nddp  = modle1.NDDP; % list of dof names
        % nnip  = modle1.NNIP; % list of basis functions
        % nnit = modle1.NNIT; % list of basis functions for element transformation
        lnu1 = numerinv1(lno1); % numeros locaux des noeuds
        lddl1 = mapddlPrim1(lnu1,:); % numeros des ddl 
        lddl1 = reshape(lddl1',prod(size(lddl1)),1); % dans l'ordre
        Ue1 = U1(lddl1,:);

        % On desassemble par element Sig1 en Sige1 
        % """"""""""""""
        % on suppose les composantes rangees dans l'ordre
%%        Sige1 = Sig1((elt1-1)*nbcomp*nbptg+1:elt1*nbcomp*nbptg,:);
        Sige1 = Sig1(jj:jj+nbcomp*nbptg-1,:);
        jj = jj+nbcomp*nbptg;

        % Boucle sur les points d'integration
        % """""""""""""""""""""""""""""""""""
        % Fie1 = zeros(nbddl,1);
        Fie1t = zeros(1,nbddl); % transpose du vecteur force interieure de l'element
        for ptg1 = 1:nbptg
            wptg  = intge1.WEIGHT(ptg1); % poids d'integration
            % disp('il faudra multiplier par abs(Jaco)')

            % On desassemble Sige1(nbcomp*nbptg,:) par point d'integration
            % (notation de Voigt, Robin, car symmetrique)
            Sigp1 = Sige1((ptg1-1)*nbcomp+1:ptg1*nbcomp,:);
            % On mappe dans une matrice Sigp1(idim,idim) symetrique
            Sigp1 = Sigp1(map) .* scale1;

%           Check symmetrie de Sigp1...
err1 = max(max(abs(Sigp1 - Sigp1')));
if (err1 > 1e-10)
    error('Sigp1 pas symetrique...')
end

            % Transformation for isoparametric elements
            dphixt1 = intge1.DPHI(:,nnit,ptg1);
            [Mjaco,Jaco] = LocalJaco2(dphixt1,xcoort1);
              clear dphixt1;
            if (ptg1 == 1)
              Jaco0 = Jaco;
            else
              if (Jaco0 * Jaco < 0.)
                Jaco0
                ptg1
                Jaco
                error('Change of sign in Jacobian')
              end
            end
            Ijaco = inv(Mjaco);
            % dphix(i,j) = Nj,i with i in idimr, j in ~nbnni for natural reference element
            dphix = intge1.DPHI(:,nnip,ptg1);
            % dphiX(i,nni) = Nnni,i valeur de la derivee par rapport a la coordonne i dans l'espace reel (de reference) de l'element, de la fonction de forme nni
            dphiX = Ijaco * dphix;
              clear dphix Ijaco Mjaco;

            % Gradients : dV_i/dX_j = G_ij * v
            % avec G_ij = dphiX(j,i:idim:end) * Ai de taille (1,nbddl)
            % dV/dX = G "*" v avec G de taille (idim,idim,nbddl)
            % N.B. CE SERAIT MIEUX DE STOCKER G A L'ENVERS
            G = zeros(idim,idim,nbddl);
            i = 1;
                for j = 1:idim
                    G(i,j,:) = dphiX(j,i:idim:end) * A1;
                end
            i = 2;
                for j = 1:idim
                    G(i,j,:) = dphiX(j,i:idim:end) * A2;
                end
            i = 3;
                for j = 1:idim
                    G(i,j,:) = dphiX(j,i:idim:end) * A3;
                end
            % Gradient du deplacement
            % Z = G "*" Ue1;
            % Z = zeros(idim,idim);
            % for i = 1:idim
            %     for j = 1:idim
            %         toto = zeros(1,nbddl);
            %         toto = squeeze(G(i,j,:));
            %         Z(i,j) = toto'*Ue1;
            %     end
            % end
            Z = tensorprod(G,Ue1,3,1); % gradient de U (idim,idim)

            gg = zeros(1,nbddl); % pb de structure de donnees... c'est nul...
            if strcmp(Option1,'Green-Lagrange')
                % P = Z' "*" G  le produit se fait sur le permier bloc
                % (idim,idim) de G
                P = tensorprod(Z',G,2,1); % ou P = pagetimes(Z',G)
                % PP = pagemtimes(Z',G);
                % [max(max(max(abs(P-PP)))) max(max(max(abs(P))))]
                pp = zeros(1,nbddl);
            end

            switch Option1
                case 'HPP',
                    % Sigp1(i,j)*0.5*(G(i,j,:)+G(j,i,:))
                    % Si Sigp1 est symetrique, c'est pareil que
                    % Sigp1(i,j)*G(i,j,:)
                    for i = 1:idim
                        for j = 1:idim
                            gg(1,:) = G(i,j,:);
                            Fie1t = Fie1t + wptg*abs(Jaco)*Sigp1(i,j)*gg;
                        end
                    end
                case 'Green-Lagrange',
                    % Sigp1(i,j)*0.5*(G(i,j,:)+G(j,i,:)+P(i,j,:)+P(j,i,:))
                    % Si Sigp1 est symetrique, c'est pareil que
                    % Sigp1(i,j)*(G(i,j,:)+P(i,j,:))
                    for i = 1:idim
                        for j = 1:idim
                            gg(1,:) = G(i,j,:);
                            pp(1,:) = P(i,j,:);
                            Fie1t = Fie1t + wptg*abs(Jaco)*Sigp1(i,j)*(gg+pp);
                        end
                    end
                case 'Material-Hencky-IsotrElas',
%%                    error('a finaliser')
                    % S:(C^-1 F^T Fpoint)_sym = tr(S CFT Fpoint)
                    % CFT = C^-1 F^T, C = F^T F, F = 1 + Z
                    % Sigp1(i,j)*CFT(i,k)*G(k,j,:)
% BUG : pas 1                    F = 1 + Z;
                    F = eye(idim) + Z;
                    C = F' * F; CFT = C \ F';
                    CFTG = tensorprod(CFT,G,2,1);
                    for i = 1:idim
                        for j = 1:idim
                            gg(1,:) = CFTG(i,j,:);
                            Fie1t = Fie1t + wptg*abs(Jaco)*Sigp1(i,j)*gg;
                        end
                    end

            end
  
            ptgt1 = ptgt1 + 1;
            if (nou1 == 1)
                W1(ptgt1,1) = wptg*abs(Jaco);
            end

%       Fin de boucle sur les points d'integration
        end

        % On assemble Fie1t en Fi1
        % """""""""""
        Fi1(lddl1,:) = Fi1(lddl1,:) + Fie1t';
    
%   Fin de boucle sur les elements de la zone
    end
    clear topo1 maile1 intge1 modle1;

% Fin de boucle sur les zones
end

if (nou1 == 1)
    varargout{1} = W1;
end