function [Epsi1] = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
                         intg1,modl1, ...
                         Option1,U1)

% Calcule les deformations associees a
% un champ de deplacement 
% (tenseur symetrique, notation de Robin... ou Voigt...
% on met racine(2) sur les termes hors diagonaux)
%
% DUREISSEIX David  LaMCoS                     le 10 / 04 / 2023
%   Element CUB8
% DUREISSEIX David  LaMCoS                     le 30 / 04 / 2023
%   Hencky strain
%
% Inputs
%   mail1           Objet maillage
%   xcoor1(nbno,idim) Coordonneees des oneuds
%   numerinv1(no)   Numerotation inverse des noeuds
%                   le noeud local no1 a pour numero global numerinv1(no1)
%   mapddlPrim1(no,idim) mapping des ddl (no noeuds locaux, idim
%                   composantes)
%   Option1         Option du calcul :
%                   'HPP' en petites perturbations
%                   'Green-Lagrange' en Lagrangien total
%                   'Material-Hencky' version Lagranangienne ou materielle
%   U1(nddl,1)      Vecteur des deplacements nodaux 
% Output
%   Epsi1(Nbcompt,1)  Vecteur des deformations

switch Option1,
    case 'HPP',
    case 'Green-Lagrange',
    case 'Material-Hencky',
    otherwise,
        Option1
        error('Bad option')
end

idim = size(xcoor1,2);
if (idim ~= 3)
    error('pas encore valide en 2D')
end

switch idim
    case 2,
        % Mapping deformation symetrique (2,2) vers vecteur (3,1) (Voigt)
        r2 = 2 ^ 0.5;
        map = [1 3; 3 2]; 
        scale1 = [1 r2; r2 1]; 
        nbco = 3;
    case 3,
        % Mapping vecteur contrainte symetrique 6,1 (Voigt) dans matrice 3x3
        r2 = 2 ^ 0.5;
        map = [1 6 5; 6 2 4; 5 4 3]; 
        scale1 = [1 r2 r2; r2 1 r2; r2 r2 1]; 
        nbco = 6;
    otherwise,
        idim
        error('Bad idim')
end

% Nombre total de composantes
Nbcompt = 0;
nzo1 = length(mail1); % nombre de zones
for zo1 = 1:nzo1
    intge1 = intg1{zo1}; % segment d'integration
    nbptg = size(intge1.DPHI,3); % Nombre de points d'integration
    maile1 = mail1{zo1};
    nbel1 = size(maile1.MAIL,1); % nombre d'elments
%   disp('Verifier element massif')
    ind1 = findoccur({maile1.TYPE},[{'TET4'} {'CUB8'} {'PRI6'}]);
    if (ind1 == 0)
        maile1.TYPE
        error('Bad element')
    end
    modle1 = modl1{zo1};
    nbcomp = length(modle1.COMP); % Nombre de composantes de la deformation (symetrique)
    if (nbcomp ~= nbco)
        error('pb avec nb de composantes')
    end
    Nbcompt = Nbcompt + nbel1*nbptg*nbcomp;
end
Epsi1 = zeros(Nbcompt,1);

% Boucle sur tous les elements
% """"""""""""""""""""""""""""
elt1 = 0;
ii = 1; % pour ranger dans Epsi1
for zo1 = 1:nzo1
    intge1 = intg1{zo1}; % segment d'integration
    % idimr : Reference dimension (reference space)
    % nbnni : Number of primal basis functions
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
    nbcomp = length(modle1.COMP); % Nombre de composantes de la deformation (symetrique)
    % checks
    if (length(modle1.COMP) ~= length(modle1.COMD))
      modle1.COMP
      modle1.COMD
      error('Different primal and dual numbers of components not implemented')
    end

    % list of node numbers
    nnop = modle1.NNOP;
    % list of dof names
    nddp = modle1.NDDP;
    % list of basis functions (one for each dof)
    nnip = modle1.NNIP;
    % list of basis functions for element transformation
    nnit = modle1.NNIT;

    maile1 = mail1{zo1}; % maillage elementaire
    topo1  = maile1.MAIL; % liste des noeuds par element
    nbel   = size(topo1,1); % nombre d'elements dans la zone
    for el1 = 1:nbel
        elt1 = elt1 + 1; % numero global de l'element

    % intge1.DPHI(ir,ni,ptg) = derivee de la fct de forme ni, par rapport a
    % la coordonnee ir, au point d'integration ptg

        % On desassemble par element U1 en Ue1
        % """"""""""""""
        lno1 = topo1(el1,:); % numeros globaux des noeuds
        xcoorel = xcoor1(lno1,:); % coordonnees des noeuds
        % coordonnees des noeuds qui participent a la transformation de l'element naturel
        xcoort1 = xcoorel(modle1.NNOT,:);
        % nnop = modle1.NNOP; % list of node numbers
        % nddp  = modle1.NDDP; % list of dof names
        % nnip  = modle1.NNIP; % list of basis functions (one for each dof)
        % nnit = modle1.NNIT; % list of basis functions for element transformation
        lnu1 = numerinv1(lno1); % numeros locaux des noeuds
        lddl1 = mapddlPrim1(lnu1,:); % numeros des ddl 
        lddl1 = reshape(lddl1',prod(size(lddl1)),1); % dans l'ordre
        Ue1 = U1(lddl1,:); % liste des dof de l'element
        % ranges autrement
        Ue2 = zeros(length(lddl1),idim);
        for i = 1:idim
            Ue2(i:idim:end,i) = Ue1(i:idim:end,:);
        end

        % Boucle sur les points d'integration
        % """""""""""""""""""""""""""""""""""
        for ptg1 = 1:nbptg

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
            % dphix(i,j) = Nj,i with i in idimr, j in nbnnip for natural reference element
            dphix = intge1.DPHI(:,nnip,ptg1);
            % dphiX(i,nnip) = Nnni,i valeur de la derivee par rapport a la coordonne i dans l'espace reel (de reference) de l'element, de la fonction de forme du ddl nni
            dphiX = Ijaco * dphix;
              clear dphix Ijaco Mjaco;
            % si Z_ij = dU_i/dX_j, alors on construit ZT = Z^T
            % ZT(i,j) = dU_j/dX_i
            ZT = dphiX*Ue2;
            switch Option1
                case 'HPP',
                    E1 = 0.5*(ZT' + ZT);
                case 'Green-Lagrange',
                    E1 = 0.5*(ZT' + ZT + ZT*ZT');
                case 'Material-Hencky',
                    % H = lnU = 0.5 lnC
                    % F = 1 + Z = Q D R^T, 
                    % C = F^T F = R D^2 R^T, H = R lnD R^T
                    F = eye(idim)+ZT';
                    [Q,D,R] = svd(F); 
                    E1 = R * diag(log(diag(D))) * R';
            end
            % on assemble pour l'element elt1, le point d'integration ptg1
            E2 = zeros(nbcomp,1); E2(map) = E1 .* scale1;
            Epsi1(ii:ii+nbcomp-1,:) = E2;
            ii = ii + nbcomp;
        end 

    end
    clear topo1 maile1 intge1 modle1;
end



