function [BZ1,varargout] = Gradient(mail1,xcoor1,intg1,


numerinv1,mapddlPrim1, ...
                         intg1,modl1, ...
                         Option1,U1)

% Calcule les matrices de calcul du gradient de la transformation
%
% DUREISSEIX David  LaMCoS                     le 18 / 06 / 2023
%
% Inputs
%   mail1           Objet maillage
%   xcoor1(nbno,idim) Coordonneees des noeuds
% Output
%   BZ1(idim,idimr,nbptgt)  Gradients
%          pour le ptg de numero global iptg, 
%          dN/d_x0 = BZ1(idim,idimr,iptg)
% Optional output
%   W1(nbptgt)   Poids d'integration


idim = size(xcoor1,2);
if (idim ~= 3)
    error('pas encore valide en 2D')
end

nbptgt = Nptg(mail1,intg1);
nou1 = nargout-1; % number of optional outputs
if (nou1 > 1)
    error('Bad number of output arguments')
end
if (nou1 == 1)
    W1 = zeros(nbptgt,1);
end

% Boucle sur tous les elements
% """"""""""""""""""""""""""""
elt1 = 0; % numero courant d'element
iptg1 = 0; % numero courant de point d'integration
nzo1 = length(mail1);
for zo1 = 1:nzo1
    intge1 = intg1{zo1}; % segment d'integration
    % idimr : Reference dimension (reference space)
    % nbnni : Number of primal basis functions
    % nbptg : Number of integration points
    [idimr nbnni nbptg] = size(intge1.DPHI);
    % checks
    if (idimr ~= 2)
        error('surface seulement')
    end

    % modle1 = modl1{zo1}; % modele
    % % checks
    % if (length(modle1.DDLP) ~= length(modle1.DDLD))
    %   modle1.DDLP
    %   modle1.DDLD
    %   error('Different primal and dual physical components not implemented')
    % end
    % if (length(modle1.NDDP) ~= length(modle1.NDDD))
    %   modle1.NDDP
    %   modle1.NDDD
    %   error('Different primal and dual ddl not implemented')
    % end
    % if ~all(modle1.NNOP == modle1.NNOD)
    %   modle1.NNOP
    %   modle1.NNOD
    %   error('Different primal and dual nodes not implemented')
    % end
    % if ~all(modle1.NNIP == modle1.NNID)
    %   modle1.NNIP
    %   modle1.NNID
    %   error('Different primal and dual basis functions not implemented')
    % end
    % 
    % % nbddl: Number of primal and dual ddl on the element
    % nbddl = length(modle1.NNIP);
    % nbcomp = length(modle1.COMP); % Nombre de composantes de la deformation (symetrique)
    % % checks
    % if (length(modle1.COMP) ~= length(modle1.COMD))
    %   modle1.COMP
    %   modle1.COMD
    %   error('Different primal and dual numbers of components not implemented')
    % end
    % 
    % % list of node numbers
    % nnop = modle1.NNOP;
    % % list of dof names
    % nddp = modle1.NDDP;
    % % list of basis functions (one for each dof)
    % nnip = modle1.NNIP;
    % % list of basis functions for element transformation
    % nnit = modle1.NNIT;

    maile1 = mail1{zo1}; % maillage elementaire
    topo1  = maile1.MAIL; % liste des noeuds par element
    nbel   = size(topo1,1); % nombre d'elements dans la zone
    for el1 = 1:nbel
        elt1 = elt1 + 1; % numero global de l'element
nbptg

    % intge1.DPHI(ir,ni,ptg) = derivee de la fct de forme ni, par rapport a
    % la coordonnee ir, au point d'integration ptg
    % BZ1(1,1,iptg1) =  intge1.DPHI(1,ni1x de chaque noeud,ptg)
    % BZ1(2,1,iptg1) =  intge1.DPHI(1,ni1y,ptg)
    % BZ1(3,1,iptg1) =  intge1.DPHI(1,ni1z,ptg)
    % BZ1(1,2,iptg1) =  intge1.DPHI(2,ni1x,ptg)
    % BZ1(2,2,iptg1) =  intge1.DPHI(2,ni1y,ptg)
    % BZ1(3,2,iptg1) =  intge1.DPHI(2,ni1z,ptg)
ou est nddl_des_noeuds dans BZ1 ?
modle1.NNIP(1:3:end) numeros des fcts de formes pour UX pour tous les noeuds
modle1.NNIP(2:3:end) numeros des fcts de formes pour UY pour tous les noeuds
modle1.NNIP(3:3:end) numeros des fcts de formes pour UZ pour tous les noeuds


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



