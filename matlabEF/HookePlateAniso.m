function [Dpl] = HookePlateAniso(D,h,mode1)
% Hooke's matrix for a plate anisotropic model, homogeneous in the
% thickness
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 04 / 01 / 2005
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 01 / 06 / 2008
%   Correction de bug, reecriture plus propre des commentaires
%
% Construction de la matrice de Hooke, pour un modele de plaque
% (Kirchhoff en dimension 3 pour l'instant...)
% anisotrope (sans excentration, constant dans l'epaisseur
% ce qui entraine le decouplage)
% On suppose les composantes rangees dans le bon ordre...
%   comp = [{'EP11'} {'EP22'} {'GA12'} {'CP11'} {'CP22'} {'CG12'}];
%   comd = [{'EF11'} {'EF22'} {'EG12'} {'MF11'} {'MF22'} {'MG12'}];
%
% Entrees
%   D(nc,nc)	matrice de Hooke du materiau en base locale a l'element
%               (la 3e coordonnee est la normale a l'element)
%   h		epaisseur
%   mode1	mode d'analyse (DKIR)
% Sorties
%   Dpl		matrice de Hooke du modele de plaque en base locale a l'element
%
% D est donnee en notation de Cowin
% epsilon = [EPS11 EPS22 EPS33 r2*EPS23 r2*EPS31 r2*EPS12]'
% sigma   = [SIG11 SIG22 SIG33 r2*SIG23 r2*SIG31 r2*SIG12]'
% sigma = D epsilon
%
% Dpl est donnee aussi en notation de Cowin pour les plaques
% e   = [EP11 EP22 r2*EP12]'
% khi = [CP11 CP22 r2*CP12]'
% N   = [EF11 EF22 r2*EF12]'
% M   = [MF11 MF22 r2*MF12]'
% ici N et M sont les integrales dans l'epaisseur (pas la moyenne <.>)
% pour que l'energie soit 0.5 * [e ; khi]' * [N ; M]
% N = h<Kcp>e + h<z.Kcp>chi
% M = h<z.Kcp>e + h<z^2.Kcp>chi
% [N;M] = Dpl [e;chi]
% Kcp relie [sigma] et [epsilon] en contraintes planes
% (ou [.] est la partie plane)
% Hcp^-1 = [D^-1] 
%
% Without excentration and with homogeneous Hooke matrix along the
% thickness, the constitutive relation for plate is decoupled
% N = h<Kcp>e       = h Kcp e
% M = h<z^2.Kcp>chi = h^3/12 Kcp chi


nc=size(D,1);
if (nc ~= 6)
  nc
  error('Bad number of components')
end

if strcmp(mode1,'DKIR')
%  Dinv = inv(D);
%  Dinvcp = Dinv([1 2 6],[1 2 6]); % plane part for plane stress
%  Dcp = inv(Dinvcp);
%  Dpl = [h*Dcp      zeros(3,3)
%         zeros(3,3) (h^3/12.)*Dcp];
% autre methode : condensation
  Ddp = D([1 2 6],[1 2 6]); % plane part for plane strain
  D12 = D([1 2 6],[3 4 5]);
  D21 = D([3 4 5],[1 2 6]);
  D22 = D([3 4 5],[3 4 5]);
  Dcp = Ddp - D12 * inv(D22) * D21; % plane part for plane stress
  Dpl = [h*Dcp      zeros(3,3)
         zeros(3,3) (h^3/12.)*Dcp];
else
  mode1
  error('Not yet implemented')
end
