function [field] = mp_V2Ch(vect,map)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 05 / 2003
%
% Transforme un vecteur(ddl,t) en champ(comp,M,t)

[nddl,nptmps] = size(vect);
[npts,ncomps] = size(map);

field = zeros(ncomps,npts,nptmps);

% Boucle sur les piquets de temps
for t = 1:nptmps
  toto = vect(:,t);
  field(:,:,t) = toto(map');
  clear toto;
end
