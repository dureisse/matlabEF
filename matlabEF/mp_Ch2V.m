function [vect] = mp_Ch2V(field,map)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 05 / 2003
%
% Transforme un champ(comp,M,t) en vecteur(ddl,t)

[npts,ncomps] = size(map);
[ncomps,npts,nptmps] = size(field);
nddl = max(max(map));

vect = zeros(nddl,nptmps);

% Boucle sur les piquets de temps
for t = 1:nptmps
  toto = zeros(nddl,1);
% DD BUG 06/07/2003  toto(map) = field(:,:,t);
  toto(map) = field(:,:,t)';
  vect(:,t) = toto;
  clear toto;
end
