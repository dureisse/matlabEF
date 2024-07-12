function [f1] = IntegrateTheta(t1,Df1,theta1,xval0)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 12 / 12 / 2002

% Integration d'une fonction Df1 sur la grille en temps t1
% par la theta-methode
% Le resultat est la fonction f1 sur la grille en temps t1
% On a besoin de la valeur initiale f1(t=0) = xval0


f1(1) = xval0;
for it1 = 2:length(t1)
  hi = t1(it1) - t1(it1-1);
  f1(it1) = f1(it1-1) + hi*(theta1*Df1(it1) + (1.-theta1)*Df1(it1-1));
end
