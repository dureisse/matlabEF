function [KE] = ElementIntegrOperator(modle1,intge1,xcoorel,mode)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 03 / 2003

% nbno : Number of primal or dual nodes
% idim : Physical dimension (real space)
[nbno idim] = size(xcoorel);

% idimr : Reference dimension (reference space)
% nbnni : Number of primal and dual basis functions
% nbptg : Number of integration points
[idimr nbnni nbptg] = size(intge1.DPHI);

% nbddl : Number of primal and dual ddl for physical field
nbddl = nbptg;

KE  = zeros(nbddl,nbddl);

% Coordinate of nodes that participate to the transformation
% of the element
xcoort1 = xcoorel(modle1.NNOT,:);

% list of basis functions for element transformation
nnit = modle1.NNIT;

% Loop on integration points
for ptg = 1:nbptg
  wptg  = intge1.WEIGHT(ptg);

% Transformation for isoparametric elements
  dphixt1 = intge1.DPHI(:,nnit,ptg);
  [Mjaco,Jaco] = LocalJaco2(dphixt1,xcoort1);
  if (ptg == 1)
    Jaco0 = Jaco;
  else
    if (Jaco0 * Jaco < 0.)
      Jaco0
      ptg
      Jaco
      error('Change of sign in Jacobian')
    end
  end

% Elemental integration operator
  KE(ptg,ptg) = wptg * abs(Jaco);
end 
