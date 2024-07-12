function [chpo1] = CoorMail(nmail1,xcoor1)

% Retourne le chpoint de coordonnees du maillage POI1 nmail1
% (sous-zone par sous-zone)

GlobalVar;

clear chpo1;
nbzone1 = TestMeshType(nmail1,'POI1');
if ~nbzone1
  error('Mesh is not POI1')
end

idim = size(xcoor1,2);

for zo1 = 1:nbzone1
  nmaile1 = nmail1{zo1};
  topo1   = nmaile1.MAIL';
  clear chpoel1;
  for comp1 = 1:idim
    xval1 = xcoor1(topo1,comp1);
    chpoel1{comp1} = struct('COMP',liste_ddlp(comp1),'UNIT','', ...
				'XVAL',xval1);
  end
  chpo1{zo1} = chpoel1;
end
