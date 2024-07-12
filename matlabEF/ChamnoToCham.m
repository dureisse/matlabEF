function [cham2] = ChamnoToCham(chamno1,mail1,intgno1,modl1,intg2);
% Transforms an element-field defined at nodes into an element-field
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 07 / 01 / 2005
%
% Transforme un champ par elements defini aux noeuds
% (chamno1,mail1,intgno1,modl1) en un champ par elements
% (cham2,mail1,intg2) par interpolation.
%
% Entrees
%  chamno1	Champ par elements aux noeuds
%  mail1	Son maillage sous tendant
%  intgno1	Son segment d'integration sous tendant
%  modl1	Son modele sous tendant
%  intg2	Nouveau segment d'integration (toujours sur mai1)
%
% Sorties
%  cham2	Champ par elements interpole

% Remarque : on suppose les fonctions de forme rangees dans l'ordre
% des noeuds dans intg2 !


nbzone1 = length(mail1);

clear cham2;

% Loop on zones
% """""""""""""
for zo1 = 1:nbzone1
  maile1   = mail1{zo1};
  chamnoe1 = chamno1{zo1};
  intgnoe1 = intgno1{zo1};
  intge2   = intg2{zo1};
  modle1   = modl1{zo1};
  clear chame2;

  [nbel,nbno] = size(maile1.MAIL);
  if (nbno ~= length(intgnoe1.WEIGHT))
    nbno
    intgnoe1
    error('Integration segment not defined at nodes')
  end
  nbptg = length(intge2.WEIGHT);

% node shape functions (not dof ones...)
  nnop = modle1.NNOP;
  nnod = modle1.NNOD;
  if ~all(nnop==nnod)
    nnop
    nnod
    error('Not same primal and dual nodes')
  end
  nnip = modle1.NNIP;
  nnid = modle1.NNID;
  if ~all(nnip==nnid)
    nnip
    nnid
    error('Not same primal and dual dofs')
  end
  [ListNod,i,j] = unique(nnop); % ListNod=nnop(i); nnop=ListNod(j);
  if length(ListNod) ~= nbno
    error('Unable to find all nodes')
  end
  ListFct = nnip(i);

% Loop on components
% """"""""""""""""""
  nbcom = length(chamnoe1);
  for icom = 1:nbcom
    nom1 = chamnoe1{icom}.COMP;
    uni1 = chamnoe1{icom}.UNIT;
    chame2{icom} = struct('COMP',nom1,'UNIT',uni1);
    nbvalt = size(chamnoe1{icom}.XVAL,2);
    if (nbvalt ~= nbno)
      nbvalt
      nbno
      error('More than 1 component per point not implemented')
    end
    if (icom ~= 1)
      if ~all(frame1==size(chamnoe1{icom}.XVAL))
        error('Complex field structure not implemented')
      end
    end
    frame1 = size(chamnoe1{icom}.XVAL);
  end

% Loop on elements
% """"""""""""""""
%  xval1 = zeros(nbcom,nbel,nbptg);
  xval1 = zeros([nbcom nbel nbptg frame1(3:end)]);
  for el1 = 1:nbel
%   values of _nodal_ shape functions at integration points
    P21 = intge2.PHI(ListFct,:)';
    for icom = 1:nbcom
%     Special treatment for special component names
      switch chamnoe1{icom}.COMP
        case 'Repere',
	  for idim = 1:frame1(4)
	    for ivec = 1:frame1(3)
              xval1(icom,el1,:,ivec,idim) = P21 * ...
	               chamnoe1{icom}.XVAL(el1,:,ivec,idim)';
	    end
	  end
	otherwise,
          xval1(icom,el1,:) = P21 * chamnoe1{icom}.XVAL(el1,:);
      end
    end
  end
  for icom = 1:nbcom
    chame2{icom}.XVAL = xval1(icom,:,:);
  end

  cham2{zo1} = chame2; 
  clear chame2 xval1;
  clear maile1 chamnoe1 intgnoe1 intge2;
end
