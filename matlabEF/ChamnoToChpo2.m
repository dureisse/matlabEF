function [chpo1,nmail3] = ChamnoToChpo2(chamno2,mail2,intg2,nmail1, ...
                                        xcoor1,mode1,varargin);
% Smooth an element-based field into a nodal-based field
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 13 / 04 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 13 / 04 / 2004
%   Correction de quelques bugs
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 04 / 2006
%   Versions geometrique ou EF : on introduit mode1 comme argument
%
% Transforme un champ par elemt defini au noeuds (chamno2,mail2,intg2)
% en un champ par point (chpo1,nmail1)
% Operation de lissage - ponderation par la mesure des elements
%
% Inputs
%   chamno2	Values of the element-based field defined at nodes
%   mail2	Mesh of the field
%   intg2	Integration informations (at nodes for geometric option)
%   nmail1	Cloud of nodes where the ouput field will be defined
%   xcoor1(nbno,idim)	Coordinates of the nodes
%   mode1	Mode of analysis
% Optional inputs
%   optio1	Optional string (default = 'FE')
%   		'FE' for (massive) finite element oriented method
%		'geometric' for a geometric approach
% Outputs
%   chpo1	Values of the node-based smoothed field
%   nmail3	Cloud of nodes where it is effectively defined

% UTILISER ChamnoToChpoOperator ! Pas avec l'option geometric

narg = nargin-6;
if (narg == 0)
  optio1 = 'FE';
elseif (narg == 1)
  optio1 = varargin{1};
else
  narg
  error('Bar number of arguments')
end

[nbnot,idim] = size(xcoor1);

switch optio1
  case 'FE',

% ======================================================================
% FE oriented approach
% ======================================================================

% On fait un assemblage de chamno2 en U2
% """"""""""""""""""""""""""""""""""""""
[listComp2,listUnit2] = ListCompCham2(chamno2);
Ncomp2 = length(listComp2);
Nptg2 = Nptg(mail2,intg2);
numer2 = [1:Nptg2*Ncomp2];
mapComp2 = MapCompCham(chamno2,mail2,intg2,numer2,listComp2);
U2 = ChamToVect3(chamno2,mail2,intg2,numer2,mapComp2,listComp2);
U2a = U2(mapComp2);


if (1 == 0)

% On construit l'operateur discret de sommation aux noeuds S
% """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
numer1 = nmail1{1}.MAIL'; % on suppose 1 seule zone
Nno1 = length(numer1);
S = sparse(Nno1,Nptg2);

numer_inv2 = InverseList(numer2,max(numer2));          
numer_inv1 = InverseList(numer1,max(numer1));          

ind2 = 0;
nbzone2 = length(mail2);
for zo2 = 1:nbzone2
  maile2 = mail2{zo2};
  topo2 = maile2.MAIL;
  [nbel2,nbno2] = size(topo2);
  list_ptg2 = [ind2+1:ind2+nbel2*nbno2];
  list_ptg2 = numer_inv2(list_ptg2);
  list_no1 = reshape(topo2',1,nbel2*nbno2);
  list_no1 = numer_inv1(list_no1);
  for i = 1:nbel2*nbno2
    S(list_no1(i),list_ptg2(i)) = 1.;
  end

  ind2 = ind2 + nbel2*nbno2;
  clear list_no1 list_ptg2 topo2 maile2;
end

% Remarque : S' est l'operateur de duplication aux ptgs
% """"""""""

end

% On construit l'operateur de duplication aux ptg
% """""""""""""""""""""""""""""""""""""""""""""""
numer1 = nmail1{1}.MAIL'; % on suppose 1 seule zone
T1 = ChpoToChamnoOperator(numer1,mail2,numer2);

% Remarque : l'operateur discret de sommation aux noeuds est T1'
% """"""""


% On construit le chamno discret de valeur aire de l'element A2
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A2 = zeros(Nptg2,1);

%% mode1 = 'COPL';
[modl2,intg2a] = ModlIntg13(mail2,'ELASTIQUE','RIGIDITE',mode1,idim,intg2);
mode1 = '';
[r1,modlr1] = IntegrOperator2(modl2,mail2,intg2a,xcoor1,mode1);

numer_inv2 = InverseList(numer2,max(numer2));          
ind2 = 0;
nbzone2 = length(mail2);
for zo2 = 1:nbzone2
  maile2 = mail2{zo2};
  topo2 = maile2.MAIL;
  [nbel2,nbno2] = size(topo2);

  re1 = r1{zo2};
  xval1 = re1.XVAL;

  for el2 = 1:nbel2
    aire2 = sum(sum(xval1(:,:,el2)));
    list_ptg2 = [ind2+1:ind2+nbno2];
    list_ptg2 = numer_inv2(list_ptg2);
    A2(list_ptg2,1) = aire2;

    ind2 = ind2 + nbno2;
    clear list_ptg2;
  end

  clear xval1 re1 topo2 maile2;
end
clear r1 modlr1;

% On construit la matrice de ponderation D
% """"""""""""""""""""""""""""""""""""""""
%D = diag(((S'*(S*A2)) .^ (-1)) .* A2);
D = diag(((T1*(T1'*A2)) .^ (-1)) .* A2);

% verif
%T1 = ChpoToChamnoOperator(numer1,mail2,numer2);
%DD = diag(((T1*S*A2) .^ (-1.)) .* A2);
%norm(norm(DD-D))/norm(norm(DD))


% On a l'operateur de lissage LI
% """"""""""""""""""""""""""""""
%LI = S * D;
LI = T1' * D;


% On l'applique
% """""""""""""
U1a = LI * U2a;

% On fait un pseudo desassemblage de U1 en chpo1
% """"""""""""""""""""""""""""""""""""""""""""""
ListComp1 = listComp2;
ListUnit1 = listUnit2;
chpo1 = ManuChpoList(nmail1,ListComp1,ListUnit1,repmat(0.,1,Ncomp2));
mapddl1 = MapddlChpo(chpo1,nmail1,numer1,ListComp1);
U1 = zeros(max(max(mapddl1)),1);
U1(mapddl1) = U1a;
[chpo1,nmail3] = VectToChpo2(U1,numer1,mapddl1,ListComp1,ListUnit1);



  case 'geometric'

% ======================================================================
% geometric approach
% ======================================================================
  intgno2 = intg2;
  [modl2,junk] = ModlIntg13(mail2,'ELASTIQUE','RIGIDITE',mode1,idim,intgno2);

% Boucle sur les zones
% """"""""""""""""""""
  clear chpo1;
  nbzo2 = length(mail2);
  for izo2 = 1:nbzo2
    maile2 = mail2{izo2};
    topo2  = maile2.MAIL;
    [nbel2,nbno2] = size(topo2);
    chamnoe2 = chamno2{izo2};
    intgnoe2 = intgno2{izo2};
    modle2   = modl2{izo2};

    switch length(nmail1)
      case 1,
        numer1 = nmail1{1}.MAIL';
      case nbzo2,
        numer1 = nmail1{izo2}.MAIL';
      otherwise,
        length(nmail1)
	nbzo2
	error('Bad number of zones for cloud of nodes')
    end

%   Boucle sur les composantes
%   """"""""""""""""""""""""""
    clear chpoe1 nmaile3;
    nbcomp2 = length(chamnoe2);
    for icomp2 = 1:nbcomp2
      xval2 = chamnoe2{icomp2}.XVAL;
      nbvalt = size(xval2,2); % Nombre total de valeurs dans la composante

      work1 = sparse(nbnot,0);
      work2 = sparse(nbnot,1);

%     Boucle sur les elements
%     """""""""""""""""""""""
      for iel2 = 1:nbel2
	lno2 = topo2(iel2,:);

%       Ponderation sur cet element
	xcoorel = xcoor1(lno2,:);
	KE = ElementIntegrOperator(modle2,intgnoe2,xcoorel,'UNKNOWN');
%        mesu2 = abs(sum(intgnoe2.WEIGHT)); pas bon : element reel pas de ref !
	mesu2 = sum(diag(KE));
	clear KE xcoorel;

%       On somme les ponderations aux noeuds (poids total)
        work2(lno2) = work2(lno2) + mesu2;

%       Assemblage des valeurs aux noeuds : boucle sur les noeuds
%       """""""""""""""""""""""""""""""""""""""""""""""""""""""""
        for ino2 = 1:nbno2
          no2 = topo2(iel2,ino2);
	  lind2 = [ino2:nbno2:nbvalt];

%         On somme les valeurs ponderees aux noeuds
	  if size(work1,2) < length(lind2)
	    work1(no2,1:length(lind2)) = 0.;
	  end
          work1(no2,1:length(lind2)) = work1(no2,1:length(lind2)) + ...
	                               mesu2*xval2(iel2,lind2);

	  clear lind2;
	end

      end

      norm2 = max(abs(work2));
      if norm2 < 1.e-8
	norm2
        error('tous les elements sont trop petits')
      end
      lno1 = find(abs(work2)/norm2 > 1.e-6);
      work1(lno1) = work1(lno1) ./ work2(lno1); % valeurs moyennees aux noeuds lno1

      xval1 = work1(numer1,:);
      comp1 = chamnoe2{icomp2}.COMP;
      unit1 = chamnoe2{icomp2}.UNIT;
      chpoe1{icomp2} = struct('COMP',comp1,'UNIT',unit1,'XVAL',xval1);
      clear comp1 unit1 xval1;
      clear work1 work2;
    end

    chpo1{izo2}  = chpoe1;
    clear chpoe1;
  end
  nmail3 = nmail1;

  otherwise
    optio1
    error('bad option')
  end
