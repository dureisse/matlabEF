function [modl2] = TireModlIntg8(modl1,intg1,motDual1,motPrimal1,mode1);
% Tire un sous-modele modl2 d'un modele modl1
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 14 / 11 / 2002
%   ajout de la distinction primal,dual
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 11 / 2002
%   modification des modeles version 1.3 passe a version 1.4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   bug correction on sorting
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajour mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 10 / 11 / 2007
%   Ajour mode TRID
%
% Tire un sous-modele modl2 correspondant a la formulation
% motDual1 (ELASTIQUE,FLUIDE) pour les quantites duales,
% motPrimal1 (ELASTIQUE,FLUIDE) pour les quantites primales,
% pour le mode d'analyse mode1 (COPL,DEPL,TRID),
% d'un modele plus general modl1.
% modl2 est associe au meme segment d'integration que modl1 
% l'etait.

GlobalVar;

nbzone1 = length(modl1);

clear modl2;
% Loop on zones
for zo1 = 1:nbzone1
  modle1 = modl1{zo1};

  if strcmp(motDual1,liste_model{1})
%   ELASTIQUE
    if (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}))
%     COPL Plane stress | DEPL Plane strain
      ddld = liste_ddld2;       % noms des ddl duaux
      comd = liste_stress2;     % noms des composantes duales
    elseif strcmp(mode1,'AXIS')
%     AXIS Axisymmetric
      ddld = liste_ddld2a;       % noms des ddl duaux
      comd = liste_stress2a;     % noms des composantes duales
    elseif strcmp(mode1,'TRID')
%     TRID Tridimensional
      ddld = liste_ddld;       % noms des ddl duaux
      comd = liste_stress;     % noms des composantes duales
    else
      motDual1
      mode1
      error('mode not implemented')
    end
  elseif strcmp(motDual1,liste_model{3})
%   FLUIDE
    if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || ...
        strcmp(mode1,'AXIS') || strcmp(mode1,'TRID'))
%     COPL Plane stress | DEPL Plane strain | AXIS Axisymmetric
%     TRID Tridimensional
      clear ddld; ddld{1} = 'FP'; % Names of dual dof
      comd = ddld;                % Names of dual components
    else
      motDual1
      mode1
      error('mode not implemented')
    end
  else
    motDual1
    error('model not implemented')
  end
  if strcmp(motPrimal1,liste_model{1})
%   ELASTIQUE
    if (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}))
%     COPL Plane stress | DEPL Plane strain
      ddlp = liste_ddlp2;       % noms des ddl primaux
      comp = liste_strain2;     % noms des composantes primales
    elseif strcmp(mode1,'AXIS')
%     AXIS Axisymmetric
      ddlp = liste_ddlp2a;       % noms des ddl primaux
      comp = liste_strain2a;     % noms des composantes primales
    elseif strcmp(mode1,'TRID')
%     TRID Tridimensional
      ddlp = liste_ddlp;       % noms des ddl primaux
      comp = liste_strain;     % noms des composantes primales
    else
      mot1
      mode1
      error('mode not implemented')
    end
  elseif strcmp(motPrimal1,liste_model{3})
%   FLUIDE
    if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || ...
        strcmp(mode1,'AXIS') || strcmp(mode1,'TRID'))
%     COPL Plane stress | DEPL Plane strain | AXIS Axisymmetric |
%     TRID Tridimensional
      clear ddlp; ddlp{1} = 'P';  % Names of primal dof
      comp = ddlp;                % Names of primal components
    else
      mot1
      mode1
      error('mode not implemented')
    end
  else
    mot1
    error('model not implemented')
  end

  modle2.NNOT = modle1.NNOT;
  modle2.NNIT = modle1.NNIT;
  modle2.DDLP = modle1.DDLP;
  modle2.DDLD = modle1.DDLD;
  modle2.COMP = modle1.COMP;
  modle2.COMD = modle1.COMD;

% Selection of DDLP
  [junk,ia,ib] = intersect(ddlp,modle1.DDLP);
% junk=ddlp(ia)=modle1.DDLP(ib)
%% DD 25 / 11 / 03
% As the result is sorted on junk, sort it on ia
  [junk,ii] = sort(ia); ib = ib(ii);
  SelectedDDLP = ib;

% Primal dof which have their names in modle1.DDLP(SelectedDDLP)
  SelectedNDDP = find(ismember(modle1.NDDP,SelectedDDLP));
  modle2.NDDP = modle1.NDDP(SelectedNDDP);
  modle2.NNOP = modle1.NNOP(SelectedNDDP);
  modle2.NNIP = modle1.NNIP(SelectedNDDP);

% Selection of COMP
  [junk,ia,ib] = intersect(comp,modle1.COMP(modle1.NCOP));
% junk=comp(ia)=modle1.COMP(modle1.NCOP(ib))
%% DD 25 / 11 / 03
% As the result is sorted on junk, sort it on ia
  [junk,ii] = sort(ia); ib = ib(ii);
  SelectedCOMP = ib;
  modle2.NCOP = modle1.NCOP(SelectedCOMP);

% Selection of DDLD
  [junk,ia,ib] = intersect(ddld,modle1.DDLD);
% junk=ddld(ia)=modle1.DDLD(ib)
%% DD 25 / 11 / 03
% As the result is sorted on junk, sort it on ia
  [junk,ii] = sort(ia); ib = ib(ii);
  SelectedDDLD = ib;

% Dual dof which have their names in modle1.DDLD(SelectedDDLD)
  SelectedNDDD = find(ismember(modle1.NDDD,SelectedDDLD));
  modle2.NDDD = modle1.NDDD(SelectedNDDD);
  modle2.NNOD = modle1.NNOD(SelectedNDDD);
  modle2.NNID = modle1.NNID(SelectedNDDD);

% Selection of COMD
  [junk,ia,ib] = intersect(comd,modle1.COMD(modle1.NCOD));
% junk=comd(ia)=modle1.COMD(modle1.NCOP(ib))
%% DD 25 / 11 / 03
% As the result is sorted on junk, sort it on ia
  [junk,ii] = sort(ia); ib = ib(ii);
  SelectedCOMD = ib;
  modle2.NCOD = modle1.NCOD(SelectedCOMD);

  modl2{zo1} = modle2;
end
