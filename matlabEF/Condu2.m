function [condu1,varargout] = Condu2(modl1,matr1,mail1,intg1, ...
                                     xcoor1,mode1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 09 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 11 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2003
%   arguments optionnels
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 10 / 2007
%   ajout anisotropie
% DUREISSEIX David  LaMCoS                           le 13 / 06 / 2014
%   ajout BARR (thermique 1D)
%
% Matrices de conductivite elementaires condu1
% Et operateur dual generalise bq1 associe
% Et operateur primal bz1 associe
% Et operateur de comportement k1 associe
%
% Entrees
%   modl1               modele
%   matr1               materiau
%   mail1               maillage
%   intg1               segment d'integration
%   xcoor1(nbno,idim)   coordonnees des noeuds
%   mode1               mode d'analyse
%  et une liste d'options parmi
%   PrimOp      operateur primal b1
%   GenDualOp   operateur dual generalise bsig1
%   ConstiOp    operateur de comportement d1
%
% Sorties
%   condu1       rigidites elementaires
%  et une liste d'operateurs ranges comme demandes dans la liste des
%  options
%
% Exemples
%   condu1 = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [condu1,bq1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                         'GenDualOp')
%   [condu1,bz1,bq1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                            'PrimOp','GenDualOp')
%   [condu1,bq1,bz1,k1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                                'GenDualOp','PrimOp','ConstiOp')
%
% On peut ainsi avoir :
% l'operateur bq1 dont le transpose donne un champ
% de flux generalise a partir d'un champ de flux aux points d'integration
% et l'operateur bz1 qui donne un champ de gradient de temperature aux
% points d'integration a partir d'un champ de temperature aux noeuds
% et l'operateur k1 qui passe d'un champ de gradient de temperature a un
% champ de flux, tous deux aux points d'integration.


% Appelle Perm6 puisque la formulation THERMIQUE est tres similaire
% a la formulation FLUIDE
% formul.   type         parametres
% THERMIQUE CONDUCTIVITE K
% FLUIDE    PERMEABILITE PERM/VISC

% REMARQUE : LE CONTRAIRE ETAIT PLUS JUDICIEUX PUISQUE SEUL
% PERM/VISC INTERVIENT COMME PARAM MATER

disp('Condu2 : ==================================================')
disp('Condu2 : ATTENTION, NE METTRE QUE LES PARAMETRES MAT UTILES')
disp('Condu2 : SOUS RISQUE DE LES ECRASER !!!')
disp('Condu2 : ==================================================')
clear matr2;
nzo2 = length(matr1);
for zo2 = 1:nzo2
  matr2e = matr1{zo2};
  ncom2 = length(matr2e);
  for i = 1:ncom2
    switch matr2e{i}.COMP
      case 'K',
        matr2e{i}.COMP = 'PERM'; [a,b] = size(matr2e{i}.XVAL);
      case 'KXX',
        matr2e{i}.COMP = 'PEXX'; [a,b] = size(matr2e{i}.XVAL);
      case 'KYY',
        matr2e{i}.COMP = 'PEYY'; [a,b] = size(matr2e{i}.XVAL);
      case 'KZZ',
        matr2e{i}.COMP = 'PEZZ'; [a,b] = size(matr2e{i}.XVAL);
      case 'KXY',
        matr2e{i}.COMP = 'PEXY'; [a,b] = size(matr2e{i}.XVAL);
      case 'KYX',
        matr2e{i}.COMP = 'PEYX'; [a,b] = size(matr2e{i}.XVAL);
      case 'KYZ',
        matr2e{i}.COMP = 'PEYZ'; [a,b] = size(matr2e{i}.XVAL);
      case 'KZY',
        matr2e{i}.COMP = 'PEZY'; [a,b] = size(matr2e{i}.XVAL);
      case 'KZX',
        matr2e{i}.COMP = 'PEZX'; [a,b] = size(matr2e{i}.XVAL);
      case 'KXZ',
        matr2e{i}.COMP = 'PEXZ'; [a,b] = size(matr2e{i}.XVAL);
      case 'KRR',
        matr2e{i}.COMP = 'PERR'; [a,b] = size(matr2e{i}.XVAL);
      case 'KRZ',
        matr2e{i}.COMP = 'PERZ'; [a,b] = size(matr2e{i}.XVAL);
      case 'KZR',
        matr2e{i}.COMP = 'PEZR'; [a,b] = size(matr2e{i}.XVAL);
    end
  end
  matr2e{ncom2+1} = struct('XVAL',ones(a,b),'COMP','VISC','UNIT','');
  matr2{zo2} = matr2e;
end
modl2 = modl1;
nzo1 = length(modl1);
for zo1 = 1:nzo1
  modl2{zo1}.DDLP = [{'P'}];
  modl2{zo1}.DDLD = [{'FP'}];
  switch mode1
    case {'DEPL','COPL'},
      modl2{zo1}.COMP = [{'PX'} {'PY'}];
      modl2{zo1}.COMD = [{'WX'} {'WY'}];
    case 'AXIS',
      modl2{zo1}.COMP = [{'PR'} {'PZ'}];
      modl2{zo1}.COMD = [{'WR'} {'WZ'}];
    case 'TRID',
      modl2{zo1}.COMP = [{'PX'} {'PY'} {'PZ'}];
      modl2{zo1}.COMD = [{'WX'} {'WY'} {'WZ'}];
    case 'BARR',
      modl2{zo1}.COMP = [{'PX'}];
      modl2{zo1}.COMD = [{'WX'}];
    otherwise,
      mode1
      error('Mode not implemented')
  end
end

nin1 = nargin-6;
nout = nargout-1;
if (nin1 ~= nout)
  nin1
  nout
  error('bad number of arguments and/or outputs')
else
  listRegularOp = [{'GenDualOp'} {'PrimOp'} {'ConstiOp'}];
  listOp = findoccur(listRegularOp,varargin);
  test = findoccur(varargin,listRegularOp);
  if ~isempty(find(test==0))
    varargin
    listRegularOp
    error('bad option')
  end
end


disp('  Condu2: Perm6 is called due to similarity in the formulation')
if nin1 == 0
  [condu1] = Perm6(modl2,matr2,mail1,intg1,xcoor1,mode1);
elseif nin1 == 1
  [condu1,varargout{1}] = ...
    Perm6(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1});
elseif nin1 == 2
  [condu1,varargout{1},varargout{2}] = ...
    Perm6(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1},varargin{2});
elseif nin1 == 3
  [condu1,varargout{1},varargout{2},varargout{3}] = ...
    Perm6(modl2,matr2,mail1,intg1,xcoor1,mode1, ...
    varargin{1},varargin{2},varargin{3});
end
