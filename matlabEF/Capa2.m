function [capa1,varargout] = Capa2(modl1,matr1,mail1,intg1, ...
                                   xcoor1,mode1,varargin)
%  Capacity (thermic) elementary matrices and associated operators
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 09 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 11 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 01 / 08 / 2006
%   Utilisation des arguments optionnels
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 10 / 2007
%   Correction d'un bug dans le nom de composantes
%
% Matrices de capacite elementaires capa1
% Et operateur Bsigma associe bsig1
% Et operateur B b1
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
%   capa1       capacites elementaires
%  et une liste d'operateurs ranges comme demandes dans la liste des
%  options
%
% Exemples
%   capa1 = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [capa1,bq1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                        'GenDualOp')
%   [capa1,bz1,bq1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                            'PrimOp','GenDualOp')
%   [capa1,bq1,bz1,c1] = Condu2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                                'GenDualOp','PrimOp','ConstiOp')
%
% Appelle Compress7 puisque la formulation THERMIQUE est tres similaire
% a la formulation FLUIDE
% formul.   type            parametres
% THERMIQUE CAPACITE        RHO*C
% FLUIDE    COMPRESSIBILITE 1/MOB

% On cherche les coefficients
clear lrho lc;
nzo2 = length(matr1);
for zo2 = 1:nzo2
  matr2e = matr1{zo2};
  ncom2 = length(matr2e);
  for i = 1:ncom2
    if strcmp(matr2e{i}.COMP,'C')
      lc(zo2) = i;
    end
    if strcmp(matr2e{i}.COMP,'RHO')
      lrho(zo2) = i;
    end
  end
  clear matr2e;
end

% On construit un pseudo-materiau
clear matr2;
for zo2 = 1:nzo2
  clear matr2e;
  xval1 = matr1{zo2}{lc(zo2)}.XVAL;
  xval1 = xval1 .* matr1{zo2}{lrho(zo2)}.XVAL;
  xval1 = xval1 .^ -1;
  matr2e{1} = struct('COMP','MOB','UNIT','','XVAL',xval1);
  matr2{zo2} = matr2e;
  clear matr2e xval1;
end

% et un pseudo-modele
modl2 = modl1;
nzo1 = length(modl1);
for zo1 = 1:nzo1
  modl2{zo1}.DDLP = [{'P'}];
  modl2{zo1}.DDLD = [{'FP'}];
  modl2{zo1}.COMP = [{'P'}];
  modl2{zo1}.COMD = [{'FP'}];
%  switch mode1
%    case {'DEPL','COPL'},
%      modl2{zo1}.COMP = [{'PX'} {'PY'}];
%      modl2{zo1}.COMD = [{'WX'} {'WY'}];
%    case 'AXIS',
%      modl2{zo1}.COMP = [{'PR'} {'PZ'}];
%      modl2{zo1}.COMD = [{'WR'} {'WZ'}];
%    case 'TRID',
%      modl2{zo1}.COMP = [{'PX'} {'PY'} {'PZ'}];
%      modl2{zo1}.COMD = [{'WX'} {'WY'} {'WZ'}];
%    otherwise,
%      mode1
%      error('Mode not yet implemented')
%  end
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

if nin1 == 0
  [capa1] = Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1);
elseif nin1 == 1
  [capa1,varargout{1}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1});
elseif nin1 == 2
  [capa1,varargout{1},varargout{2}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1},varargin{2});
elseif nin1 == 3
  [capa1,varargout{1},varargout{2},varargout{3}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1, ...
    varargin{1},varargin{2},varargin{3});
end
clear matr2 modl2;
