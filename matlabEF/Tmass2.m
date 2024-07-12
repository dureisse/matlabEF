function [mass1,varargout] = Tmass2(modl1,matr1,mail1,intg1, ...
                                   xcoor1,mode1,varargin)
%  `Mass' elementary matrices and associated operators for TRANSPORT
%  see Capa2 and Compress7
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 31 / 07 / 2007
%
% Matrices de masse elementaires mass1
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
%   mass1       masses elementaires
%  et une liste d'operateurs ranges comme demandes dans la liste des
%  options
%
% Exemples
%   mass1 = Tmass2(modl1,matr1,mail1,intg1,xcoor1,mode1)
%   [mass1,bq1] = Tmass2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                        'GenDualOp')
%   [mass1,bz1,bq1] = Tmass2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                            'PrimOp','GenDualOp')
%   [mass1,bq1,bz1,c1] = Tmass2(modl1,matr1,mail1,intg1,xcoor1,mode1, ...
%                                'GenDualOp','PrimOp','ConstiOp')
%
% Appelle Compress7 puisque la formulation TRANSPORT,MASSE est tres similaire
% a la formulation FLUIDE,COMPRESSIBILITE
% formul.   type            parametres
% TRANSPORT MASSE           RHO
% THERMIQUE CAPACITE        RHO*C
% FLUIDE    COMPRESSIBILITE 1/MOB

% On cherche le coefficient
clear lrho;
nzo2 = length(matr1);
for zo2 = 1:nzo2
  matr2e = matr1{zo2};
  ncom2 = length(matr2e);
  for i = 1:ncom2
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
  xval1 = matr1{zo2}{lrho(zo2)}.XVAL;
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
  modl2{zo1}.COMP = [{'PX'} {'PY'}];
  modl2{zo1}.COMD = [{'WX'} {'WY'}];
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
  [mass1] = Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1);
elseif nin1 == 1
  [mass1,varargout{1}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1});
elseif nin1 == 2
  [mass1,varargout{1},varargout{2}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1,varargin{1},varargin{2});
elseif nin1 == 3
  [mass1,varargout{1},varargout{2},varargout{3}] = ...
    Compress7(modl2,matr2,mail1,intg1,xcoor1,mode1, ...
    varargin{1},varargin{2},varargin{3});
end
clear matr2 modl2;
