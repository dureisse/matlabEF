function [val1] = IntegralTheta2(t1,f1,theta1,varargin)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 11 / 08 / 2004
%
% Calcul de l'integrale d'un champ f1(:,:, ... ,t) sur la grille en
% temps t1 par la theta-methode (nombre d'indices quelconque)
%
% Entrees
%  t1(ntmps)		Liste des piquets de temps
%  f1(:,:, ... ,t)	Champ a integrer
%  theta1		Parametre de la methode d'integration en temps
% Entrees optionnelles
%  weight(ntmps)	Liste des ponderations en temps eventuelle
% Sorties
%  val1(:,:, ...)	Valeur de l'integrale en temps

ntmps   = length(t1);
nin1 = nargin-3;
if (nin1 == 0)
  weight = ones(1,ntmps);
elseif (nin1 == 1)
  weight = varargin{1};
else
  nin1
  error('Bad number of arguments')
end

if (1 == 0)
% DD 15/08/2004 UN BUG ICI VOIR test_Theta
% OU EST weight ???
str = repmat(':,',1,ndims(f1)-1);
str1 = str(1:end-1);

t = 1;
% truc = f1(::,t);
  eval(['truc = f1(' str int2str(t) ');'])
  val1 = zeros(size(truc)); clear truc;

for t = 2:length(t1)
  hi = t1(t) - t1(t-1);
%  val1(::) = val1(::) + hi*(theta1*f1(::,t) + ...
%               (1.-theta1)*f1(::,t-1));
  eval(['val1(' str1 ') = val1(' str1 ...
        ') + hi*(theta1*f1(' ...
        str int2str(t) ') + (1.-theta1)*f1(' str int2str(t-1) '));'])
end
end
% VERSION PLUS BRUTALE, MAIS PLUS JUSTE...


nbindic = length(size(f1));
nbelem  = prod(size(f1));
if (ntmps ~= length(t1)) | (ntmps ~= length(weight))
  [ntmps length(t1) length(weight)]
  error('bad number of time steps')
end

% On code dans toto(:,ntmps);
toto = reshape(f1,nbelem/ntmps,ntmps);

% On integre toto dans val2(:,1)
t = 1;
  val2 = zeros(size(toto,1),1);

for t = 2:length(t1)
  hi = t1(t) - t1(t-1);
  val2(:,1) = val2(:,1) + hi*(theta1*toto(:,t)*weight(t) + ...
                              (1.-theta1)*toto(:,t-1)*weight(t-1));
end

% On decode  dans val1
dims = size(f1);
dims = dims(1:end-1);
if (length(dims) == 0)
  error('on ne doit pas passer par la...')
elseif (length(dims) == 1)
  dims = [dims 1];
end
val1 = reshape(val2,dims);
