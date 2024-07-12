function [Fnew] = GramS(F,varargin)
% Gram Schmidt orthogonalization

nin1 = nargin-1;
if (nin1 == 0)
  opti1 = 'chol';
elseif (nin1 == 1)
  opti1 = varargin{1};
else
  nin1
  error('Bad number of optional arguments')
end

M = F'*F;
switch opti1
  case 'chol',
%   For SPD matrices
    B = chol(M); % B'*B = M
    Fnew = B' \ F';
  case 'lu';
    [L,U] = lu(M); % L'*U = M
    error('option not yet understood... sorry')
  otherwise
    opti1
    errorr('Option not implemented')
end
Fnew = Fnew';

