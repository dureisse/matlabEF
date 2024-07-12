function [n1 t1] = BaseLocaleC2D(NO1,mapddlPrim1,nocon1)
% 2022/10/11   Dureisseix David             LaMCoS
%
% Base locale aux contacts n1,t1 en 2D
%
% Inputs

% Outputs

n1 = NO1(mapddlPrim1);
n1 = n1(nocon1,:);
t1 = 0*n1;
nbnc = length(nocon1); % nb de noeuds de contact potentiel
for i = 1:nbnc
    vectn = n1(i,:)';
    vectn = vectn / norm(vectn);
    n1(i,:) = vectn';
    t1(i,:) = [-vectn(2,1) vectn(1,1)];
end

end