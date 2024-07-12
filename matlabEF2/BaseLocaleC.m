function [n1 t1 t2] = BaseLocaleC(NO1,mapddlPrim1,nocon1)
% 2022/10/11   Dureisseix David             LaMCoS
%
% Base locale aux contacts n1,t1,t2
%
% Inputs

% Outputs

n1 = NO1(mapddlPrim1);
n1 = n1(nocon1,:);
t1 = 0*n1;
t2 = 0*n1;
nbnc = length(nocon1); % nb de noeuds de contact potentiel
for i = 1:nbnc
    vectn = n1(i,:)';
    vectn = vectn / norm(vectn);
    X = [vectn eye(3)];
    [Q,R] = qr(X);
    if sign(R(1,1)) < 0
        Q(:,1) = -Q(:,1);
    end
    if det(Q) < 0
        Q(:,2) = -Q(:,2);
    end
    % Q contient la base avec comme premier vecteur, vectn
    err1 = norm(vectn - Q(:,1));
    if err1 > 1.e-5
        error('wrong local basis')
    end
    n1(i,:) = vectn';
    t1(i,:) = Q(:,2)';
    t2(i,:) = Q(:,3)';
end

end