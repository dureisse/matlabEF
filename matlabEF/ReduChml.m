function [chml1] = ReduChml(chml2,mail2,mail1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2002

% Reduit un champ constant par element (chml2,mail2)
% sur un maillage mail1, en un champ constant par element
% (chml1,mail1).
% Si le champ n'est pas defini sur tout mail1, il est 
% extrapole a 0.
% (voir IntersectMesh.m)


% First, search the list listelt1 of elements of mail1
% corresponding to the list listelt2 of elements of mail2
% (global numbering, each time)
% REMARK: THIS SHOULD BE COMMON WITH IntersectMesh.m

listelt1 = []; listelt2 = [];

% For global numbering of elements
nbzone1 = length(mail1);
listnbel1 = Nelements2(mail1,[1:nbzone1]);
listnbel1 = [0 cumsum(listnbel1)];
nbzone2 = length(mail2);
listnbel2 = Nelements2(mail2,[1:nbzone2]);
listnbel2 = [0 cumsum(listnbel2)];

% List of element types of mail2
ListType2 = ListTypeMesh(mail2);

% Loop on elements of mail1
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  typ1   = maile1.TYPE;
  topo1  = maile1.MAIL;
  [nbel1,nbnn1] = size(topo1);
  listzo2 = findoccur([{typ1}],ListType2);
  if length(listzo2) ~= 1
    typ1
    ListType2
    error('Twice the same type of element not considered')
  end

  i2 = 1;
    zo2 = listzo2(i2);
    if zo2
      topo2  = mail2{zo2}.MAIL;
      topo2i = topo2(:,[nbnn1:-1:1]); % reverse order is considered
%     Keep as possible the same order of elements of mail1
      [t1,ia1,ib1]    = intersect(topo1,topo2,'rows');
      [t1i,ia1i,ib1i] = intersect(topo1,topo2i,'rows');
%%DD25/02/2003      [listel1,ia,ib] = unique([ia1 ia1i]); % result is sorted
      [listel1,ia,ib] = unique([ia1' ia1i']); % result is sorted
      listelt1 = [listelt1 (listel1 + listnbel1(zo1))];
%%DD25/02/2003      ib = [ib1 ib1i];
      ib = [ib1' ib1i'];
%%DD17/05/2005      listelt2 = [listelt2 (ib(ia) + listnbel1(zo2))];
      listelt2 = [listelt2 (ib(ia) + listnbel2(zo2))];

      clear topo2 topo2i t1 t1i ia1 ib1 ia1i ib1i ia ib listel1;
    end

  clear listzo2 topo1 maile1;
end

clear listnbel1 listnbel2 ListType2;

% Second: fill in chml1
clear chml1;
nbel = Nelements2(mail1);
nbcomp2 = length(chml2);
for ic = 1:nbcomp2
  xval2 = chml2{ic}.XVAL;
  xval1 = zeros(nbel,size(xval2,2));
  xval1(listelt1,:) = xval2(listelt2,:);
  chml1{ic} = struct('COMP',chml2{ic}.COMP,'UNIT',chml2{ic}.UNIT, ...
                     'XVAL',xval1);
  clear xval1 xval2;
end

