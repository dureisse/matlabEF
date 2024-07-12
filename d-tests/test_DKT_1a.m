clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'DKIR';
idim = 3;

% On lit les donnees
file1 = 'test_DKT_1a.inp';
xcrit1 = -1.;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
xcoor1 = xcoor1(:,1:idim);

% Maillage
mail1 = ListMesh1{1};
 
% Champ de deplacement
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};

% Materiau
matr1 = ListChml1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);

% Base locale
[chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);

% Rigidites elementaires
rigi1 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,chamno1);


%aa = [3 4 5 9 10 11 15 16 17];
% que la fleche
aa = [3 9 15];
el0 = 0;
for zo1 = 1 :length(mail1)
  maile1 = mail1{zo1};
  rigie1 = rigi1{zo1};
  topo1 = maile1.MAIL;
  nbel = size(topo1,1);
  for el1 = 1:nbel
    el0 = el0 + 1;
    youn1 = matr1{1}.XVAL(el0,1);
    nu1   = matr1{2}.XVAL(el0,1);
    rho = 1.;
    e1  = matr1{3}.XVAL(el0,1);
    xel1 = xcoor1(topo1(el1,:),1:2);
    VCORE = [xel1(1,:) xel1(2,:) xel1(3,:)];
    VPREE = [youn1 nu1 rho e1];
    KKE = DKT(VCORE,VPREE);
    KKE = KKE([1 4 7],[1 4 7]);
    KE  = rigie1.XVAL(:,:,el1);
    KE = KE(aa,aa);
    err1 = norm(KE - KKE)/norm(KE);
    if (err1 > 1e-8)
        err1
        error('pb')
    end
  end
end

disp('TEST PASSE AVEC SUCCES')
quit