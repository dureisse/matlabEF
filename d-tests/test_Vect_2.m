clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
idim = 2;

fid = fopen('test_Vect_2.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
  [xcoor2,mail2,chpo2,chml2,cara2,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');
xcoor2 = xcoor2(:,1:2);
nbno2 = size(xcoor2,1);
clear nmail2;
  nmail2{1} = struct('MAIL',[1:nbno2]','TYPE','POI1');

% On fusionne les bases de donnees
xcoor3 = [xcoor1 ; xcoor2];
  decal1 = size(xcoor1,1);
  mail2 = ShiftNodeNumber(mail2,decal1);
  nmail2 = ShiftNodeNumber(nmail2,decal1);
clear xcoor1 xcoor2;

clear nmailt1;
nmailt1{1} = nmail1{1};
nmailt1{2} = nmail2{1};
clear chpot1;
chpot1{1} = chpo1{1};
chpot1{2} = chpo2{1};

% On vectorise le Chpo
numer1 = [1:size(xcoor3,1)];
[listComp1,listUnit1] = ListCompChpo2(chpot1);
mapddl1 = MapddlChpo(chpot1,nmailt1,numer1,listComp1);
U1 = ChpoToVect3(chpot1,nmailt1,numer1,mapddl1,listComp1);

[chpop1,nmailp1] = VectToChpo2 ...
  (U1,numer1,mapddl1(:,1:2),listComp1(1:2),listUnit1(1:2)); 
[chpop2,nmailp2] = VectToChpo2 ...
  (U1,numer1,mapddl1(:,3:end),listComp1(3:end),listUnit1(3:end)); 

% Verification
U2 = ChpoToVect3(chpop1,nmailp1, ...
                 numer1,mapddl1(:,1:2),listComp1(1:2));
U3 = ChpoToVect3(chpop2,nmailp2, ...
                 numer1,mapddl1(:,3:end),listComp1(3:end));
U2a = ChpoToVect3(chpo1,nmail1, ...
                 numer1,mapddl1(:,1:2),listComp1(1:2));
U3a = ChpoToVect3(chpo2,nmail2, ...
                 numer1,mapddl1(:,3:end),listComp1(3:end));
err1 = norm(U2 - U2a)/norm(U2);
err2 = norm(U3 - U3a)/norm(U3);
if (err1 > 1.e-8 | err2 > 1.e-8)
  err1
  err2
  error('erreur')
end

% On vectorise le Cham
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE','COPL',idim);
cham1 = ChmlToCham(chml1,mail1,intg1);

% nbptg1 = 0;
% for zo1 = 1:length(cham1)
%   nbptg1 = nbptg1 + size(mail1{zo1}.MAIL,1)*size(intg1{zo1}.WEIGHT,2);
% end

nbptg1 = Nptg(mail1,intg1);
numer2 = [1:nbptg1];
[listComp2,listUnit2] = ListCompCham2(cham1);
mapcomp1 = MapCompCham(cham1,mail1,intg1,numer2,listComp2);
S1 = ChamToVect3(cham1,mail1,intg1,numer2,mapcomp1,listComp2);

% Il reste a faire un VectToCham et a le verifier

disp('TEST PASSE AVEC SUCCES')
quit

