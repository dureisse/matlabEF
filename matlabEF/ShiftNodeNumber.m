function [mail2] = ShiftNodeNumber(mail1,decal1)
 % DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 25 / 07 / 2002


 nbzone1 = size(mail1,2);
 mail2 = [];
% clear mail2;
 for zo1=1:nbzone1
   mail2{zo1} = mail1{zo1};
   topo2 = mail2{zo1}.MAIL + decal1;
   mail2{zo1}.MAIL = topo2;
 end
