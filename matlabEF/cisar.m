      function [HS4,HS5,HS6,HFX,HFY] =  cisar ...
               (ASS,BSS,AXX,AYY,BXX,BYY,RL,EXX,EYY)

%     DIMENSION ASS(6),BSS(6),AXX(3),AYY(3),BXX(3),BYY(3)                       
%     DIMENSION EXX(3),EYY(3),RL(3)                                             
      HS4 = zeros(1,9);
      HS5 = zeros(1,9);
      HS6 = zeros(1,9);
      COFA = zeros(3,3);
      COFB = zeros(3,3);
      DGX4=zeros(1,9);
      DGX5=zeros(1,9);
      DGX6=zeros(1,9);
      HFX=zeros(1,9);
      HFY=zeros(1,9);
%
      COFB(1,1)=1.D0-1.5D0*EXX(1)*ASS(4)-1.5D0*EYY(1)*BSS(4);
      COFB(2,2)=1.D0-1.5D0*EXX(2)*ASS(5)-1.5D0*EYY(2)*BSS(5);
      COFB(3,3)=1.D0-1.5D0*EXX(3)*ASS(6)-1.5D0*EYY(3)*BSS(6);
      COFB(1,2)=-1.5D0*EXX(1)*ASS(5)-1.5D0*EYY(1)*BSS(5);
      COFB(1,3)=-1.5D0*EXX(1)*ASS(6)-1.5D0*EYY(1)*BSS(6);
      COFB(2,1)=-1.5D0*EXX(2)*ASS(4)-1.5D0*EYY(2)*BSS(4);
      COFB(2,3)=-1.5D0*EXX(2)*ASS(6)-1.5D0*EYY(2)*BSS(6);
      COFB(3,1)=-1.5D0*EXX(3)*ASS(4)-1.5D0*EYY(3)*BSS(4);
      COFB(3,2)=-1.5D0*EXX(3)*ASS(5)-1.5D0*EYY(3)*BSS(5);
%     INVERSION DE LA MATRICE COFB((3,3)                                      
      COFA(1,1)= COFB(2,2)*COFB(3,3)-COFB(2,3)*COFB(3,2);
      COFA(2,2)= COFB(1,1)*COFB(3,3)-COFB(1,3)*COFB(3,1);
      COFA(3,3)= COFB(1,1)*COFB(2,2)-COFB(1,2)*COFB(2,1);
      COFA(1,2)=-COFB(1,2)*COFB(3,3)+COFB(3,2)*COFB(1,3);
      COFA(2,1)=-COFB(2,1)*COFB(3,3)+COFB(2,3)*COFB(3,1);
      COFA(1,3)= COFB(1,2)*COFB(2,3)-COFB(2,2)*COFB(1,3);
      COFA(3,1)= COFB(2,1)*COFB(3,2)-COFB(2,2)*COFB(3,1);
      COFA(2,3)=-COFB(1,1)*COFB(2,3)+COFB(2,1)*COFB(1,3);
      COFA(3,2)=-COFB(1,1)*COFB(3,2)+COFB(1,2)*COFB(3,1);
      DJAC=COFB(1,1)*COFA(1,1)+COFB(2,1)*COFA(1,2)+COFB(3,1)*COFA(1,3);
      XXXX = DJAC;
      if (DJAC ~= 0)
         XXXX=1.D0/DJAC;
      else
        error('Pb in cisar.m: DJAC=0')
      end
      for IA = 1:3
        for IB = 1:3
          COFA(IA,IB)=COFA(IA,IB)*XXXX;
        end
      end
%
      DGX4(1)=0.D0+00;
      DGX4(3)=1.5D0*(EXX(1)*AXX(1)+EYY(1)*BXX(1));
      DGX4(2)=-1.5D0*(EXX(1)*AYY(1)+EYY(1)*BYY(1));
      DGX4(6)=1.5D0*(EXX(1)*AXX(2)+EYY(1)*BXX(2))-(EXX(1)/4.D0);
      DGX4(5)=-1.5D0*(EXX(1)*AYY(2)+EYY(1)*BYY(2))+(EYY(1)/4.D0);
      DGX4(9)=1.5D0*(EXX(1)*AXX(3)+EYY(1)*BXX(3))-(EXX(1)/4.D0);
      DGX4(8)=-1.5D0*(EXX(1)*AYY(3)+EYY(1)*BYY(3))+(EYY(1)/4.D0);
      DGX4(4)=1.5D0/RL(1);
      DGX4(7)=-1.5D0/RL(1);
%                                                                               
      DGX5(1)=-1.5D0/RL(2);
      DGX5(3)=1.5D0*(EXX(2)*AXX(1)+EYY(2)*BXX(1))-EXX(2)/4.D0;
      DGX5(2)=-1.50D+00*(EXX(2)*AYY(1)+EYY(2)*BYY(1))+EYY(2)/4.D0;
      DGX5(6)=1.5D0*(EXX(2)*AXX(2)+EYY(2)*BXX(2));
      DGX5(5)=-1.5D0*(EXX(2)*AYY(2)+EYY(2)*BYY(2));
      DGX5(9)=1.5D0*(EXX(2)*AXX(3)+EYY(2)*BXX(3))-EXX(2)/4.D0;
      DGX5(8)=-1.5D0*(EXX(2)*AYY(3)+EYY(2)*BYY(3))+EYY(2)/4.D0;
      DGX5(4)=0.D0;
      DGX5(7)=1.5D0/RL(2);
%                                                                               
      DGX6(1)=1.5D0/RL(3);
      DGX6(3)=1.5D0*(EXX(3)*AXX(1)+EYY(3)*BXX(1))-EXX(3)/4.D0;
      DGX6(2)=-1.5D0*(EXX(3)*AYY(1)+EYY(3)*BYY(1))+EYY(3)/4.D0;
      DGX6(6)=1.5D0*(EXX(3)*AXX(2)+EYY(3)*BXX(2))-EXX(3)/4.D0;
      DGX6(5)=-1.5D0*(EXX(3)*AYY(2)+EYY(3)*BYY(2))+EYY(3)/4.D0;
      DGX6(9)=1.5D0*(EXX(3)*AXX(3)+EYY(3)*BXX(3));
      DGX6(8)=-1.5D0*(EXX(3)*AYY(3)+EYY(3)*BYY(3));
      DGX6(4)=-1.5D0/RL(3);
      DGX6(7)=0.D0;
      for I = 1:9
      HS4(I)=COFA(1,1)*DGX4(I)+COFA(1,2)*DGX5(I)+COFA(1,3)*DGX6(I);
      HS5(I)=COFA(2,1)*DGX4(I)+COFA(2,2)*DGX5(I)+COFA(2,3)*DGX6(I);
      HS6(I)=COFA(3,1)*DGX4(I)+COFA(3,2)*DGX5(I)+COFA(3,3)*DGX6(I);
      end
      HFX(1)=ASS(4)*HS4(1)+ASS(5)*HS5(1)+ASS(6)*HS6(1);
      HFX(4)=ASS(4)*HS4(4)+ASS(5)*HS5(4)+ASS(6)*HS6(4);
      HFX(7)=ASS(4)*HS4(7)+ASS(5)*HS5(7)+ASS(6)*HS6(7);
      HFX(2)=ASS(4)*HS4(2)+ASS(5)*HS5(2)+ASS(6)*HS6(2)-AYY(1);
      HFX(3)=ASS(4)*HS4(3)+ASS(5)*HS5(3)+ASS(6)*HS6(3)+AXX(1);
      HFX(5)=ASS(4)*HS4(5)+ASS(5)*HS5(5)+ASS(6)*HS6(5)-AYY(2);
      HFX(6)=ASS(4)*HS4(6)+ASS(5)*HS5(6)+ASS(6)*HS6(6)+AXX(2);
      HFX(8)=ASS(4)*HS4(8)+ASS(5)*HS5(8)+ASS(6)*HS6(8)-AYY(3);
      HFX(9)=ASS(4)*HS4(9)+ASS(5)*HS5(9)+ASS(6)*HS6(9)+AXX(3);
%                                                                               
      HFY(1)=BSS(4)*HS4(1)+BSS(5)*HS5(1)+BSS(6)*HS6(1);
      HFY(4)=BSS(4)*HS4(4)+BSS(5)*HS5(4)+BSS(6)*HS6(4);
      HFY(7)=BSS(4)*HS4(7)+BSS(5)*HS5(7)+BSS(6)*HS6(7);
      HFY(2)=BSS(4)*HS4(2)+BSS(5)*HS5(2)+BSS(6)*HS6(2)-BYY(1);
      HFY(3)=BSS(4)*HS4(3)+BSS(5)*HS5(3)+BSS(6)*HS6(3)+BXX(1);
      HFY(5)=BSS(4)*HS4(5)+BSS(5)*HS5(5)+BSS(6)*HS6(5)-BYY(2);
      HFY(6)=BSS(4)*HS4(6)+BSS(5)*HS5(6)+BSS(6)*HS6(6)+BXX(2);
      HFY(8)=BSS(4)*HS4(8)+BSS(5)*HS5(8)+BSS(6)*HS6(8)-BYY(3);
      HFY(9)=BSS(4)*HS4(9)+BSS(5)*HS5(9)+BSS(6)*HS6(9)+BXX(3);
