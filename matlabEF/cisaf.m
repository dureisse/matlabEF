      function [ASS,AXX,AYY,BSS,BXX,BYY] = cisaf(XE,EXX,EYY,DDHOMU,AIR)
%                                                                               
%******************* FONCTIONS DE CISAILLEMENT *****************                
%                                                                               
%          CE SOUS PROGRAMME CALCULE LES CONSTANTES as'(k),ax'(i)               
%          ay'(i),bs'(k),bx'(i),by'(i) ,ou i=1--->3     k=4---->6               
%                                                                               
%                        Y                                                      
%                        .     3                                                
%                        .     *                                                
%                        .5  *   *   4    (4,5,6,noueuds condenses)             
%                        . +        +                                           
%                        .*          *                                          
%                      1 + * * +  * *  + 2 . . . . . . .X                       
%                              6                                                
%*****************************************************************              
%                                                                               
      ASS = zeros(6,1);
      BSS = zeros(6,1);
      AXX = zeros(3,1);
      AYY = zeros(3,1);
      BYY = zeros(3,1);
      BXX = zeros(3,1);
      VV2 = zeros(3,3);
      DFON = zeros(6,3);
      VV3 = zeros(3,3);
      AAX = zeros(1,6);
      AAY = zeros(1,6);
      BBX = zeros(1,6);
      BBY = zeros(1,6);
%
      DFON = zeros(6,3);
      DFON(1,1)=4.D0;
      DFON(1,2)=4.D0;
      DFON(1,3)=4.D0;
      DFON(2,1)=4.D0;
      DFON(3,2)=4.D0;
      DFON(4,3)=4.D0;
      DFON(5,2)=-8.D0;
      DFON(6,3)=-4.D0;
      DFON(5,3)=-4.D0;
      DFON(6,1)=-8.D0;
      DPSIX=(XE(2,3)-XE(2,1))/(2.D0*AIR);
      DPSIY=-(XE(1,3)-XE(1,1))/(2.D0*AIR);
      DETAX=(XE(2,1)-XE(2,2))/(2.D0*AIR);
      DETAY=-(XE(1,1)-XE(1,2))/(2.D0*AIR);
%
      VV3 = zeros(3,3);
%     calcul de la matrice jacobienne du second ordre                           
      VV3(1,1)=DPSIX*DPSIX;
      VV3(1,2)=DETAX*DETAX;
      VV3(2,1)=DPSIY*DPSIY;
      VV3(2,2)=DETAY*DETAY;
      VV3(1,3)=2.D0*DPSIX*DETAX;
      VV3(2,3)=2.D0*DPSIY*DETAY;
      VV3(3,2)=DETAX*DETAY;
      VV3(3,1)=DPSIY*DPSIX;
      VV3(3,3)=DPSIY*DETAX+DETAY*DPSIX;
%
      RAF=DDHOMU(4,4);
      RAC1=DDHOMU(7,7);
      RAC2=DDHOMU(8,8);
      RAFC1= RAF/RAC1;
      RAFC2= RAF/RAC2;
      PXY=DDHOMU(4,5)/DDHOMU(4,4);
      RAPU=DDHOMU(6,6)/DDHOMU(4,4);
      RAFF=DDHOMU(5,5)/DDHOMU(4,4);
      for I = 1:6
      AAX(I)=VV3(1,1)*DFON(I,1)+VV3(1,2)*DFON(I,2)+VV3(1,3)*DFON(I,3) ...
      +RAPU*(VV3(2,1)*DFON(I,1)+VV3(2,2)*DFON(I,2)+VV3(2,3)*DFON(I,3));
      AAY(I)=(VV3(3,1)*DFON(I,1)+VV3(3,2)*DFON(I,2)+VV3(3,3)*DFON(I,3)) ...
      *PXY ...
      +RAPU*(VV3(3,1)*DFON(I,1)+VV3(3,2)*DFON(I,2)+VV3(3,3)*DFON(I,3));
      BBX(I)=(VV3(3,1)*DFON(I,1)+VV3(3,2)*DFON(I,2)+VV3(3,3)*DFON(I,3)) ...
      *PXY ...
      +RAPU*(VV3(3,1)*DFON(I,1)+VV3(3,2)*DFON(I,2)+VV3(3,3)*DFON(I,3));
      BBY(I)=(VV3(2,1)*DFON(I,1)+VV3(2,2)*DFON(I,2)+VV3(2,3)*DFON(I,3)) ...
      *RAFF ...
      +RAPU*(VV3(1,1)*DFON(I,1)+VV3(1,2)*DFON(I,2)+VV3(1,3)*DFON(I,3));
       AAX(I)=RAFC1*AAX(I);
       AAY(I)=RAFC1*AAY(I);
       BBX(I)=RAFC2*BBX(I);
       BBY(I)=RAFC2*BBY(I);
       end
%
       ASS(4)=AAX(4)*EXX(1)+AAY(4)*EYY(1);
       ASS(5)=AAX(5)*EXX(2)+AAY(5)*EYY(2);
       ASS(6)=AAX(6)*EXX(3)+AAY(6)*EYY(3);
       BSS(4)=BBX(4)*EXX(1)+BBY(4)*EYY(1);
       BSS(5)=BBX(5)*EXX(2)+BBY(5)*EYY(2);
       BSS(6)=BBX(6)*EXX(3)+BBY(6)*EYY(3);
       AXX(1)=AAX(1)+0.5D0*EYY(2)*(AAX(5)*EYY(2)-AAY(5)*EXX(2)) ...
      +0.5D0*EYY(3)*(AAX(6)*EYY(3)-AAY(6)*EXX(3));
       AXX(2)=AAX(2)+0.5D0*EYY(1)*(AAX(4)*EYY(1)-AAY(4)*EXX(1)) ...
      +0.5D0*EYY(3)*(AAX(6)*EYY(3)-AAY(6)*EXX(3));
       AXX(3)=AAX(3)+0.5D0*EYY(1)*(AAX(4)*EYY(1)-AAY(4)*EXX(1)) ...
      +0.5D0*EYY(2)*(AAX(5)*EYY(2)-AAY(5)*EXX(2));
       AYY(1)=AAY(1)-0.5D0*EXX(2)*(AAX(5)*EYY(2)-AAY(5)*EXX(2)) ...
      -0.5D0*EXX(3)*(AAX(6)*EYY(3)-AAY(6)*EXX(3));
       AYY(2)=AAY(2)-0.5D0*EXX(1)*(AAX(4)*EYY(1)-AAY(4)*EXX(1)) ...
      -0.5D0*EXX(3)*(AAX(6)*EYY(3)-AAY(6)*EXX(3));
       AYY(3)=AAY(3)-0.5D0*EXX(1)*(AAX(4)*EYY(1)-AAY(4)*EXX(1)) ...
      -0.5D0*EXX(2)*(AAX(5)*EYY(2)-AAY(5)*EXX(2));
       BXX(1)=BBX(1)+0.5D0*EYY(2)*(BBX(5)*EYY(2)-BBY(5)*EXX(2)) ...
      +0.5D0*EYY(3)*(BBX(6)*EYY(3)-BBY(6)*EXX(3));
       BXX(2)=BBX(2)+0.5D0*EYY(1)*(BBX(4)*EYY(1)-BBY(4)*EXX(1)) ...
      +0.5D0*EYY(3)*(BBX(6)*EYY(3)-BBY(6)*EXX(3));
       BXX(3)=BBX(3)+0.5D0*EYY(1)*(BBX(4)*EYY(1)-BBY(4)*EXX(1)) ...
      +0.5D0*EYY(2)*(BBX(5)*EYY(2)-BBY(5)*EXX(2));
       BYY(1)=BBY(1)-0.5D0*EXX(2)*(BBX(5)*EYY(2)-BBY(5)*EXX(2)) ...
      -0.5D0*EXX(3)*(BBX(6)*EYY(3)-BBY(6)*EXX(3));
       BYY(2)=BBY(2)-0.5D0*EXX(1)*(BBX(4)*EYY(1)-BBY(4)*EXX(1)) ...
      -0.5D00*EXX(3)*(BBX(6)*EYY(3)-BBY(6)*EXX(3));
       BYY(3)=BBY(3)-0.5D0*EXX(1)*(BBX(4)*EYY(1)-BBY(4)*EXX(1)) ...
      -0.5D0*EXX(2)*(BBX(5)*EYY(2)-BBY(5)*EXX(2));

