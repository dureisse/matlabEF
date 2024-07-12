      function [B] = bmato(SX,SY,HXABS,HYABS,HXORD,HYORD)
% *********************************************************                     
% ce sous programme calcule la matrice [B] qui lie les                          
% courbures aux d.d.l. de l element                                             
% *********************************************************                     
      B = zeros(3,9);
      B11 = zeros(1,9);
      B12 = zeros(1,9);
      B21 = zeros(1,9);
      B22 = zeros(1,9);
      B31 = zeros(1,9);
      B32 = zeros(1,9);
      B33 = zeros(1,9);
      B34 = zeros(1,9);

      for I = 1:9
      B11(I)=SY(2)*HXABS(I);
      B12(I)=SY(3)*HXORD(I);
      B21(I)=-SX(2)*HYABS(I);
      B22(I)=-SX(3)*HYORD(I);
      B31(I)=-SX(2)*HXABS(I);
      B32(I)=-SX(3)*HXORD(I);
      B33(I)=SY(2)*HYABS(I);
      B34(I)=SY(3)*HYORD(I);
      end
%
      AIR2=abs(SX(2)*SY(3)-SX(3)*SY(2));
      for I = 1:9
      B(1,I)=(B11(I)+B12(I))/AIR2;
      B(2,I)=(B21(I)+B22(I))/AIR2;
      B(3,I)=(B31(I)+B32(I)+B33(I)+B34(I))/AIR2;
      end
