      function [BGENE,DJAC] = bmfdst(IGAU,XE,QSI,ETA,SHPTOT,HS4,HS5,HS6,BGENE)
%
%            CALCUL LA MATRICE B  RELATIVE AUX EFFETS DE                        
%            MEMBRANE ET DE FLEXION                                             
%
      SX = zeros(1,3);
      SY = zeros(1,3);
      RL = zeros(1,3);
      EXX = zeros(1,3);
      EYY = zeros(1,3);
      HXABS = zeros(1,9);
      HXORD = zeros(1,9);
      HYABS = zeros(1,9);
      HYORD = zeros(1,9);
      B = zeros(3,9);
%                                                                               
%      MATRICE B RELATIVE A L'EFFET DE MEMBRANE                                 
%                                                                               
      SHP = zeros(6,3);
      for NPOI = 1:3
      SHP(1,NPOI)=SHPTOT(1,NPOI,IGAU);
      SHP(2,NPOI)=SHPTOT(2,NPOI,IGAU);
      SHP(3,NPOI)=SHPTOT(3,NPOI,IGAU);
      end
      [SHP,DJAC] = jacobidst(XE,SHP,2,3);
      K=1;
      for NPOI = 1:3
      BGENE(1,K  )=SHP(2,NPOI);
      BGENE(1,K+1)=0.D0;
      BGENE(2,K  )=0.D0;
      BGENE(2,K+1)=SHP(3,NPOI);
      BGENE(3,K  )=SHP(3,NPOI);
      BGENE(3,K+1)=SHP(2,NPOI);
      K=K+6;
      end
%
%     MATRICE B RELATIVE A L'EFFET DE FLEXION                                   
%                                                                               
       for K = 4:6
        if (K == 4)
         IJ=1;
         I=2;
         J=3;
        elseif (K == 5)
         IJ=2;
          I=3;
          J=1;
        else
         IJ=3;
          I=1;
          J=2;
        end
        SX(IJ)=XE(1,I)-XE(1,J);
        SY(IJ)=XE(2,I)-XE(2,J);
        RL(IJ)=sqrt(SX(IJ)*SX(IJ)+SY(IJ)*SY(IJ));
        EXX(IJ)=-SX(IJ)/RL(IJ);
        EYY(IJ)=-SY(IJ)/RL(IJ);
       end
       AIR=abs(0.5D0*(SX(1)*SY(2)-SX(2)*SY(1)));
%
       [HXABS,HXORD,HYABS,HYORD] = derivo ...
         (SX,SY,RL,QSI(IGAU),ETA(IGAU),HS4,HS5,HS6,EXX,EYY);

       [B] = bmato(SX,SY,HXABS,HYABS,HXORD,HYORD);
%
      K=2;
      KK=0;
      for NPOI = 1:3                                                            
      for IX = 1:3                                                              
      for IY = 1:3                                                              
        BGENE(3+IX,K+IY)=B(IX,IY+KK);
      end
      end
      KK=KK+3;
      K=K+6;
      end
