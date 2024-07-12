function [HS4,HS5,HS6,BGENE] = rcdst(XE,DDHOMU)

      NSTRS = 8;
      LRE = 18; %% ?
      BGENE = zeros(NSTRS,LRE);
      BCISA = zeros(2,9);

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

      [ASS,AXX,AYY,BSS,BXX,BYY] = cisaf(XE,EXX,EYY,DDHOMU,AIR);
      [HS4,HS5,HS6,HFX,HFY] = cisar(ASS,BSS,AXX,AYY,BXX,BYY,RL,EXX,EYY);

        for IC=1:9
          BCISA(1,IC)=HFX(IC);
          BCISA(2,IC)=HFY(IC);
	end

        for IL = 1:2
          IL1=IL+6;
          for IC = 1:9
             IC1=floor((IC-1)/3);
             IC2=IC1*3+IC+2;
             BGENE(IL1,IC2)=BCISA(IL,IC);
          end
        end
