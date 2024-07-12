      function [SHP,DJac] = jacobidst(XEL,SHP,IDIM,NBNO)

      XZero=0.;
      XUn=1.;
      dJInv=XZero;

% ===========================
%  1 - Cas de la DIMENSION 1
% ===========================
      if (IDIM == 1) 
        DJac=XZero;
        for i = 1:NBNO
          DJac=DJac+SHP(2,i)*XEL(1,i);
        end
        if (DJac ~= 0) 
          dJInv=XUn/DJac;
        else
          error('Pb in jacobidst.m: DJac=0')
        end
        for i = 1:NBNO
          SHP(2,i)=SHP(2,i)*dJInv;
        end
% ===========================
%  2 - Cas de la DIMENSION 2
% ===========================
      elseif (IDIM == 2) 
        dXdQsi=XZero;
        dXdEta=XZero;
        dYdQsi=XZero;
        dYdEta=XZero;
        for i = 1:NBNO
          dXdQsi=dXdQsi+SHP(2,i)*XEL(1,i);
          dXdEta=dXdEta+SHP(3,i)*XEL(1,i);
          dYdQsi=dYdQsi+SHP(2,i)*XEL(2,i);
          dYdEta=dYdEta+SHP(3,i)*XEL(2,i);
        end
        DJac=dXdQsi*dYdEta-dXdEta*dYdQsi;
        if (DJac ~= 0) 
          dJInv=XUn/DJac;
        else
          error('Pb in jacobidst.m: DJac=0 bis')
        end
        dQsidX= dYdEta*dJInv;
        dQsidY=-dXdEta*dJInv;
        dEtadX=-dYdQsi*dJInv;
        dEtadY= dXdQsi*dJInv;
        for i = 1:NBNO
          V1=SHP(2,i)*dQsidX+SHP(3,i)*dEtadX;
          SHP(3,i)=SHP(2,i)*dQsidY+SHP(3,i)*dEtadY;
          SHP(2,i)=V1;
        end
% ===========================
%  3 - Cas de la DIMENSION 3
% ===========================
      elseif (IDIM == 3) 
        D11=XZero;
        D21=XZero;
        D31=XZero;
        D12=XZero;
        D22=XZero;
        D32=XZero;
        D13=XZero;
        D23=XZero;
        D33=XZero;
        for i = 1:NBNO
          D11=D11+SHP(2,i)*XEL(1,i);
          D21=D21+SHP(3,i)*XEL(1,i);
          D31=D31+SHP(4,i)*XEL(1,i);
          D12=D12+SHP(2,i)*XEL(2,i);
          D22=D22+SHP(3,i)*XEL(2,i);
          D32=D32+SHP(4,i)*XEL(2,i);
          D13=D13+SHP(2,i)*XEL(3,i);
          D23=D23+SHP(3,i)*XEL(3,i);
          D33=D33+SHP(4,i)*XEL(3,i);
        end
        DInv11=D22*D33-D23*D32;
        DInv12=D32*D13-D12*D33;
        DInv13=D12*D23-D22*D13;
        DJac=D11*DInv11+D21*DInv12+D31*DInv13;
        if (DJac ~= 0) 
          dJInv=XUn/DJac;
        else
          error('Pb in jacobidst.m: DJac=0 ter')
        end
        DInv11=DInv11*dJInv;
        DInv12=DInv12*dJInv;
        DInv13=DInv13*dJInv;
        DInv21=(D23*D31-D21*D33)*dJInv;
        DInv22=(D11*D33-D13*D31)*dJInv;
        DInv23=(D21*D13-D11*D23)*dJInv;
        DInv31=(D21*D32-D22*D31)*dJInv;
        DInv32=(D12*D31-D11*D32)*dJInv;
        DInv33=(D11*D22-D12*D21)*dJInv;
        for i = 1:NBNO
          V1=DInv11*SHP(2,i)+DInv12*SHP(3,i)+DInv13*SHP(4,i);
          V2=DInv21*SHP(2,i)+DInv22*SHP(3,i)+DInv23*SHP(4,i);
          V3=DInv31*SHP(2,i)+DInv32*SHP(3,i)+DInv33*SHP(4,i);
          SHP(2,i)=V1;
          SHP(3,i)=V2;
          SHP(4,i)=V3;
        end
% =======================================
%  4 - Cas particulier de l'element JOI3
% =======================================
      elseif (IDIM == 86) 
        dXdQsi=XZero;
        dYdQsi=XZero;
        for i = 1:NBNO
          dXdQsi=dXdQsi+SHP(2,i)*XEL(1,i);
          dYdQsi=dYdQsi+SHP(2,i)*XEL(2,i);
        end
        DJac=SQRT(dXdQsi*dXdQsi+dYdQsi*dYdQsi);
      else
        error('Pb in jacobidst.m: bad IDIM')
      end
