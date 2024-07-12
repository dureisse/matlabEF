function  [VB]=BDKT(VL,VCOS,VJI,KSI,ETA)
KI=[3 1 2];  % construire la matrice VB(3,9)
% N function d'interpolation              N1=1-Ksi-Eta , N2=Ksi , N3=Eta
VN=[-1.0 1.0 0.0; -1.0 0.0 1.0]; % VN derives en Ksi et Eta de N linear VN(2*3)
% VP derives en Ksi et Eta de P linear VP(2*3)
C3=1.0-KSI-ETA;
VP(1,1)=4.0*(C3-KSI);
VP(2,1)=-4.0*KSI;
VP(1,2)=4.0*ETA;
VP(2,2)=4.0*KSI;
VP(1,3)=-4.0*ETA;
VP(2,3)=4.0*(C3-ETA);
for J=1:3
    CL=1.5/VL(J);
    AC=VCOS(1,J);
    AS=VCOS(2,J);
    for I=1:2
        VPC(I,J)=CL*VP(I,J)*AC;
        VPS(I,J)=CL*VP(I,J)*AS;
        VPCC(I,J)=0.75*VP(I,J)*AC*AC;
        VPSS(I,J)=0.75*VP(I,J)*AS*AS;
        VPCS(I,J)=0.75*VP(I,J)*AC*AS;
    end
end
% VB(3,9) : BETAX,x   BETAY,y   BETAX,y+BETAY,x
CJ11=VJI(1,1);
CJ12=VJI(1,2);
CJ21=VJI(2,1);
CJ22=VJI(2,2);
% pour chaque noeud -->>  3ddl  w , dw/dx , dw/dy
JJ=1;
for J=1:3
    JM=KI(J);
%   BETAX,x=VJ(1,1)*BETAX,KSI+VJ(1,2)*BETAX,ETA
    VB(1,JJ)=CJ11*(VPC(1,J)-VPC(1,JM))+ ...
             CJ12*(VPC(2,J)-VPC(2,JM));
    VB(1,JJ+1)=-CJ11*(-VN(1,J)+VPCC(1,J)+VPCC(1,JM))- ...
                CJ12*(-VN(2,J)+VPCC(2,J)+VPCC(2,JM));
    VB(1,JJ+2)=-CJ11*(VPCS(1,J)+VPCS(1,JM))- ...
                CJ12*(VPCS(2,J)+VPCS(2,JM));
%   BETAY,y=VJ(2,1)*BETAY,KSI+VJ(2,2)*BETAY,ETA
    VB(2,JJ)=CJ21*(VPS(1,J)-VPS(1,JM))+ ...
             CJ22*(VPS(2,J)-VPS(2,JM));
    VB(2,JJ+1)=-CJ21*(VPCS(1,J)+VPCS(1,JM))- ...
               CJ22*(VPCS(2,J)+VPCS(2,JM));
    VB(2,JJ+2)=-CJ21*(-VN(1,J)+VPSS(1,J)+VPSS(1,JM))- ...
               CJ22*(-VN(2,J)+VPSS(2,J)+VPSS(2,JM));
%   BETAX,y+BETAY,x
    VB(3,JJ)=CJ21*(VPC(1,J)-VPC(1,JM))+ ...
             CJ22*(VPC(2,J)-VPC(2,JM))+ ...
             CJ11*(VPS(1,J)-VPS(1,JM))+ ...
             CJ12*(VPS(2,J)-VPS(2,JM));
    VB(3,JJ+1)=-CJ21*(-VN(1,J)+VPCC(1,J)+VPCC(1,JM))- ...
               CJ22*(-VN(2,J)+VPCC(2,J)+VPCC(2,JM))- ...
               CJ11*(VPCS(1,J)+VPCS(1,JM))- ...
               CJ12*(VPCS(2,J)+VPCS(2,JM));
    VB(3,JJ+2)=-CJ21*(VPCS(1,J)+VPCS(1,JM))- ...
               CJ22*(VPCS(2,J)+VPCS(2,JM))- ...
               CJ11*(-VN(1,J)+VPSS(1,J)+VPSS(1,JM))- ...
               CJ12*(-VN(2,J)+VPSS(2,J)+VPSS(2,JM));
    JJ=JJ+3;
end
