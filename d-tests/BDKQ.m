function  [VB]=BDKT(VL,VCOS,JI,KSI,ETA)
%VN(1,i) function d'interpolation-VN(2,i) deriver par rapport de ksi
%VN(3,i) dériver par rapport de eta
VN(1,1)=0.25*(1-KSI)*(1-ETA);    VN(1,2)=0.25*(1+KSI)*(1-ETA);
VN(1,3)=0.25*(1+KSI)*(1+ETA);    VN(1,4)=0.25*(1-KSI)*(1+ETA);

VN(2,1)=-0.25*(1-ETA);           VN(2,2)=0.25*(1-ETA) ;
VN(2,3)=0.25*(1+ETA) ;           VN(2,4)=-0.25*(1+ETA);

VN(3,1)=-0.25*(1-KSI);           VN(3,2)=-0.25*(1+KSI);
VN(3,3)=0.25*(1+KSI) ;           VN(3,4)=0.25*(1-KSI) ;

% VP derives en Ksi et Eta de P linear VP(2*3)
%VP(1,i) -VP(2,i) deriver par rapport de ksi- VP(3,i) dériver par rapport de eta

VP(1,1)=0.5*(1-KSI^2)*(1-ETA);   VP(1,2)=0.5*(1-ETA^2)*(1+KSI);
VP(1,3)=0.5*(1-KSI^2)*(1+ETA);   VP(1,4)=0.5*(1-ETA^2)*(1-KSI);

VP(2,1)=-KSI*(1-ETA);            VP(2,2)=0.5*(1-ETA^2);
VP(2,3)=-KSI*(1+ETA);            VP(2,4)=-0.5*(1-ETA^2);

VP(3,1)=-0.5*(1-KSI^2);          VP(3,2)=-ETA*(1+KSI);
VP(3,3)=0.5*(1-KSI^2);           VP(3,4)=-ETA*(1-KSI);

nn=[1 4 2 1 3 2 4 3];  % Vecteur de numerotation de VN , VP
ii=0;
j=1;
% Nx    ;   Nx,ksi    ;   Nx,eta
for i=1:2:8
    VNij(1,1+ii)=1.5/VL(nn(i))*VP(1,nn(i))*VCOS(1,nn(i))-1.5/VL(nn(i+1))...
        *VP(1,nn(i+1))*VCOS(1,nn(i+1));
    VNij(1,2+ii)=VN(1,j)-0.75*VP(1,nn(i))*VCOS(1,nn(i))^2-0.75*VP(1,nn(i+1))...
        *VCOS(1,nn(i+1))^2;
    VNij(1,3+ii)=-0.75*VP(1,nn(i))*VCOS(1,nn(i))*VCOS(2,nn(i))-0.75*VP(1,nn(i+1))...
        *VCOS(1,nn(i+1))*VCOS(2,nn(i+1));
%--------------
    VNijksi(1,1+ii)=1.5/VL(nn(i))*VP(2,nn(i))*VCOS(1,nn(i))-1.5/VL(nn(i+1))...
        *VP(2,nn(i+1))*VCOS(1,nn(i+1));
    VNijksi(1,2+ii)=VN(2,j)-0.75*VP(2,nn(i))*VCOS(1,nn(i))^2-0.75*VP(2,nn(i+1))...
        *VCOS(1,nn(i+1))^2;
    VNijksi(1,3+ii)=-0.75*VP(2,nn(i))*VCOS(1,nn(i))*VCOS(2,nn(i))-0.75*VP(2,nn(i+1))...
        *VCOS(1,nn(i+1))*VCOS(2,nn(i+1));
%--------------
    VNijeta(1,1+ii)=1.5/VL(nn(i))*VP(3,nn(i))*VCOS(1,nn(i))-1.5/VL(nn(i+1))...
        *VP(3,nn(i+1))*VCOS(1,nn(i+1));
    VNijeta(1,2+ii)=VN(3,j)-0.75*VP(3,nn(i))*VCOS(1,nn(i))^2-0.75*VP(3,nn(i+1))...
        *VCOS(1,nn(i+1))^2;
    VNijeta(1,3+ii)=-0.75*VP(3,nn(i))*VCOS(1,nn(i))*VCOS(2,nn(i))-0.75*VP(3,nn(i+1))...
        *VCOS(1,nn(i+1))*VCOS(2,nn(i+1));
%--------------    
    ii=ii+3;
    j=j+1;
end
%-------------------------------------------------------------------
ii=0;
j=1;
% Ny    ;   Ny,ksi    ;   Ny,eta
for i=1:2:8
    VNij(2,1+ii)=1.5/VL(nn(i))*VP(1,nn(i))*VCOS(2,nn(i))-1.5/VL(nn(i+1))*VP(1,nn(i+1))...
        *VCOS(2,nn(i+1));
    VNij(2,3+ii)=VN(1,j)-0.75*VP(1,nn(i))*VCOS(2,nn(i))^2-0.75*VP(1,nn(i+1))...
        *VCOS(2,nn(i+1))^2;
    VNij(2,2+ii)=VNij(1,2+ii+1);
%--------------
    VNijksi(2,1+ii)=1.5/VL(nn(i))*VP(2,nn(i))*VCOS(2,nn(i))-1.5/VL(nn(i+1))*VP(2,nn(i+1))...
        *VCOS(2,nn(i+1));
    VNijksi(2,3+ii)=VN(2,j)-0.75*VP(2,nn(i))*VCOS(2,nn(i))^2-0.75*VP(2,nn(i+1))*...
        VCOS(2,nn(i+1))^2;
    VNijksi(2,2+ii)=VNijksi(1,2+ii+1);
%--------------
    VNijeta(2,1+ii)=1.5/VL(nn(i))*VP(3,nn(i))*VCOS(2,nn(i))-1.5/VL(nn(i+1))*VP(3,nn(i+1))...
        *VCOS(2,nn(i+1));
    VNijeta(2,3+ii)=VN(3,j)-0.75*VP(3,nn(i))*VCOS(2,nn(i))^2-0.75*VP(3,nn(i+1))...
        *VCOS(2,nn(i+1))^2;
    VNijeta(2,2+ii)=VNijeta(1,2+ii+1);
%--------------    
    j=j+1;
    ii=ii+3;
end
VB=zeros(3,12);
DETJ=det(JI);
for i=1:12
    VB(1,i)=VNijksi(1,i)*JI(2,2)-VNijeta(1,i)*JI(1,2);
    VB(2,i)=-VNijksi(2,i)*JI(2,1)+VNijeta(2,i)*JI(1,1);
    VB(3,i)=-VNijksi(1,i)*JI(2,1)+VNijeta(1,i)*JI(1,1)+VNijksi(2,i)*JI(2,2)...
        -VNijeta(2,i)*JI(1,2);
end
VB=VB/DETJ;
