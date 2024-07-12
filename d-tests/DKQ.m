function [KE]=DKQ(VCORE,VPREE)
%   Ahmad Mobasher Amini     GéM Ecole Central de Nantes  (ECN) 20/06/2005
%  4 noeuds et 3 ddl pour chaque noeude   w, - dw/dx, - dw/dy
%                       l * * * * * * * * k
%                         *             *
%                         *             *
%                         *             *
%                         *             *
%                       i * * * * * * * * j
% Function de Calule regidite d'élémentair 
% VCORE is the vector coordinate of four nodes of the quadrangle  [X1 Y1 X2 Y2 X3 Y3 X4 Y4]
% VPROPE is the vector property of the quadrangle element  % [module young,nu,ro,epais]
%clear all
%DD VCORE=[0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0];
%DD VPROPE=[10.0 0.3 7800.0 1.0];
PG1 =0.774596669241483;
PG2= 0.000000000000000;
WG1=0.555555555555556;
WG2=0.888888888888889;

% Coordonnées KSI et ETA de points d'intégration
VPG=[-PG1  -PG1  -PG1   PG2  PG2  PG2   PG1  PG1  PG1; ...
     -PG1   PG2   PG1  -PG1  PG2  PG1  -PG1  PG2  PG1];
% Poids associés aux points d'integration
VWG=[WG1  WG1  WG1  WG1  WG2  WG1  WG1  WG1  WG1; ...
     WG1  WG2  WG1  WG2  WG2  WG2  WG1  WG2  WG1];
IKE=12;  % Dimension of stiffness matrix
IPG=4;  % Points integration gusse
IH=3;   % Dimension matrice de comportement de flexion HF
CO1=VCORE(1);   CO2=VCORE(2);
CO3=VCORE(3);   CO4=VCORE(4);
CO5=VCORE(5);   CO6=VCORE(6);
CO7=VCORE(7);   CO8=VCORE(8);
% Calcule des cos direction et longueur de 4 cotes 1-2   2-3   3-4   4-1 
% VCOS(2,4)     cos IJ    sin  IJ  , IJ=12   23   34   41
%          VL            Longueur des cotes  12  23  34  41
KI=[1 2 3 4 5 6 7 8 1 2];
II=1;
for I=1:4
    CX=VCORE(KI(II+2))-VCORE(KI(II));
    CY=VCORE(KI(II+3))-VCORE(KI(II+1));
    CL=sqrt(CX^2+CY^2);
    VL(I)=CL;
    VCOS(1,I)=CX/CL; % Cosinous
    VCOS(2,I)=CY/CL; % Sinous
    II=II+2;
end
clear KI II CL CX CY
KE=zeros(IKE,IKE);
% integration numerique :   Gauss 2 x 2
j=1;
for IG=1:IH
    for JG=1:IH
    % Matrice B pour la Flexion
    [J,DETJ]=jacobDKQ(VCORE,VPG(1,j),VPG(2,j));
    [VB]=BDKQ(VL,VCOS,J,VPG(1,j),VPG(2,j));  %construire la matrice VB(3,12) 
% construire la matrice proprietes elastiques
% VHf    matrice des coefficients de flexion
    VHf=zeros(3,3);
    %Module d'elastisité , Poison coefficion , epaisseur de plaque
    E=VPREE(1);    NU=VPREE(2);    t=VPREE(4);
%cas isotrope
    C=E*t^3/(12.0*(1.0-NU^2));
    VHf(1,1)=C;               VHf(2,2)=C;
    VHf(3,3)=C*(1-NU)/2.0;    VHf(1,2)=C*NU;
    VHf(2,1)=VHf(1,2);
% multiplier VH, par POIDS*DETJ
    C=VWG(1,j)*VWG(2,j)*DETJ;
    KE=VB'*(VHf*C)*VB+KE;
    j=j+1;
    end   % JG
end  % IG
