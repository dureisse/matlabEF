function [KE]=DKT(VCORE,VPROPE)
%   Ahmad Mobasher Amini     G=E9M Ecole Central de Nantes  (ECN) 10/05/2005
%   3 ddl pour chaque noeude   w, - dw/dx, - dw/dy
%
%                                            * j
%                                          *  *
%                                        *     *
%                                    k *        *
%                                          *     *
%                                              *  *
%                                                  * i
%
% Function de Calule regidite d'=E9l=E9mentair
% VCORE is the vector coordinate of three nodes of the triangle  [X1 Y1 X2= Y2 X3 Y3]
% VPROPE is the vector property of the triangle element  % [module= young,nu,ro,epais]
%DD ????????????????? clear all
%VCORE=[0.0 0.0 1.0 0.0 0.0 1.0];
%VPROPE=[200.0 0.3 7800.0 0.01];
PG =0.16666666666;
PG1=0.66666666666;
% coordonn=E9es KSI et ETA de trois points d'int=E9gration
VPG=[ PG  PG1  PG ; PG  PG  PG1];
%poids associ=E9s aux trios points d'integration
VWG=[PG PG PG];

IKE=9;  % Dimension of stiffness matrix
IPG=3;  % Points integration gusse
IH=3;   % Dimension matrice de comportement de flexion HF

CO1=VCORE(1);   CO2=VCORE(2);
CO3=VCORE(3);   CO4=VCORE(4);
CO5=VCORE(5);   CO6=VCORE(6);

% Calcule des cos direction et longueur de 3 cotes 1-2   2-3   3-1
%Entre   : VCORE (6)
%Sorties : VCOS(2,3)     cos IJ    sin  IJ  , IJ=12   23   31
%          VL            Longueur des cotes  12  23  31
KI=[1 2 3 4 5 6 1 2];
II=1;
for I=1:3
    CX=VCORE(KI(II+2))-VCORE(KI(II));
    CY=VCORE(KI(II+3))-VCORE(KI(II+1));
    II=II+2;
    CL=sqrt(CX^2+CY^2);
    VL(I)=CL;
    VCOS(1,I)=CX/CL; % Cosinous
    VCOS(2,I)=CY/CL; % Sinous
end
clear KI,II,CL,CX,CY;

%Calcule de la matrice Jacobien invers

J11=CO3-CO1;
J12=CO4-CO2;
J21=CO5-CO1;
J22=CO6-CO2;
J=[J11 J12;J21 J22];
DETJ=det(J);
if (DETJ<0.0)
    disp('Erreur Det Negative verifier donner entre')
    return
end
VJI=zeros(2,2);
VJI(2,2)= J11/DETJ;
VJI(2,1)=-J21/DETJ;
VJI(1,2)=-J12/DETJ;
VJI(1,1)= J22/DETJ;
clear J11,J12,J21,J22;

KE=zeros(IKE,IKE);
% integration numerique :   3 points
for IG=1:IPG
    % Matric B pour la Flexion
    [VB]=BDKT(VL,VCOS,VJI,VPG(1,IG),VPG(2,IG));  %construire la matrice= VB(3,9)=20
% construire la matrice proprietes elastiques
% VH    matrice des coefficients de flexion
    VH=zeros(3,3);
    E=VPROPE(1);   %Module d'elastisit=E9
    NU=VPROPE(2);  % Poison coefficion
    t=VPROPE(4);   % epaisseur de plaque=20

%cas isotrope
    C=E*t^3/(12.0*(1.0-NU^2));
    VH(1,1)=C;
    VH(2,2)=C;
    VH(3,3)=C*(1-NU)/2.0;
    VH(1,2)=C*NU;
    VH(2,1)=VH(1,2);
% multiplier VH, par POIDS*DETJ
    C=DETJ*VWG(1,IG);
    for I=1:IH
        for J=1:IH
            VH(I,J)=VH(I,J)*C;
        end % J
    end     % I
    KE=VB'*VH*VB+KE;
end   % IPG
disp(' Calcule fini !!!')
