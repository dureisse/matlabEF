function [JI,DETJ]=jacobKDQ(VCORE,ksi,eta)
% Calcule de matrice de Jacobienne et déterminante 
CO1=VCORE(1);   CO2=VCORE(2);
CO3=VCORE(3);   CO4=VCORE(4);
CO5=VCORE(5);   CO6=VCORE(6);
CO7=VCORE(7);   CO8=VCORE(8);

x21=CO3-CO1;    x34=CO5-CO7;
y21=CO4-CO2;    y34=CO6-CO8;

x32=CO5-CO3;    x41=CO7-CO1;
y32=CO6-CO4;    y41=CO8-CO2;

J11=0.25*(x21+x34+eta*(x34-x21));
J12=0.25*(y21+y34+eta*(y34-y21));
J21=0.25*(x32+x41+ksi*(x32-x41));
J22=0.25*(y32+y41+ksi*(y32-y41));
JI=[J11 J12;J21 J22];
DETJ=det(JI);
if (DETJ<0.0)
    disp('Erreur Det Negative verifier donner entre')
    return
end
clear J11 J12 J21 J22
