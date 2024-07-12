% Geometric element names from AVS
global list_type_AVS;
list_type_AVS{1} = 'quad';
list_type_AVS{2} = 'tri';
list_type_AVS{3} = 'line';
list_type_AVS{4} = 'pt';
list_type_AVS{5} = 'hex';
list_type_AVS{6} = 'tet';
list_type_AVS{7} = 'prism';

% Geometric element names from Cast3M
% ("linear" and "quadratic" elements)
global list_type_C3M1 list_type_C3M2;
QUA4 = 'QUA4';
list_type_C3M1{1} = QUA4; list_type_C3M2{1} = 'QUA8';
list_type_C3M1{2} = 'TRI3'; list_type_C3M2{2} = 'TRI6';
list_type_C3M1{3} = 'SEG2'; list_type_C3M2{3} = 'SEG3';
list_type_C3M1{4} = 'POI1'; list_type_C3M2{4} = 'POI1';
list_type_C3M1{5} = 'CUB8'; list_type_C3M2{5} = 'CU20';
list_type_C3M1{6} = 'TET4'; list_type_C3M2{6} = 'TE10';
list_type_C3M1{7} = 'PRI6'; list_type_C3M2{7} = 'PR15';

% Geometric element names from gmsh
list_type_GMSH{1} = 'SEG2';
list_type_GMSH{2} = 'TRI3';
list_type_GMSH{3} = 'QUA4';
list_type_GMSH{4} = 'TET4';
list_type_GMSH{5} = 'CUB8';
list_type_GMSH{6} = 'PRI6';
list_type_GMSH{7} = 'PYR5';
list_type_GMSH{8} = 'SEG3';
list_type_GMSH{9} = 'TRI6';
list_type_GMSH{10} = 'QUA9';
list_type_GMSH{11} = 'TE10';
list_type_GMSH{12} = 'CU27';
list_type_GMSH{13} = 'PR18';
list_type_GMSH{14} = 'PY14';
list_type_GMSH{15} = 'POI1';

% Reading format for AVS
global list_format1 list_size1 list_size2;
list_format1{1} = '%6i %5i %5i %5i\n'; list_size1{1}=4; list_size2{1}=8;
list_format1{2} = '%6i %5i %5i\n';     list_size1{2}=3; list_size2{2}=6;
list_format1{3} = '%6i %5i\n';         list_size1{3}=2; list_size2{3}=3;
list_format1{4} = '%6i\n';             list_size1{4}=1; list_size2{4}=1;
list_format1{5} = '%6i %5i %5i %5i %5i %5i %5i %5i\n';
                                       list_size1{5}=8; list_size2{5}=20;
list_format1{6} = '%6i %5i %5i %5i\n'; list_size1{6}=4; list_size2{6}=10;
list_format1{7} = '%6i %5i %5i %5i %5i %5i\n'; list_size1{7}=6; list_size2{7}=15;

list_format2{1} = '%6i %5i %5i %5i %5i %5i %5i %5i\n';
list_format2{2} = '%6i %5i %5i %5i %5i %5i\n';
list_format2{3} = '%6i %5i %5i\n';
list_format2{4} = '%6i\n';
list_format2{5} = '%6i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n';
list_format2{6} = '%6i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n';
list_format2{7} = '%6i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n';

% Different types of analyse
global liste_mode;
TRID = 'TRID';
liste_mode{1} = TRID; % Tridimensional
liste_mode{2} = 'COPL'; % Plane stress
liste_mode{3} = 'DEPL'; % Plane strain
liste_mode{4} = 'BARR'; % Barre droite en traction-compression
liste_mode{5} = 'TIMO'; % Poutre de Timoshenko, courbure negligee
liste_mode{6} = 'POUT'; % Poutre d'Euler-Bernoulli, courbure negligee
liste_mode{7} = 'DKIR'; % Plaque de Kirchhoff discret, supposee plane
liste_mode{8} = 'JOIN'; % Modele de joint (epaisseur nulle)
liste_mode{9} = 'DSHE'; % Plaque avec cisaillement discret, supposee plane
% Deformations planes generalisees, contraintes planes generalisees,
% modes de Fourier...

% Different models
global liste_model;
liste_model{1} = 'ELASTIQUE';
liste_model{2} = 'POREUX';    % Couple fluide-structure
liste_model{3} = 'FLUIDE';
liste_model{4} = 'THERMIQUE'; % En pratique, quasi-identique a FLUIDE
liste_model{5} = 'THERMELAS'; % Couple thermo-elastique,
                              % en pratique, quasi-identique a POREUX
liste_model{6} = 'TRANSPORT'; % Transport Galerkin discontinu / Volumes finis

% Different integrations
global liste_intg;
% Contraintes==RIGIDITE
liste_intg{1} = 'RIGIDITE';
liste_intg{2} = 'MASSE';
liste_intg{3} = 'COMPRESSIBILITE';
liste_intg{4} = 'PERMEABILITE';
liste_intg{5} = 'NOEUDS';

liste_intg{6} = 'CAPACITE';      % cf COMPRESSIBILITE
liste_intg{7} = 'CONDUCTIVITE';  % cf PERMEABILITE

liste_intg{8} = 'ADVECTION';
liste_intg{9} = 'FLUX';


% ddl names for massive elements
global liste_ddlp2 liste_ddld2 liste_ddld2a liste_ddlp liste_ddld;
  liste_ddlp2{1} = 'UX'; liste_ddld2{1} = 'FX';
  liste_ddlp2{2} = 'UY'; liste_ddld2{2} = 'FY';
  liste_ddlp = [liste_ddlp2 {'UZ'}]; liste_ddld = [liste_ddld2 {'FZ'}];
  liste_ddlp2a{1} = 'UR'; liste_ddld2a{1} = 'FR';
  liste_ddlp2a{2} = 'UZ'; liste_ddld2a{2} = 'FZ';
global liste_local_ddlp2 liste_local_ddld2;
global liste_local_ddlp liste_local_ddld;
  liste_local_ddlp2{1} = 'U1'; liste_local_ddld2{1} = 'F1';
  liste_local_ddlp2{2} = 'U2'; liste_local_ddld2{2} = 'F2';
  liste_local_ddlp = [liste_local_ddlp2 {'U3'}];
  liste_local_ddld = [liste_local_ddld2 {'F3'}];

% ddl names for structural elements
global liste_ddlpr2 liste_ddldr2;
global liste_ddlpr liste_ddldr;
  liste_ddlpr2{1} = 'UX'; liste_ddldr2{1} = 'FX';
  liste_ddlpr2{2} = 'UY'; liste_ddldr2{2} = 'FY';
  liste_ddlpr2{3} = 'RZ'; liste_ddldr2{3} = 'MZ';
  liste_ddlpr{1} = 'UX'; liste_ddldr{1} = 'FX';
  liste_ddlpr{2} = 'UY'; liste_ddldr{2} = 'FY';
  liste_ddlpr{3} = 'UZ'; liste_ddldr{3} = 'FZ';
  liste_ddlpr{4} = 'RX'; liste_ddldr{4} = 'MX';
  liste_ddlpr{5} = 'RY'; liste_ddldr{5} = 'MY';
  liste_ddlpr{6} = 'RZ'; liste_ddldr{6} = 'MZ';
global liste_local_ddlpr2 liste_local_ddldr2;
global liste_local_ddlpr liste_local_ddldr;
  liste_local_ddlpr2{1} = 'U1'; liste_local_ddldr2{1} = 'F1';
  liste_local_ddlpr2{2} = 'U2'; liste_local_ddldr2{2} = 'F2';
  liste_local_ddlpr2{3} = 'RZ'; liste_local_ddldr2{3} = 'MZ';
  liste_local_ddlpr{1} = 'U1';  liste_local_ddldr{1} = 'F1';
  liste_local_ddlpr{2} = 'U2';  liste_local_ddldr{2} = 'F2';
  liste_local_ddlpr{3} = 'U3';  liste_local_ddldr{3} = 'F3';
  liste_local_ddlpr{4} = 'R1';  liste_local_ddldr{4} = 'M1';
  liste_local_ddlpr{5} = 'R2';  liste_local_ddldr{5} = 'M2';
  liste_local_ddlpr{6} = 'R3';  liste_local_ddldr{6} = 'M3';
  

% Material coefficient names
global liste_rigi_isotropic;            % Rigidity
  liste_rigi_isotropic{1} = 'YOUN';     %  Young's modulus
  liste_rigi_isotropic{2} = 'NU';       %  Poisson's ratio
global liste_barr;                      % Structural element
  liste_barr{1} = 'SECT';               %  Section
global liste_pout;                      % Structural element
  liste_pout{1} = 'SECT';               %  Section
  liste_pout{2} = 'INRZ';               %  Inertial surface moment
%  liste_pout{3} = 'INRY';               %  Inert. surf. momt (3D only)
%  liste_pout{4} = 'TORS';               %  Inert. polar surf. momt (3D)
global liste_mass_isotropic;            % Mass
  liste_mass_isotropic{1} = 'RHO';      %  Specific mass
global liste_compr_isotropic;           % Diffusion
  liste_compr_isotropic{1} = 'MOB';     %  Biot's modulus
global liste_rigicompr_isotropic;       % Coupled Rigidity-Compress
  liste_rigicompr_isotropic{1} = 'COB'; %  Biot's coefficient
global liste_perm_isotropic;            % Permeability
  liste_perm_isotropic{1} = 'PERM';     %  Intrinsic permeability
  liste_perm_isotropic{2} = 'VISC';     %  Dynamic viscosity


% Classical strain and stress names for massive elements
% (Voigt notation: GAxx = 1/sqrt(2) * EPxx, TAxx= 1/sqrt(2) * SMxx)
global liste_strain liste_stress;
global liste_strain1 liste_stress1;
global liste_strain2 liste_stress2 liste_strain2a liste_stress2a;
  liste_strain = [{'EPXX'} {'EPYY'} {'EPZZ'} {'GAYZ'} {'GAZX'} {'GAXY'}];
  liste_stress = [{'SMXX'} {'SMYY'} {'SMZZ'} {'TAYZ'} {'TAZX'} {'TAXY'}];
  liste_strain1 = [{'EP11'} {'EP22'} {'EP33'} {'GA23'} {'GA31'} {'GA12'}];
  liste_stress1 = [{'SM11'} {'SM22'} {'SM33'} {'TA23'} {'TA31'} {'TA12'}];
  liste_strain2 = [{'EPXX'} {'EPYY'} {'GAXY'} {'EPZZ'}];
  liste_stress2 = [{'SMXX'} {'SMYY'} {'TAXY'} {'SMZZ'}];
  liste_strain2a = [{'EPRR'} {'EPZZ'} {'GARZ'} {'EPTT'}];
  liste_stress2a = [{'SMRR'} {'SMZZ'} {'TARZ'} {'SMTT'}];
% Classical strain and stress names for structural elements
global liste_barr_strain liste_barr_stress;
  liste_barr_strain = [{'EPS1'}];
  liste_barr_stress = [{'EFF1'}];
global list_pout_strain list_pout_stress;
global list_pout_strain2 list_pout_stress2;
list_pout_strain=[{'EPS1'} {'GA12'} {'GA13'} {'C1'}   {'C2'}   {'C3'}];
list_pout_stress=[{'EFF1'} {'EFF2'} {'EFF3'} {'MOM1'} {'MOM2'} {'MOM3'}];
list_pout_strain2 = [{'EPS1'} {'GA12'} {'C3'}];
list_pout_stress2 = [{'EFF1'} {'EFF2'} {'MOM3'}];
global list_plaq_strain list_plaq_stress;
list_plaq_strain=[{'EP11'} {'EP22'} {'GA12'} {'CP11'} {'CP22'} {'CG12'}];
list_plaq_stress=[{'EF11'} {'EF22'} {'EG12'} {'MF11'} {'MF22'} {'MG12'}];
global list_plaq_strain_cis list_plaq_stress_cis;
list_plaq_strain_cis = ...
  [{'EP11'} {'EP22'} {'GA12'} {'CP11'} {'CP22'} {'CG12'} {'G1'} {'G2'}];
list_plaq_stress_cis = ...
  [{'EF11'} {'EF22'} {'EG12'} {'MF11'} {'MF22'} {'MG12'} {'T1'} {'T2'}];
