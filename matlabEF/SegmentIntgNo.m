function [intg1] = SegmentIntgNo(mail1,varargin);
% Give a (pseudo) integration segment at element nodes or centroid
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 09 / 04 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 21 / 07 / 2003
%   Ajout SEG2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 05 / 2005
%   Ajout RAC2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 09 / 2005
%   Ajout RAC3
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 28 / 02 / 2006
%   Ajout centre de gravite
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 04 / 06 / 2006
%   Ajout CUB8, TET4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 12 / 2006
%   Ajout PRI6, POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 03 / 03 / 2007
%   Ajout CUB8, SEG3, CU20, PR15
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 21 / 10 / 2007
%   Ajout nombre de points de Gauss variable
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 11 / 2007
%   Factorisation de la definition des points d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2007
%   Ajout TE10
%
% Entrees
%  mail1        Maillage
% Entrees optionnelles	
%  type		'nodes', 'centroid' ou le nombre de pts d'integration
% Sorties
%  intg1        Son segment d'integration

nin1 = nargin - 1;
if nin1 == 0
  type = 'nodes'; % default value
elseif nin1 == 1
  type = varargin{1};
else
  nin1
  error('bad number of optional arguments')
end

clear intg1;
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  [nbel,nbnod] = size(maile1.MAIL);
  type1 = maile1.TYPE;

  if ischar(type)
%   Chaine de caracteres
  switch type
    case 'nodes',
      nbptg = nbnod;
      COOR = EF_CoorRefNod(type1);
    case 'centroid',
      nbptg = 1;
      COOR = EF_CoorRefNod(type1);
      COOR = sum(COOR,1) / nbnod;
    otherwise,
      type
      error('bad location type')
  end

  switch type1
    case 'POI1',
      WEIGHT = 1.; % par convention
    case {'SEG2','SEG3'},
      WEIGHT = 2. / nbptg * ones(1,nbptg);
    case {'TRI3','TRI6'},
      WEIGHT = 0.5 / nbptg * ones(1,nbptg);
    case {'QUA4','QUA8'},
      WEIGHT = 4. / nbptg * ones(1,nbptg);
    case {'RAC2','RAC3'},
      nbptg = nbptg / 2;
      WEIGHT = 1. / nbptg * ones(1,nbptg);
    case {'CUB8','CU20'},
      WEIGHT = 8. / nbptg * ones(1,nbptg);
    case {'TET4','TE10'},
      WEIGHT = (1./6.) / nbptg * ones(1,nbptg);
    case {'PRI6','PR15'},
      WEIGHT = 1. / nbptg * ones(1,nbptg);
    otherwise,
      type1
      error('Type of element not yet implemented for integration at nodes')
  end
%
  elseif isnumeric(type)
%   Nombre de points d'integration

  nbptg = type;
  switch type1
    case 'POI1',
      if type~=1
        type1
        type
        error('For POI1, only 1 integration point')
      end
      WEIGHT = 1.; % par convention
    case 'SEG2',
%      r2 = sqrt(2.)/2.;
%      r3 = sqrt(3.)/3.;
      O577=0.577350269189626;
      switch type
        case 1,
          COOR   = [0.];
          WEIGHT = [2.];
        case 2,
          COOR   = [-O577
                     O577];
          WEIGHT = [1. 1.]; 
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
    case 'SEG3',
      O577=0.577350269189626;
      X774=0.774596669241483;
      P555=5./9.;
      P888=8./9.;
      switch type
        case 1,
          COOR   = [0.];
          WEIGHT = [2.];
        case 2,
          COOR   = [-O577
                     O577];
          WEIGHT = [1. 1.]; 
        case 3,
          COOR   = [-X774
                     0.
                     X774];
          WEIGHT = [P555 P888 P555]; 
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
    case 'TRI3',
      nbptg = type;
      [COOR,WEIGHT] = TRI_Integr(nbptg);
%%%     INTEGRATION 7 POINTS TRIANGLE
%%      X101=.1012865073235D0;
%%      X797=.7974269853531D0;
%%      X470=.4701420641051D0;
%%      X059=.0597158717898D0;
%%      P125=.1259391805448D0;
%%      P132=.1323941527885D0;
%%      switch type
%%        case 1,
%%          COOR   = [1./3. 1./3.];
%%          WEIGHT = [0.5];
%%        case 3,
%%          COOR   = [1./6. 1./6.
%%                    2./3. 1./6.
%%                    1./6. 2./3.];
%%          WEIGHT = [1./6. 1./6. 1./6.];
%%        case 4,
%%error('Incoherence dans les poids... a revoir')
%%          COOR   = [0.33333333  0.33333333
%%                    0.73333333  0.13333333
%%                    0.13333333  0.73333333
%%                    0.13333333  0.13333333];
%%          WEIGHT = [0.28125000  0.26041666  0.26041666  0.26041666]*0.5;
%%        case 7,
%%          COOR   = [X101 X101
%%                    X797 X101
%%                    X101 X797
%%                    X470 X059
%%                    X470 X470
%%                    X059 X470
%%                    1./3. 1./3.];
%%          WEIGHT = [P125*0.5 P125*0.5 P125*0.5 P132*0.5 ...
%%                    P132*0.5 P132*0.5 .1125D0];
%%        otherwise
%%          type
%%          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
%%        end
    case 'QUA4',
      r2 = sqrt(2.)/2.;
      r3 = sqrt(3.)/3.;
      r35 = sqrt(3./5.);
      switch type
        case 1,
          COOR   = [0. 0.];
          WEIGHT = [4];
        case 4,
          COOR   = [-r3 -r3
                     r3 -r3
                     r3  r3
                    -r3  r3];
          WEIGHT = [1. 1. 1. 1.]; 
        case 9,
          COOR   = [-r35 -r35
                     0.  -r35
                     r35 -r35
                    -r35  0.
                     0.   0.
                     r35  0.
                    -r35  r35
                     0.   r35
                     r35  r35];
          WEIGHT = [25./81. 40./81. 25./81. ...
                    40./81. 64./81  40./81. ...
                    25./81. 40./81. 25./81.]; 
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
    case 'QUA8',
      X774 = .774596669241483D0;
      P555 = .555555555555555555D0;
      P888 = .888888888888888888D0;
      switch type
        case 1,
          COOR   = [0. 0.];
          WEIGHT = [4];
        case 9,
          COOR   = [-X774 -X774
	             0.   -X774
		     X774 -X774
                    -X774 0.
		     0.   0.
		     X774 0.
                    -X774 X774
		     0.   X774
		     X774 X774];
          WEIGHT = [P555*P555
                    P888*P555
                    P555*P555
		    P888*P555
		    P888*P888
		    P888*P555
		    P555*P555
		    P888*P555
		    P555*P555]'; 
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
    case 'TRI6',
      nbptg = type;
      [COOR,WEIGHT] = TRI_Integr(nbptg);
%%      switch type
%%        case 3,
%%          COOR   = [1./6. 1./6.
%%                    2./3. 1./6.
%%                    1./6. 2./3.];
%%          WEIGHT = [1./6. 1./6. 1./6.];
%%        case 6,
%%          A = 0.445948490915965;
%%          B = 0.091576213509771;
%%          C = 0.111690794839005;
%%          D = 0.054975871827661;
%%          COOR = [ ...
%%            A          A
%%            (1. -2.*A) A
%%            A          (1.-2.*A)
%%            B          B
%%            (1.-2.*B)  B
%%            B          (1.-2.*B)];
%%          WEIGHT = [C C C D D D];
%%        case 7,
%%          COOR = [...
%%            1./3.                 1./3.
%%            (9.-2.*sqrt(15.))/21. (6.+sqrt(15.))/21.
%%            (6.+sqrt(15.))/21.    (9.-2.*sqrt(15.))/21.
%%            (6.+sqrt(15.))/21.    (6.+sqrt(15.))/21.
%%            (9.+2.*sqrt(15.))/21. (6.-sqrt(15.))/21.
%%            (6.-sqrt(15.))/21.    (9.+2.*sqrt(15.))/21.
%%            (6.-sqrt(15.))/21.    (6.-sqrt(15.))/21.];
%%          WEIGHT = [ 9./80. ...
%%            (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. ...
%%            (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400.];
%%        otherwise
%%          type
%%          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
%%        end
    case 'CUB8',
      X577 = 0.577350269189626;
      switch type
        case 8,
          COOR   = [-X577 -X577 -X577
                     X577 -X577 -X577
                     X577  X577 -X577
                    -X577  X577 -X577
                    -X577 -X577  X577
                     X577 -X577  X577
                     X577  X577  X577
                    -X577  X577  X577];
          WEIGHT = [1. 1. 1. 1. 1. 1. 1. 1.];
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
%    case 'CU20',
    case 'TET4',
      X138 = 0.1381966011250105D0;
      X585 = 0.5854101966249684D0;
      switch type
        case 1,
          COOR   = [0.25 0.25 0.25];
          WEIGHT = [1./6.];
        case 4,
          COOR   = [X138 X138 X138
                    X585 X138 X138
                    X138 X585 X138
                    X138 X138 X585];
          WEIGHT = [1./24. 1./24. 1./24. 1./24.];
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
    case 'PRI6',
      X138 = 0.1381966011250105D0;
      X585 = 0.5854101966249684D0;
      X577 = 0.577350269189626D0;
      switch type
        case 6,
          COOR   = [1./6. 1./6. -X577
                    2./3. 1./6. -X577
                    1./6. 2./3. -X577
                    1./6. 1./6. X577
                    2./3. 1./6. X577
                    1./6. 2./3. X577];
          WEIGHT = [1./6. 1./6. 1./6. 1./6. 1./6. 1./6.];
        otherwise
          type
          error(['For ' type1 ', ' int2str(type) ' points not implemented'])
        end
%    case 'PR15',
    otherwise,
      type1
      error('Type of element not yet implemented for integration at nodes')
  end

%
  else
    type
    error('Type not implemented')
  end

  nbnni = nbnod;
  idimr = size(COOR,2);
  PHI    = zeros(nbnni,nbptg);
  DPHI   = zeros(idimr,nbnni,nbptg);
  for ptg = 1:nbptg
    PHI(:,ptg)    = EF_Shape(type1,COOR(ptg,:))'; % attention au '
    DPHI(:,:,ptg) = EF_Dshape(type1,COOR(ptg,:));
  end
  intge1 = struct('PHI',PHI,'DPHI',DPHI,'COOR',COOR,'WEIGHT',WEIGHT);
  intg1{zo1} = intge1;
  clear intge1 COOR maile1;
end
