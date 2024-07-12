function [rigi1,Bt1,bt1] = RigiCondu7(modl1,matr1,mail1,...
                                      intg1,xcoor1,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 12 / 2003
%
% Matrices elementaires du systeme couple rigidite-conductivite
% Et matrices B elementaire associees
%
% Appelle RigiCompress7 puisque la formulation THERMELAS est tres
% similaire a la formulation POREUX

% On calcule les coefficients materiaux pour la formulation POREUX
% a partir de ceux de la formulation THERMELAS :
% POREUX    THERMELAS
% YOUN   =  YOUN
% NU     =  NU
% IMOB   =  0
% COB    =  ALPH*BULK

% Material coefficients
% """""""""""""""""""""
  [listComp1,listUnit1] = ListCompCham2(matr1,zo1);

% ALPH: dilatability coefficient
  [junk,ia,ib_ALPH] = intersect([{'ALPH'}],listComp1);
  if ~ib_ALPH
    listComp1
    error('No dilatability coefficient')
  end

  [junk,ia,ib] = intersect([{'MOB'}],listComp1);
  [ia,junk] = sort(ia); ib_MOB=ib(junk);
  if xor((length(ib_MOB) == 1),(length(ib_IMOB) == 1))
    isotropic1 = 1;
  else
    listComp1
    error('The correct material coefficient has not been found')
  end

% BULK: bulk modulus, or
% LAMBDA,MU: Lame coefficient, or
% YOUN,NU: Young and Poisson coefficients
% BULK = 3*LAMBDA + 2*MU = YOUN / (1 - 2*NU)
  [junk,ia,ib_BULK] = intersect([{'BULK'}],listComp1);
  [junk,ia,ib_YN] = intersect([{'YOUN'} {'NU'}],listComp1);
  [ia,junk] = sort(ia); ib_YN=ib(junk);
  [junk,ia,ib_LAME] = intersect([{'LAMBDA'} {'MU'}],listComp1);
  [ia,junk] = sort(ia); ib_LAME=ib(junk);
  if ~xor((length(ib_BULK) == 1),(length(ib_YN) == 2), ...
          (length(ib_LAME) == 2))
    listComp1
    error('The correct material coefficient has not been found')
  end
  
clear matr2;
nzo2 = length(matr1);
for zo2 = 1:nzo2
  matr2e = matr1{zo1};
  ncom2 = length(matr2e);
  i = 1;
    siz1 = size(matr2e{i}.XVAL);
  i = ncom2 + 1;
    matr2e{i} = struct('COMP','IMOB', ...
                       'UNIT','', ...
                       'XVAL',zeros(siz1));
  i = ncom2 + 2;
    xval = matr2e{ib_ALPH}.XVAL;
    if ib_BULK
      xval = xval .* matr2e{ib_BULK}.XVAL;
    elseif ib_YN
      xval = xval .* (matr2e{ib_YN(1)}.XVAL ./ ...
                      (1. - 2.*matr2e{ib_YN(2)}.XVAL));
    else
      xval = xval .* (3.*matr2e{ib_LAME(1)}.XVAL .+ ...
                      2.*matr2e{ib_LAME(2)}.XVAL);
    end
    matr2e{i} = struct('COMP','COB', ...
                       'UNIT','', ...
                       'XVAL',xval);
    clear xval;

  matr2{zo2} = matr2e;
  clear matr2e;
end


% Modele
% """"""
modl2 = modl1;
nzo1 = length(modl1);
for zo1 = 1:nzo1
  modl2{zo1}.DDLP = [{'P'}];
  modl2{zo1}.DDLD = [{'FP'}];
  modl2{zo1}.COMP = [{'PX'} {'PY'}];
  modl2{zo1}.COMD = [{'WX'} {'WY'}];
end


% On appelle la routine de la formulation POREUX
% """"""""""""""""""""""""""""""""""""""""""""""
[rigi1,bsig1,b1] = RigiCompress7(modl2,matr2,mail1,intg1,xcoor1,mode1);
