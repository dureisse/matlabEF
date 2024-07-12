function [chamno3] = FusChamno(chamno1,chamno2);
%  Fusion of 2 element-based fields defined at nodes, into one
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 03 / 2007
%
% Inputs
%   chamno1		First field
%   chamno2		Second field
% Outputs
%   chamno3		Resulting field
%
% The two fields should have the same number of zones

nbzo1 = length(chamno1);
if ~nbzo1
  chamno3 = chamno2;
  return
end

nbzo2 = length(chamno2);
if ~nbzo2
  chamno3 = chamno1;
  return
end

if (nbzo1 ~= nbzo2)
  nbzo1
  nbzo2
  error('Different numbers of zones')
end
nbzo3 = max(nbzo1,nbzo2);

clear chamno3;
for zo1 = 1:nbzo3
  chamnoe1 = chamno1{zo1};
  nbcomp1 = length(chamnoe1);
  chamnoe2 = chamno2{zo1};
  nbcomp2 = length(chamnoe2);

  clear chamnoe3;
  icomp3 = 0;
  for icomp1 = 1:nbcomp1
    icomp3 = icomp3 + 1;
    chamnoe3{icomp3} = chamnoe1{icomp1}; 
  end
  for icomp2 = 1:nbcomp2
    icomp3 = icomp3 + 1;
    chamnoe3{icomp3} = chamnoe2{icomp2}; 
  end
  chamno3{zo1} = chamnoe3;
  clear chamnoe3;

end
