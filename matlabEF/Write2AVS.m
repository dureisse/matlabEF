function [] = Write2AVS(fid,size_comps,comp_label,units_label,xval)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002

% Ecriture d'un champ d'un fichier AVS UCD

%   Components
    num_comp = size(size_comps,2);
    line = zeros(1,1);
    line(1,1) = num_comp; % number of components
%    format = '%i'; size1 = 1;
%DD 28/07/2002    format = '%4i'; size1 = 1;
    format = '%5i'; size1 = 1;
    count = fprintf(fid,format,line);
    line = size_comps'; % size of each component
    format = '%5i'; size1 = num_comp;
    count = fprintf(fid,format,line);
    count = fprintf(fid,'\n');
    for icomp=1:num_comp
%      format = '%s'; size1 = 1;
%DD 28/07/2002      format = '%3s'; size1 = 1;
      format = ' %-4s'; size1 = 1;
      line = comp_label{icomp};
      count = fprintf(fid,format,line);
%DD 28/07/2002      format = '%3s'; size1 = 1;
      format = '%s'; size1 = 1;
      line = ',';
      count = fprintf(fid,format,line);
%      format = '%s'; size1 = 1;
      format = '%-4s\n'; size1 = 1;
      line = units_label{icomp};
      count = fprintf(fid,format,line);
    end

%   Data
    [num_n num_data] = size(xval);
    for ino=1:num_n
%      format = '%i'; size1 = 1;
      format = '%6i'; size1 = 1;
      id = ino; % id
      line = zeros(1,1);
      line(1,1) = id;
      count = fprintf(fid,format,line);
%      format = '%e'; size1 = num_data;
      format = '  %12.7E'; size1 = num_data;
      line = xval(id,:)';
      count = fprintf(fid,format,line);
      count = fprintf(fid,'\n');
    end
