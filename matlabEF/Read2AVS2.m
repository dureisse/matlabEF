function [size_comps,comp_label,units_label,xval] = ...
	     Read2AVS2(fid,num_data,num_n,echo)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2002
%  ajout echo

% Lecture d'un champ d'un fichier AVS UCD

%   Components
    format = '%i'; size1 = 1;
    [line,count] = fscanf(fid,format,size1);
    num_comp = line(1,1); % number of components
    format = '%i'; size1 = num_comp;
    [line2,count] = fscanf(fid,format,size1);
    size_comps = line2'; % size of each component
    line1 = fgetl(fid); % up to the end of line
if echo
  if size(line1,1)
    disp([int2str(line') ' ' int2str(line2) ' ' line1])
  end
end
    for icomp=1:num_comp
      line1 = fgetl(fid);
if echo
disp(line1)
end
      oldindex = 1;
%DD 29/07/2002      format = '%s'; size1 = 1;
      format = ' %[^, ]'; size1 = 1;
      [line,count,errmsg,nextindex] = ...
	 sscanf(line1(1,oldindex:end),format,size1);
      if errmsg
        error1 = errmsg;
        return;
      end
      comp_label{icomp} = line(1,:); % component label
      oldindex = oldindex + nextindex;
%DD 29/07/2002      format = '%s'; size1 = 1;
      format = ' %s '; size1 = 1;
      [line,count,errmsg,nextindex] = ...
	 sscanf(line1(1,oldindex:end),format,size1); % ,
      if errmsg
        error1 = errmsg;
        return;
      end
      oldindex = oldindex + nextindex;
%DD 29/07/2002      format = '%s'; size1 = 1;
      format = ' %[^, ]'; size1 = 1;
      [line,count,errmsg,nextindex] = ...
	 sscanf(line1(1,oldindex:end),format,size1);
      if errmsg
        error1 = errmsg;
        return;
      end
      if count
        units_label{icomp} = line(1,:); % unit label
      else
        units_label{icomp} = [''];
      end
    end

%   Data
    xval = zeros(0,num_data);
    for ino=1:num_n
      format = '%i'; size1 = 1;
      [line,count] = fscanf(fid,format,size1);
      id          = line(1,1); % id
      format = '%e'; size1 = num_data;
      [line2,count] = fscanf(fid,format,size1);
if echo
disp([int2str(line') ' ' num2str(line2')])
end
      xval(id,:) = line2';
    end
