function [ind] = findoccur(val,list_val)
% Find occurence of each value in a list of values
%
% ind = findoccur(val,list_val);
% works on the differents cases :
% - val is a string, and list_val a cell of strings, e.g. 
%   val=[{'UX'}] list_val=[{'UX'} {'UY'}]
% - val is a number, and list_val a list of numbers; e.g.
%   val=1 list_val=[ 1 2 3 4 ]
%
% inputs
%   val(nval)		values that occurencies are required
%   list_val(nlval)	list of values
% outputs
%   ind(nval)		occurencies
%
% val(i) = list_val(ind(i)) if ind(i)~=0
% if ind(i)==0, then val(i) is not in list_val


dim1 = length(val);
list_dim1 = length(list_val);

if list_dim1 == 0
  ind = zeros(1,dim1);
else

  ind = [];
%%  if iscell(val)
%%%   on ne cherche l'egalite que dans les cellules
%%    if length(val) ~= 1
%%      val
%%      list_val
%%      error('dans le cas des cellules, 1 seule entree prevue')
%%    end
%%    i = 1;
%%      ind1 = 0;
%%      for j = 1:list_dim1
%%        if iscell(list_val(j))
%%          if length(list_val(j)) == length(val(i))
%%            list_val1 = list_val(j);
%%            val1 = val(i);
%%	    ind1 = j;
%%            for k = 1, length(val1)
%%              if val1(k) ~= list_val1(k)
%%                ind1 = 0;
%%                break
%%              end
%%            end
%%            if ind1
%%              break
%%            end
%%          end
%%        end
%%      end
%%      ind(i) = ind1;
%%  else
  if iscell(list_val)
    for i = 1:dim1
      ind1 = 0;
      for j = 1:list_dim1
        if strcmp(val(i),list_val(j))
	   ind1 = j;
	   break
        end
      end
      ind(i) = ind1;
    end
  else
    for i = 1:dim1
      ind1 = find(list_val == val(i));
      if (length(ind1) == 1)
        ind(i) = ind1;
      elseif (length(ind1) == 0)
        ind(i) = 0;
      else
        error('multiple occurencies')
      end
    end
  end
%%  end

end
