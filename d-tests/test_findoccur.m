clear all
close all
path(path,'../matlabUtils')

toto = [{'UX'} {'P'} {'UX'} {'T'}];
titi = [{'UX'} {'UY'} {'T'}];
ind1 = findoccur(toto,titi);
if ~all(ind1 == [1 0 1 3])
  error('Probleme')
end

ind1 = findoccur(toto(1),titi);
if ~all(ind1 == [1])
  error('Probleme')
end

ind1 = findoccur(toto([1 4 3]),titi);
if ~all(ind1 == [1 3 1])
  error('Probleme')
end

toto = [{'UX'}];
ind1 = findoccur(toto,titi);
if ~all(ind1 == [1])
  error('Probleme')
end

ind1 = findoccur({'UX'},titi);
if ~all(ind1 == [1])
  error('Probleme')
end

clear toto;
toto{1} = 'UX';
toto{2} = 'P';
toto{3} = 'UX';
toto{4} = 'T';
ind1 = findoccur(toto,titi);
if ~all(ind1 == [1 0 1 3])
  error('Probleme')
end

ind1 = findoccur(toto(1),titi);
if ~all(ind1 == [1])
  error('Probleme')
end

clear toto;
toto{1} = 'UX';
ind1 = findoccur(toto,titi);
if ~all(ind1 == [1])
  error('Probleme')
end

clear toto;
toto{1} = 'UX';
ind1 = findoccur(toto,[]);
if ~all(ind1 == [0])
  error('Probleme')
end
clear toto;
toto = [1 2 3];
ind1 = findoccur(toto,[]);
if ~all(ind1 == [0 0 0])
  error('Probleme')
end



disp('TEST PASSE AVEC SUCCES')
quit
