clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

error1 = AVStoGMSH('test_AVStoGMSH_2.inp',1,'test_AVStoGMSH_2');
if error1
  error('TEST FAILED')
end

%% error1 = AVStoGMSH('test_AVStoGMSH_1a.inp',1,'test_AVStoGMSH_1a');
%% if error1
%%    error('TEST FAILED')
%% end

%% error1 = AVStoGMSH('test_AVStoGMSH_1b.inp',1,'test_AVStoGMSH_1b');
%% if error1
%%  error('TEST FAILED')
%% end

disp('TEST PASSE AVEC SUCCES')
disp('verifier les fichiers gmsh')

quit
