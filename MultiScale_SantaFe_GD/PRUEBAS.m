clc
clear all
% load('SAVE.mat','y')
q = ones(10,1);
% for i=1:10    
    
    save('SAVE.mat','q','-append') ;
    q2=q*2 ;
    save('SAVE.mat','q2','-append') ;
% end
q
clear all
load('SAVE.mat','q')
q

% delete('/mytests/*.mat' %borra todos los archivos con extension .mat