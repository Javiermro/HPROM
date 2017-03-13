
% ******************
% MULTISCALE PROGRAM
%     MAIN FILE
% ******************

% clear all

% ********* TESTS **********

% TEST_DATA(1).path_file= '/home/mcaicedo/Desktop/MultiFailure2/Examples/ME';
TEST_DATA(1).path_file= 'C:/Users/jlmroginski/Downloads/MultiScale_SantaFe/MS_Henky';
TEST_DATA(1).file = 'Macro.mfl';
TEST_DATA(1).nLab = 1;
isMICRO(1).MICRO =0; % For macro & multiscale models
 
for iTEST = 1:length(TEST_DATA)
%   try
        analysis(TEST_DATA(iTEST).path_file,TEST_DATA(iTEST).file,TEST_DATA(iTEST).nLab,isMICRO);
%        matlabpool close
%   catch
%       warning('El proceso ha parado su ejecucion por falta de convergencia!');
%       matlabpool close
%   end
%
end

