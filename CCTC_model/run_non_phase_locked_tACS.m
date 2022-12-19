clear;clc;

createConnections;clear;

% Stimulation amplitude (in nA)
ampparam = 0.005;

% Specify paths for storing simulation files
formatSpec = '%s \n';
opt_fileID = fopen('trialname_non_pl.txt','w');
f_PCap = strcat('recordings_non_pl/PCap.txt');
f_DCN = strcat('recordings_non_pl/v_DCN.txt');
f_Vim = strcat('recordings_non_pl/v_Vim.txt');
f_ION = strcat('recordings_non_pl/v_ION.txt');
f_sin = strcat('recordings_non_pl/i_sin.txt');
fprintf(opt_fileID,formatSpec,f_PCap);
fprintf(opt_fileID,formatSpec,f_DCN);
fprintf(opt_fileID,formatSpec,f_Vim);
fprintf(opt_fileID,formatSpec,f_ION);
fprintf(opt_fileID,formatSpec,f_sin);
fclose(opt_fileID);
	
formatSpec = '%.8f \n';

opt_fileID = fopen('sinparams_non_pl.txt','w');
fprintf(opt_fileID,formatSpec,ampparam);
fclose(opt_fileID);

% Setup random seeds
formatSpec = '%.0f \n';
opt_fileID = fopen('rngSeeds_non_pl.txt','w');
fprintf(opt_fileID,formatSpec,ceil(rand(16,1)*300 + 1000));
fclose(opt_fileID);

% Start NEURON simulation
[status,cmdout] = system('nrniv -nogui ET_tACS_non_pl.hoc');
	
disp('Done');