clear;clc;
% Initialize synaptic connections between neurons
createConnections;clear;

ampparam = 0.004; % Stimulation amplitude (in nA)
phaseparam = 10*pi/12; % Stimulation phase

t_stop = 11500; % Simulation period (ms)
t_step = 0.25; % Time step (ms)

%% Set up params for ecHT
% Compute the phase of a reference sinusoidal sequence at each time point
intrinsicPhase = rem(2*pi*7*(t_step:t_step:t_stop)'/1000,2*pi);
intrinsicPhase(intrinsicPhase>pi) = intrinsicPhase(intrinsicPhase>pi) - 2*pi;

% Start iterating through each parameter pairs
currentAmp_all = [zeros(1500/t_step,1);ampparam*ones(10000/t_step,1)];

% Specify paths for storing simulation files
formatSpec = '%s \n';
opt_fileID = fopen('trialname.txt','w');
f_PCap = strcat('recordings_full/PCap.txt');
f_ION = strcat('recordings_full/ION.txt');
f_DCN = strcat('recordings_full/DCN.txt');
f_Vim = strcat('recordings_full/Vim.txt');
f_PCv = strcat('recordings_full/v_PC.txt');
% The instantaneous activity of Vim needs to be analyzed in real-time
% Only the 5/25 Vim TC neurons will be analyzed. See ET_tACS.hoc.
f_Vim1 = strcat('temporary/v_Vim1.txt');
f_Vim2 = strcat('temporary/v_Vim2.txt');
f_Vim3 = strcat('temporary/v_Vim3.txt');
f_Vim4 = strcat('temporary/v_Vim4.txt');
f_Vim5 = strcat('temporary/v_Vim5.txt');
f_sin = strcat('recordings_full/i_sin.txt');
fprintf(opt_fileID,formatSpec,f_PCap);
fprintf(opt_fileID,formatSpec,f_ION);
fprintf(opt_fileID,formatSpec,f_DCN);
fprintf(opt_fileID,formatSpec,f_Vim);
fprintf(opt_fileID,formatSpec,f_PCv);
fprintf(opt_fileID,formatSpec,f_Vim1);
fprintf(opt_fileID,formatSpec,f_Vim2);
fprintf(opt_fileID,formatSpec,f_Vim3);
fprintf(opt_fileID,formatSpec,f_Vim4);
fprintf(opt_fileID,formatSpec,f_Vim5);
fprintf(opt_fileID,formatSpec,f_sin);
fclose(opt_fileID);

formatSpec = '%.8f \n';

opt_fileID = fopen('sinparams.txt','w');
fprintf(opt_fileID,formatSpec,0); % Start with no-stim (till 1500 ms)
fprintf(opt_fileID,formatSpec,phaseparam);
fclose(opt_fileID);

% Update RNG seeds that determine synaptic/membrane noises
formatSpec = '%.0f \n';
opt_fileID = fopen('rngSeeds.txt','w');
fprintf(opt_fileID,formatSpec,ceil(rand(16,1)*300 + 1000));
fclose(opt_fileID);

Vim_v_all = zeros(t_stop/t_step,5);
Vim_name_all = {f_Vim1,f_Vim2,f_Vim3,f_Vim4,f_Vim5};

% Call NEURON from the operating system; save state for the first
% timestep; initialize recording files
[status,cmdout] = system('nrniv -nogui INIT_ET_tACS.hoc');

% Read instantaneous membrane potentials of Vim
formatSpec = '%f';
Vim_fileID1 = fopen(f_Vim1,'r');
Vim_v_all(1,1) = fscanf(Vim_fileID1,formatSpec);
fclose(Vim_fileID1);
Vim_fileID2 = fopen(f_Vim2,'r');
Vim_v_all(1,2) = fscanf(Vim_fileID2,formatSpec);
fclose(Vim_fileID2);
Vim_fileID3 = fopen(f_Vim3,'r');
Vim_v_all(1,3) = fscanf(Vim_fileID3,formatSpec);
fclose(Vim_fileID3);
Vim_fileID4 = fopen(f_Vim4,'r');
Vim_v_all(1,4) = fscanf(Vim_fileID4,formatSpec);
fclose(Vim_fileID4);
Vim_fileID5 = fopen(f_Vim5,'r');
Vim_v_all(1,5) = fscanf(Vim_fileID5,formatSpec);
fclose(Vim_fileID5);

% Start with no-stim until 1500 ms (Tremor triggered at 1000 ms)
for i = 2:1500/t_step % Initialize for 1500 ms
	if mod(i,150/t_step)==1
		disp(i);
	end
	
	formatSpec = '%.0f \n';
	opt_fileID = fopen('rngSeeds.txt','w');
	fprintf(opt_fileID,formatSpec,ceil(rand(16,1)*300 + 1000));
	fclose(opt_fileID);
	
	% Load state -> proceed with a new time step -> save state
	[status,cmdout] = system('nrniv -nogui ET_tACS.hoc');
	
	formatSpec = '%f';
	for n = 1:5
		Vim_fileID = fopen(Vim_name_all{n},'r');
		Vim_v_all(i,n) = fscanf(Vim_fileID,formatSpec);
		fclose(Vim_fileID1);
		% Fix peak membrane potentials at +30 mV for each spike (to
		% compensate for the low temporal resolution)
		if Vim_v_all(i,n)-Vim_v_all(i-1,n) > 30
			Vim_v_all(i,n) = 30;
		end
	end
end

%% Real-Time Phase Adaptation
% Initialize phase estimation-related variables
newphase_all = zeros(1000,1);
newphase = 0;
tseries = 1:1e3;

% Improve the accuracy of phase tracking with a time lag (in ms; optional);
% however, it reduces the responsiveness of phase-locked stimulation
phalag = 150;

% Proceed from 1500 ms
for i = (1500/t_step+1):(t_stop/t_step) % 10 seconds
	if mod(i,150/t_step)==1
		disp(i);
	end
	formatSpec = '%f'; % float-point numbers
	
	for n = 1:5
		Vim_fileID = fopen(Vim_name_all{n},'r');
		Vim_v_all(i,n) = fscanf(Vim_fileID,formatSpec);
		fclose(Vim_fileID1);
		if Vim_v_all(i,n)-Vim_v_all(i-1,n) > 30
			Vim_v_all(i,n) = 30;
		end
	end

	% Real-time update of phase-tracking every 10ms
	if mod(i,10/t_step)==0
		spike_mat_tACS = zeros(length(tseries),5);
		for n = 1:5
			[~,tmpVimspk] = spike_times(Vim_v_all(i-(1000/t_step-1):i,n),-30);
			tmpVimspk = round(tmpVimspk*t_step);
			for j = 1:length(tseries)-1
				for k = 1:length(tmpVimspk)
					tmpspike = tmpVimspk(k);
					if tseries(j)<=tmpspike && tseries(j+1)>tmpspike
						spike_mat_tACS(j,n) = 1;
					end
				end
			end
		end
		% Real-time estimation of instantaneous oscillation phase in
		% the Vim through ecHT
		test_tACS_ec = angle(echt(sum(spike_mat_tACS,2), 6, 10, 1e3));
		% Update new phase
		newphase = test_tACS_ec(end-phalag) + phaseparam - intrinsicPhase(i);
		newphase_all(i/(10/t_step)) = newphase;
	end
	
	formatSpec = '%.8f \n';
	opt_fileID = fopen('sinparams.txt','w');
	fprintf(opt_fileID,formatSpec,currentAmp_all(i)); % Update instantaneous current amplitude
	fprintf(opt_fileID,formatSpec,newphase); % Update instantaneous stimulation phase
	fclose(opt_fileID);
	
	formatSpec = '%.0f \n';
	opt_fileID = fopen('rngSeeds.txt','w');
	fprintf(opt_fileID,formatSpec,ceil(rand(16,1)*300 + 1000));
	fclose(opt_fileID);
	
	% Execute NEURON
	[status,cmdout] = system('nrniv -nogui ET_tACS.hoc');

end

% Save all Vim recordings
f_Vim_full1 = strcat('recordings_full/v_Vim1.txt');
f_Vim_full2 = strcat('recordings_full/v_Vim2.txt');
f_Vim_full3 = strcat('recordings_full/v_Vim3.txt');
f_Vim_full4 = strcat('recordings_full/v_Vim4.txt');
f_Vim_full5 = strcat('recordings_full/v_Vim5.txt');
f_newphase_full = strcat('recordings_full/newphase.txt');

formatSpec = '%f \n';
opt_fileID = fopen(f_Vim_full1,'w');
fprintf(opt_fileID,formatSpec,Vim_v_all(:,1));
fclose(opt_fileID);
opt_fileID = fopen(f_Vim_full2,'w');
fprintf(opt_fileID,formatSpec,Vim_v_all(:,2));
fclose(opt_fileID);
opt_fileID = fopen(f_Vim_full3,'w');
fprintf(opt_fileID,formatSpec,Vim_v_all(:,3));
fclose(opt_fileID);
opt_fileID = fopen(f_Vim_full4,'w');
fprintf(opt_fileID,formatSpec,Vim_v_all(:,4));
fclose(opt_fileID);
opt_fileID = fopen(f_Vim_full5,'w');
fprintf(opt_fileID,formatSpec,Vim_v_all(:,5));
fclose(opt_fileID);
opt_fileID = fopen(f_newphase_full,'w');
fprintf(opt_fileID,formatSpec,newphase_all);
fclose(opt_fileID);

disp('Done');