clear;close all;clc;
% Load all data
% E = effective phase (0°); I = ineffective phase (210°)
PCap_E = dlmread('./Effective_phase/PCap.txt');
PCap_I = dlmread('./Ineffective_phase/PCap.txt');
i_tACS_E = dlmread('./Effective_phase/i_sin.txt');
i_tACS_I = dlmread('./Ineffective_phase/i_sin.txt');
ION_E = dlmread('./Effective_phase/ION.txt');
ION_I = dlmread('./Ineffective_phase/ION.txt');
vDCN1_E = dlmread('./Effective_phase/v_DCN1.txt');
vDCN2_E = dlmread('./Effective_phase/v_DCN2.txt');
vDCN3_E = dlmread('./Effective_phase/v_DCN3.txt');
vDCN4_E = dlmread('./Effective_phase/v_DCN4.txt');
vDCN5_E = dlmread('./Effective_phase/v_DCN5.txt');
vDCN1_I = dlmread('./Ineffective_phase/v_DCN1.txt');
vDCN2_I = dlmread('./Ineffective_phase/v_DCN2.txt');
vDCN3_I = dlmread('./Ineffective_phase/v_DCN3.txt');
vDCN4_I = dlmread('./Ineffective_phase/v_DCN4.txt');
vDCN5_I = dlmread('./Ineffective_phase/v_DCN5.txt');
vVim1_E = dlmread('./Effective_phase/v_Vim1.txt');
vVim2_E = dlmread('./Effective_phase/v_Vim2.txt');
vVim3_E = dlmread('./Effective_phase/v_Vim3.txt');
vVim4_E = dlmread('./Effective_phase/v_Vim4.txt');
vVim5_E = dlmread('./Effective_phase/v_Vim5.txt');
vVim1_I = dlmread('./Ineffective_phase/v_Vim1.txt');
vVim2_I = dlmread('./Ineffective_phase/v_Vim2.txt');
vVim3_I = dlmread('./Ineffective_phase/v_Vim3.txt');
vVim4_I = dlmread('./Ineffective_phase/v_Vim4.txt');
vVim5_I = dlmread('./Ineffective_phase/v_Vim5.txt');

% 1) ION raster plots - effective phase
figure(1);
subplot(4,1,1);
hold on;
for i = 1:40
    tmpind = find(ION_E(:,1)==i-1);
    tmpION = ION_E(tmpind,2);
    yyy = diff(tmpION);
    mmm = find(abs(yyy)<1)+1;
    tmpION(mmm)=[];
    plot(tmpION,ones(length(tmpION),1)*i,'k.');
end
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('ION number');
title('ION Raster Plot: Effective Phase');

% 1) ION raster plots - ineffective phase
figure(2);
subplot(4,1,1);
hold on;
for i = 1:40
    tmpind = find(ION_I(:,1)==i-1);
    tmpION = ION_I(tmpind,2);
    yyy = diff(tmpION);
    mmm = find(abs(yyy)<1)+1;
    tmpION(mmm)=[];
    plot(tmpION,ones(length(tmpION),1)*i,'k.');
end
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('ION number');
title('ION Raster Plot: Ineffective Phase');

% 2) PC raster plots - effective phase
figure(1);
subplot(4,1,2);
hold on;
for i = 1:200
    tmpind = find(PCap_E(:,1)==i-1);
    tmpPC = PCap_E(tmpind,2);
    yyy = diff(tmpPC);
    mmm = find(abs(yyy)<1)+1;
    tmpPC(mmm)=[];
    plot(tmpPC,ones(length(tmpPC),1)*i,'k.','MarkerSize',2);
end
plot(0.25:0.25:11500,i_tACS_E*1.25e4+100,'r','LineWidth',2);
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('PC number');
title('PC Raster Plot: Effective Phase');

% 2) PC raster plots - ineffective phase
figure(2);
subplot(4,1,2);
hold on;
for i = 1:200
    tmpind = find(PCap_I(:,1)==i-1);
    tmpPC = PCap_I(tmpind,2);
    yyy = diff(tmpPC);
    mmm = find(abs(yyy)<1)+1;
    tmpPC(mmm)=[];
    plot(tmpPC,ones(length(tmpPC),1)*i,'k.','MarkerSize',2);
end
plot(0.25:0.25:11500,i_tACS_I*1.25e4+100,'r','LineWidth',2);
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('PC number');
title('PC Raster Plot: Ineffective Phase');

% 3) DCN raster plots - effective phase
figure(1);
subplot(4,1,3);
hold on;
for i = 1:5
    [~,tmpDCN] = spike_times(eval(strcat('vDCN',num2str(i),'_E')),-30);
    tmpDCN = tmpDCN*0.25;
    for j = 1:length(tmpDCN)
        plot([tmpDCN(j) tmpDCN(j)],[i-0.4 i+0.4],'k');
    end
end
xlim([1000 2800]);ylim([0.5 5.5]);
xlabel('Time (ms)');
ylabel('DCN number');
title('DCN Raster Plot: Effective Phase');
set(gca,'ytick',1:5);

% 3) DCN raster plots - ineffective phase
figure(2);
subplot(4,1,3);
hold on;
for i = 1:5
    [~,tmpDCN] = spike_times(eval(strcat('vDCN',num2str(i),'_I')),-30);
    tmpDCN = tmpDCN*0.25;
    for j = 1:length(tmpDCN)
        plot([tmpDCN(j) tmpDCN(j)],[i-0.4 i+0.4],'k');
    end
end
xlim([1000 2800]);ylim([0.5 5.5]);
xlabel('Time (ms)');
ylabel('DCN number');
title('DCN Raster Plot: Ineffective Phase');
set(gca,'ytick',1:5);

% 4) Vim spectrogram - effective phase (Wavelet Toolbox required)
tfin = 11500; % Because of downsampling (*20 originally)
tbin = 1;
tseries = tbin:tbin:tfin;

figure(1);
subplot(4,1,4);
spike_mat_all = zeros(5,5000);
tfin = 5000;
% Extract APs from each Vim neuron
for n = 1:5
    [~,tmpVim] = spike_times(eval(strcat('vVim',num2str(i),'_E')),-30);
    tmpVim = tmpVim*0.25;
    yyy = diff(tmpVim);
    mmm = find(abs(yyy)<1)+1;
    tmpVim(mmm)=[];
    
    tbin = 1;
    tseries = tbin:tbin:tfin;
    spike_mat_tmp = zeros(1,length(tseries));
    for j = 1:length(tseries)-1
        for k = 1:length(tmpVim)
            tmpspike = tmpVim(k);
            if tseries(j)<=tmpspike && tseries(j+1)>tmpspike
                spike_mat_tmp(j) = 1;
            end
        end
    end
    spike_mat_all(n,:) = spike_mat_tmp;
    
end
% Temporal summation of APs across the 5 Vim neurons
spike_mat = sum(spike_mat_all,1);
spike_mat(spike_mat>1)=1;
% Compute spectrogram
[wt2,f] = cwt(spike_mat,1e3);
wt2 = wt2(end:-1:1,:);
f = f(end:-1:1);
wt2 = wt2(6:55,:);
f = f(6:55);
pcolor(1:5000,f,abs(wt2));shading interp
ylim([1 14]);
h = colorbar;
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
caxis([0 0.025]);
title('Vim Spectrogram');

% 4) Vim spectrogram - ineffective phase (Wavelet Toolbox required)
figure(2);
subplot(4,1,4);
spike_mat_all = zeros(5,5000);
tfin = 5000;
% Extract APs from each Vim neuron
for n = 1:5
    [~,tmpVim] = spike_times(eval(strcat('vVim',num2str(i),'_I')),-30);
    tmpVim = tmpVim*0.25;
    yyy = diff(tmpVim);
    mmm = find(abs(yyy)<1)+1;
    tmpVim(mmm)=[];
    
    tbin = 1;
    tseries = tbin:tbin:tfin;
    spike_mat_tmp = zeros(1,length(tseries));
    for j = 1:length(tseries)-1
        for k = 1:length(tmpVim)
            tmpspike = tmpVim(k);
            if tseries(j)<=tmpspike && tseries(j+1)>tmpspike
                spike_mat_tmp(j) = 1;
            end
        end
    end
    spike_mat_all(n,:) = spike_mat_tmp;
    
end
% Temporal summation of APs across the 5 Vim neurons
spike_mat = sum(spike_mat_all,1);
spike_mat(spike_mat>1)=1;
% Compute spectrogram
[wt2,f] = cwt(spike_mat,1e3);
wt2 = wt2(end:-1:1,:);
f = f(end:-1:1);
% Extract the 4~12 Hz range
wt2 = wt2(6:55,:);
f = f(6:55);
pcolor(1:5000,f,abs(wt2));shading interp
ylim([1 14]);
h = colorbar;
xlim([1000 2800]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
caxis([0 0.025]);
title('Vim Spectrogram');
