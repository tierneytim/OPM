%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)
%% read Sensor Level data
S =[];
S.data = 'meg.bin';
D = spm_opm_create(S);

%% Compute PSD

S=[];
S.D =D;
S.trialength=10000;
S.bc=0;
S.channels='MEG';
S.plot=1;

[pxx,f]=spm_eeg_psd(S);
xlim([0,100])
ylim([1,1e6])
delete(D)