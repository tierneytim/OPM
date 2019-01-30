%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)

%% read data
S =[];
S.data = 'meg.bin';
S.sMRI='T1w.nii';
D = spm_opm_create(S);

%% filter the data
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [1 80];
S.dir = 'twopass';
S.order = 2;
fD = spm_eeg_filter(S);
delete(D)
%% denoising
S=[];
S.D=fD;
S.confounds={'REF'};
dD = spm_opm_synth_gradiometer(S);
delete(fD)
%% epoch the data
S =[];
S.D=dD;
S.timewin=[-100 300];
eD= spm_opm_epoch_trigger(S);
delete(dD)
%% Average
S =[];
S.D=eD;
mD =spm_eeg_average(S);
delete(eD)
%% baseline correct
S=[];
S.D=mD;
S.timewin=[-100 -20];
bD = spm_eeg_bc(S);
delete(mD)

