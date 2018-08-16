%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)

%% Reading labview
S = [];
S.filename= 'QZFM_6.zip';
S.nchannels=81;
S.trigThresh=4;
S.decimalTriggerInds=74:81;
S.binaryTriggerInds=[];
lbv = spm_opm_read_lvm(S);

%% raw data with labels, triggers and sensor positions
S =[];
S.data = lbv.B';
S.trig = lbv.decimalTrigs';
S.fs =1200;
S.scale = 1e6/2.7;
S.pinout= 'OPMpinout_20171018';
S.sensorsUsed='OPM2cast_MedianNerve';
S.pos = 'SEF_coarse';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% filter the data
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [1 80];
S.dir = 'twopass';
S.order = 2;
D = spm_eeg_filter(S);

%% epoch the data
S =[];
S.D=D;
S.timewin=[-100 300];
D= spm_opm_epoch_trigger(S);

%% denoising
S=[];
S.D=D;
S.confounds={'REF'};
S.gs=0;
S.derivative=1;
D = spm_opm_synth_gradiometer(S);
%% Detecting outlier Trials
S=[];
S.D=D;
S.thresh=3;
D = spm_opm_removeOutlierTrials(S);
%% Average
S =[];
S.D=D;
D =spm_eeg_average(S);
%% baseline correct
S=[];
S.D=D;
S.timewin=[-100 -20];
D = spm_eeg_bc(S);
