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
S.precision = 'double';
D = spm_opm_create(S);

%% filter the data
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [2 100];
S.dir = 'twopass';
S.order = 5;
fD = spm_eeg_filter(S);
delete(D)
%% denoising
chans = {'(CL|CM|CO|CP)','(CL|CM|CO|CP)','(CL|CM|CO|CP)'};
lp =  [25, 34, 53];
hp =  [15, 30, 47];

S=[];
 S.D=fD;
 S.lp=lp;
 S.hp=hp;
 S.confounds=chans;
 dD = spm_opm_synth_gradiometer(S);
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
%% plot 
inds = selectchannels(bD,'MEG');
figure()
plot(bD.time(),bD(inds,:,:)')
xlabel('Time(ms)')
ylabel('Field(fT)')
ax = gca; % current axes
ax.FontSize = 13;
ax.TickLength = [0.02 0.02];
fig= gcf;
fig.Color=[1,1,1];
maxmag = round(max(max(bD(:,:,:))));

if(maxmag<970 &&  maxmag>950)
    display('integartion test passed')
else
    display('integartion test failed')
end
