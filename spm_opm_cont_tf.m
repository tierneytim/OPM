function [pD] = spm_opm_cont_tf(S)
% computed a sliding window of power in  particular frquency band 
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D           - SPM MEEG object                   - Default: no Default
%   S.win         - size of window in ms              - Default: no Default              
%   S.band        - n X 2 vector of band of interest  - Default: no Default              
% Output:
%   pD               - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging
%
% a = [3:99]';
% b = [4:100]';
% 
% S=[];
% S.band=[a,b];
% S.win=700;
% S.D= D;
% S.chunkSize=100;
% S.fsnew=240;    % note downsampling to prevent data size explosion
% pD = spm_opm_bandpow(S);
% Tim Tierney
% $Id$

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),      error('an MEEG object must be supplied.'); end
if ~isfield(S, 'win'),    error('a time window(ms) must be supplied.'); end
if ~isfield(S, 'band'),   error('a frequency band must be supplied'); end
if ~isfield(S, 'chunkSize'),   S.chunkSize=200; end
if ~isfield(S , 'fsnew')
        downsample =0;
else
    downsample=1;
end

if size(S.D,3)>1
    error('epoched data is currently not supported');
end
%- downsample data (so output is  not enormous)
%--------------------------------------------------------------------------
if(downsample)
args = [];
args.D = S.D;
args.fsample_new=S.fsnew;
args.trigNN = 1;
D = spm_eeg_downsample(args);

else
    D= S.D;
end


%- Create output dataset
%--------------------------------------------------------------------------
ind = indchantype(D,'MEG');
otherind = setdiff(1:size(D,1),ind);

fprintf('Creating output dataset\n');
outname = fullfile(path(D),['pow_' fname(D)]);
pD = clone(S.D,outname,[size(D,1), size(S.band,1) size(D,2) 1], 0, 1);
pD.save();


pD = pD.frequencies(':', mean(S.band,2));
pD = timeonset(pD, D.time(1));
pD= fsample(pD, D.fsample);
pD = transformtype(pD, 'TF');

for j = 1:size(S.band,1)
%-Filter to specific band
%--------------------------------------------------------------------------
args = [];
args.D = D;
args.type = 'butterworth';
args.band = 'bandpass';
args.freq = S.band(j,:);
args.dir = 'twopass';
args.order = 2;
fD = spm_eeg_ffilter(args);

%-get number of samples 
%--------------------------------------------------------------------------
nsamples = fD.fsample*S.win/1000;
if mod(nsamples,2)==0
    nsamples = nsamples +1;
end


%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(fD,1))*1e6);
begs=1:chunkSamples:size(fD,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(fD,2))
    ends(end)= size(fD,2);
end

%-initial value for variance calcultion 
%--------------------------------------------------------------------------

w = (nsamples-1)/2;
begwin = begs-w;
endwin = ends+w;
prevpad= zeros(size(fD,1),w);
postpad= zeros(size(fD,1),w);



%-loop over chunks 
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Computing Windowed variance',spm('time'));

for i =1:length(begs)
    
    
display(['completed chunk ' num2str(i) ' of ' num2str(length(begs))]);
if  ((begwin(i))<0 && (endwin(i)>size(fD,2)))
    input = [prevpad squeeze(fD(:,begs(i):ends(i),:)) postpad]';
elseif (endwin(i))>size(fD,2)
    input = [squeeze(fD(:,begwin(i):ends(i),:)) postpad]'; 
elseif (begwin(i))<0
    input = [prevpad squeeze(fD(:,begs(i):endwin(i),:))]'; 
else
    input = squeeze(fD(:,begwin(i):endwin(i),:))';
end

a = input((w+1):(end-w),:);
output = zeros(size(a,2),1,size(a,1));
%non MEG chnnels should be unchanged
hmm = sqrt(movvar(input,w*2+1));
hmm(:,otherind)=input(:,otherind);

output(:,1,:) = hmm((w+1):(end-w),:)';
pD(:,j,begs(i):ends(i),:)= output;


end
pD.save();
end
end 



