function [pD] = spm_opm_bandpow(S)
% computed a sliding window of power in  particular frquency band 
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D           - SPM MEEG object                   - Default: no Default
%   S.win         - size of window in ms              - Default: no Default              
%   S.band        - 2 X 1 vector of band of interest  - Default: no Default              
% Output:
%   pD               - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging
%
% S=[];
% S.band=[15,30];
% S.win=700;
% S.D= D;
% pD = spm_opm_bandpow(S);
%
% Tim Tierney
% $Id$

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),      error('an MEEG object must be supplied.'); end
if ~isfield(S, 'win'),    error('a time window(ms) must be supplied.'); end
if ~isfield(S, 'band'),   error('a frequency band must be supplied'); end
if ~isfield(S, 'chunkSize'),   S.chunkSize=200; end

if size(S.D,3)>1
    error('epoched data is currently not supported');
end
%-Filter to specific band
%--------------------------------------------------------------------------
args = [];
args.D = S.D;
args.type = 'butterworth';
args.band = 'bandpass';
args.freq = S.band;
args.dir = 'twopass';
args.order = 2;
fD = spm_eeg_ffilter(args);

%-get number of samples 
%--------------------------------------------------------------------------
nsamples = S.D.fsample*S.win/1000;
if mod(nsamples,2)==0
    nsamples = nsamples +1;
end

%- Create output dataset
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),['pow_' fname(S.D)]);
pD = clone(S.D,outname);
pD.save();
ind = indchantype(fD,'MEG');
otherind = setdiff(1:size(fD,1),ind);

%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
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
%https://nestedsoftware.com/2019/09/26/incremental-average-and-standard-deviation-with-sliding-window-470k.176143.html
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


output = input((w+1):(end-w),:);
%non MEG chnnels should be unchanged
hmm = sqrt(movvar(input,w*2+1));
output(:,ind) = hmm((w+1):(end-w),ind);
pD(:,begs(i):ends(i),:)= output';


end



pD.save();


end 



