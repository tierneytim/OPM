function [pD] = spm_opm_bandpow(S)
% computed a sliding window of power in  particular frquency band 
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                 - Default: no Default
%   S.win           - size of window in ms            - Default: no Default              
%   S.band          - band to compute power in        - Default: no Default              
% Output:
%   pD               - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id$

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),      error('an MEEG object must be supplied.'); end
if ~isfield(S, 'win'),    error('a time window(ms) must be supplied.'); end
if ~isfield(S, 'band'),   error('a frequency band must be supplied'); end
if ~isfield(S, 'chunkSize'),   S.chunkSize=200; end


%-Filter to specific band
%--------------------------------------------------------------------------
args = [];
args.D = S.D;
args.type = 'butterworth';
args.band = 'bandpass';
args.freq = S.band;
args.dir = 'twopass';
args.order = 2;
fD = spm_eeg_filter(args);


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

ssqprev = 0;
muprev = 0;
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

if(begwin(i))<0
    input = [prevpad squeeze(fD(:,begs(i):endwin(i),:))]';
elseif (endwin(i))>size(fD,2)
    input = [squeeze(fD(:,begwin(i):ends(i),:)) postpad]'; 
else
    input = squeeze(fD(:,begwin(i):endwin(i),:))';
end

output = input((w+1):(end-w),:);

ssqprev =  (var(input(1:(2*w+1),:))*(nsamples-1));
muprev =  mean(input(1:(2*w+1),:));
output(1,:) = ssqprev/(nsamples-1);
x0ind = 1;
xnind = x0ind+2*w+1;

for j =1:(size(output,1)-1)
    xo = input(x0ind,:);
    xn = input(xnind,:);
    xnind=xnind+1;
    x0ind=x0ind+1;
    mu = muprev + (xn-xo)/nsamples;
    ssq = ssqprev +(xn-xo).*(xn+xo-mu-muprev);
    output(j+1,:) = ssq/(nsamples-1);
    ssqprev = ssq;
    muprev = mu;
end


pD(:,begs(i):ends(i),:)= sqrt(output');
%non MEG chnnels should be unchanged
pD(otherind,begs(i):ends(i),:)=S.D(otherind,begs(i):ends(i),:);

end
pD.save();

end 
