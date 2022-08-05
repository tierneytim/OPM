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


ssqprev =  (var(input(1:(2*w+1),ind))*(nsamples-1));
muprev =  mean(input(1:(2*w+1),ind));
output(1,ind) = ssqprev/(nsamples-1);
x0ind = 1;
xnind = x0ind+2*w+1;

for j =1:(size(output,1)-1)
    xo = input(x0ind,ind);
    xn = input(xnind,ind);
    xnind=xnind+1;
    x0ind=x0ind+1;
    xd = xn-xo;
    mu = muprev + (xd)/nsamples;
    ssq = ssqprev +(xd).*(xn+xo-mu-muprev);
    output(j+1,ind) = ssq/(nsamples-1);
    ssqprev = ssq;
    muprev = mu;
end


output(:,ind) = sqrt(output(:,ind));
pD(:,begs(i):ends(i),:)= output';

%non MEG chnnels should be unchanged
%pD(otherind,begs(i):ends(i),:)=S.D(otherind,begs(i):ends(i),:);

end


pD.save();


end 
function Dnew = spm_opm_filter(S)
% Filter M/EEG data
% FORMAT D = spm_eeg_filter(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%   S.band    - filterband [low|high|bandpass|stop]
%   S.freq    - cutoff frequency(-ies) [Hz]
%
%  Optional fields:
%   S.type    - filter type [default: 'butterworth']
%                 'butterworth': Butterworth IIR filter
%                 'fir':         FIR filter (using MATLAB fir1 function)
%   S.order   - filter order [default: 5 for Butterworth]
%   S.dir     - filter direction [default: 'twopass']
%                 'onepass':         forward filter only
%                 'onepass-reverse': reverse filter only, i.e. backward in time
%                 'twopass':         zero-phase forward and reverse filter
%   S.prefix  - prefix for the output file [default: 'f']
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_filter.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

if ~isfield(S, 'dir'),    S.dir    = 'onepass';     end
if ~isfield(S, 'chunkSize'),    S.chunkSize    = 200;     end
if ~isfield(S, 'prefix'), S.prefix = 'f';           end
if ~isfield(S, 'order'),  S.order=5; end

if ~isfield(S, 'band') 
    error('A frequency band must be supplied')
end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

%-Check band
%--------------------------------------------------------------------------
switch lower(S.band)
    
    case {'low','high'}
        if numel(S.freq)~=1
            error('Cutoff frequency should be a single number.');
        end
        
        if S.freq < 0 || S.freq > D.fsample/2
            error('Cutoff must be > 0 & < half sample rate.');
        end
        
    case {'bandpass','stop'}
        if S.freq(1) < 0 || S.freq(2) > D.fsample/2 || S.freq(1) > S.freq(2)
            error('Incorrect frequency band specification.');
        end
        
    otherwise
        error('Incorrect filter band.')
end

fprintf('%-40s: %30s\n',...
    ['Filter ' S.band ' (' 'butterworth' ', ' S.dir ')'],...
    ['[' num2str(S.freq) '] Hz']);                                      %-#

%-Filter
%==========================================================================

%-Generate new meeg object with new filenames
Dnew = copy(D, [S.prefix fname(D)]);

%-Determine channels for filtering
Fchannels = D.indchantype('Filtered');

if isempty(Fchannels)
    warning('No channels suitable for filterning found. Please check your channel type specification.');
end
    

%- filter coeffciients (just lowpass for now)
%--------------------------------------------------------------------------

unstable =1;
while (unstable)
    [B, A] = butter(S.order,S.freq/(S.D.fsample/2),S.band);
    unstable = any(abs(roots(A))>=1);
    if(unstable)
        ft_warning('instability detected - reducing the %dth order filter to an %dth order filter', S.order, S.order-1);
        S.order=S.order-1;
        if(S.order<1)
            error('Cannot stabilise filter');
        end
    end
end
   

%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
end

%-Forward direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (Forward)',spm('time'));
zf = zeros(max(length(A),length(B))-1,length(Fchannels));
for i =1:length(begs)
    %display(['completed chunk ' num2str(i) ' of ' num2str(length(begs))]);
    input = squeeze(S.D(:,begs(i):ends(i),:))';
    output = input;
    [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
    Dnew(:,begs(i):ends(i),:)= output';
end
Dnew.save();

%-backward  direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (backward)',spm('time'));
switch S.dir
    case 'twopass'
        zf = zeros(max(length(A),length(B))-1,length(Fchannels));
        ends = flip(ends);
        begs = flip(begs);
        for i =1:length(ends)
     %       display(['completed chunk ' num2str(i) ' of ' num2str(length(begs))]);
            input = squeeze(Dnew(:,begs(i):ends(i),:))';
            input = flip(input);
            output = input;
            [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
            Dnew(:,begs(i):ends(i),:)= flip(output)';
        end
    case 'onepass'
        
    otherwise
       error('supplied direction not supported');
end


Dnew.save();
%     
%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#

end




