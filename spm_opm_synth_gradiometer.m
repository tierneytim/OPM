function D = spm_opm_synth_gradiometer(S)
% Denoise OPM data
% FORMAT D = spm_opm_synth_gradiometer(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                       - Default: no Default
%   S.confounds     - n x 1 cell array containing           - Default: REF
%                     channel types used for denoising.
%   S.derivative    - flag to denoise using derivatives     - Default: 1
%   S.gs            - flag to denoise using global signal   - Default: 0
%   S.prefix        - string prefix for output MEEG object  - Default 'd'
%   S.lp            - 1xn matrix containing lowpass cutoff  - Default: no filter applied to references
%   S.hp            - 1xn matrix containing highpass cutoff - Default: no filter applied to references
%   S.filtord       - IIR filter order(butterworth)         - Default: 5
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_synth_gradiometer.m 7414 2018-09-07 11:00:29Z spm $


%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),          error(errorMsg); end
if ~isfield(S, 'confounds'),  S.confounds = {'REF'}; end
if ~isfield(S, 'derivative'), S.derivative = 1; end
if ~isfield(S, 'gs'),         S.gs = 0; end
if ~isfield(S, 'prefix'),     S.prefix = 'd'; end
if ~isfield(S, 'Y'),          S.Y = {'MEG'}; end
if ~isfield(S, 'lp'),         S.lp=0; end
if ~isfield(S, 'hp'),         S.hp=0; end
if ~isfield(S, 'filtord'),    S.filtord=5; end

%- Determine X and Y
%--------------------------------------------------------------------------
filtlp = [];
lpord=[];
hpord=[];
filthp = [];
refInd = [];

for j = 1:length(S.confounds)

% Check if any confounds refer to channel types.
tempReftype = S.D.indchantype(S.confounds{j});
% now check if any confounds refer to channel labels
tempRefChan = S.D.indchannel(S.confounds{j});
tempRefInd = [tempReftype tempRefChan];

nRef = length(tempRefInd);

if nRef>0
    filtlp = [filtlp repmat(S.lp(j),1,nRef)];
    filthp = [filthp repmat(S.hp(j),1,nRef)];
else
    filtlp= filtlp;
    filthp= filthp;
end
refInd= [refInd tempRefInd];
end


% now select channeltypes to denoise
megind=S.D.indchantype(S.Y);
% now check if denoising refers to specific channels
megind=[megind S.D.indchannel(S.Y)];

%- Loop and Regress
%--------------------------------------------------------------------------
megres = zeros(size(S.D(megind,:,:)));
ref = S.D(refInd,:,:);

% need to loop over Ntrials of size winSize
Ntrials = size(S.D,3);
winSize = size(S.D,2);

% loop (inefficient due to continued memory reallocation but ...)
for i =1:Ntrials
     % add a mean column to the reference regressors
     intercept = ones(winSize,1);
     reference = ref(:,:,i)';
    
    
    filt = reference';
    for j =1:size(filt,1)
        
        if filthp(j)>0
            filt(j,:) = ft_preproc_highpassfilter(filt(j,:),S.D.fsample,filthp(j),S.filtord,'but','twopass','reduce');
        end
        %  optional low pass filter for references
        if filtlp(j)>0
            filt(j,:) = ft_preproc_lowpassfilter(filt(j,:),S.D.fsample,filtlp(j),S.filtord,'but','twopass','reduce');
        end
    end
    reference = filt'; 
    
       % optionally add derivatives
     if(S.derivative)
        drefer = diff(reference);
        drefer = [drefer(1,:); drefer];
        reference = [drefer reference];
     end
     % optionally add global signal 
    if(S.gs)
        trial     = S.D(megind,:,i)';
        gsrefer   = mean(trial,2);
        reference = [gsrefer reference];
    end
    % reference = detrend(reference,'constant');
    reference= [reference ones(size(reference,1),1)];
    
    % reference is column major so transpose sensors
    beta = pinv(reference)*S.D(megind,:,i)';
    megres(:,:,i) = (S.D(megind,:,i)'- reference*beta)';
end


%- Return Outputs
%--------------------------------------------------------------------------
res = S.D(:,:,:); % new file will only have selected sensors changed
res(megind,:,:) = megres;

% make sure output has the singleton dimension if matrix supplied
if ((length(size(res)))==2)
    outsize = [size(res) 1];
else
    outsize = size(res);
end


inFile  = fnamedat(S.D);
[a,b]   = fileparts(inFile);
outfile = fullfile(a,[S.prefix b '.dat']);
D = clone(S.D,outfile,outsize);
D(:,:,:) = res;
