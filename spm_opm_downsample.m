function [data_out]=spm_opm_downsample(S)
% Filter raw data off of the FIL's OPM scanner and then downsample it via
% decimation.
%
% Note that non-MEG channels (i.e. trigger and flux channels)
% will not be filtered, but just downsampled via decimation.
%
% Note also that changes to the D object are *NOT* saved to disk.
%
% FORMAT D = spm_opm_downsample(S)
%   S               - input structure
% SAMPLING FREQUENCIES
%   S.fs_orig           - Original Sampling frequency (Hz) - Default: 6e3
%   S.fs_new            - New Sampling frequency (Hz) - Default: 2e3
% DATA
%   S.D                 - SPM D object currently loaded in RAM (no default)
%   OR
%   S.filepath          - Path to raw .bin file (no default). Note that you
%                         have to be in the same BIDS directory as this file  
%
% TO DO:
% OPTIONS
%   S.channelwise       - If set to 1, the filter is applied
%                         channel-by-channel. This can be useful if the
%                         dataset which is being downsampled is
%                         particularly large and may not fit into RAM.
%                         Defaults to 0, i.e. the filter is applied all in
%                         one hit. See spm_eeg_downsample.m for inspiration
%
% Example usage:
% S=[];
% S.fs_orig=6e3;
% S.fs_new=2e3;
% S.D=D;
% D_ds=spm_opm_downsample(S)
%
% Filters designed by Tim Tierney and Ryan Timms, 2021.
%
% Requires the MATLAB signal processing toolbox.


SVNrev = '$Rev: 7125 $';
spm('FnBanner', mfilename, SVNrev);

% Check Inputs:
if ~isfield(S, 'fs_orig'),      S.fs_orig=6e3; end
if ~isfield(S, 'fs_new'),       S.fs_new=2e3; end
if ~isfield(S, 'D'),            S.D = []; end
if ~isfield(S, 'filepath'),     S.filepath= []; end
if ~isfield(S, 'channelwise'),  S.channelwise= 0; end

if isempty(S.D) && isempty(S.filepath)
    error('At least the D object or the path to the raw binary file must be specified.')
end

assert(mod(S.fs_orig,S.fs_new)==0,'The new sampling frequency must be an integer factor of the original sampling frequency.')

% Define filter parameters
fs = S.fs_orig;  % Original F sample
fsNew = S.fs_new;   % Target F sample
betaNew = 7;    % Kaiser window. higher values makes filter more gaussian (worse transition but better attenuation and ripple)
transitionStart= .4; % Transition start
transitionEnd= .8;  % Transition end
orderMult=20;   % Increasing sharpens transition but also increases ripple
D=S.D;

% Calculate the integer we need to downsample the data
step = fs/fsNew;

% Get filter order:
L1 = round(orderMult*(fs)/fsNew + 1);
if(mod(L1,2)==0)
    L1=L1+1; % Make sure it is odd
end

% Define transition regions
fc = transitionStart*fsNew/(fs);
fc2 = transitionEnd*fsNew/(fs);

% Define the filter:
h12 = firls( L1-1, [0 fc fc2 1], [1 1 0 0]).*kaiser(L1,betaNew)' ;
h12 = h12/sum(h12); % Normalise the coefficients

% Do the filter:
if isempty(S.D)==0; % If we have a D object
    
    meginds=D.indchantype('MEG'); % Get the MEG channels
    X = D(meginds,:,:); % And their data...
    fprintf('\nFiltering...\n');
    X = (filtfilt(h12,1,X'))'; % Do the filter
    fprintf('Downsampling...\n');
    X = X(:,1:step:end); % Downsample
    
    % Now downsample the non-MEG data
    tmp=1:size(D,1);
    tmp(D.indchantype('MEG'))=[];
    Y = D(tmp,:,:);
    Y = Y(:,1:step:end);
    fprintf('Saving...\n');
    
    % Populate the new D object
    data_out=zeros(size(D,1),size(Y,2),1);
    data_out(tmp,:,:)=Y;
    data_out(meginds,:,:)=X;
    
    D_new = clone(D, ['ds_',D.fname], [D.nchannels length(data_out) D.ntrials]);
    D_new(:,:,:)=data_out;
    D_new = fsample(D_new, fsNew);
    save(D_new);
    data_out=D_new;
    fprintf('Done.\n');
    
else % We're dealing with a .bin file
    
    fname=S.filepath;
    frawBinaries=fopen(fname);
    data = fread(frawBinaries,Inf,'single',0,'b');
    fclose(frawBinaries);
    
    base = strsplit(fname,'_meg');
    chanFile= [base{1},'_channels.tsv'];
    megfile = [base{1},'_meg.json'];
    channels = spm_load(chanFile);
    megjson=spm_load(megfile);
    
    % format data into  nchannels x nsamples matrix
    nc = size(channels.name,1);
    data = reshape(data,nc,numel(data)/nc);

    % identify trigger channels
    trigInds = strmatch('TRIG',channels.type);
    dataInds = setdiff(1:size(data,1),trigInds);
    tmpTrigs = data(trigInds,:);
    tmpData = data(dataInds,:);
    
    % Filter and downsample
    fprintf('\nFiltering...\n');
    rData = (filtfilt(h12,1,tmpData'))'; % Do the filter
    fprintf('Downsampling...\n');
    rData = rData(:,1:step:end); % Downsample
    rTrigs = tmpTrigs(:,1:step:end); % Just downsample the periph channels.
    fprintf('Saving...\n');
    
    % create ouput dataset
    output = zeros(nc,length(rData));
    output(dataInds,:)=rData;
    output(trigInds,:)=rTrigs;
    
    % For now, we will just save a copy of the downsampled data. This may
    % change in the future
    frawBinaries = fopen(['ds_',fname],'w');
    fwrite(frawBinaries,output,'single',0,'b');
    fclose(frawBinaries);
    
    % Update meg.json
    megjson.SamplingFrequency=fsNew;
    spm_save(['ds_',megfile],megjson);
    data_out=[];
    fprintf('Done.\n');

end

end
