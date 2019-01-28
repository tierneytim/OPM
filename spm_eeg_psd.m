function [po,freq] = spm_eeg_psd(S)
% Denoise OPM data
% FORMAT D = spm_opm_synth_gradiometer(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                       - Default: no Default
%   S.triallength   - window size (ms)                      - Default: 1000
%   S.bc            - boolean to dc correct                 - Default: 0
%   S.channels      - channels to estimate PSD from         - Default: 'ALL'
%   S.plot          - boolean to plot or nor                - Default: 0
% Output:
%   psd             - power spectral density
%   f               - frequencies psd is sampled at
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id$

%-epoch dataset
%--------------------------------------------------------------------------
args =[];
args.D= S.D;
args.trialength=S.trialength;
args.bc=S.dc;
eD = spm_eeg_epochs(S);

%- set window
%--------------------------------------------------------------------------
fs = eD.fsample();
N =size(eD,2);
Nf= ceil((N+1)/2);
nepochs=size(eD,3);
pow = zeros(Nf,size(eD,1),nepochs);
wind  = window(@hanning,size(eD,2));
wind = repmat(wind,1,size(eD,1));

%- create PSD
%--------------------------------------------------------------------------

for j = 1:nepochs
    Btemp=eD(:,:,j)';
    Btemp = Btemp.*wind;
    mu=mean(Btemp);
    zf = bsxfun(@minus,Btemp,mu);
    nchan = size(zf,2);
    if(S.dc)
        fzf = zf;
    else
        fzf=Btemp;
    end
    
    N= length(fzf);
    xdft = fft(fzf);
    xdft=xdft(1:floor(N/2+1),:);
    psdx = abs(xdft)./sqrt(N*fs);
    freq = 0:fs/size(fzf,1):fs/2;
    odd=mod(size(fzf,1),2)==1;
    if(odd)
        psdx(2:end) = 2*psdx(2:end);
    else
        psdx(2:end-1) = 2*psdx(2:end-1);
    end
    pow(:,:,j) =psdx;
end

%- plot
%--------------------------------------------------------------------------
labs = S.channels;
chans = [indchannel(eD,labs) indchantype(eD,labs)];
labs= chanlabels(eD,chans);
po = median(pow(:,chans,:),3);

if(S.plot)
    figure()
   
    p1= semilogy(freq,po,'LineWidth',2);
    hold on
    p2 =plot(0:round(freq(end)),ones(1,round(freq(end))+1)*15,'--k');
    p2.LineWidth=2;
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (fT/rtHz)')
    grid on
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.TickLength = [0.02 0.02];
    fig= gcf;
    fig.Color=[1,1,1];
    % title(meg.Instructions)
    leg = labs;
    leg{end+1}= '15fT';
    legend(leg);
    
end
end
