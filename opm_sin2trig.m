function Dout = opm_mic2trig(S)
% Helper function for working with microphone data (not for SPM)
% FORMAT D = opm_mic2trig(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%  Optional fields:
%   S.speakerThreshold   - 
%   S.ISI   - 
%   S.experiment   - 
%   S.mic  - 
%   S.speaker   - 
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging


if ~isfield(S, 'speakerThreshold'),   S.speakerThreshold =.1; end  
if ~isfield(S, 'ISI'),                S.ISI =0.5; end  
if ~isfield(S, 'experiment'),         S.experiment = [-Inf Inf]; end  
if ~isfield(S, 'channel'),            S.channel='Trigger 6 [Z]'; end  
if ~isfield(S, 'plot'),               S.plot=1; end  
if ~isfield(S, 'threshold'),          S.threshold=.02; end  


%- variables
%--------------------------------------------------------------------------
Dout=copy(S.D,[path(S.D),'/t', fname(S.D)]);
speaker = Dout(selectchannels(Dout,S.channel),:,:)';
speakerThreshold = S.threshold;
ISI=S.ISI;
experiment =S.experiment;

%- binary threshold audio data
%--------------------------------------------------------------------------
thresh = abs(speaker)>speakerThreshold;
samps = find(thresh);
sampsEdit = samps;
inds=[];
%- assume next audio trigger is at least ISI*.9 away and get first index
%- above threshold, removing all others
%--------------------------------------------------------------------------
while(~isempty(sampsEdit))
    inds=[inds sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+Dout.fsample*ISI*.9))=[];
end

%- Do same thing on reversed time series to get audio offset
%--------------------------------------------------------------------------
thresh = abs(speaker)>speakerThreshold;
thresh = flipud(thresh);
samps = find(thresh);
sampsEdit = find(thresh);

indsend =[];
while(~isempty(sampsEdit))
    indsend=[indsend sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+Dout.fsample*ISI*.9))=[];
end

%- use audio onset and offset to create trigger
%--------------------------------------------------------------------------
trEnd= fliplr(size(Dout,2)-indsend)';
trchan = zeros(size(thresh,1),1);
for i =1:length(inds)
    trchan((inds(i):trEnd(i)))=1;
end



%- plots
%--------------------------------------------------------------------------
Dout(selectchannels(Dout,S.channel),:,:) = trchan';
Dout.save();

if(S.plot)

 figure()
 plot(Dout.time(),trchan)
 hold on 
  plot(S.D.time(),speaker/max(speaker))
 xlabel('Time (s)')
end

end
