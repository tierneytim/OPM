function D = opm_mic2trig(S)
% Helper function for working with microphone data (not for SPM)
% FORMAT D = opm_mic2trig(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%  Optional fields:
%   S.speakerThreshold   - 
%   S.micThreshold   - 
%   S.ISI   - 
%   S.experiment   - 
%   S.mic  - 
%   S.speaker   - 
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging


if ~isfield(S, 'speakerThreshold'),   S.speakerThreshold =.05; end  
if ~isfield(S, 'micThreshold'),       S.micThreshold =.2; end  
if ~isfield(S, 'ISI'),                S.ISI =6; end  
if ~isfield(S, 'experiment'),         S.experiment = [-Inf Inf]; end  
if ~isfield(S, 'mic'),                S.mic='Trigger 8 [Z]'; end  
if ~isfield(S, 'speaker'),            S.speaker='Trigger 6 [Z]'; end  
if ~isfield(S, 'plot'),               S.plot=0; end  


%- variables
%--------------------------------------------------------------------------
D=S.D;
mic = D(selectchannels(D,S.mic),:,:)';
speaker = D(selectchannels(D,S.speaker),:,:)';
speakerThreshold = S.speakerThreshold;
micThreshold =  S.micThreshold;
ISI=S.ISI;
micWindow = [1,S.ISI/2+1];
experiment =S.experiment;
%- high pas filter data 
%--------------------------------------------------------------------------
[b,a] = butter(4,60/(D.fsample/2),'high');
speaker = filtfilt(b,a,speaker);
mic = filtfilt(b,a,mic);
figure()
plot(D.time(),speaker)
hold on 
figure()
plot(D.time(),mic)

%- binary threshold audio data
%--------------------------------------------------------------------------
thresh = abs(speaker)>speakerThreshold;

inds=[];
figure()
plot(D.time(),abs(speaker))
hold on 
hmm = medfilt1(thresh*1,3);
hmm(hmm>.3)=1;
plot(D.time(),hmm)
samps = find(hmm);
sampsEdit = samps;

%- assume next audio trigger is at least ISI*.9 away and get first index
%- above threshold, removing all others
%--------------------------------------------------------------------------
while(~isempty(sampsEdit))
    inds=[inds sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+D.fsample*ISI*.9))=[];
end

%- Do same thing on reversed time series to get audio offset
%--------------------------------------------------------------------------
thresh = abs(speaker)>speakerThreshold;
thresh = flipud(thresh);
hmm = medfilt1(thresh*1,3);
hmm(hmm>.3)=1;
samps = find(hmm);
sampsEdit = samps;

indsend =[];
while(~isempty(sampsEdit))
    indsend=[indsend sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+D.fsample*ISI*.9))=[];
end

%- use audio onset and offset to create trigger
%--------------------------------------------------------------------------
trEnd= fliplr(size(D,2)-indsend)';
trchan = zeros(size(thresh,1),1);
for i =1:length(inds)
    trchan((inds(i):trEnd(i)))=1;
end

%- look for speech in window after speaker trigger
%--------------------------------------------------------------------------
test =  inds+micWindow(1)*D.fsample;
testend = inds +micWindow(2)*D.fsample;
st=[];

for i = 1:length(test)
   hmm = test(i):testend(i);
   thresh = find(abs(mic(hmm))>micThreshold );
   if ~isempty(thresh) 
    st = [st hmm(thresh(1))];
   end
end

%- create a termporary trigger of mic speecj
%--------------------------------------------------------------------------
trchan2 = zeros(size(trchan));
st = (st);

for i =1:length(st)
    trchan2((st(i):(st(i)+D.fsample*.05)))=1;
end

%- create a termporary trigger of mic speecj
%--------------------------------------------------------------------------
indsend =[];
test = st+1.5*D.fsample;
ends = zeros(length(st),1);

for i = 1:length(st)
    indices = st(i):test(i);
    window = find(abs(mic(indices))>micThreshold);
    ends(i) = indices(window(end));
end

%- create final trigger
%--------------------------------------------------------------------------
trchan2 = zeros(size(trchan));
for i =1:length(st)
    trchan2((st(i):(ends(i))))=1;
end

%- window it
%--------------------------------------------------------------------------
trchan2(D.time()<experiment(1) | D.time()>experiment(2))=0;
react=[];
 for i = 1:length(st)
 h = st(i)-trEnd;
 react = [react min(h(h>0))];
 end
 RTs =react/D.fsample;

%- plots
%--------------------------------------------------------------------------
D(selectchannels(D,S.mic),:,:) = trchan2';
D(selectchannels(D,S.speaker),:,:) = trchan';
D.save();

if(S.plot)
 figure()
 histogram(RTs)
 xlabel('Reaction Times (s)')
 figure()
 plot(D.time(),trchan)
 hold on 
 plot(D.time(),abs(mic))
 plot(D.time(),trchan2)
 xlabel('Time (s)')
 
end

end
