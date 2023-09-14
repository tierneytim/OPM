function D = opm_pd2trig(S)
% Helper function for working with photodiode data (not for SPM)
% FORMAT D = opm_mic2trig(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%  Optional fields:
%   S.pdThreshold   - 
%   S.ISI   - 
%   S.experiment   - 
%   S.pd  - 
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging


if ~isfield(S, 'pdThreshold'),        S.pdThreshold =3.1; end  
if ~isfield(S, 'ISI'),                S.ISI =.25; end  
if ~isfield(S, 'experiment'),         S.experiment = [-Inf Inf]; end  
if ~isfield(S, 'pd'),                 S.pd = 'Trigger 6 Z'; end 
if ~isfield(S, 'where'),              S.pd = 'BOTH'; end 
if ~isfield(S, 'plot'),               S.plot =0; end  

%- variables
%--------------------------------------------------------------------------
D=S.D;
pd = D(selectchannels(D,S.pd),:,:)';
pdThreshold =  S.pdThreshold;
ISI=S.ISI*1.5;

%- binary threshold PD data
%--------------------------------------------------------------------------
figure()
a = double(diff(pd)>.5);


%thresh = hampel(a,8);
winsamps = S.experiment*S.D.fsample;
window  = zeros(size(a));
window(winsamps(1):winsamps(2))=1;
thresh = a.*window;
plot(D.time(),[0;thresh])
hold on 
plot(D.time(),[pd])
plot(D.time(),D(1,:,:)')
inds=[];
samps = find(thresh);
sampsEdit = samps;

%- assume next audio trigger is at least ISI*.9 away and get first index
%- above threshold, removing all others
%--------------------------------------------------------------------------

if(strcmp(S.where,'RISE')||strcmp(S.where,'BOTH'))
while(~isempty(sampsEdit))
    inds=[inds sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+D.fsample*ISI))=[];
end
end
%- Do same thing on reversed time series to get audio offset

if(strcmp(S.where,'FALL')||strcmp(S.where,'BOTH'))

samps = find(flipud(thresh));
sampsEdit = samps;
indsend =[];
while(~isempty(sampsEdit))
    indsend=[indsend sampsEdit(1)];
    sampsEdit(sampsEdit<(sampsEdit(1)+D.fsample*ISI))=[];
end

trEnd= fliplr(size(D,2)-indsend)';
end

%- use audio onset and offset to create trigger
%--------------------------------------------------------------------------

trchan = zeros(size(thresh,1),1);

if(strcmp(S.where,'BOTH'))
evsamps = [inds trEnd'];
end

if(strcmp(S.where,'RISE'))
evsamps = [inds];
end


if(strcmp(S.where,'FALL'))
evsamps = [trEnd'];
end


for i =1:length(evsamps)
    trchan((evsamps(i):(evsamps(i)+30)))=3.1;
end
trchan= trchan(1:length(thresh));

disp(['Identified ' num2str(length(evsamps)) 'events' ])
%- plots
%--------------------------------------------------------------------------
D(selectchannels(D,S.pd),:,:) = trchan';

if(S.plot)
 figure()
 hold on 
 plot(D.time(),abs(pd),'r')
 plot(D.time(),[0 ;trchan],'g')
 %plot(D.time(),D(1,:,:),'b')
 xlabel('Time (s)')
 legend('PD','trig','trig station')
end
D.save()

end
