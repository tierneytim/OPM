function [opti] = spm_opm_read_optitrack(S)
% FORMAT [opti] = spm_opm_read_lvm(S)
%   S               - input structure
% Optional fields of S:
%   S.fname         - filepath to csv file             -Default: no Default                      
%   S.headerlength  - integer specifying how many      -Default: 8
%                     lines of file are header
%   S.avPos         - flag to average pos              -Deafult: 1
%   S.optiTrigger   - optiTrack Trigger                -Deafult: no default
%   S.cleanupThresh - threshold to clean optitrack     -Default: don't clean     
%
% Output: opto - output Structure
%  Fields of opto:
%   opti.motion     - motion data
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'fname'),         error('fname needs to be provided');end
if ~isfield(S, 'headerlength'),  S.headerlength = 8; end
if ~isfield(S, 'avPos'),         S.avPos = 1; end
if ~isfield(S, 'cleanupThresh'), S.avPos = []; end

if  isfield(S, 'optiTrigger')  
    interpolate = 1; 
else 
    interpolate = 0;
end


%-find start of data
%--------------------------------------------------------------------------

data = dlmread(S.fname, ',',8,0);
time= data(:,2);
fs = 1/mean(diff(time));
rot=data(:,3:5);
posi=data(:,6:end);
nposChans=size(posi,2);

if(S.avPos)
    x= mean(posi(:,1:3:(nposChans-2)),2);
    y= mean(posi(:,2:3:(nposChans-1)),2);
    z= mean(posi(:,3:3:(nposChans)),2);
    pos=[x,y,z];
else
    pos=posi;
end
motion=[pos,rot];
%- Interpolate to other space
%--------------------------------------------------------------------------

if(interpolate)
x=time;                           % time in seconds(optiTrack)
y=motion;                         % motion 

trigger= S.optiTrigger;
dtrig= diff(trigger);
beginSamp= find(dtrig==-1)+1;     % begining of optiTrack in other space

msg= 'Time varaible for trigger is mandatory if interpolating';
if  ~isfield(S, 'triggerTime'),  error(msg); end

xq=S.triggerTime(beginSamp:end);  % query points in trigger space(seconds)
xq=xq-xq(1)+x(1);                 % same time zero for opti and trigger
xq(xq>max(x)) = max(x);           % Do not extrapolate outside max(x)

n = size(S.triggerTime,1);        % number of output timepoints
p = size(motion,2);               % number of output channels
optiInterp= zeros(n,p);           % empty matrix to hold output

for i =1:p        % loop over channels and interpolate at query points(xq)
optiInterp(beginSamp:end,i)= spline(x,y(:,i),xq);
end
                  % fill samples before begin with constant
fill=repmat(optiInterp(beginSamp,:),beginSamp-1,1);
motion = optiInterp;
motion(1:beginSamp-1,:)=fill;
fs = 1/mean(diff(S.triggerTime));
time= S.triggerTime;
end

% Cleanup
%--------------------------------------------------------------------------
% if

%-Output Struct
%--------------------------------------------------------------------------

opti.motion= motion;
opti.fs= fs;
opti.time=time;
end