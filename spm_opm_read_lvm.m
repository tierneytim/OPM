function [lbv] = spm_opm_read_lvm(S)
% FORMAT [lbv] = spm_opm_read_lvm(S)
%   S               - input structure
% Optional fields of S:
%   S.filename      - filepath to LVM file             -Default: no Default                      
%   S.nchannels     - integer specifying how many      -Default: 81
%                     channels
%   S.headerlength  - integer specifying how many      -Default: 23
%                     lines of file are header
%   S.timeind       - integer specifying which         -Default: 1
%                     column is time variable
%   S.triggerInds   - Indices of trigger Channels      -Default: 73:80
%   S.trigThresh    - Value to threshold triggers at   -Default: Auto
%   S.trigAsBinary  - Read Triggers as Binary          -Deafult: 1
%   S.nbits         - number of bits to use            -nTriggers 
% Output: lbv - output Structure
%  Fields of lbv:
%   lbv.B           - MEG data
%   lbv.Time        - Time variable
%   lbv.trigs       - Trigger Channels
%   lbv.pinout      - pinout of lbv file(coming soon)
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%-Set default values
%--------------------------------------------------------------------------
msg= 'filename needs to be provided';
if ~isfield(S, 'filename'),            error(msg); end
if ~isfield(S, 'nchannels'),           S.nchannels = 81; end
if ~isfield(S, 'headerlength'),        S.headerlength = 23; end
if ~isfield(S, 'timeind'),             S.timeind = 1; end
if ~isfield(S, 'decimalTriggerInds'),  S.decimalTriggerInds = 74:81; end
if ~isfield(S, 'binaryTriggerInds'),   S.binaryTriggerInds = []; end


%-Check for zipped files
%--------------------------------------------------------------------------
[~,~,c]= fileparts(S.filename);
zipped = strmatch(c,'.zip');

if(zipped)
cellFile=unzip(S.filename);
    S.filename= cellFile{1};
end
%-find start of data
%--------------------------------------------------------------------------

fid = fopen(S.filename);
for i =1:S.headerlength
    fgetl(fid);
end
data_start_position = ftell(fid);
fseek(fid, data_start_position, 'bof');

%-read and return data
%--------------------------------------------------------------------------
raw_data = fread(fid,'uint8=>char');
data = sscanf(raw_data, '%f',[S.nchannels,Inf])';
fclose(fid);

chans = 1:S.nchannels;
time = data(:,S.timeind);
Bind = setdiff(chans,S.timeind);
B = data(:,Bind);

%-Binarise triggers
%--------------------------------------------------------------------------
decTrigs = data(:,S.decimalTriggerInds);
binTrigs = data(:,S.binaryTriggerInds);

if(isfield(S,'trigThresh'))
    tTrigsDecimal= decTrigs>S.trigThresh;
    tTrigsBinary= binTrigs>S.trigThresh;
else
    
    thresh= repmat(std(decTrigs)*2,size(decTrigs,1),1);
    tTrigsDecimal = decTrigs>thresh;
    
    thresh= repmat(std(binTrigs)*2,size(binTrigs,1),1);
    tTrigsBinary = binTrigs>thresh;
    
end

msg='Triggers are not bitwise consistent. Please confirm trigger accuracy';

    for i = 1:size(tTrigsDecimal,2)
        for j= 1:size(tTrigsDecimal,2)
            di =tTrigsDecimal(:,i)-tTrigsDecimal(:,j);
            prLoc = [find(di==1);find(di==-1)];
            if(all(abs(diff(prLoc))>1))
                display(msg)
                tTrigsDecimal(prLoc,i)=1;
                tTrigsDecimal(prLoc,j)=1;
            end
            
        end
    end

%-convert triggers
%--------------------------------------------------------------------------
S.nbits = length(S.decimalTriggerInds);
    cFactor = repmat([2.^(0:S.nbits-1)],size(tTrigsDecimal,1),1);
    aTrigs = sum(cFactor.*tTrigsDecimal,2);
    uTrigs = unique(aTrigs);
    nTrigs = length(uTrigs)-1;
    
    outTrigs = zeros(length(aTrigs),nTrigs);
    
    for i =1:nTrigs
        outTrigs(:,i)= (aTrigs==uTrigs(i+1))*uTrigs(i+1);
    end

%-Output Struct
%--------------------------------------------------------------------------
lbv=[];
lbv.B= B;
lbv.time= time;
lbv.decimalTrigs= outTrigs;
lbv.binaryTrigs= tTrigsBinary;

end