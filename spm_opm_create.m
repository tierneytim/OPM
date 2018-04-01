function [D,L] = spm_opm_create(S)
%--------------------------------------------------------------------------
% Create an MEG object with a valid forward model for a 2 or 3 dimensional
% array. If no data is provided an empty dataset is generated with a valid 
% forward model for user provided  sensor positions and orientations in MRI
% world space. If these positions are not provided  but data is provided the 
% function will create a dataset suitable for sensor level anaysis. If 
% neither data nor postional information is provided the function will 
% attempt to place sensors evenly accross the scalp surface of either an 
% individual structural MRI(if provided) or accross the scalp surface of 
% the MNI template brain. One can also provide custom meshes in order to
% target structures not included in the inverse normalised canonical
% template. Lastly lead fields for different forward models can be computed
% for sensors sampled at different densities and offsets from the scalp.
% the examples below should highlight most of this functions uses.
%--------------------------------------------------------------------------
%
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
%   S.data          - 2/3 dimensaional array       - Default: empty dataset 
%   S.pinout        - filepath to pinout file      - Default: no labels
%   S.fname         - filename for  dataset        - Default: 'simOPM.dat'
%   S.fs            - Sampling rate of dataset     - Default: 1000
%   S.scale         - scale factor (to fT)         - Default: 1
%   S.trig          - trigger matrix               - Default: no triggers 
%   S.pinout        - Mapping of labels to S.data  - Default: generate labels 
%   S.sensorsUsed   - organisation of S.data       - Default: all of S.data is considered MEG data 
%   S.sMRI          - Filepath to  MRI file        - Default: uses template
%   S.cortex        - Custom cortical mesh         - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh            - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh      - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh      - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type  - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)       - Default: 1
%   S.fid           - 3 x 3 Fiducial  matrix       - Default: [0 0 0; -1 0 0; 1 0 0]
%   S.wholehead     - whole head coverage flag     - Deafult: 0
%   S.space         - space between sensors(mm)    - Default: 25
%   S.offset        - scalp to sensor distance(mm) - Default: 6.5
%   S.nSamples      - number of samples            - Default: 1000
%   S.lead          - flag to compute lead field   - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Sensor Level Example - data supplied 
% S =[];
% S.data = ones(10,2)*5;
% S.fs =1000;
% D = spm_opm_create(S);
%
% %% Sensor Level Example - data supplied with labels
% S =[];
% S.data = ones(80,2)*5;
% S.fs =1000;
% S.pinout = 'OPMpinout_20171018';
% S.sensorsUsed = 'opm2cast_130218_1';
% D = spm_opm_create(S);
%
% %% Sensor Level Example - data supplied with labels and triggers
% S =[];
% S.data = ones(80,2)*5;
% S.trig = ones(3,2)*1;
% S.fs =1000;
% S.pinout = 'OPMpinout_20171018';
% S.sensorsUsed = 'opm2cast_130218_1';
% D = spm_opm_create(S);
%
% %% Source Level Example - labels, triggers, positions and MRI file
% S =[];
% S.data = ones(80,2)*5;
% S.trig = ones(3,2)*1;
% S.fs =1000;
% S.pinout = 'OPMpinout_20171018';
% S.sensorsUsed = 'opm2cast_130218_1';
% S.pos = 'SB_headcast';
% S.sMRI= 'rsMQ02589-0008-00001-000224-01.nii';
% D = spm_opm_create(S);
%
% %% Simulation Example - no data, but positions and MRI file supplied
% S =[];
% S.pos = 'SB_headcast';
% S.sMRI= 'rsMQ02589-0008-00001-000224-01.nii';
% D = spm_opm_create(S);
%
% %% Simulation Example - Individual MRI, custom spacing between sensors
% S =[];
% S.space = 25;
% S.sMRI= 'rsMQ02589-0008-00001-000224-01.nii';
% D = spm_opm_create(S);
%
% %% Simulation Example - Individual MRI, custom mesh file
% S =[];
% S.space = 25;
% S.cortex='testCustom.gii';
% S.sMRI= 'rsMQ02589-0008-00001-000224-01.nii';
% D = spm_opm_create(S);
%
% %% Simulation Example - Individual MRI, wholehead coverage
% S =[];
% S.space = 15;
% S.wholehead = 1;
% S.sMRI= 'rsMQ02589-0008-00001-000224-01.nii';
% D = spm_opm_create(S);
%
% %% Simulation Example - Template MRI file, custom  sensor spacing
% S =[];
% S.space = 15;
% D = spm_opm_create(S);
%
% %% Simulation Example - Template MRI file, custom mesh
% S =[];
% S.cortex = 'testCustom.gii';
% D = spm_opm_create(S);
%
% %% Simulation Example - Template MRI file, custom offsets (SQUID vs OPM) 
% %results consistent with Iivanainen et al (2017) and  Boto et al (2017) 
%
% S =[];
% S.space = 30;
% S.offset = 1; 
% S.lead = 1
% [D,L1] = spm_opm_create(S);
% sc1 = D.inv{1}.forward.scale;
% OPM1 = mean(max(abs(L1)))/sc1;
% 
% S =[];
% S.space = 30;
% S.offset = 20; 
% S.lead = 1
% [D,L2] = spm_opm_create(S);
% sc2 = D.inv{1}.forward.scale;
% squid1 = mean(max(abs(L2)))/sc2;
% 
% 
% S =[];
% S.space = 30;
% S.offset = 6; 
% S.lead = 1
% [D,L1] = spm_opm_create(S);
% sc1 = D.inv{1}.forward.scale;
% OPM2 = mean(max(abs(L1)))/sc1;
% 
% S =[];
% S.space = 30;
% S.offset = 50; 
% S.lead = 1
% [D,L2] = spm_opm_create(S);
% sc2 = D.inv{1}.forward.scale;
% squid2 = mean(max(abs(L2)))/sc2;
% 
% %Iivanainen et al(2017) = ~2.75, Boto et al(2017) = ~5
% [OPM1/squid1,OPM2/squid2]
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%-Initialise
%--------------------------------------------------------------------------

% sensor level data with(out) forward model
if (~isfield(S, 'sMRI') && isfield(S, 'data'))
    if isfield(S, 'pos')
        error('if S.pos is not empty a structural  MRI is required')
    end
    forward = 0;
else
    forward = 1;
end

% template based simulation
if ~isfield(S, 'sMRI') 
    S.sMRI   = 1;
    template = 1;
else
    template = 0;
end

% labeled data
if (isfield(S,'pinout') && isfield(S,'data'))
    
    % read the pinout and opm2cast file
    pinout = readtable(S.pinout,...
        'ReadVariableNames',false,'Delimiter','\t');
    sensorsUsed = readtable(S.sensorsUsed,...
        'ReadVariableNames',false,'Delimiter','\t');
  
    labeledData = 1;
else
    labeledData = 0;
end


%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'fname'),       S.fname   = 'OPM.dat'; end
if ~isfield(S, 'fs'),          S.fs   = 1000; end
if ~isfield(S, 'nSamples'),    S.nSamples   = 1000; end
if ~isfield(S, 'space'),       S.space  = 25; end
if ~isfield(S, 'offset'),      S.offset  = 6.5; end
if ~isfield(S, 'fid'),         S.fid  = [0 0 0; -1 0 0; 1 0 0]; end
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell'; end
if ~isfield(S, 'meshres'),     S.meshres = 1; end
if ~isfield(S, 'wholehead'),   S.wholehead = 0; end
if ~isfield(S, 'scale'),       S.scale = 1; end
if ~isfield(S, 'data'),        S.data = zeros(1,S.nSamples); end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'lead'),        S.lead = 0; end

%-File Management
%--------------------------------------------------------------------------
[a,b]= fileparts(S.fname);
ma = [a,b,'.mat'];
da = [a,b,'.dat'];

if exist(ma,'file')==2
    delete(ma);
end
if exist(da,'file')==2
    delete(da);
end

%- Select Channels
%--------------------------------------------------------------------------
% labelled data of different types(MEG, REF)
if labeledData
    used = zeros(size(sensorsUsed,1),1);
    for i = 1:length(used)
        used(i) = find(strcmp(sensorsUsed.Var1{i},pinout.Var2));
    end
    matPos= pinout.Var1(used);
    labs = pinout.Var2(used);
    refInd = find(strcmp('REF',sensorsUsed.Var2));
    nMeg = length(used)-length(refInd);
    S.data = S.data(matPos,:);
    megInd= setdiff(1:size(used,1),refInd)';
else
    refInd=[];
    nMeg = size(S.data,1);
    megInd = 1:nMeg;
end
%- Account for triggers
%--------------------------------------------------------------------------
if isfield(S ,'trig')
    binTrig = zeros(size(S.trig));
    S.data = [S.data;binTrig];
    for j =1:size(S.trig,1)
        binTrig(j,:)= S.trig(j,:);
    end
    St = [];
    St.base='TRIG';
    St.n=size(binTrig,1);
    triglabs = createLabels(St);
else
    binTrig = [];
    triglabs={};
end
%- Account for PHYS type
%--------------------------------------------------------------------------
if isfield(S ,'other')
    other = zeros(size(S.other));
    S.data = [S.data;other];
    St = [];
    St.base='PHYS';
    St.n=size(S.other,1);
    physlabs = createLabels(St);
else
    physlabs={};
end
%- Account for different channel types
%--------------------------------------------------------------------------
% initially say its all MEG data
cType = {};
cType{1} = 'MEG';
cType = repmat(cType,size(S.data,1),1);

% then add REF labels if they exist
for i = 1:length(refInd)
    cType{refInd(i),:} = 'REF';
end
% then add TRIG type if they exist
for i = 0:(size(binTrig,1)-1)
    cType{(end-i)} ='TRIG';
    trigInd = find(strcmp('TRIG',cType));
end

% then add PHYS type if they exist
if(isfield(S,'other'))
stPhys=size(cType,1)-size(binTrig,1)-size(S.other,1)+1;
endPhys=size(cType,1)-size(binTrig,1);
for i = stPhys:endPhys
    cType{i,:} ='PHYS';
end
physInd = find(strcmp('PHYS',cType));
end
%-Sensor Level Data
%--------------------------------------------------------------------------
D = spm_opm_convert(S.data,S.fname,S.fs,S.scale);

% initially D is filled with zeros in trigger channel. This prevents
% scaling of the triggers
if isfield(S ,'trig')
    D(trigInd,:,:)= binTrig;
    save(D)
end

if isfield(S ,'other')
    D(physInd,:,:)= S.other;
    save(D)
end
%-Scalp Extraction
%--------------------------------------------------------------------------
if forward
    %initially used inverse normalised meshes
    D = spm_eeg_inv_mesh_ui(D,1,S.sMRI,S.meshres);
    save(D);
    % then fill in custom meshes(if they exist)
    args = [];
    args.D = D;
    args.scalp = S.scalp;
    args.cortex = S.cortex;
    args.iskull = S.iskull;
    args.oskull = S.oskull;
    args.template = template;
    D = spm_opm_customMeshes(args);
    save(D);
end

%- Create the Sensor Array
%--------------------------------------------------------------------------
if forward
    % if user supplies postions and orientations
    if isfield(S, 'pos')
        posOri = readtable(S.pos,'Delimiter','\t');
        
        % if its labelled data subset the postions and orientaitons
        if(labeledData)
            megInd = find(strcmp('MEG',cType));
            slots = {sensorsUsed.Var2{megInd}}';
            posSlots = posOri.Var7;
            slotind = zeros(length(slots),1);
            for j = 1:length(slots)
                slotind(j)= find(strcmp(slots{j},posSlots));
            end
            px = posOri.Var1(slotind);
            py = posOri.Var2(slotind);
            pz = posOri.Var3(slotind);
            pos = [px,py,pz];
            ox = posOri.Var4(slotind);
            oy = posOri.Var5(slotind);
            oz = posOri.Var6(slotind);
            ori = [ox,oy,oz];
            nSensors = size(pos,1);
            
            % if not labeled data then don't subset
        else
            pos = [posOri.Var1,posOri.Var2,posOri.Var3];
            ori = [posOri.Var4,posOri.Var5,posOri.Var6];
            nSensors = length(pos);
            D = clone(D,S.fname,[nSensors,S.nSamples,1],1);
            megInd = 1:nSensors;
        end
        
        % if no postions and orientations provided then create them
    else
        args = [];
        args.D =D;
        args.offset = S.offset;
        args.space = S.space;
        args.wholehead = S.wholehead;
        [pos,ori] = spm_opm_createSensorArray(args);
        nSensors = length(pos);
        D = clone(D,S.fname,[nSensors,S.nSamples,1],1);
        megInd = 1:nSensors;
    end
end

%-units, labels, and types
%--------------------------------------------------------------------------
nSensors = size(D,1);

% update channel type (sensor number changes when simulating)
if length(cType)<nSensors
    cType = repmat(cType,nSensors,1);
end

% if we don't have labelled data then create labels
if ~labeledData
    St= [];
    St.base ='Sensor';
    St.n = nSensors;
    labs = createLabels(St);
end

% now we know how many sensors we have we can set units labels and types
labs = [labs;physlabs';triglabs'];
D = units(D,[megInd;refInd],'fT');
D = chantype(D,1:nSensors,cType);
D = chanlabels(D,1:size(D,1),labs');
save(D);

%-Place Sensors  in object
%--------------------------------------------------------------------------
if forward
    grad= [];
    grad.label = {labs{megInd}}';
    grad.coilpos = pos;
    grad.coilori = ori;
    grad.tra = eye(numel(grad.label));
    grad.chanunit = repmat({'T'}, numel(grad.label), 1);
    grad.chantype= 'MEG';
    grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
    D = sensors(D, 'MEG', grad);
    save(D);
end
%- fiducial settings
%--------------------------------------------------------------------------
% because we're simulating we already know where things are so no need to
% supply realistic values unless required for compatibility with other
% code. If using a template fiducials are adjusted to be the MNI
% fiducials. MUST UPDATE THIS TO PUT IN VAGUELY SENSIBLE FIDUCIALS BY
% TRANSFORMING MNI FIDUCIALS INTO THE INDIVIDUAL SPACE...

if forward
    fid.fid.label = {'nas', 'lpa', 'rpa'}';
    if(template)
        fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
    else
        fid.fid.pnt = S.fid;
    end
    D = fiducials(D, fid);
    save(D);
end
%- Coregister
%--------------------------------------------------------------------------
% We already have the meshes and sensors in same space but we need to add
% the datareg field so SPM won't get confused at the forward model stage
if forward
    f = fiducials(D);
    f.pnt =zeros(0,3);
    M = f;
    M.pnt = D.inv{1}.mesh.fid.pnt;
    D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
end
%- Foward  model specification
%--------------------------------------------------------------------------
if forward
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    nverts = length(D.inv{1}.forward.mesh.vert);
    if(S.lead)
    [L,D] = spm_eeg_lgainmat(D,1:nverts);
    end
    spm_eeg_inv_checkforward(D,1,1)
end
save(D)
end