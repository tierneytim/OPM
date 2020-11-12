function [D,L] = spm_opm_create(S)
% Read magnetometer data and optionally set up forward model
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.data          - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.channels      - channels.tsv file                        - Default: REQUIRED
%   S.fs            - Sampling frequency (Hz)                  - Default: REQUIRED if S.meg is empty
%   S.meg           - meg.json file                            - Default: REQUIRED if S.fs is empty
%   S.precision     - 'single' or 'double'                     - Default: 'single'
% SOURCE LEVEL INFO
%   S.coordsystem   - coordsystem.json file                    - Default:
%   S.positions     - positions.tsv file                       - Default:
%   S.sMRI          - Filepath to  MRI file                    - Default: uses template
%   S.cortex        - Custom cortical mesh                     - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh                        - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh                  - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh                  - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type              - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)                   - Default: 1
%   S.lead          - flag to compute lead field   - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_create.m 7778 2020-02-05 13:52:28Z tim $
spm('FnBanner', mfilename);

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell'; end
if ~isfield(S, 'meshres'),     S.meshres = 1; end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'fname'),       S.fname = 'sim_opm'; end
if ~isfield(S, 'precision'),   S.precision = 'single'; end
if ~isfield(S, 'lead'),        S.lead = 0; end


%- Read Binary File
%----------------------------------------------------------------------
try % to read data
    [direc, dataFile] = fileparts(S.data);
    dat = fopen(S.data);
    S.data = fread(dat,Inf,S.precision,0,'b');
    fclose(dat);
    binData=1;
catch % if not readable check if it is numeric
    if ~isa(S.data,'numeric') % if not numeric throw error
        error('A valid dataest or file was not supplied')
    end
    binData=0;
    direc = pwd();
    dataFile=S.fname;
    if ~isfield(S, 'channels')
       error('A channels.tsv file must be supplied'); 
    end

end
%- identify potential BIDS Files
%----------------------------------------------------------------------
base = strsplit(dataFile,'meg');
chanFile= fullfile(direc,[base{1},'channels.tsv']);
megFile= fullfile(direc,[base{1},'meg.json']);
posFile= spm_select('FPList',direc,[base{1},'positions.tsv']);
coordFile= fullfile(direc,[base{1},'coordsystem.json']);

%- Check for channel Info
%--------------------------------------------------------------------------
try % to load a channels file
    channels = spm_load(S.channels);
catch
    try % to load a BIDS channel file
        channels = spm_load(chanFile);
    catch
        try  % use channel struct if supplied
            channels = S.channels;
        catch % create channel struct
            warning('A valid channels.tsv file was not found. Setting default type as MEG and unit as fT');
            args=[];
            args.base='Chan';
            args.n= size(S.data,1);
            labs= spm_create_labels(args);
            channels = [];
            channels.name=labs;
            channels.type=repmat({'MEG'},size(S.data,1),1);
            channels.units= repmat({'fT'},size(S.data,1),1);
        end
    end
end

%- reformat data according to channel info
%--------------------------------------------------------------------------
nc = size(channels.name,1);

if binData
    S.data = reshape(S.data,nc,numel(S.data)/nc);
elseif nc~=size(S.data,1)
    error('numer of channels in S.data different to S.channels')
else
    S.data =S.data;
end

%- Check for MEG Info
%--------------------------------------------------------------------------
try % to load a meg file
    meg = spm_load(S.meg);
catch
    try % to load a BIDS meg file
        meg = spm_load(megFile);
        
    catch
        try % to use meg struct
            meg = S.meg;
        catch
            try % to use S.fs argument to get sampling frequency
                meg =[];
                meg.SamplingFrequency=S.fs;
            catch
                error('A valid meg.json file is required if S.fs is empty');
            end
        end
    end
end

%- Position File check
%----------------------------------------------------------------------
try % to load a channels file
    posOri = spm_load(S.positions);
    positions =1;
catch
    try % to load a BIDS channel file
        posOri = spm_load(posFile);
        positions =1;
    catch
        try % to assign a BIDS struct of positions
            if (ismatrix(S.positions))
                
            end
        catch
            warning('No position information found')
            positions=0;
        end
    end
end
%- Forward model Check
%----------------------------------------------------------------------
subjectSource = positions & isfield(S,'sMRI');
if subjectSource
    forward =1;
    template =0;
    % elseif positions
    %     forward =1;
    %     template =1;
    %     S.sMRI=1;
    %     error('template methods not implemented yet') % need to get back to this
else
    forward =0;
    template =0;
end


%- Create SPM object of simulated or real data
%--------------------------------------------------------------------------
Dtemp = meeg(size(S.data,1),size(S.data,2),size(S.data,3));
Dtemp = fsample(Dtemp,meg.SamplingFrequency);
Dtemp = fname(Dtemp,[dataFile,'.mat']);
Dtemp = path(Dtemp,direc);
Dtemp = chanlabels(Dtemp,1:size(Dtemp,1),channels.name);
Dtemp = units(Dtemp,1:size(Dtemp,1),channels.units);
Dtemp = chantype(Dtemp,1:size(Dtemp,1),channels.type);

%- Overwrite and Save
%--------------------------------------------------------------------------
ma = fullfile(direc,[dataFile,'.mat']);
da = fullfile(direc,[dataFile,'.dat']);

ae = exist(fname(Dtemp),'file')==2;
if(ae)
    delete(ma);
    delete(da);
end
Dtemp.save();
% create data file and insert data
D= blank(Dtemp,[dataFile,'.dat']);
dim=size(D);
D(1:dim(1),1:dim(2),1:dim(3)) = S.data;
D.save();

%- Create Meshes
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
    D = opm_customMeshes(args);
    save(D);
end


%-Place Sensors  in object
%--------------------------------------------------------------------------
if positions
    pos = [posOri.Px,posOri.Py,posOri.Pz];
    ori = [posOri.Ox,posOri.Oy,posOri.Oz];
    cl = posOri.name;
        
    grad= [];
    grad.label = cl;
    grad.coilpos = pos;
    grad.coilori = ori;
    grad.tra = eye(numel(grad.label));
    grad.chanunit = repmat({'T'}, numel(grad.label), 1);
    grad.chantype= 'MEG';
    grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
    D = sensors(D, 'MEG', grad);
    save(D);

%- 2D view based on mean orientation of sensors
%--------------------------------------------------------------------------

    n1=mean(grad.coilori); n1= n1./sqrt(dot(n1,n1));
    t1=cross(n1,[0 0 1]);
    t2=cross(t1,n1);
    pos2d =zeros(size(grad.coilpos,1),2);
    for i=1:size(grad.coilpos,1)
        pos2d(i,1)=dot(grad.coilpos(i,:),t1);
        pos2d(i,2)=dot(grad.coilpos(i,:),t2);
    end
    
    nMEG = length(indchantype(D,'MEG'));
    if nMEG~=size(pos2d,1)
        m1 = '2D positions could not be set as there are ';
        m2 =num2str(nMEG);
        m3 = ' channels but only ';
        m4 = num2str(size(pos2d,1));
        m5 =  ' channels with position information.';
        message = [m1,m2,m3,m4,m5];
        warning(message);
    else
        args=[];
        args.D=D;
        args.xy= pos2d';
        args.label=grad.label;
        args.task='setcoor2d';
        D=spm_eeg_prep(args);
        D.save;
    end
end


%- fiducial settings
%--------------------------------------------------------------------------
if subjectSource
    miMat = zeros(3,3);
    fiMat = zeros(3,3);
    fid=[];
    try % to read the coordsystem.json
        coord = spm_load(S.coordsystem);
        fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
        fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
        fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
        miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
        miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
        miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
        fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
        fid.fid.pnt =fiMat;
        fid.pos= []; % headshape field that is left blank (GRB)
        M = fid;
        M.fid.pnt=miMat;
        M.pnt = D.inv{1}.mesh.fid.pnt;
    catch
        try % to find the BIDS coordsystem.json
            coord = spm_load(coordFile);
            fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
            fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
            fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
            miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
            miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
            miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
            fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
            fid.fid.pnt =fiMat;
            fid.pos= []; % headshape field that is left blank (GRB)
            M = fid;
            M.fid.pnt=miMat;
            M.pnt = D.inv{1}.mesh.fid.pnt;
        catch % DEFAULT: transform between fiducials and anatomy is identity
            fid.fid.label = {'nas', 'lpa', 'rpa'}';
            fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:); % George O'Neill
            fid.pos= [];
            M = fid;
            M.pnt = D.inv{1}.mesh.fid.pnt; % George O'Neill
        end
    end
end

if(template) %make
    fid.fid.label = {'nas', 'lpa', 'rpa'}';
    fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
    fid.pos= []; % headshape field that is left blank (GRB)
    M = fid;
    M.pnt = D.inv{1}.mesh.fid.pnt;
end

%- Coregistration
%--------------------------------------------------------------------------
if(subjectSource||template)
    D = fiducials(D, fid);
    save(D);
    f=fiducials(D);
    f.pnt =zeros(0,3);
    D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
end

%- Foward  model specification
%--------------------------------------------------------------------------
if forward
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    spm_eeg_inv_checkforward(D,1,1);
end
save(D);

%- Foward  model specification
%--------------------------------------------------------------------------

%- Create lead fields
%--------------------------------------------------------------------------
D.inv{1}.forward.voltype = S.voltype;
D = spm_eeg_inv_forward(D);
nverts = length(D.inv{1}.forward.mesh.vert);
if(S.lead)
    [L,D] = spm_eeg_lgainmat(D,1:nverts);
end
spm_eeg_inv_checkforward(D,1,1);
save(D);
fprintf('%-40s: %30s\n','Completed',spm('time'));


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Custom Meshes                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = opm_customMeshes(S)
% wrapper for adding custom meshes to MEG object
% FORMAT D = spm_opm_customMeshes(S)
%   S               - input structure
% Fields of S:
%   S.D             - Valid MEG object         - Default:
%   S.cortex        - Cortical mesh file       - Default: Use inverse normalised cortical mesh
%   S.scalp         - Scalp mesh file          - Default: Use inverse normalised scalp mesh
%   S.oskull        - Outer skull mesh file    - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Inner skull mesh file    - Default: Use inverse normalised inner skull mesh
%   S.template      - is mesh in MNI space?    - Default: 0
% Output:
%  D           - MEEG object
%--------------------------------------------------------------------------

%- Default values & argument check
%--------------------------------------------------------------------------
if ~isfield(S, 'D'),           error('MEG object needs to be supplied'); end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'template'),    S.template = 0; end

D = S.D;
if ~isfield(D.inv{1}.mesh,'sMRI')
    error('MEG object needs to be contain inverse normalised meshes already')
end

%- add custom scalp and skull meshes if supplied
%--------------------------------------------------------------------------
if ~isempty(S.scalp)
    D.inv{1}.mesh.tess_scalp = S.scalp;
end

if ~isempty(S.oskull)
    D.inv{1}.mesh.tess_oskull = S.oskull;
end

if ~isempty(S.iskull)
    D.inv{1}.mesh.tess_iskull = S.iskull;
end

%- add custom cortex and replace MNI cortex with warped cortex
%--------------------------------------------------------------------------
if ~isempty(S.cortex)
    D.inv{1}.mesh.tess_ctx = S.cortex;
    if(S.template)
        D.inv{1}.mesh.tess_mni = S.cortex;
    else
        defs.comp{1}.inv.comp{1}.def = {D.inv{1}.mesh.def};
        defs.comp{1}.inv.space = {D.inv{1}.mesh.sMRI};
        defs.out{1}.surf.surface = {D.inv{1}.mesh.tess_ctx};
        defs.out{1}.surf.savedir.savesrc = 1;
        out = spm_deformations(defs);
        D.inv{1}.mesh.tess_mni  = export(gifti(out.surf{1}), 'spm');
    end
end
save(D);

end
