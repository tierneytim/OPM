function [tHelm] = spm_affine_opreg(S)
% Read magnetometer data and optionally set up forward model
% FORMAT tHelm = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
%   S.helmetfile     - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headfile     - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.helmetref    - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headhelmetref - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headfid
%   S.headhelmetfid
%   S.affine  
%   S.frontonly
% Output:
%  tHelm       - transformed helmet object
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id$

if ~isfield(S, 'sensreg'),               S.sensreg = 0; end
if ~isfield(S, 'affine'),                S.affine = 1; end

%   S.frontonly
%-read and convert meshes to mm
%--------------------------------------------------------------------------
Native= gifti(S.headfile);

inmetres = sqrt(sum((max(Native.vertices)-min(Native.vertices)).^2))<10;
if inmetres
Native.vertices=Native.vertices*1000;
end


helmetInfo = load(S.helmetfile);
helmet = gifti(helmetInfo.Helmet_info.helmet_surface);

inmetres = sqrt(sum((max(helmet.vertices)-min(helmet.vertices)).^2))<10;
if inmetres
    helmet.vertices=helmet.vertices*1000;
end

%%-  helmet --> head with helmet
%--------------------------------------------------------------------------
helm2headhelm = spm_eeg_inv_rigidreg(S.headhelmetref',S.helmetref');


%-  head with helmet  --> head
%--------------------------------------------------------------------------
headhelm2head = spm_eeg_inv_rigidreg(S.headfid',S.headhelmetfid');

%- templates 
%--------------------------------------------------------------------------
scalp = gifti(fullfile(spm('dir'), 'canonical','scalp_2562.surf.gii'));
brain = gifti(fullfile(spm('dir'), 'canonical','cortex_5124.surf.gii'));

%- head space --> template (6 paramater)
%--------------------------------------------------------------------------
fid_template = spm_eeg_fixpnt(ft_read_headshape(fullfile(spm('dir'), 'EEGtemplates','fiducials.sfp')));
head2templatescalp = spm_eeg_inv_rigidreg(fid_template.fid.pnt(1:3,:)',S.headfid');
scalpTemplate6 = spm_mesh_transform(Native,head2templatescalp);

%- sensors --> template (6 paramater)
%--------------------------------------------------------------------------
sens = helmetInfo.Helmet_info.sens_pos;
inmetres = sqrt(sum((max(sens)-min(sens)).^2))<10;
if inmetres
   sens=sens*1000;
end
senstemplate6 = spm_mesh_transform(sens,head2templatescalp*headhelm2head*helm2headhelm);

%- cut everything below the ears off of scanned mesh
%--------------------------------------------------------------------------
p = double(scalpTemplate6.vertices);
cut = min(fid_template.fid.pnt(1:3,3));
p(p(:,3)<cut,:) =[];

%- cut everything  left and right of scalp mesh
%--------------------------------------------------------------------------
ml = min(scalp.vertices);
mr = max(scalp.vertices);
cut = (p(:,1)<ml(1)) | (p(:,1)>mr(1)) ;
p(cut,:) =[];

%- cut everything  behind the mesh
%--------------------------------------------------------------------------
ml = min(scalp.vertices);
mr = max(scalp.vertices);
cut = (p(:,2)<ml(2)) | (p(:,2)>mr(2)) ;
p(cut,:) =[];


%- affine registration (12 paramater) 

%--------------------------------------------------------------------------
if S.sensreg
hmm= spm_eeg_inv_icp(double(scalp.vertices'),senstemplate6',[],[],[],[],S.affine);
scalpTemplate12 = spm_mesh_transform(scalpTemplate6,hmm);
else
hmm= spm_eeg_inv_icp(double(scalp.vertices'),p',[],[],[],[],S.affine);
scalpTemplate12 = spm_mesh_transform(scalpTemplate6,hmm);
end


%-  sensors to template 
%--------------------------------------------------------------------------
sens2temp= hmm*head2templatescalp*headhelm2head*helm2headhelm;
sens=[];
sens.vertices =  helmetInfo.Helmet_info.sens_pos*1000;
senstemplate = spm_mesh_transform(sens,sens2temp);
helmettemplate = spm_mesh_transform(helmet,sens2temp);


%-  check registration
%--------------------------------------------------------------------------

p1 = [];
p1.faces = scalpTemplate12.faces;
p1.vertices = scalpTemplate12.vertices;
p2 = [];
p2.faces = scalp.faces;
p2.vertices = scalp.vertices;
p3 = [];
p3.vertices = helmettemplate.vertices; 
p3.faces = helmettemplate.faces; 

p4=[];
p4.faces = brain.faces; 
p4.vertice = brain.vertices; 


f = figure()

patch(p1,'FaceColor','red','FaceAlpha',1,'EdgeColor','none')
hold on 
patch(p2,'FaceColor','blue','FaceAlpha',.3,'EdgeColor','none')
patch(p3,'FaceColor','green','FaceAlpha',.15,'EdgeColor','none')
patch(p4,'FaceColor','black','FaceAlpha',1,'EdgeColor','none')

ax = gca();
camlight;
camlight(-80,-10);
hold on 
 sensp = senstemplate.vertices;
plot3(sensp(:,1),sensp(:,2),sensp(:,3),'.y','MarkerSize',20)
%plot3(ps(:,1),ps(:,2),ps(:,3),'.b','MarkerSize',20)

daspect([1,1,1])

%-  save new sensor positions
%--------------------------------------------------------------------------
tHelm = helmetInfo.Helmet_info;
tHelm.sens_pos= senstemplate.vertices;
Helmet_info = tHelm;
[a,b,c] = fileparts(S.helmetfile);
outfile = fullfile(a,['trans_',b,c]);
save(outfile,'Helmet_info');



end