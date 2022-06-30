%% Housekeeping

clearvars

% if isempty(which('spm'))
%     addpath(spmpath)
%     spm('defaults', 'eeg')
% end

testdir = fileparts(which('testCerca'));
if isempty(testdir)
    error('Please move to directory where testCerca.m exists!')
end

addpath(fullfile(testdir,'..'))
datadir = fullfile(testdir,'..','testData','cerca');

cd(datadir)

%% Unzip cMEG data 

zipped = fullfile(fullfile(datadir,'20220623_095714.tar.gz'));
tmptar = [tempname '.tar'];
tmpgz = [tmptar '.gz'];

copyfile(zipped,tmpgz);
gunzip(tmpgz);
untar(tmptar,datadir);
delete(tmptar);
delete(tmpgz);


%% Simple import, data with no known sensor locations

S = [];
S.data = fullfile(datadir,'20220623_095714.cMEG');
S.path = datadir;

D = spm_opm_create(S);

%% tidyup - purge all SPM meeg objects saved in cerca folder

d = spm_select('List',datadir,'^*.dat');

for ii = 1:size(d)
    dat = fullfile(datadir,deblank(d(ii,:)));
    tmp = spm_eeg_load(dat);
    delete(tmp);
end

clear tmp
