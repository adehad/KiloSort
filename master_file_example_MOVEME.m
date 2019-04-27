% Customised for NBits Lab Imperial College London
% Usage requires configfile, masterfile and binary data (see formatBin_wrapper)
% files in the same folder

addpath(genpath('D:\CODE\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('D:\CODE\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = ''; % put masterfile in same folder as config
                           % ops.root in configfile should be the folder 
                           % where the masterfile is
run(fullfile(pathToYourConfigFile, 'StandardConfig_MOVEME.m'))

tic; % start timer

if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% A formatBin_wrapper.m is used for our data, to allow channel padding
% if strcmp(ops.datatype , 'openEphys')
%    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
% end

[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


% save matlab results 
saveOutput(rez,ops,'pre')


%% AUTOMERGES
% Kilosort's AUTO merges should not be confused with the "best" merges 

% rez = merge_posthoc2(rez);

% saveOutput(rez,ops,'post')

% remove temporary file
% delete(ops.fproc);

disp('Done')

% Custom save function
function saveOutput(rez,ops,preOrPost)

if strcmpi(preOrPost,'pre')
    runType = 'preAutoMerge';
    rezName = 'rez.mat';
elseif strcmpi(preOrPost,'post')
    runType = 'postAutoMerge';
    rezName = 'rez_post.mat';
else
    error('Must be ''pre'' or ''post'' ')
end

% save matlab results file
save(fullfile(ops.root,  rezName), 'rez', '-v7.3');

% save python results file for Phy, in a new folder
mkdir( [ops.root, runType] ) 
rezToPhy(rez, [ops.root, filesep, runType, filesep]);


% To prevent conflicts that arise due to the .phy folder that is
% created in the same location as the binary file
% We use the 'link' method, which prevents taking up actual hard drive space
if ismac        % Code to run on Mac platform
    fprintf(['Mac currently not supported. To use Phy: \n', ...
    'Please move the binary file to the folder \n', ...
    'that you want to visualise \n'])

elseif isunix    % Code to run on Linux platform
    command = ['ln',' ',ops.fbinary,' ',runType, filesep, ops.fbinary];
    [status,cmdout] = system(command);
    
elseif ispc      % Code to run on Windows platform
    command = ['mklink /H',' ', runType, filesep,ops.fbinary,' ',ops.fbinary];
    [status,cmdout] = system(command);
    
else
    status = -1;
    cmdout = 'Platform not supported';
    
end
fprintf('Link Status: %i. \n %s',status,cmdout);

end