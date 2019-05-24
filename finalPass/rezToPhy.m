function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% spikeTimes will be in samples, not seconds


outputs = {'amplitudes.npy', 'channel_map.npy', 'channel_positions.npy', 'pc_features.npy', ...
           'pc_feature_ind.npy', 'similar_templates.npy', 'spike_clusters.npy', 'spike_templates.npy', ...
           'spike_times.npy', 'templates.npy', 'templates_ind.npy', 'template_features.npy', ...
           'template_feature_ind.npy', 'whitening_mat.npy', 'whitening_mat_inv.npy', 'params.py'};

fs = dir(fullfile(savePath, '*.*py'));
for i = 1:length(fs)
    fname = fs(i).name;
    % don't delete .npy files which have nothing to do with us
    if find(strcmp(fname, outputs))
        delete(fullfile(savePath, fname));
    end
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

% clean up rez duplicates
[~,uniqueIdxs,~] = unique(rez.st3,'rows');
rez.st3 = rez.st3(uniqueIdxs,:);
rez.cProj = rez.cProj(uniqueIdxs,:);
rez.cProjPC = rez.cProjPC(uniqueIdxs,:,:);

spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;

% try
%     load(rez.ops.chanMap);
% catch
%    chanMap0ind  = [0:Nchan-1]';
%    connected    = ones(Nchan, 1);
%    xcoords      = ones(Nchan, 1);
%    ycoords      = (1:Nchan)';
% end
% chanMap0 = chanMap(connected>1e-6);

connected   = rez.connected(:);
xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

% for i = 1:length(chanMap0)
%     chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
% end
% [~, invchanMap0] = sort(chanMap0);

templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
for iNN = 1:rez.ops.Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

if ~isempty(savePath)
    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(int32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(int32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
%     Fs = rez.ops.fs;
    conn        = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savePath, 'channel_map.npy'));
    %writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savePath, 'channel_positions.npy'));
    
    writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    whiteningMatrix = rez.Wrot/200;
    whiteningMatrixInv = whiteningMatrix^-1;
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end
    
    
    % Make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        % include relative path elements in dat_path
        if rez.linkBinary
            dat_path = rez.ops.fbinary;
            rez = rmfield(rez,'linkBinary'); % remove temporary variable
        else
            try
                dat_path = getRelativePath( GetFullPath(rez.ops.fbinary),GetFullPath(savePath) );
            catch
                warning('Could not specify a relative path, you may need to manually change params.py ')
                [dat_path, fname, ext] = fileparts(rez.ops.fbinary);
                dat_path = split(dat_path, filesep);
                dat_path = fullfile( dat_path{find(strcmp(dat_path,'..'),1):end}, [fname ext]);
            end
        end
        
        dat_path = replace(dat_path, '\', '/'); % / is the preferred file separator (stops frprintf conflicts on Windows)
        
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        fprintf(fid,['dat_path = ''',dat_path '''\n']);  % locations of binary file, can be relative - e.g. ../file.bin
        fprintf(fid,['dir_path = ''.',filesep '''\n']);  % location of spike output, '.' means current folder, keeps phy from wandering everywhere
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False');
        fclose(fid);
    end

else
    warning(['Empty path specified, skipping save.', newline, ...
            'To save in current folder provide ''.'' as the path', newline, ...
            'e.g.: rezToPhy(rez, [''.''])'])
end



function rel_path = getRelativePath(tgt_path, act_path)
%RELATIVEPATH  returns the relative path from an actual path to the target path.
%   Both arguments must be strings with absolute paths.
%   The actual path is optional, if omitted the current dir is used instead.
%   In case the volume drive letters don't match, an absolute path will be returned.
%   If a relative path is returned, it always starts with './' or '../'
%
%   Syntax:
%      rel_path = getRelativePath(target_path, actual_path)
%   
%   Parameters:
%      target_path        - Path which is targetted
%      actual_path        - Start for relative path (optional, default = current dir)
%
%   Examples:
%      relativepath('/local/data/matlab', '/local') = './data/matlab/'
%      relativepath('/MyProject/', '/local')        = '/myproject\'
%
%      relativepath('/local/data/matlab', pwd) is the same as
%      relativepath('/local/data/matlab')
%
%   See also:  ABSOLUTEPATH PATH
%   Jochen Lenz
%   Modified by Ken Chatfield to make case sensitive
%            and to also support filenames for the tgt_path parameter
% https://uk.mathworks.com/matlabcentral/fileexchange/41253-improved-relativepath-m
% 2nd parameter is optional:
if nargin < 2
    act_path = pwd;
end
% Predefine return string:
rel_path = '';
% Ensure act_path ends with a filesep character
if isempty(act_path) || ~isequal(act_path(end), filesep)
    act_path = [act_path filesep];
end
% If there is a file with an extension, save it for later
[tgt_path, tgt_fname, tgt_ext] = fileparts(tgt_path);
if isempty(tgt_ext)
    % treat extensionless files as part of the path
    tgt_path = fullfile(tgt_path, tgt_fname);
    tgt_fname = '';
else
    tgt_fname = [tgt_fname tgt_ext];
end
% Ensure tgt_path ends with a filesep character
if isempty(tgt_path) || ~isequal(tgt_path(end), filesep)
    tgt_path = [tgt_path filesep];
end
% Create a cell-array containing the directory levels
act_path_cell = pathparts(act_path);
tgt_path_cell = pathparts(tgt_path);
% If volumes are different, return absolute path on Windows
process_paths = true;
if ispc
   if ~isequal(act_path_cell{1}, tgt_path_cell{1})
       rel_path = tgt_path;
       process_paths = false;
   end
end
if process_paths
    % Remove level by level, as long as both are equal
    while ~isempty(act_path_cell)  && ~isempty(tgt_path_cell)
        if isequal(act_path_cell{1}, tgt_path_cell{1})
            act_path_cell(1) = [];
            tgt_path_cell(1) = [];
        else
            break
        end
    end
    % As much levels down ('../') as levels are remaining in "act_path"
    rel_path = [repmat(['..' filesep],1,length(act_path_cell)), rel_path];
    % Relative directory levels to target directory:
    rel_dirs_cell = cellfun(@(x) [x filesep], tgt_path_cell, 'UniformOutput', false);
    rel_dirs = [rel_dirs_cell{:}];
    rel_path = [rel_path rel_dirs];
    % Start with '.' or '..' :
    if isempty(rel_path)
        rel_path = ['.' filesep];
    elseif ~isequal(rel_path(1),'.')
        rel_path = ['.' filesep rel_path];
    end
end
% add on the original filename if appropriate
rel_path = fullfile(rel_path, tgt_fname);
function path_cell = pathparts(path_str)
    path_str = [filesep path_str filesep];
    path_cell = regexp(path_str, filesep, 'split');
    path_cell(strcmp(path_cell, '')) = [];
end

end

function File = GetFullPath(File, Style)
% GetFullPath - Get absolute canonical path of a file or folder
% Absolute path names are safer than relative paths, when e.g. a GUI or TIMER
% callback changes the current directory. Only canonical paths without "." and
% ".." can be recognized uniquely.
% Long path names (>259 characters) require a magic initial key "\\?\" to be
% handled by Windows API functions, e.g. for Matlab's FOPEN, DIR and EXIST.
%
% FullName = GetFullPath(Name, Style)
% INPUT:
%   Name:  String or cell string, absolute or relative name of a file or
%          folder. The path need not exist. Unicode strings, UNC paths and long
%          names are supported.
%   Style: Style of the output as string, optional, default: 'auto'.
%          'auto': Add '\\?\' or '\\?\UNC\' for long names on demand.
%          'lean': Magic string is not added.
%          'fat':  Magic string is added for short names also.
%          The Style is ignored when not running under Windows.
%
% OUTPUT:
%   FullName: Absolute canonical path name as string or cell string.
%          For empty strings the current directory is replied.
%          '\\?\' or '\\?\UNC' is added on demand.
%
% NOTE: The M- and the MEX-version create the same results, the faster MEX
%   function works under Windows only.
%   Some functions of the Windows-API still do not support long file names.
%   E.g. the Recycler and the Windows Explorer fail even with the magic '\\?\'
%   prefix. Some functions of Matlab accept 260 characters (value of MAX_PATH),
%   some at 259 already. Don't blame me.
%   The 'fat' style is useful e.g. when Matlab's DIR command is called for a
%   folder with les than 260 characters, but together with the file name this
%   limit is exceeded. Then "dir(GetFullPath([folder, '\*.*], 'fat'))" helps.
%
% EXAMPLES:
%   cd(tempdir);                    % Assumed as 'C:\Temp' here
%   GetFullPath('File.Ext')         % 'C:\Temp\File.Ext'
%   GetFullPath('..\File.Ext')      % 'C:\File.Ext'
%   GetFullPath('..\..\File.Ext')   % 'C:\File.Ext'
%   GetFullPath('.\File.Ext')       % 'C:\Temp\File.Ext'
%   GetFullPath('*.txt')            % 'C:\Temp\*.txt'
%   GetFullPath('..')               % 'C:\'
%   GetFullPath('..\..\..')         % 'C:\'
%   GetFullPath('Folder\')          % 'C:\Temp\Folder\'
%   GetFullPath('D:\A\..\B')        % 'D:\B'
%   GetFullPath('\\Server\Folder\Sub\..\File.ext')
%                                   % '\\Server\Folder\File.ext'
%   GetFullPath({'..', 'new'})      % {'C:\', 'C:\Temp\new'}
%   GetFullPath('.', 'fat')         % '\\?\C:\Temp\File.Ext'
%
% COMPILE:
%   Automatic: InstallMex GetFullPath.c uTest_GetFullPath
%   Manual:    mex -O GetFullPath.c
%   Download:  http://www.n-simon.de/mex
% Run the unit-test uTest_GetFullPath after compiling.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
% Assumed Compatibility: higher Matlab versions
% Author: Jan Simon, Heidelberg, (C) 2009-2016 matlab.2010(a)n(MINUS)simon.de
%
% See also: CD, FULLFILE, FILEPARTS.
% $JRev: R-G V:032 Sum:zBDFj0/m8a0f Date:15-Jan-2013 01:06:12 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_GetFullPath $
% $File: Tools\GLFile\GetFullPath.m $
% History:
% 001: 20-Apr-2010 22:28, Successor of Rel2AbsPath.
% 010: 27-Jul-2008 21:59, Consider leading separator in M-version also.
% 011: 24-Jan-2011 12:11, Cell strings, '~File' under linux.
%      Check of input types in the M-version.
% 015: 31-Mar-2011 10:48, BUGFIX: Accept [] as input as in the Mex version.
%      Thanks to Jiro Doke, who found this bug by running the test function for
%      the M-version.
% 020: 18-Oct-2011 00:57, BUGFIX: Linux version created bad results.
%      Thanks to Daniel.
% 024: 10-Dec-2011 14:00, Care for long names under Windows in M-version.
%      Improved the unittest function for Linux. Thanks to Paul Sexton.
% 025: 09-Aug-2012 14:00, In MEX: Paths starting with "\\" can be non-UNC.
%      The former version treated "\\?\C:\<longpath>\file" as UNC path and
%      replied "\\?\UNC\?\C:\<longpath>\file".
% 032: 12-Jan-2013 21:16, 'auto', 'lean' and 'fat' style.
% Initialize: ==================================================================
% Do the work: =================================================================
% #############################################
% ### USE THE MUCH FASTER MEX ON WINDOWS!!! ###
% #############################################
% Difference between M- and Mex-version:
% - Mex does not work under MacOS/Unix.
% - Mex calls Windows API function GetFullPath.
% - Mex is much faster.
% Magix prefix for long Windows names:
if nargin < 2
   Style = 'auto';
end
% Handle cell strings:
% NOTE: It is faster to create a function @cell\GetFullPath.m under Linux, but
% under Windows this would shadow the fast C-Mex.
if isa(File, 'cell')
   for iC = 1:numel(File)
      File{iC} = GetFullPath(File{iC}, Style);
   end
   return;
end
% Check this once only:
isWIN    = strncmpi(computer, 'PC', 2);
MAX_PATH = 260;
% Warn once per session (disable this under Linux/MacOS):
persistent hasDataRead
if isempty(hasDataRead)
   % Test this once only - there is no relation to the existence of DATAREAD!
   %if isWIN
   %   Show a warning, if the slower Matlab version is used - commented, because
   %   this is not a problem and it might be even useful when the MEX-folder is
   %   not inlcuded in the path yet.
   %   warning('JSimon:GetFullPath:NoMex', ...
   %      ['GetFullPath: Using slow Matlab-version instead of fast Mex.', ...
   %       char(10), 'Compile: InstallMex GetFullPath.c']);
   %end
   
   % DATAREAD is deprecated in 2011b, but still available. In Matlab 6.5, REGEXP
   % does not know the 'split' command, therefore DATAREAD is preferred:
   hasDataRead = ~isempty(which('dataread'));
end
if isempty(File)  % Accept empty matrix as input:
   if ischar(File) || isnumeric(File)
      File = cd;
      return;
   else
      error(['JSimon:', mfilename, ':BadTypeInput1'], ...
         ['*** ', mfilename, ': Input must be a string or cell string']);
   end
end
if ischar(File) == 0  % Non-empty inputs must be strings
   error(['JSimon:', mfilename, ':BadTypeInput1'], ...
      ['*** ', mfilename, ': Input must be a string or cell string']);
end
if isWIN  % Windows: --------------------------------------------------------
   FSep = '\';
   File = strrep(File, '/', FSep);
   
   % Remove the magic key on demand, it is appended finally again:
   if strncmp(File, '\\?\', 4)
      if strncmpi(File, '\\?\UNC\', 8)
         File = ['\', File(7:length(File))];  % Two leading backslashes!
      else
         File = File(5:length(File));
      end
   end
   
   isUNC   = strncmp(File, '\\', 2);
   FileLen = length(File);
   if isUNC == 0                        % File is not a UNC path
      % Leading file separator means relative to current drive or base folder:
      ThePath = cd;
      if File(1) == FSep
         if strncmp(ThePath, '\\', 2)   % Current directory is a UNC path
            sepInd  = strfind(ThePath, '\');
            ThePath = ThePath(1:sepInd(4));
         else
            ThePath = ThePath(1:3);     % Drive letter only
         end
      end
      
      if FileLen < 2 || File(2) ~= ':'  % Does not start with drive letter
         if ThePath(length(ThePath)) ~= FSep
            if File(1) ~= FSep
               File = [ThePath, FSep, File];
            else                        % File starts with separator:
               File = [ThePath, File];
            end
         else                           % Current path ends with separator:
            if File(1) ~= FSep
               File = [ThePath, File];
            else                        % File starts with separator:
               ThePath(length(ThePath)) = [];
               File = [ThePath, File];
            end
         end
         
      elseif FileLen == 2 && File(2) == ':'   % "C:" current directory on C!
         % "C:" is the current directory on the C-disk, even if the current
         % directory is on another disk! This was ignored in Matlab 6.5, but
         % modern versions considers this strange behaviour.
         if strncmpi(ThePath, File, 2)
            File = ThePath;
         else
            try
               File = cd(cd(File));
            catch    % No MException to support Matlab6.5...
               if exist(File, 'dir')  % No idea what could cause an error then!
                  rethrow(lasterror);
               else  % Reply "K:\" for not existing disk:
                  File = [File, FSep];
               end
            end
         end
      end
   end
   
else         % Linux, MacOS: ---------------------------------------------------
   FSep = '/';
   File = strrep(File, '\', FSep);
   
   if strcmp(File, '~') || strncmp(File, '~/', 2)  % Home directory:
      HomeDir = getenv('HOME');
      if ~isempty(HomeDir)
         File(1) = [];
         File    = [HomeDir, File];
      end
      
   elseif strncmpi(File, FSep, 1) == 0
      % Append relative path to current folder:
      ThePath = cd;
      if ThePath(length(ThePath)) == FSep
         File = [ThePath, File];
      else
         File = [ThePath, FSep, File];
      end
   end
end
% Care for "\." and "\.." - no efficient algorithm, but the fast Mex is
% recommended at all!
if ~isempty(strfind(File, [FSep, '.']))
   if isWIN
      if strncmp(File, '\\', 2)  % UNC path
         index = strfind(File, '\');
         if length(index) < 4    % UNC path without separator after the folder:
            return;
         end
         Drive            = File(1:index(4));
         File(1:index(4)) = [];
      else
         Drive     = File(1:3);
         File(1:3) = [];
      end
   else  % Unix, MacOS:
      isUNC   = false;
      Drive   = FSep;
      File(1) = [];
   end
   
   hasTrailFSep = (File(length(File)) == FSep);
   if hasTrailFSep
      File(length(File)) = [];
   end
   
   if hasDataRead
      if isWIN  % Need "\\" as separator:
         C = dataread('string', File, '%s', 'delimiter', '\\');  %#ok<REMFF1>
      else
         C = dataread('string', File, '%s', 'delimiter', FSep);  %#ok<REMFF1>
      end
   else  % Use the slower REGEXP, when DATAREAD is not available anymore:
      C = regexp(File, FSep, 'split');
   end
   
   % Remove '\.\' directly without side effects:
   C(strcmp(C, '.')) = [];
   
   % Remove '\..' with the parent recursively:
   R = 1:length(C);
   for dd = reshape(find(strcmp(C, '..')), 1, [])
      index    = find(R == dd);
      R(index) = [];
      if index > 1
         R(index - 1) = [];
      end
   end
   
   if isempty(R)
      File = Drive;
      if isUNC && ~hasTrailFSep
         File(length(File)) = [];
      end
      
   elseif isWIN
      % If you have CStr2String, use the faster:
      %   File = CStr2String(C(R), FSep, hasTrailFSep);
      File = sprintf('%s\\', C{R});
      if hasTrailFSep
         File = [Drive, File];
      else
         File = [Drive, File(1:length(File) - 1)];
      end
      
   else  % Unix:
      File = [Drive, sprintf('%s/', C{R})];
      if ~hasTrailFSep
         File(length(File)) = [];
      end
   end
end
% "Very" long names under Windows:
if isWIN
   if ~ischar(Style)
      error(['JSimon:', mfilename, ':BadTypeInput2'], ...
         ['*** ', mfilename, ': Input must be a string or cell string']);
   end
   
   if (strncmpi(Style, 'a', 1) && length(File) >= MAX_PATH) || ...
         strncmpi(Style, 'f', 1)
      % Do not use [isUNC] here, because this concerns the input, which can
      % '.\File', while the current directory is an UNC path.
      if strncmp(File, '\\', 2)  % UNC path
         File = ['\\?\UNC', File(2:end)];
      else
         File = ['\\?\', File];
      end
   end
end
% return;
end

end