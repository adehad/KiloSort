ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)		
ops.parfor              = 1; % whether to use parfor to accelerate some parts of the algorithm		
ops.verbose             = 1; % whether to print command line progress		
ops.showfigures         = 0; % whether to plot figures during optimization		
		
ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'		
ops.fbinary             = 'C:\DATA\Spikes\Piroska\piroska_example_short.dat'; % will be created for 'openEphys'		
ops.fproc               = 'C:\DATA\Spikes\Piroska\temp_wh.dat'; % residual from RAM of preprocessed data		
ops.root                = [pathToYourConfigFile,filesep]; % we use a custom 'openEphys' converter
% ops.OEroot              = 'C:\DATA\Spikes\Piroska'; % 'openEphys' only: where raw files are		
		
% ops.fs                  = 25000;        % sampling rate		(omit if already in chanMap file)
% ops.NchanTOT            = 32;           % total number of channels (omit if already in chanMap file)
% ops.Nchan               = 32;           % number of active channels (omit if already in chanMap file)
ops.Nfilt               = 12;             % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 2; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
			% must be less than number of channels
ops.nNeigh              = 3; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
			% Must be less than Nfilt 

% options for channel whitening		
ops.whitening           = 'full'; 	% type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 1; 		% compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 1; 		% how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
% define the channel map as a filename (string) or simply an array		
ops.chanMap             = 'C:\DATA\Spikes\Piroska\chanMap.mat'; % make this file using createChannelMapFile.m		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file		
		
% other options for controlling the model and optimization		
ops.Nrank               = 3;    	% matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    	% number of complete passes through data during optimization (6)		
ops.maxFR               = 20000;  	% maximum number of spikes to extract per batch (20000)		
ops.fshigh              = 300;   	% frequency for high pass filtering		
% ops.fslow             = 2000;   	% frequency for low pass filtering (optional)
ops.ntbuff              = 64;    	% samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   	% int16 scaling of whitened data		
ops.NT                  = 32*1024+ ops.ntbuff; % this is the batch size (try decreasing if out of memory) 		
			% for GPU should be multiple of 32 + ntbuff		

ops.nt0 = 21; % set samples to take on either side of the detected peak (61)
		
	
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
ops.Th               = [4 10 10];   % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam              = [5 20 20];   % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;           % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 400]; % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;           % allow merges and splits during optimization (1)	
ops.freqUpdate       = 200;         % after how many batches, should try split and merge (400)	
ops.mergeT           = 1e-3;        % upper threshold for merging (.1)		
ops.splitT           = 1e-3;        % lower threshold for splitting (.1)		
ops.muTh             = 15;          % minimum mu ("variance") required per cluster (10)
ops.minSpks          = 200;         % minimum number of spikes allowed per cluster (200)
ops.forgetFac        = 0.9975;      % dbins forgetful factor (0.9975)
		
% options for initializing spike templates from data		
ops.initialize      = 'fromData'; 	%'fromData' or 'no'		
ops.spkTh           = -5;      		% spike threshold in standard deviations, KiloSort looks for a minimum peak, so must be negative (-4)		
ops.loc_range       = [10  1];  	% ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; 		% ranges to detect isolated peaks ([30 6])	
ops.maxSpkPeaks     = 2;       		% in the loc_range around a peak, how many other threshold crossings are allowed to consider this an isolated spike (1)	
ops.maskMaxChannels = 5;       		% how many channels to mask up/down ([5])		
ops.crit            = .65;     		% upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10e3;   		% maximum "unique" spikes to consider (10000)		
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = load('PCspikes2.mat'); % Entire PC output, ideally compute for your own data	('PCspikes2.mat')	
ops.centrePC        = 21;                    % (element of spike peak)+1 (21) 		
ops.wPCA            = dd.Wi(:,1:7);   		 % Take first 7 PCs 	
ops.saveInitTemps   = 0;                     % if saveInitTemps         ( 0 )
                                                % true: will try and store initial templates
                                                % false: load from template initFile ('WUinit.mat')
ops.initFilePath    = ops.root;              % if saveInitTemps =       ( [] )
                                                % true: [path to store initial templates, 'WUinit.mat'] to store initial templates
                                                % false: loads initial templates from [path to store initial templates, 'WUinit.mat']
% ops.nt0min = ops.centrePC-1;               %  comment out if you have not computed your PC, used in template centering, ideally should be ops.centrePC-1 ()
		
% options for posthoc merges (under construction)		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		
		
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
