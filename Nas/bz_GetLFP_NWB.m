function [lfp] = bz_GetLFP_NWB(varargin)
% bz_GetLFP - Get local field potentials.
%
%  Load local field potentials from disk. No longer dependent on
%  FMAT/SetCurrentSession.
%
%  USAGE
%
%    [lfp] = bz_GetLFP(channels,<options>)
%
%  INPUTS
%
%    channels(required) -must be first input, numeric  
%                        list of channels to load (use keyword 'all' for all)
%                        channID is 0-indexing, a la neuroscope
%  Name-value paired inputs:
%    'basepath'           - folder in which .lfp file will be found (default
%                           is pwd)
%                           folder should follow buzcode standard:
%                           whateverPath/baseName
%                           and contain file baseName.lfp
%    'basename'           -base file name to load
%    'intervals'          -list of time intervals [0 10; 20 30] to read from 
%                           the LFP file (default is [0 inf])
%    'downsample'         -factor to downsample the LFP (i.e. 'downsample',5
%                           will load a 1250Hz .lfp file at 250Hz)
%    'noPrompts'          -logical (default) to supress any user prompts
%
%  OUTPUT
%
%    lfp             struct of lfp data. Can be a single struct or an array
%                    of structs for different intervals.  lfp(1), lfp(2),
%                    etc for intervals(1,:), intervals(2,:), etc
%    .data           [Nt x Nd] matrix of the LFP data
%    .timestamps     [Nt x 1] vector of timestamps to match LFP data
%    .interval       [1 x 2] vector of start/stop times of LFP interval
%    .channels       [Nd X 1] vector of channel ID's
%    .samplingRate   LFP sampling rate [default = 1250]
%    .duration       duration, in seconds, of LFP interval
%
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    lfp = bz_GetLFP(5,'restrict',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    lfp = bz_GetLFP(5,'restrict',[0 120;240.2 265.23]);
%    % multiple channels
%    lfp = bz_GetLFP([1 2 3 4 10 17],'restrict',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    lfp = bz_GetLFP(3,'restrict',[0 120],'select','number');

% Copyright (C) 2004-2011 by Michaël Zugaro
% editied by David Tingley, 2017
%
% NOTES
% -'select' option has been removed, it allowed switching between 0 and 1
%   indexing.  This should no longer be necessary with .lfp.mat structs
%
% TODO
% add saveMat input 
% expand channel selection options (i.e. region or spikegroup)
% add forcereload



% Added support for NWB: Konstantinos Nasiotis 2019


%% Parse the inputs!

channelsValidation = @(x) isnumeric(x) || strcmp(x,'all');

% parse args
p = inputParser;
addRequired(p,'channels',channelsValidation)
addParameter(p,'basename','',@isstr)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'downsample',1,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);

parse(p,varargin{:})
basename = p.Results.basename;
channels = p.Results.channels;
downsamplefactor = p.Results.downsample;
basepath = p.Results.basepath;
noPrompts = p.Results.noPrompts;

% doing this so you can use either 'intervals' or 'restrict' as parameters to do the same thing
intervals = p.Results.intervals;
restrict = p.Results.restrict;
if isempty(intervals) && isempty(restrict) % both empty
    intervals = [0 Inf];
elseif isempty(intervals) && ~isempty(restrict) % intervals empty, restrict isn't
    intervals = restrict;
end

%% let's check that there is an appropriate NWB file
if isempty(basename)
   %disp('No basename given, so we look for a *nwb file...')
   d = dir([basepath filesep '*nwb']);
   if length(d) > 1 % we assume one .nwb file or this should break
       error('there is more than one .nwb file in this directory?');
   elseif length(d) == 0
       d = dir([basepath filesep '*nwb']);
       if isempty(d)
           error('could not find an nwb file..')
       end
   end
   lfp.Filename = d.name;
   basename = strsplit(lfp.Filename,'.');
   if length(basename) > 2
       base = [];
       for i=1:length(basename)-1
          base = [base basename{i} '.'];
       end
       basename = base(1:end-1);  % this is an fugly hack to make things work with Kenji's naming system...
   else
       basename = basename{1};
   end
   
else
   lfp.Filename = [basename '.lfp'];
end

%% things we can parse from sessionInfo or xml file

sessionInfo = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);

try
    samplingRate = sessionInfo.lfpSampleRate;
catch
    samplingRate = sessionInfo.rates.lfp; % old ugliness we need to get rid of
end
samplingRateLFP = samplingRate./downsamplefactor;

if mod(samplingRateLFP,1)~=0
    error('samplingRate/downsamplefactor must be an integer')
end
%% Channel load options
%Right now this assumes that all means channels 0:nunchannels-1 (neuroscope
%indexing), we could also add options for this to be select region or spike
%group from the xml...
if strcmp(channels,'all')
    channels = sessionInfo.channels;
end

%% get the data
disp('loading LFP file...')
nIntervals = size(intervals,1);


nwb2 = nwbRead([basepath filesep basename '.nwb']);


try
    % Check if the data is in LFP format
    all_lfp_keys = keys(nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries);

    for iKey = 1:length(all_lfp_keys)
        if ismember(all_lfp_keys{iKey}, {'lfp','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
            iLFPDataKey = iKey;
            LFPDataPresent = 1;
            break % Once you find the data don't look for other keys/trouble
        else
            LFPDataPresent = 0;
        end
    end
catch
    LFPDataPresent = 0;
end

if ~LFPDataPresent
    error ('No LFP signals exist in this .nwb file')
end





% returns lfp/bz format
for i = 1:nIntervals
    lfp(i).duration = (intervals(i,2)-intervals(i,1));
    lfp(i).interval = [intervals(i,1) intervals(i,2)];

    % Load data and put into struct
    % we assume 0-indexing like neuroscope, but bz_LoadBinary uses 1-indexing to
    % load....
              
    if intervals(i,2) == Inf
        use_this_bracket = sessionInfo.samples_NWB/sessionInfo.lfpSampleRate;
        data_temp = zeros(length(channels),ceil((use_this_bracket - intervals(i,1))*samplingRateLFP));
    else
        use_this_bracket = intervals(i,2);
        data_temp = zeros(length(channels),ceil((intervals(i,2) - intervals(i,1))*samplingRateLFP));
    end
    
    % This is not optimized yet
    
    for iChannel = 1:length(channels)
        disp(['Now concatenating channel: ' num2str(channels(iChannel))])
        data_temp(iChannel,:) = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.load([channels(iChannel)+1, intervals(i,1)*samplingRate+1],[1, downsamplefactor],[channels(iChannel)+1, use_this_bracket*samplingRate]);
    end 
      
    lfp(i).data = data_temp'; clear data_temp
    lfp(i).timestamps = [lfp(i).interval(1):(1/samplingRateLFP):...
                        (lfp(i).interval(1)+(length(lfp(i).data)-1)/...
                        samplingRateLFP)]';
    lfp(i).channels = channels;
    lfp(i).samplingRate = samplingRateLFP;
    % check if duration is inf, and reset to actual duration...
    if lfp(i).interval(2) == inf
        lfp(i).interval(2) = length(lfp(i).timestamps)/lfp(i).samplingRate;
        lfp(i).duration = (lfp(i).interval(i,2)-lfp(i).interval(i,1));
    end
    
    if isfield(sessionInfo,'region') && isfield(sessionInfo,'channels')
        [~,~,regionidx] = intersect(lfp(i).channels,sessionInfo.channels,'stable');
        lfp(i).region = sessionInfo.region(regionidx); % match region order to channel order..
    end
end
