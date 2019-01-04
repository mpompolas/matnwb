function spikes = bz_GetSpikes_NWB_clean(varargin)
% bz_getSpikes - Get spike timestamps.
%
% USAGE
%
%    spikes = bz_getSpikes(varargin)
% 
% INPUTS
%
%    spikeGroups     -vector subset of shank IDs to load (Default: all)
%    region          -string region ID to load neurons from specific region
%                     (requires sessionInfo file or units->structures in xml)
%    UID             -vector subset of UID's to load 
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    getWaveforms    -logical (default=true) to load mean of raw waveform data
%    forceReload     -logical (default=false) to force loading from
%                     res/clu/spk files
%    saveMat         -logical (default=false) to save in buzcode format
%    noPrompts       -logical (default=false) to supress any user prompts
%    
% OUTPUTS
%
%    spikes - cellinfo struct with the following fields
%          .sessionName    -name of recording file
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%          .region         -region ID for each neuron (especially important large scale, high density probes)
%          .shankID        -shank ID that each neuron was recorded on
%          .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%          .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%          .cluID          -cluster ID, NOT UNIQUE ACROSS SHANKS
%           
% NOTES
%
% This function can be used in several ways to load spiking data.
% Specifically, it loads spiketimes for individual neurons and other
% sessionInfodata that describes each neuron.  Spiketimes can be loaded using the
% UID(1-N), the shank the neuron was on, or the region it was recorded in.
% The default behavior is to load all spikes in a recording. The .shankID
% and .cluID fields can be used to reconstruct the 'units' variable often
% used in FMAToolbox.
% units = [spikes.shankID spikes.cluID];
% 
% 
% first usage recommendation:
% 
%   spikes = bz_getSpikes('saveMat',true); Loads and saves all spiking data
%                                          into buzcode format .cellinfo. struct
% other examples:
%
%   spikes = bz_getSpikes('spikeGroups',1:5); first five shanks
%
%   spikes = bz_getSpikes('region','CA1'); cells tagged as recorded in CA1
%
%   spikes = bz_getSpikes('UID',[1:20]); first twenty neurons
%
%
% written by David Tingley, 2017



% Added support for NWB: Konstantinos Nasiotis 2019
 
%% Deal With Inputs 
spikeGroupsValidation = @(x) assert(isnumeric(x) || strcmp(x,'all'),...
    'spikeGroups must be numeric or "all"');

p = inputParser;
addParameter(p,'spikeGroups','all',spikeGroupsValidation);
addParameter(p,'region','',@isstr); % won't work without sessionInfodata 
addParameter(p,'UID',[],@isvector);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);

parse(p,varargin{:})

spikeGroups = p.Results.spikeGroups;
region = p.Results.region;
UID = p.Results.UID;
basepath = p.Results.basepath;
getWaveforms = p.Results.getWaveforms;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;
noPrompts = p.Results.noPrompts;


[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);


spikes.samplingRate = sessionInfo.rates.wideband;
nChannels = sessionInfo.nChannels;


%% Load the nwb
% Get the name of the nwb based on the basepath set
[base,name] = fileparts(basepath);

nwb2 = nwbRead(fullfile(base,[name,'.nwb']));


%% Convert the spikes saved in NWB to Buzsaki spikes format:

spikes.UID = 0 % 1x101 double % 1:101
spikes.times =0 % 1x101 cell
spikes.shankID = 0% 1x101 double
spikes.cluID  = 0% 1x101 double
spikes.rawWaveform = 0% 1x101 cell
spikes.maxWaveformCh = 0% 1x101 double  %%%%% THIS IS IMPORTANT - TALK TO BEN
spikes.sessionName = 0% 'm120_25secs'
spikes.numcells = 0% 101
spikes.spindices  =0 % 4514x2 double









%% save to buzcode format (before exclusions)
if saveMat
    save([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes')
end


%% filter by spikeGroups input
if ~strcmp(spikeGroups,'all')
    [toRemove] = ~ismember(spikes.shankID,spikeGroups);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.region = removeEmptyCells(spikes.region);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
    end
end
%% filter by region input
if ~isempty(region)
    if ~isfield(spikes,'region') %if no region information in metadata
        error(['You selected to load cells from region "',region,...
            '", but there is no region information in your sessionInfo'])
    end
    
  toRemove = ~ismember(spikes.region,region);
    if sum(toRemove)==length(spikes.UID) %if no cells from selected region
        warning(['You selected to load cells from region "',region,...
            '", but none of your cells are from that region'])
    end
  
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.region = removeEmptyCells(spikes.region);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    if isfield(spikes,'rawWaveform')
        for r = 1:length(toRemove)
            if toRemove(r) == 1
             spikes.rawWaveform{r} = [];
            end
        end
        spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
        spikes.maxWaveformCh(toRemove) = [];
    end
    end
end
%% filter by UID input
if ~isempty(UID)
        [toRemove] = ~ismember(spikes.UID,UID);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.region = removeEmptyCells(spikes.region);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
    end
end

%% Generate spindices matrics
spikes.numcells = length(spikes.UID);
for cc = 1:spikes.numcells
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if spikes.numcells>0
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

%% Check if any cells made it through selection
if isempty(spikes.times) | spikes.numcells == 0
    spikes = [];
end




