function spikes = bz_GetSpikes_bypass_clu_fet_spk_NWB(varargin)
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

% Konstantinos Nasiotis 2019



%% TODO
% 1. The template waveform should be filled by: nwb2.units.waveform_mean      The example dataset didn't have any. I added a template
% 2. getWaveforms is not supported. Not sure nwb holds the spiking waveforms - Haven't seen a field





%% Deal With Inputs 
spikeGroupsValidation = @(x) assert(isnumeric(x) || strcmp(x,'all'),...
    'spikeGroups must be numeric or "all"');

p = inputParser;
addParameter(p,'nwb_file','',@isstr);
addParameter(p,'spikeGroups','all',spikeGroupsValidation);
addParameter(p,'region','',@isstr); % won't work without sessionInfodata 
addParameter(p,'UID',[],@isvector);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);

parse(p,varargin{:})

%% Adding support only for spikeGroups, region and UID for now

nwb_file      = p.Results.nwb_file;
spikeGroups   = p.Results.spikeGroups;
region        = p.Results.region;
UID           = p.Results.UID;
% basepath     = p.Results.basepath;
% getWaveforms = p.Results.getWaveforms;
% forceReload  = p.Results.forceReload;
% saveMat      = p.Results.saveMat;
% noPrompts    = p.Results.noPrompts;

nwb2 = nwbRead(nwb_file);
[the_path, name,~] = fileparts(nwb_file);
new_path_for_files = [the_path filesep name];
load([new_path_for_files filesep name '.sessionInfo.mat'])



    nNeurons = length(nwb2.units.id.data.load);

    spikes = struct;
    spikes.samplingRate = sessionInfo.rates.wideband;
    
    % Assign neuron to region based on the region that its Shank belongs to
    shank_that_neurons_belong_to = nwb2.units.vectordata.get('shank').data.load;
    
    for iNeuron = 1:nNeurons
        all_region{iNeuron} = sessionInfo.region{sessionInfo.SpkGrps(shank_that_neurons_belong_to(iNeuron)).Channels(1)+1}; % The channels are 0 indexed
    end
    
    
    
    %% Make the check here of what will be loaded
    
    selected_neurons_UID         = false(1,nNeurons);
    selected_neurons_spikeGroups = false(1,nNeurons);
    selected_neurons_region      = false(1,nNeurons);
    
    % Check which neurons were selected
    if ~isempty(UID)
        selected_neurons_UID(UID) = true;
    end
    %Check which Shanks were selected
    if ~isempty(spikeGroups)
        for iShank = spikeGroups
            selected_neurons_spikeGroups(find(shank_that_neurons_belong_to==iShank)) = true;
        end
    end
    %Check which regions were selected
    if ~isempty(region)
        selected_neurons_region(find(contains(all_region,region))) = true;
    end
    
    UID = find(selected_neurons_UID | selected_neurons_spikeGroups | selected_neurons_region);
    
    if isempty(UID)
        warning('No preference for which Neurons to load was selected. Loading everything..')
        UID = 1:nNeurons;
    end
        
        
        
    %% Get Spikes
    spikes.UID = UID;

    % The template waveform should be filled by: nwb2.units.waveform_mean
    
    template_Waveform = [21.2963012297339,20.2585005424487,21.2206998551635,21.3478476214865,21.5609060407305,...      % This is from a template I used on Kilosort.  
                         22.8392565561944,24.6983630854040,28.2241362812803,30.1107342194246,29.1416620544762,...      % I assign the same on every neuron.
                         25.2447548379813,23.9595314702837,24.8908029479470,24.3822118826549,23.4681225355758,...      % Check if I can get this from nwb
                         23.8598751128954,21.7980194427923,18.0866792366068,11.9973321575690,-2.96486715514583,...
                         -44.0439049558331,-110.074832790885,-186.628097395696,-240.085142069235,-255.067959938651,...
                         -244.823973684355,-220.122942756520,-183.947685024562,-143.562805299476,-104.992358564081,...
                         -69.3806747152833,-35.0576506603005,-3.29132763624548,22.3478476214865,44.7121087898714,...
                         59.7499094771566,73.1141706455415,79.3615933259538,82.9595314702837,82.9182943568817,...
                         82.0591878276721,78.3134833603181,75.0076414359195,70.1691534634109,65.2413184118645,...
                         59.5643424668473,55.0969885149573,51.9698407486342,45.6296345630672,41.3822118826549,...
                         37.1519713328267,31.5334146317958];

    times       = cell(1,length(UID));
    rawWaveform = cell(1,length(UID));
    spindices   = [];
    
    entry = 0;
    for iNeuron = UID
        entry = entry+1;
        if iNeuron == 1
            times_temp = nwb2.units.spike_times.data.load(1:sum(nwb2.units.spike_times_index.data.load(iNeuron)));
        else
            times_temp = nwb2.units.spike_times.data.load(sum(nwb2.units.spike_times_index.data.load(iNeuron-1))+1:sum(nwb2.units.spike_times_index.data.load(iNeuron)));
        end
        times{entry} = times_temp(times_temp~=0);

        rawWaveform{entry} = template_Waveform;     % The template waveform should be filled by: nwb2.units.waveform_mean     The example dataset didn't have any
        spindices = [spindices ; times{entry} ones(length(times{entry}),1)*iNeuron];
    end

    % Spindices have to be sorted according to when each spike occured
    [~,sortedIndices] = sort(spindices(:,1));
    spindices = spindices(sortedIndices,:);
    
    shankID_of_selected_Neurons = nwb2.units.vectordata.get('shank').data.load';
    shankID_of_selected_Neurons = shankID_of_selected_Neurons(UID);

    spikes.times        = times;
    spikes.shankID      = shankID_of_selected_Neurons;
    spikes.cluID        = ones(1,length(UID))*(-1);   % THESE ARE THE SPIKING TEMPLATES. THEY ARE FILLED FROM KILOSORT. I ADD A NEGATIVE VALUE TO SEE IF IT CAUSES AN ERROR SOMEWHERE
    spikes.rawWaveform  = rawWaveform;
    spikes.maxWaveformCh = ones(1,length(UID))*(-1);     % THESE ASSIGN THE MAXIMUM WAVEFORM TO A ACHANNEL. CHECK HOW TO ADD THIS. I ADD A NEGATIVE VALUE TO SEE IF IT CAUSES AN ERROR SOMEWHERE
    spikes.sessionName  = sessionInfo.FileName;
    spikes.numcells     = length(UID);
    spikes.spindices    = spindices;            % This holds the timing of each spike, sorted, and the neuron it belongs to.

    

    for iNeuron = 1:length(UID)
        spikes.region{iNeuron} = sessionInfo.region{sessionInfo.SpkGrps(shank_that_neurons_belong_to(UID(iNeuron))).Channels(1)+1}; % The channels are 0 indexed
    end
    


end








