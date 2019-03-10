function spikes = bz_GetSpikes_NWB(varargin)
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
disp(' 1. The template waveform should be filled by: nwb2.units.waveform_mean      The example dataset didnt have any. I added a template')
disp(' 2. The maxWaveformCh should be added on the example dataset')
disp(' 3. CluID is probably filled by Kilosort - shouldbe added on the example dataset')
% 4. getWaveforms is not supported. Not sure nwb holds the spiking waveforms - Haven't seen a field





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
saveMat       = p.Results.saveMat;
% noPrompts    = p.Results.noPrompts;

nwb2 = nwbRead(nwb_file);
[the_path, name,~] = fileparts(nwb_file);
new_path_for_files = [the_path filesep name];
load([new_path_for_files filesep name '.sessionInfo.mat'])



    nNeurons = length(nwb2.units.id.data.load);

    
    
    
% % % % % % % % % % % % % % % % % % % % % %     
% CHEK THIS     nwb2.units.electrode_group.data(1).refresh(nwb2)
% % % % % % % % % % % % % % % % % % % %     
    

    % Assign neuron to region based on the region that its Shank belongs to

     % Get a single electrode that belongs in that shank
    electrodes2Shank  = nwb2.general_extracellular_ephys_electrodes.vectordata.get('group_name').data;    
    electrodes2Region = nwb2.general_extracellular_ephys_electrodes.vectordata.get('location').data;
    
    neurons2ShankNames           = cell(1,nNeurons);
    all_region                   = cell(1,nNeurons);
    
    for iNeuron = 1:length(nwb2.units.electrode_group.data)
        path = nwb2.units.electrode_group.data(iNeuron).path;
        inds = strfind(path, '/');
        neurons2ShankNames{iNeuron} = path(inds(end)+1:end);
        
        ElectrodesInThatShank = find(strcmp(electrodes2Shank,neurons2ShankNames{iNeuron}));        
        all_region{iNeuron}   = electrodes2Region{ElectrodesInThatShank(1)};
        
    end

    
    uniqueShanks = unique(neurons2ShankNames,'legacy');
    shankID_of_selected_Neurons = zeros(1, nNeurons);
    for iNeuron = 1:nNeurons
        shankID_of_selected_Neurons(iNeuron) = find(strcmp(uniqueShanks, neurons2ShankNames{iNeuron}));
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
            selected_neurons_spikeGroups(find(shankID_of_selected_Neurons==iShank)) = true;
        end
    end
    %Check which regions were selected
    if ~isempty(region)
        selected_neurons_region(find(contains(all_region,region))) = true;
    end
    
    UID = find(selected_neurons_UID | selected_neurons_spikeGroups | selected_neurons_region);
    
    if isempty(UID)
        warning(['Warning: You selected to load cells from region "' region '", but none of your cells are from that region'])
        spikes = [];
        return
    end
        
    
    
    %% The first time that this function is loaded, make sure to save a .spikes.cellInfo.mat that contains information from all neurons
    
    if saveMat
        UID_forSaveMAt = 1:nNeurons;
        temp = get_the_spikes_from_selected_UIDs(UID_forSaveMAt, nwb2, sessionInfo, saveMat, all_region, neurons2ShankNames); clear temp       
    end
        
    spikes = get_the_spikes_from_selected_UIDs(UID, nwb2, sessionInfo, 0, all_region, neurons2ShankNames);
    
    
    
    
    
    
    
    

end




function spikes = get_the_spikes_from_selected_UIDs(UID, nwb2, sessionInfo, saveMat, all_region, neurons2ShankNames)
        
    %% Get Spikes
    spikes = struct;
    spikes.samplingRate = sessionInfo.rates.wideband;
    spikes.UID = UID;

    % The template waveform should be filled by: nwb2.units.waveform_mean

    template_Waveform = [21.2963012297339,21.2206998551635,21.5609060407305,...      % This is from a template I used on Kilosort.  
                         22.8392565561944,28.2241362812803,30.1107342194246,...      % I assign the same on every neuron.
                         23.9595314702837,24.3822118826549,23.4681225355758,...      % Check if I can get this from nwb
                         18.0866792366068,11.9973321575690,-2.96486715514583,...
                         -110.074832790885,-186.628097395696,-240.085142069235,...
                         -183.947685024562,-143.562805299476,-104.992358564081,...
                         -3.29132763624548,22.3478476214865,44.7121087898714,...
                         73.1141706455415,79.3615933259538,82.9182943568817,...
                         75.0076414359195,70.1691534634109,65.2413184118645,...
                         51.9698407486342,45.6296345630672,41.3822118826549,...
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

    
    uniqueShanks = unique(neurons2ShankNames,'legacy');
    shankID_of_selected_Neurons = zeros(1, length(UID));
    for iNeuron = 1:length(UID)
        shankID_of_selected_Neurons(iNeuron) = find(strcmp(uniqueShanks, neurons2ShankNames{UID(iNeuron)}));
    end
    

    spikes.times        = times;
    spikes.shankID      = shankID_of_selected_Neurons;
    spikes.cluID        = ones(1,length(UID))*(-1);      % THESE ARE THE SPIKING TEMPLATES. THEY ARE FILLED FROM KILOSORT. I ADD A NEGATIVE VALUE TO SEE IF IT CAUSES AN ERROR SOMEWHERE
    spikes.rawWaveform  = rawWaveform;
    spikes.maxWaveformCh = ones(1,length(UID))*(-1);     % THESE ASSIGN THE MAXIMUM WAVEFORM TO A ACHANNEL. CHECK HOW TO ADD THIS. I ADD A NEGATIVE VALUE TO SEE IF IT CAUSES AN ERROR SOMEWHERE
    spikes.sessionName  = sessionInfo.FileName;
    spikes.numcells     = length(UID);
    spikes.spindices    = spindices;                     % This holds the timing of each spike, sorted, and the neuron it belongs to.

    spikes.region = all_region(UID);

    if saveMat
        save([sessionInfo.FileName filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes')
    end
end





