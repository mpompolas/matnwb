classdef Neuroscope2NWB
 
% This function wraps all electrophysiological, behavioral and analyzed
% data into a single NWB file.

% It was tested on the YutaMouse41-150903 dataset:
% https://buzsakilab.nyumc.org/datasets/SenzaiY/YutaMouse41/YutaMouse41-150903/

% Only the folder needs to be specified. It assumes that all Neuroscope,
% .eeg and .mat files exist within the specified folder.


% Konstantinos Nasiotis 2019

    
    methods(Static)
       
        function xml = GetXMLInfo(folder_path)
            %% Enter the folder and run everything from in there (easier paths).
            %  When the converter is done, go back to the initial folder

            [previous_paths, name] = fileparts(folder_path);

%             current_folder = pwd;
%             cd (folder_path)

            %% This uses an .xml importer downloaded from MathWorks - File Exchange
            %  https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
            %  The loadxml from the Buzcode repo gave errors

            all_files_in_folder = dir(folder_path);

            iXML = [];
            for iFile = 1:length(all_files_in_folder)
                if strfind(all_files_in_folder(iFile).name,'.xml')
                    iXML = [iXML iFile];
                end
            end
            if isempty(iXML)
                error 'There are no .xml files in this folder'
            elseif length(iXML)>1
                error 'There is more than one .xml in this folder'
            end

            xml               = xml2struct([folder_path filesep all_files_in_folder(iXML).name]);
            xml               = xml.parameters;
            xml.folder_path   = folder_path;
            xml.name          = name;
        end
        
        
        
        function nwb = GeneralInfo(xml)
          
          
            %% General Info
            nwb_version = '2.0b';

            session_start_time = datetime(xml.generalInfo.date.Text, ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
                'TimeZone', 'local');
            timestamps_reference_time = datetime(xml.generalInfo.date.Text, ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
                'TimeZone', 'local');

            file_create_date = datetime(datestr(clock), ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
                'TimeZone', 'local');


            nwb = nwbfile( ...
                'session_description'          , 'Mouse in open exploration and theta maze', ...
                'identifier'                   , xml.name, ...
                'session_start_time'           , session_start_time,...
                'file_create_date'             , file_create_date,...
                'general_experimenter'         , xml.generalInfo.experimenters.Text,...
                'general_session_id'           , xml.name,...
                'general_institution'          , 'NYU'  ,...
                'general_lab'                  , 'Buzsaki',...
                'subject'                      , 'YutaMouse',...
                'general_related_publications' , 'DOI:10.1016/j.neuron.2016.12.011',...
                'timestamps_reference_time'    , session_start_time);

            nwb.general_subject = types.core.Subject( ...
                'description', 'mouse 5', 'genotype', 'POMC-Cre::Arch', 'age', '9 months', ...
                'sex', 'M', 'subject_id', xml.name, 'species', 'Mus musculus');
        end
        
        
        function nwb = getElectrodeInfo (xml,nwb)
            %% Get the electrodes' info

            nShanks = length(xml.spikeDetection.channelGroups.group);

            groups = xml.spikeDetection.channelGroups.group; % Use this for simplicity

            all_shank_channels = cell(nShanks,1); % This will hold the channel numbers that belong in each shank

            % Initialize variables
            x                 = [];
            y                 = [];
            z                 = [];
            imp               = [];
            location          = [];
            shank             = [];
            group_name        = [];
            group_object_view = [];
            filtering         = [];
            shank_channel     = [];
            amp_channel_id    = [];

            device_name = 'implant';
            nwb.general_devices.set(device_name, types.core.Device());
            device_link = types.untyped.SoftLink(['/general/devices/' device_name]);

            for iGroup = 1:nShanks
                for iChannel = 1:length(groups{iGroup}.channels.channel)
                    all_shank_channels{iGroup} = [all_shank_channels{iGroup} str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
                    shank_channel     = [shank_channel; iChannel-1];
                    amp_channel_id    = [amp_channel_id; str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
                    shank             = [shank; iGroup];
                    group_name        = [group_name; ['shank' num2str(iGroup)]];
                    group_object_view = [group_object_view; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(iGroup)]])];

                    if ~isfield(groups{iGroup}.channels.channel{iChannel},'position')
                        x = [x; NaN];
                        y = [y; NaN];
                        z = [z; NaN];
                    end
                    if ~isfield(groups{iGroup}.channels.channel{iChannel},'imp')
                        imp = [imp; NaN];
                    end  
                    if ~isfield(groups{iGroup}.channels.channel{iChannel},'location')
                        location{end+1,1} = 'unknown';
                    end  
                    if ~isfield(groups{iGroup}.channels.channel{iChannel},'filtering')
                        filtering = [filtering; NaN];
                    end      
                end
                nwb.general_extracellular_ephys.set(['shank' num2str(iGroup)], ...
                    types.core.ElectrodeGroup( ...
                    'description', ['electrode group for shank' num2str(iGroup)], ...
                    'location', 'unknown', ...
                    'device', device_link));
            end

            variables = {'x'; 'y'; 'z'; 'imp'; 'location'; 'filtering'; 'group'; 'group_name'; 'shank'; 'shank_channel'; 'amp_channel'};

            % In order to insert string to a table, they need to be converted to a cell
            % first (e.g. location(iElectrode))
            for iElectrode = 1:length(x)
                if iElectrode == 1
                    tbl = table(x(iElectrode),y(iElectrode),z(iElectrode),imp(iElectrode),{location{iElectrode}},filtering(iElectrode),group_object_view(iElectrode),{group_name(iElectrode,:)},shank(iElectrode),shank_channel(iElectrode),amp_channel_id(iElectrode),...
                               'VariableNames', variables);
                else
                    tbl = [tbl; {x(iElectrode),y(iElectrode),z(iElectrode),imp(iElectrode),{location{iElectrode}},filtering(iElectrode),group_object_view(iElectrode),{group_name(iElectrode,:)},shank(iElectrode),shank_channel(iElectrode),amp_channel_id(iElectrode)}];
                end
            end

            % add the |DynamicTable| object to the NWB file in
            % /general/extracellular_ephys/electrodes
            electrode_table = util.table2nwb(tbl, 'metadata about extracellular electrodes');
            nwb.general_extracellular_ephys_electrodes = electrode_table;

        end
        
        
        
        
        function nwb = getUnitsInfo(xml, nwb)
            
            %% Add the units info (copied from bz_GetSpikes)

            getWaveforms = 1; % Set this to true if you want to add waveforms on the NWB file


            spikes.samplingRate = str2double(xml.acquisitionSystem.samplingRate.Text);


            disp('loading spikes from clu/res/spk files..')
            % find res/clu/fet/spk files here
            cluFiles = dir([xml.folder_path filesep '*.clu*']);  
            resFiles = dir([xml.folder_path filesep '*.res*']);
            if any(getWaveforms)
                spkFiles = dir([xml.folder_path filesep '*.spk*']);
            end

            % remove *temp*, *autosave*, and *.clu.str files/directories
            tempFiles = zeros(length(cluFiles),1);
            for i = 1:length(cluFiles) 
                dummy = strsplit(cluFiles(i).name, '.'); % Check whether the component after the last dot is a number or not. If not, exclude the file/dir. 
                if ~isempty(findstr('temp',cluFiles(i).name)) | ~isempty(findstr('autosave',cluFiles(i).name)) | isempty(str2num(dummy{length(dummy)})) | find(contains(dummy, 'clu')) ~= length(dummy)-1  
                    tempFiles(i) = 1;
                end
            end
            cluFiles(tempFiles==1)=[];
            tempFiles = zeros(length(resFiles),1);
            for i = 1:length(resFiles)
                if ~isempty(findstr('temp',resFiles(i).name)) | ~isempty(findstr('autosave',resFiles(i).name))
                    tempFiles(i) = 1;
                end
            end
            if any(getWaveforms)
                resFiles(tempFiles==1)=[];
                tempFiles = zeros(length(spkFiles),1);
                for i = 1:length(spkFiles)
                    if ~isempty(findstr('temp',spkFiles(i).name)) | ~isempty(findstr('autosave',spkFiles(i).name))
                        tempFiles(i) = 1;
                    end
                end
                spkFiles(tempFiles==1)=[];
            end

            if isempty(cluFiles)
                disp('no clu files found...')
                spikes = [];
                return
            end


            % ensures we load in sequential order (forces compatibility with FMAT
            % ordering)
            for i = 1:length(cluFiles)
                temp = strsplit(cluFiles(i).name,'.');
                shanks(i) = str2num(temp{length(temp)});
            end
            [shanks ind] = sort(shanks);
            cluFiles = cluFiles(ind); %Bug here if there are any files x.clu.x that are not your desired clus
            resFiles = resFiles(ind);
            if any(getWaveforms)
                spkFiles = spkFiles(ind);
            end

            % check if there are matching #'s of files
            if length(cluFiles) ~= length(resFiles) && length(cluFiles) ~= length(spkFiles)
                error('found an incorrect number of res/clu/spk files...')
            end

            % use the .clu files to get spike ID's and generate UID and spikeGroup
            % use the .res files to get spike times
            count = 1;

            ecephys = types.core.ProcessingModule;

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is copied from the ElectrodesInfo
            nShanks = length(xml.spikeDetection.channelGroups.group);
            groups = xml.spikeDetection.channelGroups.group; % Use this for simplicity
            all_shank_channels = cell(nShanks,1); % This will hold the channel numbers that belong in each shank
            shank = [];
            group_object_view = [];
            
            for iGroup = 1:nShanks
            % Get all_shank_channls again  for iGroup = 1:nShanks
                for iChannel = 1:length(groups{iGroup}.channels.channel)
                    all_shank_channels{iGroup} = [all_shank_channels{iGroup} str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
                    shank = [shank iGroup];
                    group_object_view = [group_object_view; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(iGroup)]])];
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

            for iShank=1:length(cluFiles) 
                disp(['working on ' cluFiles(iShank).name])

                temp = strsplit(cluFiles(iShank).name,'.');
                shankID = str2num(temp{length(temp)}); %shankID is the spikegroup number
                clu = load(fullfile(xml.folder_path,cluFiles(iShank).name));
                clu = clu(2:end); % toss the first sample to match res/spk files
                res = load(fullfile(xml.folder_path,resFiles(iShank).name));
                spkGrpChans = all_shank_channels{iShank};

                if any(getWaveforms) && sum(clu)>0 %bug fix if no clusters 
                    nSamples = str2double(xml.spikeDetection.channelGroups.group{iShank}.nSamples.Text);
                    % load waveforms
                    chansPerSpikeGrp = length(all_shank_channels{iShank});
                    fid = fopen(fullfile(xml.folder_path,spkFiles(iShank).name),'r');
                    wav = fread(fid,[1 inf],'int16=>int16');
                    try %bug in some spk files... wrong number of samples?
                        wav = reshape(wav,chansPerSpikeGrp,nSamples,[]);
                    catch
                        if strcmp(getWaveforms,'force')
                            wav = nan(chansPerSpikeGrp,nSamples,length(clu));
                            display([spkFiles(iShank).name,' error.'])
                        else
                        error(['something is wrong with ',spkFiles(iShank).name,...
                            ' Use ''getWaveforms'', false to skip waveforms or ',...
                            '''getWaveforms'', ''force'' to write nans on bad shanks.'])
                        end
                    end
                    wav = permute(wav,[3 1 2]);




                    %% Get the DynamicTableRegion field for each shank
                    
                    % First check if the electrodes field has been filled
                    if isempty(nwb.general_extracellular_ephys_electrodes)
                        nwb = Neuroscope2NWB.getElectrodeInfo(xml, nwb);
                    end

                    electrodes_field = types.core.DynamicTableRegion('table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),'description',['shank' num2str(iShank) ' region'],'data',nwb.general_extracellular_ephys_electrodes.id.data(find(shank == iShank)'));
                    SpikeEventSeries = types.core.SpikeEventSeries('data', wav, 'electrodes', electrodes_field, 'timestamps', res./ spikes.samplingRate);

                    %% This section assigns the spike-waveforms in the .NWB
                    ecephys.nwbdatainterface.set(['SpikeEventSeries' num2str(iShank)],SpikeEventSeries);

                end


                cells  = unique(clu);
                % remove MUA and NOISE clusters...
                cells(cells==0) = [];
                cells(cells==1) = [];  % consider adding MUA as another input argument...?


                for c = 1:length(cells)
                   spikes.UID(count) = count; % this only works if all shanks are loaded... how do we optimize this?
                   ind = find(clu == cells(c));
                   spikes.times{count} = res(ind) ./ spikes.samplingRate;
                   spikes.shankID(count) = shankID;
                   spikes.cluID(count) = cells(c);

                   %Waveforms    
                   if any(getWaveforms)
                       wvforms = squeeze(mean(wav(ind,:,:)))-mean(mean(mean(wav(ind,:,:)))); % mean subtract to account for slower (theta) trends
                       if prod(size(wvforms))==length(wvforms)%in single-channel groups wvforms will squeeze too much and will have amplitude on D1 rather than D2
                           wvforms = wvforms';%fix here
                       end
                       for t = 1:size(wvforms,1)
                          [a(t) b(t)] = max(abs(wvforms(t,:))); 
                       end
                       [aa bb] = max(a,[],2);
                       spikes.rawWaveform{count} = wvforms(bb,:);
                       spikes.maxWaveformCh(count) = spkGrpChans(bb);  % Use this in Brainstorm
            %            %Regions (needs waveform peak)
            %            if isfield(xml,'region') %if there is regions field in your metadata
            %                 spikes.region{count} = 'unknown';
            %            elseif isfield(xml,'units') %if no regions, but unit region from xml via Loadparamteres
            %                 %Find the xml Unit that matches group/cluster
            %                 unitnum = cellfun(@(X,Y) X==spikes.shankID(count) && Y==spikes.cluID(count),...
            %                     {sessionInfo.Units(:).spikegroup},{sessionInfo.Units(:).cluster});
            %                 if sum(unitnum) == 0
            %                     display(['xml Missing Unit - spikegroup: ',...
            %                         num2str(spikes.shankID(count)),' cluster: ',...
            %                         num2str(spikes.cluID(count))])
            %                     spikes.region{count} = 'missingxml';
            %                 else %possible future bug: two xml units with same group/clu...              
            %                     spikes.region{count} = sessionInfo.Units(unitnum).structure;
            %                 end
            %            end
                       clear a aa b bb
                   end

                   count = count + 1;

                end

                ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
                nwb.processing.set('ecephys', ecephys);
            end


            % Serialize spiketimes and cluIDs
            spike_times       = [];
            spike_times_index = [];

            current_index = 0;
            for iNeuron = 1:length(spikes.UID)
                spike_times = [spike_times ; spikes.times{iNeuron}];
                spike_times_index = [spike_times_index; int64(length(spikes.times{iNeuron})+current_index)];
                current_index = spike_times_index(end);
            end


            % electrode_group - Assigns the group_object_view that was defined above at
            % the electrodes, to specific neurons - I need to find how each neuron is
            % assigned to a shank
            electrode_group = [];
            shank_that_neurons_belongs_to = zeros(length(spikes.UID),1);
            for iNeuron = 1:length(spikes.UID)
                shank_that_neurons_belongs_to(iNeuron) = str2double(xml.units.unit{iNeuron}.group.Text);
                first_electrode_in_shank = find(shank == shank_that_neurons_belongs_to(iNeuron));
                first_electrode_in_shank = first_electrode_in_shank(1);
                electrode_group = [electrode_group; group_object_view(first_electrode_in_shank)];
            end

            electrode_group = types.core.VectorData('data', electrode_group, 'description','the electrode group that each spike unit came from');

            % Initialize the fields needed
            spike_times       = types.core.VectorData        ('data', spike_times, 'description', 'the spike times for each unit');
            spike_times_index = types.core.VectorIndex       ('data', spike_times_index, 'target', types.untyped.ObjectView('/units/spike_times')); % The ObjectView links the indices to the spike times
            waveform_mean     = types.core.VectorData        ('data', spikes.rawWaveform, 'description', 'The mean Waveform for each unit');
            id                = types.core.ElementIdentifiers('data', [0:length(xml.units.unit)-1]');

            
            
            %% Extra Unit Info
            % FOR THE VECTORDATA, IDEALLY I NEED FILE: DG_all_6__UnitFeatureSummary_add (ACCORDING TO BEN'S CONVERTER - THIS HOLDS INFO ABOUT THE TYPE OF THE NEURONS ETC.)

            all_extra_info = types.untyped.Set();
            all_extra_info.set('cluID', types.core.VectorData('description', 'cluster ID', 'data', spikes.cluID));
            all_extra_info.set('maxWaveformCh', types.core.VectorData('description', 'The electrode where each unit showed maximum Waveform', 'data', spikes.maxWaveformCh));
            
            vectordata_units = types.core.VectorData('description','additional unit information', 'data', all_extra_info);
            

            %% 
            nwb.units = types.core.Units( ...
                'electrode_group', electrode_group, 'electrodes', [], 'electrodes_index', [], 'obs_intervals', [], 'obs_intervals_index', [], ...
                'spike_times', spike_times, 'spike_times_index', spike_times_index, 'waveform_mean', waveform_mean, 'waveform_sd', [], ...
                'colnames', {'shank_id'; 'spike_times'; 'electrode_group'; 'cell_type'; 'global_id'; 'max_electrode'}, ...
                'description', 'Generated from Neuroscope2NWB', 'id', id, 'vectordata', vectordata_units, 'vectorindex', []);

        end
        
        
        
        function nwb = getEvents(xml,nwb)
            %% Add events: nwb2.stimulus_presentation
            eventFiles = dir([xml.folder_path filesep '*.evt']);

            for iFile = 1:length(eventFiles)
                events = LoadEvents(fullfile(eventFiles(iFile).folder,eventFiles(iFile).name)); % LoadEvents is a function from the Buzcode - THIS DOESN'T LOAD ANYTHING FOR 'PulseStim_0V_10021ms_LD0' - maybe because it has a single entrty???

                if ~isempty(events.time) && ~isempty(events.description)
                    AnnotationSeries = types.core.AnnotationSeries('data',events.description,'timestamps',events.time);
                    nwb.stimulus_presentation.set(events.description{1}, AnnotationSeries);
                end
            end
        end
        
        
        
        function nwb = getBehavior(xml,nwb)
            
            %% Add behavioral data: nwb2.processing.get('behavior').nwbdatainterface
            behavior_NWB = types.core.ProcessingModule;

            
            
            % Just for the test, delete after
            xml.folder_path = 'C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge'
            
            
            
            behavioralFiles = dir([xml.folder_path filesep '*behavior.mat']);


            
            
            
            for iFile = 1:length(behavioralFiles)

                % The label of the behavior
                behavioral_Label = erase(behavioralFiles(iFile).name,'.behavior.mat');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                behavioral_Label = strsplit(behavioral_Label,'.'); % '20170505_396um_0um_merge'    'track'
                behavioral_Label = behavioral_Label{2};            % This section is hardcoded, maybe improve.
                                                                   % I assumed here that the standardized way of saving behavior files is: experimentName.BehaviorLabel.behavior.mat
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Load the file
                behaviorstruct = load(behavioralFiles(iFile).name); % The example I have has the variable "behavior" saved
                behavior       = behaviorstruct.behavior;

                
                %% This is based on the Buzcode tutorial Behavior file: 20170505_396um_0um_merge.track.behavior.mat
                %  and the Buzcode wiki: https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards#behavior
                behavioral_timestamps = behavior.timestamps;

                
                behavioral_signals_NWB = types.core.Position;
                behavior_field_names = fieldnames(behaviorstruct.behavior)';
                
                %% First add the fields that contain channels
                iChannelTypeFields = find(ismember(behavior_field_names,{'position','orientation','pupil'}));
                
                for iField = iChannelTypeFields

                    channel_fields      = fieldnames(behavior.(behavior_field_names{iField}));
                    behavioral_channels = zeros(length(behavior.(behavior_field_names{iField}).(channel_fields{1})),length(channel_fields));
                    
                    for iChannel = 1:length(channel_fields)
                        behavioral_channels(:,iChannel) = behavior.(behavior_field_names{iField}).(channel_fields{iChannel});
                    end
                        
                    spatial_series = types.core.SpatialSeries('data', behavioral_channels, 'timestamps', behavioral_timestamps,'data_unit', behavior.units);
                    behavioral_signals_NWB.spatialseries.set(behavior_field_names{iField}, spatial_series);
                end
                behavior_NWB.nwbdatainterface.set(behavioral_Label,behavioral_signals_NWB);

                
                %% Add behaviorinfo field
                behaviorinfoField = find(ismember(behavior_field_names,{'behaviorinfo'}));
                
                behaviorInfo_fields = fieldnames(behavior.(behavior_field_names{behaviorinfoField}));
                
                ErrorPerMarkerSignal = behavior.(behavior_field_names{behaviorinfoField}).errorPerMarker;
                
                
                spatial_series = types.core.SpatialSeries('data', ErrorPerMarkerSignal, 'timestamps', behavioral_timestamps,'data_unit', behavior.units, ...
                                                          'description', behavior.(behavior_field_names{behaviorinfoField}).description, 'comments', behavior.(behavior_field_names{behaviorinfoField}).acquisitionsystem);
                behavioral_signals_NWB.spatialseries.set(behavior_field_names{behaviorinfoField}, spatial_series);
                behavior_NWB.nwbdatainterface.set(behavioral_Label,behavioral_signals_NWB);
                
                
            end

            behavior_NWB.description = 'Behavioral signals';
            nwb.processing.set('behavior', behavior_NWB);

        end
        
        
        
        function nwb = getElectrophysiology(xml,nwb)
            
            lfpFile = dir([xml.folder_path filesep '*.eeg']);

            if length(lfpFile)>1
                error('More than one .eeg files are present here. Weird')
            end


            % Get the samples number, based on the size of the file
            % Check for the precision that samples are saved first

            hdr.nBits     = str2double(xml.acquisitionSystem.nBits.Text);
            hdr.nChannels = str2double(xml.acquisitionSystem.nChannels.Text);
            hdr.sRateOrig = str2double(xml.acquisitionSystem.samplingRate.Text);
            hdr.Gain      = str2double(xml.acquisitionSystem.amplification.Text);
            hdr.sRateLfp  = str2double(xml.fieldPotentials.lfpSamplingRate.Text);


            % Get data type
            switch lower(hdr.nBits)
                case 16;
                    hdr.byteSize   = 2;
                    hdr.byteFormat = 'int16';
                case 32;
                    hdr.byteSize   = 4;
                    hdr.byteFormat = 'int32';
            end
            % Guess the number of time points based on the file size
            dirInfo = dir(fullfile(lfpFile.folder,lfpFile.name));
            hdr.nSamples = floor(dirInfo.bytes ./ (hdr.nChannels * hdr.byteSize));


            
            
            
            lfp_data = bz_LoadBinary(fullfile(lfpFile.folder,lfpFile.name), 'duration',Inf, 'frequency',hdr.sRateLfp,'nchannels',hdr.nChannels, 'channels', [1:hdr.nChannels]); % nSamples x 64

            
            % If the electrode Information has not already been filled, 
            % do it now
            if isempty(nwb.general_extracellular_ephys_electrodes)
                nwb = Neuroscope2NWB.getElectrodeInfo(xml,nwb);
            end
            
            
            electrodes_field = types.core.DynamicTableRegion('table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),'description','electrode table reference','data',nwb.general_extracellular_ephys_electrodes.id.data);

            lfp = types.core.ElectricalSeries('data', lfp_data', 'electrodes',electrodes_field, 'description', 'lfp signal for all shank electrodes', 'starting_time', 0, 'starting_time_rate', hdr.sRateLfp);
            % I TRANSPOSED THE MATRIX HERE TO BE COMPATIBLE WITH THE FUNCTION BZ_GET_LFP

            LFP = types.core.LFP;
            LFP.electricalseries.set('lfp',lfp); 
            
            % Check if the ecephys field is already created
            if ~ismember(keys(nwb.processing),'ecephys')
                ecephys = types.core.ProcessingModule('description', '');
                ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
                nwb.processing.set('ecephys', ecephys);
            end
            
            nwb.processing.get('ecephys').nwbdatainterface.set('LFP', LFP);
            
        end
        
        
        function nwb = getEpochs(xml,nwb)
            % Buz function
            
            
            %% Add Epochs
            % The epochs are based on separate behavioral files
            % Each separate behavioral file = 1 epoch
            
            behavioralFiles = dir([xml.folder_path filesep '*Position*']);
            
            intervals_epochs = types.core.TimeIntervals();
            
            id_epochs  = types.core.ElementIdentifiers('data',1:length(behavioralFiles));


            start_time_epochs = zeros(length(behavioralFiles),1); % Start time
            stop_time_epochs  = zeros(length(behavioralFiles),1); % Stop time
            colnames          = cell(length(behavioralFiles),1); % Labels
            
            for iFile = 1:length(behavioralFiles)

                % The label of the behavior
                behavioral_Label = strsplit(behavioralFiles(iFile).name,'__');
                behavioral_Label = erase(behavioral_Label{2},'.mat');

                position_signals = load(fullfile(behavioralFiles(iFile).folder,behavioralFiles(iFile).name));

                % Some behavioral signal files might have more than one
                % type of signals in them (e.g. twhl_linearized, twhl_norm)
                field_names = fieldnames(position_signals);

                the_position_field_NWB = types.core.Position;
                
                for iField = 1:length(field_names)
                    the_position_field_NWB.spatialseries.set(field_names{iField}, types.core.SpatialSeries('description', [field_names{iField} ' position signals from ' behavioral_Label ' epoch'], 'data', position_signals.(field_names{iField})(:,2:end)));
                end
                
                start_time_epochs(iFile) = position_signals.(field_names{iField})(1,1);
                stop_time_epochs(iFile)  = position_signals.(field_names{iField})(end,1);
                colnames{iFile}              = behavioral_Label;

                intervals_epochs.vectordata.set(behavioral_Label,types.core.VectorData('description', ['Position signals from ' behavioral_Label ' epoch'], 'data', the_position_field_NWB));
                
            end
            
            % Transform doubles to VectorData
            start_time_epochs = types.core.VectorData('description','Starting timestamp of epochs', 'data', start_time_epochs);
            stop_time_epochs  = types.core.VectorData('description','Stopping timestamp of epochs', 'data', stop_time_epochs);
                   
            intervals_epochs.start_time  = start_time_epochs;
            intervals_epochs.stop_time   = stop_time_epochs;
            intervals_epochs.description = 'Experimental epochs';
            intervals_epochs.id          = id_epochs;
            intervals_epochs.colnames    = colnames;
            
            nwb.intervals_epochs = intervals_epochs;

        
        end
        
        
        
        
        function nwb = getTrials(xml,nwb)
            
            
            
            
%             FIGURE OUT HOW TO COMBINE BOTH CASES
%             YUTA AND NON YUTA TO LOAD FROM BZ_LOADBEHAVIOR
%             
%             
%             TRY TO AVOID STORING THE .X, .Y, .Z FIELDS, AND ONLY DEAL WITH THE CONTINUOUS SIGNAL AND THE TRIAL START AND TRIAL STOP
%             
            
            
            
            
            
            
            % Just Buzcode
            
%             % THIS SECTION LOADS THE INFO FROM *RUN.MAT - YUTA NON STANDARD FORMAT
%             % Add Trials
% 
%             % This file holds a matrix with the trial info
%             trialsFile = dir([xml.folder_path filesep '*Run.mat']);
%             trialsInfo = load(fullfile(trialsFile.folder, trialsFile.name));
%             the_field = fieldnames(trialsInfo);
%             trialsInfo = trialsInfo.(the_field{1});
%             
%             colname = the_field;
% 
%             % This file holds a cell array with the labels for the matrix above (...)
%             runFile = dir([xml.folder_path filesep '*RunInfo.mat']);
%             runInfo = load(fullfile(runFile.folder,runFile.name));
%             the_field = fieldnames(runInfo);
%             runInfo = runInfo.(the_field{1});
% 
% 
%             start_time_trials = types.core.VectorData('description','Starting timepoint of Each Trial','data',{trialsInfo(:,1)});
%             stop_time_trials  = types.core.VectorData('description','Ending timepoint of Each Trial','data',{trialsInfo(:,2)});
% 
%             id_trials = types.core.ElementIdentifiers('data',int64(0:size(trialsInfo,1)-1)');
% 
%             conditions_trials = cell(size(trialsInfo,1),1);
%             for iTrial = 1:size(trialsInfo,1)
%                 conditions_trials{iTrial} = runInfo{find(trialsInfo(iTrial,3:4))+2};
%             end
% 
%             vectordata_trials = types.untyped.Set('both_visit', trialsInfo(:,7), 'condition', conditions_trials, 'error_run', trialsInfo(:,5), 'stim_run',trialsInfo(:,6));
% 
%             
%             
%             
%             intervals_trials = types.core.TimeIntervals('start_time',start_time_trials,'stop_time',stop_time_trials,...
%                                                        'colnames',{colname},...
%                                                        'description','experimental trials','id',id_trials, 'vectordata',vectordata_trials);
% 
%             nwb.intervals_trials = intervals_trials;
        
            
            
            
            
            
            
            
            
            %% THIS SECTION LOADS THE TRIAL INFO FROM THE .BEHAVIOR.MAT FILES
            % Since the behavior.mat files can be more than one, I don't
            % fill the fields within nwb.intervals_trials
            % I use nwb.processing

            % First check if the 'behavior' field has already been created
            % (this can be done from the getBehavior function)
            if ~(sum(ismember(keys(nwb.processing), 'behavior')) == 1)%% Add behavioral data: nwb2.processing.get('behavior').nwbdatainterface
                behavior_NWB = types.core.ProcessingModule;
                behavior_NWB.description = 'Behavioral signals';
                nwb.processing.set('behavior', behavior_NWB);
            end

            
            % Just for the test, delete after
            xml.folder_path = 'C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge';
            
            
            
            behavioralFiles = dir([xml.folder_path filesep '*behavior.mat']);
            
            all_start_times = cell(length(behavioralFiles),1);
            all_stop_times  = cell(length(behavioralFiles),1);
            colnames        = cell(length(behavioralFiles),1);

    
            intervals_trials = types.core.TimeIntervals();

            for iFile = 1:length(behavioralFiles)

                % The label of the behavior
                behavioral_Label = erase(behavioralFiles(iFile).name,'.behavior.mat');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                behavioral_Label = strsplit(behavioral_Label,'.'); % '20170505_396um_0um_merge'    'track'
                behavioral_Label = behavioral_Label{2};            % This section is hardcoded, maybe improve.
                                                                   % I assumed here that the standardized way of saving behavior files is: experimentName.BehaviorLabel.behavior.mat
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                colnames{iFile} = behavioral_Label;
                
                % Load the file
                behaviorstruct = load(behavioralFiles(iFile).name); % The example I have has the variable "behavior" saved
                behavior       = behaviorstruct.behavior;

                
                
                %% This is based on the Buzcode tutorial Behavior file: 20170505_396um_0um_merge.track.behavior.mat
                %  and the Buzcode wiki: https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards#behavior

                
                
                %% Add optional events field
                if isfield(behavior, 'events')
                    
                    trials                 = behavior.events.trials;
                    map                    = behavior.events.map;
                    all_start_times{iFile} = behavior.events.trialIntervals(:,1); %types.core.VectorData('description','Starting timepoint of Each Trial','data',behavior.events.trialIntervals(:,1));
                    all_stop_times{iFile}  = behavior.events.trialIntervals(:,2); %types.core.VectorData('description','Ending timepoint of Each Trial','data',behavior.events.trialIntervals(:,2));
                  
                    all_trials = cell(length(trials),1);
                    % Assign the trials to an intervals_trials field
                    for iTrial = 1:length(trials)
                        
                        all_trial_fields = fieldnames(trials{iTrial});
                        presentFields = all_trial_fields(ismember(all_trial_fields, {'x', 'y', 'z', 'errorPerMarker', 'mapping','orientation', 'timestamps', 'direction', 'type'})');
                        
                        data_trial = types.untyped.Set();

                        for iField = 1:length(presentFields)
                            
                            % The orientation fields has subfields (other channels)
                            if strcmp(presentFields{iField}, 'orientation')
                                all_Orientation_fields = fieldnames(trials{iTrial}.(presentFields{iField}));
                                
                                data_orientation = types.untyped.Set();
                                for iOrientationField = 1:length(all_Orientation_fields)
                                    presentOrientationFields = all_Orientation_fields(ismember(all_Orientation_fields, {'x', 'y', 'z', 'w'})');

                                    data_orientation_field = types.core.VectorData('description', [all_Orientation_fields{iOrientationField} ' field from the orientation field in a single trial'], 'data', trials{iTrial}.(presentFields{iField}).(presentOrientationFields{iOrientationField}));
                                    data_orientation.set(presentOrientationFields{iOrientationField}, data_orientation_field);
                                    
                                end
                                data_field = types.core.VectorData('description', [presentFields{iField} ' field from a single trial'], 'data', data_orientation);
                                
                            else
                                data_field = types.core.VectorData('description', [presentFields{iField} ' field from a single trial'], 'data', trials{iTrial}.(presentFields{iField}));
                            end
                            
                            data_trial.set(presentFields{iField}, data_field);
                        end
                        
                        all_trials{iTrial} = data_trial;
                    end
                    
                    trials_vectordata = types.core.VectorData('description', ['Trials field from ' behavioral_Label '.behavior.mat file'], 'data', all_trials);
                    map_vectordata    = types.core.VectorData('description', ['Map field from ' behavioral_Label '.behavior.mat file'], 'data', map');
                    
                    % Combine both fields
                    trials_map_field = types.untyped.Set('trials', trials_vectordata, 'map', map_vectordata);
                    
                    intervals_trials.vectordata.set(behavioral_Label,types.core.VectorData('description', ['Trials and map field from ' behavioral_Label '.behavior.mat file'], 'data', trials_map_field));
                end
                start_time_trials = types.core.VectorData('description','Starting timepoint of Each Trial for every behavioral condition','data',all_start_times);
                stop_time_trials  = types.core.VectorData('description','Starting timepoint of Each Trial for every behavioral condition','data',all_stop_times);
                    
            end
            intervals_trials.start_time  = start_time_trials;
            intervals_trials.stop_time   = stop_time_trials;   
            intervals_trials.colnames    = colnames;
            intervals_trials.description = ['Trial info from ' num2str(length(behavioralFiles)) ' behavioral files'];
            
            nwb.intervals_trials = intervals_trials;

        end
        
        
        
        
        function nwb = special_YutaMouse_recordings(xml,nwb)
        
            %% Add raw recordings

            % value taken from Yuta's spreadsheet

            % HOW ABOUT POSITION0 - POSITION1 CHANNELS???

            hdr.nChannels = str2double(xml.acquisitionSystem.nChannels.Text);
            hdr.sRateLfp  = str2double(xml.fieldPotentials.lfpSamplingRate.Text);
            
            lfpFile = dir([xml.folder_path filesep '*.eeg']);
            
            special_electrode_labels  = {'ch_wait','ch_arm','ch_solL','ch_solR','ch_dig1','ch_dig2','ch_entL','ch_entR','ch_SsolL','ch_SsolR'};
            special_electrode_indices = [79,78,76,77,65,68,72,71,73,70]; 

            acquisition = types.untyped.Set;
            for iSpecialElectrode = 1:length(special_electrode_labels)

                special_Electrode_data = bz_LoadBinary(fullfile(lfpFile.folder,lfpFile.name), 'duration',Inf, 'frequency',hdr.sRateLfp,'nchannels',hdr.nChannels, 'channels', special_electrode_indices(iSpecialElectrode));
                single_Electrode = types.core.TimeSeries('description','environmental electrode recorded inline with neural data','data',special_Electrode_data,'starting_time', 0, 'starting_time_rate', hdr.sRateLfp, 'data_unit','V');
                acquisition.set(special_electrode_labels{iSpecialElectrode}, single_Electrode);
            end

            nwb.acquisition = acquisition;
        
        end
        
        
    end
        
end
