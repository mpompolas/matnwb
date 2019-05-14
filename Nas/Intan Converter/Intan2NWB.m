classdef Intan2NWB
 
% This function wraps all electrophysiological, behavioral and analyzed
% data into a single NWB file.

% It was tested on the A1407_190416_113437 dataset:

% Only the folder needs to be specified. It assumes that all Intan,
% .dat and .rhd files exist within the specified folder.

% Intan has 3 ways of saving acquisition files:
% 1. One file that saves everything .rhd
% 2. One separate file per channel
% 3. One file for each channel TYPE

% The following is tested on type (3).

% Konstantinos Nasiotis 2019

    
methods(Static)
       
        function rhd = GetRHDInfo(folder_path)
            
            %% This uses an .rhd importer downloaded from:
            % http://intantech.com/files/RHD2000_MATLAB_functions_v2_01.zip

            all_files_in_folder = dir(folder_path);

            iRHD = [];
            for iFile = 1:length(all_files_in_folder)
                if strfind(all_files_in_folder(iFile).name,'.rhd')
                    iRHD = [iRHD iFile];
                end
            end
            if isempty(iRHD)
                error 'There are no .rhd files in this folder'
            elseif length(iRHD)>1
                error 'There is more than one .rhd in this folder'
            end

            rhd = read_Intan_RHD2000_file([folder_path filesep all_files_in_folder(iRHD).name]);
            
            rhd.folder_path = folder_path;
            
        end
        
        
        function nwb = GeneralInfo(rhd)
          
            %% General Info
            nwb_version = '2.0b';

            session_start_time = datetime('2019-01-01', ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
                'TimeZone', 'local');
            timestamps_reference_time = datetime('2019-01-01', ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
                'TimeZone', 'local');

%             file_create_date = datetime(datestr(clock), ...
%                 'Format', 'yyyy-MM-dd''T''HH:mm:ssZZ', ...
%                 'TimeZone', 'local');

            file_create_date = [];

            nwb = nwbfile( ...
                'session_description'          , 'Mouse in open exploration', ...
                'identifier'                   , 'PoorMouse', ...
                'session_start_time'           , session_start_time,...
                'file_create_date'             , file_create_date,...
                'general_experimenter'         , 'John',...
                'general_session_id'           , 'A1407_190416_113437',...
                'general_institution'          , 'McGill'  ,...
                'general_lab'                  , 'Peyrache',...
                'subject'                      , 'PoorMouse',...
                'general_related_publications' , '',...
                'timestamps_reference_time'    , timestamps_reference_time);

            nwb.general_subject = types.core.Subject( ...
                'description', 'mouse 5', 'genotype', 'Control', 'age', '9 months', ...
                'sex', 'M', 'subject_id', 'PoorMouse', 'species', 'Mus musculus');
        end
        
        
        
        function nwb = getElectrodeInfo (rhd,nwb)
            
            %% Get the electrodes' info
            nShanks = length(unique({rhd.amplifier_channels.port_prefix}));

            groups = unique({rhd.amplifier_channels.port_prefix}); % 'A'

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

            for iShank = 1:nShanks
                for iChannel = 1:rhd.num_amplifier_channels
                    
                    if strcmp(groups{iShank}, rhd.amplifier_channels(iChannel).port_prefix)
                        all_shank_channels{iShank} = [all_shank_channels{iShank} iShank];
                        shank_channel     = [shank_channel; iChannel-1];
                        amp_channel_id    = [amp_channel_id; rhd.amplifier_channels.native_order];
                        shank             = [shank; iShank];
                        group_name        = [group_name; ['shank' num2str(iShank)]];
                        group_object_view = [group_object_view; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(iShank)]])];

                        if ~isfield(rhd.amplifier_channels,'position') % Don't know were to get the positions from
                            x = [x; NaN];
                            y = [y; NaN];
                            z = [z; NaN];
                        end
                        if ~isfield(rhd.amplifier_channels,'imp')
                            imp = [imp; NaN];
                        end  
                        if ~isfield(rhd.amplifier_channels,'location')
                            location{end+1,1} = 'unknown';
                        end  
                        if ~isfield(rhd.amplifier_channels,'filtering')
                            filtering = [filtering; NaN];
                        end      
                    end
                end
                nwb.general_extracellular_ephys.set(['shank' num2str(iShank)], ...
                    types.core.ElectrodeGroup( ...
                    'description', ['electrode group for shank' num2str(iShank)], ...
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

            disp('Electrode info added..')
        end
        
        
        
        function nwb = getElectrophysiology(rhd,nwb)
        
            %% Load the Electrophysiological data
            nChannels = rhd.num_amplifier_channels;
            nSamples  = Inf;

% % % % % % %             % Get data
% % % % % % %             fid = fopen(fullfile(rhd.folder_path, 'amplifier.dat'), 'r');
% % % % % % %             
% % % % % % %             % CHECK IF DATA CAN BE LOADED WITH SMALLER PRECISION HERE TO SAVE SPACE
% % % % % % %             data = fread(fid, [nChannels, nSamples], 'int16');
% % % % % % %             data = data * 0.195; % Convert to microvolts
% % % % % % %             fclose(fid);
            
            % Get timestamps
            fid = fopen(fullfile(rhd.folder_path, 'time.dat'), 'r');
            time = fread(fid, 'int32');
            fclose(fid);
            
            % If the electrode Information has not already been filled, 
            % do it now
            if isempty(nwb.general_extracellular_ephys_electrodes)
                nwb = Intan2NWB.getElectrodeInfo(rhd,nwb);
            end
            
            electrodes_field = types.core.DynamicTableRegion('table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),'description','electrode table reference','data',nwb.general_extracellular_ephys_electrodes.id.data);

            
            
            
            
            
            
            
            
            
            data = [zeros(32,100)];
            
            
            
            
            
            
            
            
            % CHECK IF DATA CAN BE LOADED WITH SMALLER PRECISION HERE TO SAVE SPACE
            lfp = types.core.ElectricalSeries('data', data, 'electrodes',electrodes_field, 'description', 'raw signal for all shank electrodes', 'starting_time', 0, 'starting_time_rate', rhd.frequency_parameters.amplifier_sample_rate);

            LFP = types.core.LFP;
            LFP.electricalseries.set('lfp',lfp); 
            
            % Check if the ecephys field is already created           
            if isempty(keys(nwb.processing))
                ecephys = types.core.ProcessingModule('description', '');
                ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
                nwb.processing.set('ecephys', ecephys);
            else
                if ~ismember(keys(nwb.processing),'ecephys')
                    ecephys = types.core.ProcessingModule('description', '');
                    ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
                    nwb.processing.set('ecephys', ecephys);
                end
            end
            
            nwb.processing.get('ecephys').nwbdatainterface.set('LFP', LFP);
            disp('Electrophysiological signals added..')
            
            disp('Electrophysiology data added..')

        end
        
        
        
        function nwb = getAuxiliary(rhd,nwb)
            
            %% Load the Auxiliary data
            nChannels = length(rhd.aux_input_channels);
            nSamples  = Inf;

            fid = fopen(fullfile(rhd.folder_path, 'auxiliary.dat'), 'r');
            data = fread(fid, [nChannels, nSamples], 'int16');
            fclose(fid);
            
            data_object = types.core.SpatialSeries('data', data, 'description', 'auxiliary_signals', 'starting_time', 0, 'starting_time_rate', rhd.frequency_parameters.aux_input_sample_rate);

            nwb.acquisition.set('auxiliary', data_object);
            disp('Auxiliary data added..')

        end
        
        
        function nwb = getAnalogIn(rhd,nwb)
            
             %% Load the Auxiliary data
            nChannels = length(rhd.board_adc_channels);
            nSamples  = Inf;

            fid = fopen(fullfile(rhd.folder_path, 'analogin.dat'), 'r');
            data = fread(fid, [nChannels, nSamples], 'int16');
            fclose(fid);
            
            data_object = types.core.SpatialSeries('data', data, 'description', 'analogin signals', 'starting_time', 0, 'starting_time_rate', rhd.frequency_parameters.board_adc_sample_rate);

            nwb.acquisition.set('analogin', data_object);
            disp('AnalogIn data added..')

        end
        
        
end


end