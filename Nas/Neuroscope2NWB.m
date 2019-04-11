% Convert Neuroscope to NWB

% Konstantinos Nasiotis 2019


folder_path = 'F:\NWBtoBuzcode\YutaMouse41-150903';


%% Enter the folder and run everything from in there (easier paths).
%  When the converter is done, go back to the initial folder

[previous_paths, name] = fileparts(folder_path);



current_folder = pwd;
cd (folder_path)

%% This uses an .xml importer downloaded from MathWorks - File Exchange
%  The loadxml from the Buzcode repo doesn't work on my datasets
%  https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

all_files_in_folder = dir(folder_path);

iXML = [];
for iFile = 1:length(all_files_in_folder)
    if strfind(all_files_in_folder(iFile).name,'.xml');
        iXML = [iXML iFile];
    end
end
if isempty(iXML)
    error 'There are no .xml files in this folder'
elseif length(iXML)>1
    error 'There is more than one .xml in this folder'
end


xml = xml2struct(all_files_in_folder(iXML).name);
xml = xml.parameters;

%% Get the celltypes
code = [0 1 2 3 4 5 6 8 9 10];
type = {'unknown', ...
        'granule cells (DG) or pyramidal cells (CA3)  (need to use region info. see below.)', ...
        'mossy cell', ...
        'narrow waveform cell', ...
        'optogenetically tagged SST cell', ...
        'wide waveform cell (narrower, exclude opto tagged SST cell)', ...
        'wide waveform cell (wider)', ...
        'positive waveform unit (non-bursty)', ...
        'positive waveform unit (bursty)', ...
        'positive negative waveform unit'
       };
   
for iCellType = 1:length(code)
    celltype_dict(iCellType).Type = type(iCellType);
    celltype_dict(iCellType).Code = code(iCellType);
end

%% Values taken from Yuta's spreadsheet                      
                      
code = [79,78,76,77,65,68,72,71,73,70];
type = {'ch_wait', ...
        'ch_arm', ...
        'ch_solL cell', ...
        'ch_solR', ...
        'ch_dig1', ...
        'ch_dig2', ...
        'ch_entL', ...
        'ch_entR', ...
        'ch_SsolL', ...
        'ch_SsolR'
       };
    
for iSpecialElectrode = 1:length(code)
    SpecialElectrode_dict(iSpecialElectrode).Type = type(iSpecialElectrode);
    SpecialElectrode_dict(iSpecialElectrode).Code = code(iSpecialElectrode);
end

%% Conversion values taken from methods of Senzai, Buzsaki 2017

conversion = [[], 0.46, 0.46, 0.46, 0.46, 0.46, 0.65/2];
type = {'OpenFieldPosition_ExtraLarge', ...
        'OpenFieldPosition_New_Curtain', ...
        'OpenFieldPosition_New', ...
        'OpenFieldPosition_Old_Curtain', ...
        'OpenFieldPosition_Old', ...
        'OpenFieldPosition_Oldlast', ...
        'EightMazePosition'
       };
    
for iType = 1:length(conversion)
    task_types(iType).Type       = type(iType);
    task_types(iType).conversion = conversion(iType);
end



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
    'session_description'  , 'Mouse in open exploration and theta maze', ...
    'identifier'           , name, ...
    'session_start_time'   , session_start_time,...
    'file_create_date'     , file_create_date,...
    'general_experimenter' , xml.generalInfo.experimenters.Text,...
    'general_session_id'   , name,...
    'institution'          , 'NYU'  ,...
    'lab'                  , 'Buzsaki',...
    'subject'              , 'YutaMouse',...
    'related_publications' , 'DOI:10.1016/j.neuron.2016.12.011');





nwb.general_subject = types.core.Subject( ...
    'description', 'mouse 5', 'age', '9 months', ...
    'sex', 'M', 'species', 'Mus musculus');


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

variables = {'x'; 'y'; 'z'; 'imp'; 'location'; 'filtering'; 'group'; 'group_name'; 'shank'; 'shank_channel'; 'amp_channel_id'};

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



%% Add the units info (copied from bz_GetSpikes)

getWaveforms = 0; % Set this to true if you want to add waveforms on the NWB file


spikes.samplingRate = str2double(xml.acquisitionSystem.samplingRate.Text);


disp('loading spikes from clu/res/spk files..')
% find res/clu/fet/spk files here
cluFiles = dir([folder_path filesep '*.clu*']);  
resFiles = dir([folder_path filesep '*.res*']);
if any(getWaveforms)
    spkFiles = dir([folder_path filesep '*.spk*']);
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

for iShank=1:length(cluFiles) 
    disp(['working on ' cluFiles(iShank).name])
    
    temp = strsplit(cluFiles(iShank).name,'.');
    shankID = str2num(temp{length(temp)}); %shankID is the spikegroup number
    clu = load(fullfile(folder_path,cluFiles(iShank).name));
    clu = clu(2:end); % toss the first sample to match res/spk files
    res = load(fullfile(folder_path,resFiles(iShank).name));
    spkGrpChans = all_shank_channels{iShank};
    
    if any(getWaveforms) && sum(clu)>0 %bug fix if no clusters 
        nSamples = str2double(xml.spikeDetection.channelGroups.group{iShank}.nSamples.Text);
        % load waveforms
        chansPerSpikeGrp = length(all_shank_channels{iShank});
        fid = fopen(fullfile(folder_path,spkFiles(iShank).name),'r');
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
           spikes.maxWaveformCh(count) = spkGrpChans(bb);  
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
end

    

% Serialize spiketimes and cluIDs
spike_times       = [];
spike_times_index = [];

current_index = 0;
for iNeuron = 1:length(spikes.UID)
    spike_times = [spike_times ; spikes.times{iNeuron}];
    spike_times_index = [spike_times_index; length(spikes.times{iNeuron})+current_index];
    current_index = spike_times_index + spike_times_index(end);
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
spike_times_index = types.core.VectorIndex       ('data', spike_times_index);
id                = types.core.ElementIdentifiers('data', [0:length(xml.units.unit)-1]');
% electrode_group   = types.core.VectorData        ('data',electrode_group, 'description', 'the electrode group that each spike unit came from'); % THIS NEEDS TO BE AN OBJECTVIEW OBJECT



% First, instantiate the table, listing all of the columns that will be
% added and us the |'id'| argument to indicate the number of rows. Ifa
% value is indexed, only the column name is included, not the index. For
% instance, |'spike_times_index'| is not added to the array.
nwb.units = types.core.Units( ...
    'electrode_group', electrode_group, 'electrodes', [], 'electrodes_index', [], 'obs_intervals', [], 'obs_intervals_index', [], ...
    'spike_times', spike_times, 'spike_times_index', spike_times_index, 'waveform_mean', [], 'waveform_sd', [], ...
    'colnames', {'shank_id'; 'spike_times'; 'electrode_group'; 'cell_type'; 'global_id'; 'max_electrode'}, ...
    'description', 'Generated from Neuroscope2NWB', 'id', id, 'vectordata', [], 'vectorindex', []);

% Then you can add the data column-by-column:
waveform_mean = types.core.VectorData('data', ones(30, 3), ...
    'description', 'mean of waveform');
nwb.units.waveform_mean = waveform_mean;

quality = types.core.VectorData('data', [.9, .1, .2],...
    'description', 'sorting quality score out of 1');
nwb.units.vectordata.set('quality', quality);

spike_times_cells = {[0.1, 0.21, 0.5, 0.61], [0.34, 0.36, 0.66, 0.69], [0.4, 0.43]};
[spike_times_vector, spike_times_index] = util.create_indexed_column( ...
    spike_times_cells, '/units/spike_times');
nwb.units.spike_times = spike_times_vector;
nwb.units.spike_times_index = spike_times_index;

[electrodes, electrodes_index] = util.create_indexed_column( ...
    {[0,1], [1,2], [2,3]}, '/units/electrodes', [], [], ...
    electrodes_object_view);
nwb.units.electrodes = electrodes;
nwb.units.electrodes_index = electrodes_index;






    
    
    






%% Export to nwb

cd (current_folder)

nwbExport(nwb, 'YutaMouse41_converted.nwb')

nwb3 = nwbRead('YutaMouse41_converted.nwb');

