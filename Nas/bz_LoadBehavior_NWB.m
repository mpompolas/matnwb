function behavior = bz_LoadBehavior_NWB( nwb2,behaviorName )


%  nwb2_fileName = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
%  nwb2 = nwbRead(nwb2_fileName);

    % Check if behavior fields exists in the dataset
    try
        all_raw_keys = keys(nwb2.acquisition);

        allBehaviorKeys = all_raw_keys(~ismember(all_raw_keys, {'ECoG','bla bla bla'})); % Here I exclude all the standard Raw signal keys - Make Sure this is standardized
     
        behavior_exist_here = ~isempty(allBehaviorKeys);
        if ~behavior_exist_here
            disp('No behavior in this .nwb file')
            return
            
        else
            disp(' ')
            disp('The following behavior types are present in this dataset')
            disp('------------------------------------------------')
            for iBehavior = 1:length(allBehaviorKeys)
                disp(allBehaviorKeys{iBehavior})
            end
            disp(' ')
        end
    catch
        disp('No behavior in this .nwb file')
        return
    end


    
    % Check if a specific behavior was called to be loaded. If not, display a
    % pop-up list with the available behaviors for selection
    if ~exist('behaviorName','var')
        [iBehavior, ~] = listdlg('PromptString','Which behavior type would you like to load?',...
                                 'ListString',allBehaviorKeys,'SelectionMode','single');
    else
        if sum(ismember(all_event_keys, behaviorName))>0
            iBehavior = find(ismember(allBehaviorKeys, behaviorName));
        end
    end
    
    
    
    % Fill the behavior 
    behavior = struct; % Initialize

    behavior.samplingRate       = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).starting_time_rate;

    Info = allBehaviorKeys(iBehavior);
    behavior.description        = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).description;

    % Check if timestamps and data have values. If not something weird
    % is going on
    if ~isempty(nwb2.acquisition.get(allBehaviorKeys(iBehavior)).timestamps)
        behavior.timestamps     = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).timestamps.load;
    else
        behavior.timestamps     = [];
        warning(['Behavior: ' Info{1} ' --- Timestamps are empty: weird'])
    end
    if ~isempty(nwb2.acquisition.get(allBehaviorKeys(iBehavior)).data)
        datasubstruct.data     = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).data.load;
    else
        datasubstruct.data     = [];
        warning(['Behavior: ' Info{1} ' --- Data is empty: weird'])
    end


    % position: .x, .y, and .z
    % units: millimeters, centimeters, meters[default]
    % orientation: .x, .y, .z, and .w
    % rotationType: euler or quaternion[default]
    % pupil: .x, .y, .diameter

    behaviorinfo.description          = 'Fill Me';
    behaviorinfo.acquisitionsystem    = 'Fill Me';
    behaviorinfo.substructnames       = {'data', 'reference_frame', 'starting_time_unit', 'timestamps_interval', 'comments', 'control' 'control_description', 'data_resolution', 'units', 'starting_time', 'help'};
    behavior.behaviorinfo         = behaviorinfo;

%       datasubstruct.reference_frame     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).reference_frame; This seems to be present only on the position_sensor channels
    datasubstruct.units               = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).data_unit;
    datasubstruct.starting_time_unit  = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).starting_time_unit;
    datasubstruct.timestamps_interval = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).timestamps_interval;
    datasubstruct.comments            = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).comments;
    datasubstruct.control             = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).control;
    datasubstruct.control_description = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).control_description;
    datasubstruct.data_resolution     = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).data_resolution;
    datasubstruct.starting_time       = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).starting_time;
    datasubstruct.help                = nwb2.acquisition.get(allBehaviorKeys(iBehavior)).help;

    behavior.datasubstruct        = datasubstruct;
    
    
    
    
    %Check that the behavior structure meets buzcode standards
    [isBehavior] = bz_isBehavior(behavior);
    switch isBehavior
        case false
            warning('Your behavior structure does not meet buzcode standards. Sad.')
    end
    
    
    
    
end

