function createBehaviorFiles_NWB(nwb2_fileName)


    % Creates the ~.behavior.mat files that can be loaded later from bz_LoadBehavior
    % straight from the .nwb files
    
    % Konstantinos Nasiotis - 2019

    nwb2_fileName = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
    [the_path, name,~] = fileparts(nwb2_fileName);

    nwb2 = nwbRead(nwb2_fileName);
    
    
    all_raw_keys = keys(nwb2.acquisition);

    
    allBehaviorKeys = [];
    for iKey = 1:length(all_raw_keys)
        singleBehaviorKey = all_raw_keys(~ismember(all_raw_keys{iKey}, {'ECoG','bla bla bla'})); % Here I exclude all the standard Raw signal keys - Make Sure this is standardized
        allBehaviorKeys = [allBehaviorKeys iKey];
    end


    for iBehavior = 1:length(allBehaviorKeys)
        
        behaviorName = struct; % Initialize

        behaviorName.samplingRate       = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).starting_time_rate;
        
        Info = all_raw_keys(allBehaviorKeys(iBehavior));
        behaviorName.description        = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).description;
                
        % Check if timestamps and data have values. If not something weird
        % is going on
        if ~isempty(nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).timestamps)
            behaviorName.timestamps     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).timestamps.load;
        else
            behaviorName.timestamps     = [];
            warning(['Behavior: ' Info{1} ' --- Timestamps are empty: weird'])
        end
        if ~isempty(nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).data)
            datasubstruct.data     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).data.load;
        else
            datasubstruct.data     = [];
            warning(['Behavior: ' Info{1} ' --- Data is empty: weird'])
        end

        
        % position: .x, .y, and .z
        % units: millimeters, centimeters, meters[default]
        % orientation: .x, .y, .z, and .w
        % rotationType: euler or quaternion[default]
        % pupil: .x, .y, .diameter
        
        
        
         %% Exclude special characters from the behavior name
        
        clean_string = regexprep(Info{1},'[^a-zA-Z0-9_]','');
        if ~strcmp(clean_string, Info{1})
            disp(['The variable name (' Info{1} ') of the Behavior was changed to exclude special characters'])
            Info{1} = clean_string;
        end
        
        behaviorinfo.description          = 'Fill Me';
        behaviorinfo.acquisitionsystem    = 'Fill Me';
        behaviorinfo.substructnames       = {'data', 'reference_frame', 'starting_time_unit', 'timestamps_interval', 'comments', 'control' 'control_description', 'data_resolution', 'units', 'starting_time', 'help'};
        behaviorName.behaviorinfo         = behaviorinfo;
        
%       datasubstruct.reference_frame     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).reference_frame; This seems to be present only on the position_sensor channels
        datasubstruct.units               = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).data_unit;
        datasubstruct.starting_time_unit  = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).starting_time_unit;
        datasubstruct.timestamps_interval = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).timestamps_interval;
        datasubstruct.comments            = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).comments;
        datasubstruct.control             = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).control;
        datasubstruct.control_description = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).control_description;
        datasubstruct.data_resolution     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).data_resolution;
        datasubstruct.starting_time       = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).starting_time;
        datasubstruct.help                = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).help;
        
        behaviorName.datasubstruct        = datasubstruct;
        
       

        
        % THIS LINE HERE NAMES A NEW VARIABLE AS INFO{1}, AND COPIES THE
        % VALUES FROM behaviorName in it        
        feval(@()assignin('caller',Info{1}, behaviorName)); % feval is needed since this is inside a function

        
        
        %% Save a separate .mat file for each behavior
        save([the_path filesep name filesep name '.' Info{1} '.behavior.mat'], Info{1})
        
    end



end