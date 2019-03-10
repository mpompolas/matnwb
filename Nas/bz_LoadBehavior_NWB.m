function behavior = bz_LoadBehavior_NWB( nwb2,behaviorName )


%  nwb2_fileName = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
%  nwb2 = nwbRead(nwb2_fileName);

    % Check if behavior fields exists in the dataset
    try
        allBehaviorKeys = keys(nwb2.processing.get('behavior').nwbdatainterface);
        
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

    
    % Select the key that is within the Behavior selection
    allBehaviorSUBKeys = keys(nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries);    
    if length(allBehaviorSUBKeys)>1
        [iSubKeyBehavior, ~] = listdlg('PromptString','Multiple Behavior channel-selections within the selected behavior',...
                                     'ListString',allBehaviorSUBKeys,'SelectionMode','single');
    else
        iSubKeyBehavior = 1;
    end
        
    selected_behavior = nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries.get(allBehaviorSUBKeys{iSubKeyBehavior}); % This is for easier reference through the code

    
    warning('THE SAMPLING RATE IS NOT REGULAR HERE - IMPROVE')
    behavior.samplingRate = 1000 ./ nanmean(diff(selected_behavior.timestamps.load))./1000;
    behavior.units        = selected_behavior.data_unit;

    Info = allBehaviorKeys{iBehavior};
%     behavior.description  = selected_behavior.description;

    % Check if timestamps and data have values. If not something weird
    % is going on
    if ~isempty(selected_behavior.timestamps)
        behavior.timestamps = selected_behavior.timestamps.load;
    else
        behavior.timestamps     = [];
        warning(['Behavior: ' Info ' --- Timestamps are empty: weird'])
    end
    if ~isempty(selected_behavior.data)
        
        % Specifically for position signals, seperate to x and y
        if ~isempty(strfind(allBehaviorKeys(iBehavior),'position'))
            data = selected_behavior.data.load;
            datasubstruct.x = data(1,:)';
            if size(data,1)>1
                datasubstruct.y = data(2,:)';
            end
            if size(data,1)>2
                datasubstruct.z = data(3,:)';
            end
        else
            datasubstruct.data = selected_behavior.data.load;
            
        end
        
        
    else
        datasubstruct.data = [];
        warning(['Behavior: ' Info ' --- Data is empty: weird'])
    end


    % position: .x, .y, and .z
    % units: millimeters, centimeters, meters[default]
    % orientation: .x, .y, .z, and .w
    % rotationType: euler or quaternion[default]
    % pupil: .x, .y, .diameter


%       datasubstruct.reference_frame     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).reference_frame; This seems to be present only on the position_sensor channels
   
    % Conside removing these
    datasubstruct.starting_time_unit  = selected_behavior.starting_time_unit;
    datasubstruct.timestamps_interval = selected_behavior.timestamps_interval;
    datasubstruct.comments            = selected_behavior.comments;
    datasubstruct.control             = selected_behavior.control;
    datasubstruct.control_description = selected_behavior.control_description;
    datasubstruct.data_resolution     = selected_behavior.data_resolution;
    datasubstruct.starting_time       = selected_behavior.starting_time;
    datasubstruct.help                = selected_behavior.help;

    
    
    %% Exclude special characters from the behavior name
    BehaviorLabel = allBehaviorKeys{iBehavior};    
    
    clean_string = regexprep(BehaviorLabel,'[^a-zA-Z0-9_]','');
    if ~strcmp(clean_string, BehaviorLabel)
        disp(['The variable name (' BehaviorLabel ') of the Behavior was changed to exclude special characters'])
        BehaviorLabel = clean_string;
    end
    
    %% If this is a position signal, rename to position to have compatibility with the Buzsaki functions
    if ~isempty(strfind(allBehaviorKeys(iBehavior),'position'))
        BehaviorLabel = 'position';
    end
    
    
    behavior.(BehaviorLabel)    = datasubstruct;
    
    % BehaviorInfo
    behaviorinfo.description       = selected_behavior.description;
    behaviorinfo.acquisitionsystem = 'Fill Me';
    behaviorinfo.substructnames    = {BehaviorLabel};
    behavior.behaviorinfo          = behaviorinfo;


    
    %% Fill the events substructure
    
    nTrials          = length(nwb2.intervals_trials.start_time.data.load);
    uniqueConditions = unique(nwb2.intervals_trials.vectordata.get('condition').data);
    
    
    all_conditions = nwb2.intervals_trials.vectordata.get('condition').data;
    for iCondition = 1:length(all_conditions)
        
        events.trialConditions(iCondition) = find(strcmp(all_conditions{iCondition}, uniqueConditions));%1x221 double     which unique condition each trialCondition belongs to - INDEX

    end
    events.trialIntervals  = [nwb2.intervals_trials.start_time.data.load nwb2.intervals_trials.stop_time.data.load];       % 221x2 double   start-end of trial in seconds
    events.conditionType   = unique(nwb2.intervals_trials.vectordata.get('condition').data)';     % 1x10 cell   (central, wheel ...)
    
    
  
    
%     nwb2.intervals_trials.id.data.load
%     nwb2.intervals_trials.colnames
%     nwb2.intervals_trials.vectordata
%     nwb2.intervals_trials.vectordata.get('both_visit').data.load
%     nwb2.intervals_trials.vectordata.get('condition').data
%     nwb2.intervals_trials.vectordata.get('error_run').data.load
%     nwb2.intervals_trials.vectordata.get('stim_run').data.load
%     nwb2.intervals_trials.vectorindex
    

    for iTrial = 1:nTrials
        
        %% Load only the samples from each trial
        % Find the indices of the timestamps that are selected
        
        position_timestamps =  selected_behavior.timestamps.load;
        
        trialTimestampBounds = [nwb2.intervals_trials.start_time.data.load nwb2.intervals_trials.stop_time.data.load];
        [~, iPositionTimestamps] = histc(trialTimestampBounds(iTrial,:), position_timestamps);
        
        if iPositionTimestamps(1)~=0 && iPositionTimestamps(2)==0
            iPositionTimestamps(2) = length(position_timestamps);
        end

        position = selected_behavior.data.load([1, iPositionTimestamps(1)], [Inf, iPositionTimestamps(2)]);
%         clr = rand(1,3);
%         plot(position(1,:),position(2,:),'.','color',clr);
%         drawnow
                
        trials{1,iTrial}.x = position(1,:)';% 608 x 1
        
        if size(position,1)>1
            trials{1,iTrial}.y = position(2,:)';% 608 x 1
        end
        if size(position,1)>2
            trials{1,iTrial}.z = position(3,:)';% 608 x 1
        end
        
%         trials{1,iTrial}.z = 1% 608 x 1
        trials{1,iTrial}.timestamps     = position_timestamps(iPositionTimestamps(1):iPositionTimestamps(2));% 608 x 1
    %     trials{1,iTrial}.errorPerMarker = 608 x 1
        trials{1,iTrial}.direction      = 'FILL ME'; % 'clockwise' 'counterclockwise'
        trials{1,iTrial}.type           = nwb2.intervals_trials.vectordata.get('condition').data{iTrial}; % 'central alternation'
    %     trials{1,iTrial}.orientation    = % 1x1 struct
    
        % Map the trials to a 1x201 vector
        
        warning('MAPPING FIELD IS WRONG (?) - FIX')
        
        templateVector = 1:201;
        target = size(position,2);
        n_ent = length(templateVector);
        trials{1,iTrial}.mapping = round(interp1( 1:n_ent, templateVector, linspace(1, n_ent, target) ))';
    
    
    end
    
    events.trials =  trials; %   (1x221 cell)
    

    %% Create the events.map field template
    
    % FOR NOW I WILL USE THE MEDIAN OF ALL TRIALS FOR EACH POSITION
    % LINE 357 in getBehav_events
    
    map = cell(1, length(uniqueConditions));

    for iCondition = 1:length(uniqueConditions)
        iTrialsInCondition = find(strcmp(nwb2.intervals_trials.vectordata.get('condition').data, uniqueConditions{iCondition}))';

        % Initialize template vectors
        map{iCondition}.x = zeros(length(templateVector),1); % 201x1
        x = cell(length(templateVector),1);
        
        if size(position,1)>1
            map{iCondition}.y = zeros(length(templateVector),1); % 201x1
            y = cell(length(templateVector),1);
        end
        if size(position,1)>2
            map{iCondition}.z = zeros(length(templateVector),1); % 201x1
            z = cell(length(templateVector),1);
        end
        
        % Concatenate all coordinates that are assigned to each template
        % bin
        for iTrial = iTrialsInCondition
            for iCoordinate = 1:length(templateVector) % 1:201
                
                x{iCoordinate} = [x{iCoordinate} trials{1,iTrial}.x(trials{1,iTrial}.mapping == iCoordinate)'];
                if size(position,1)>1
                    y{iCoordinate} = [y{iCoordinate} trials{1,iTrial}.y(trials{1,iTrial}.mapping == iCoordinate)'];
                end
                if size(position,1)>2
                    z{iCoordinate} = [z{iCoordinate} trials{1,iTrial}.z(trials{1,iTrial}.mapping == iCoordinate)'];
                end
            end
        end
        
        % Take the median (?) for each bin
        for iCoordinate = 1:length(templateVector) % 1:201
            map{iCondition}.x(iCoordinate) = median(x{iCoordinate});
            if size(position,1)>1
                map{iCondition}.y(iCoordinate) = median(y{iCoordinate});
            end
            if size(position,1)>2
                map{iCondition}.z(iCoordinate) = median(z{iCoordinate});
            end
        end
    end
    
    events.map = map;
    behavior.events = events;
    
    %% Check that the behavior structure meets buzcode standards
    [isBehavior] = bz_isBehavior(behavior);
    switch isBehavior
        case false
            warning('Your behavior structure does not meet buzcode standards. Sad.')
    end
end