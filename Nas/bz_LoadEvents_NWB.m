function events = bz_LoadEvents_NWB(nwb2,eventsName)


%  nwb2_fileName = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
%  nwb2 = nwbRead(nwb2_fileName);

    % Check if an events field exists in the dataset
    try
        events_exist_here = ~isempty(nwb2.processing.get('events').nwbdatainterface);
        if ~events_exist_here
            disp('No events in this .nwb file')
            return
            
        else
            all_event_keys = keys(nwb2.processing.get('events').nwbdatainterface);
            disp(' ')
            disp('The following event types are present in this dataset')
            disp('------------------------------------------------')
            for iEvent = 1:length(all_event_keys)
                disp(all_event_keys{iEvent})
            end
            disp(' ')
        end
    catch
        disp('No events in this .nwb file')
        return
    end
    
    
    % Check if a specific event was called to be loaded. If not, display a
    % pop-up list with the available events for selection
    if ~exist('eventsName','var')
        [iEvent, ~] = listdlg('PromptString','Which event type would you like to load?',...
                                 'ListString',all_event_keys,'SelectionMode','single');
    else
        if sum(ismember(all_event_keys, eventsName))>0
            iEvent = find(ismember(all_event_keys, eventsName));
        end
    end
    
    events = struct;

    events.timestamps                      = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).timestamps.load;% neuroscope compatible matrix with 1-2 columns - [starts stops] (in seconds)  
    events.detectorinfo.detectorname       = 'N/A';% substructure with information about the detection method (fields below)
    events.detectorinfo.detectionparms     = 'N/A';% parameters used for detection  
    events.detectorinfo.detectiondate      = 'N/A';% date of detection
    events.detectorinfo.detectionintervals = 'N/A';% [start stop] pairs of intervals used for detection (optional)
%         events.detectorinfo.detectionchannel   = ;% channel used for detection (optional)  

%         % Optional
%         events.amplitudes  = 1;% [Nx1 matrix]
%         events.frequencies = 1;% [Nx1 matrix]
%         events.durations   = 1;% [Nx1 matrix]


    % I add here the rest of the parameters that are saved on the nwb file
    events.starting_time_unit  = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).starting_time_unit;
    events.timestamps_interval = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).timestamps_interval;
    events.timestamps_unit     = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).timestamps_unit;
    events.comments            = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).comments;
    events.control             = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).control;
    events.control_description = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).control_description;
    events.data                = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).data;
    events.data_conversion     = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).data_conversion;
    events.data_resolution     = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).data_resolution;
    events.data_unit           = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).data_unit;
    events.description         = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).description;
    events.starting_time       = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).starting_time;
    events.starting_time_rate  = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).starting_time_rate;
    events.help                = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}).help;
    
    
    
    

end
