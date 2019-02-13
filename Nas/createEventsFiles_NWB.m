function createEventsFiles_NWB(nwb2_fileName)



    % Creates .mat events files from compatible to Buzsaki format



%     nwb2_fileName = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
    [the_path, name,~] = fileparts(nwb2_fileName);

    nwb2 = nwbRead(nwb2_fileName);
    
    
    all_event_keys = keys(nwb2.processing.get('events').nwbdatainterface);




    for iEvent = 1:length(all_event_keys)
        
        eventName = all_event_keys{iEvent};

        
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
      
        
        
        %% Exclude special characters from the behavior name
        
        clean_string = regexprep(eventName,'[^a-zA-Z0-9_]','');
        if ~strcmp(clean_string, eventName)
            disp(['The variable name (' eventName ') of the Behavior was changed to exclude special characters'])
            eventName = clean_string;
        end
        
        
        % THIS LINE HERE NAMES A NEW VARIABLE AS INFO{1}, AND COPIES THE
        % VALUES FROM behaviorName in it
        feval(@()assignin('caller',eventName, events));  % feval is needed since this is inside a function
        
        
        %% Save a separate .mat file for each behavior
        save([the_path filesep name filesep name '.' eventName '.events.mat'], eventName)
        
    end


    
end