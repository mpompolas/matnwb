function nwb = addEvents(xml,nwb)
    %% Add events: nwb2.stimulus_presentation
    
    % Info about IntervalSeries here: https://nwb-schema.readthedocs.io/en/latest/format.html#sec-intervalseries
    
    % All the timestamps enter in a vectorized way (even for 2 dimensional
    % events with start and stop).
    % the start and stop flags are indicated in the .data field by a
    % positive (start) or negative (stop) value
    
    eventFiles = dir([xml.folder_path filesep '*.events.mat']);

    if ~isempty(eventFiles)
        for iFile = 1:length(eventFiles)
            events = load(fullfile(eventFiles(iFile).folder,eventFiles(iFile).name));
            
            name = fields(events); name = name{1}; % Event name (e.g. ripples)
            
            % Get the timestamps
            timestamps = events.(name).timestamps;
       
            if size(timestamps,2) == 2 % Start-stop stimulation - The events need to be saved as IntervalSeries
                IntervalSeries = types.core.IntervalSeries();
                IntervalSeries.timestamps = [timestamps(:,1) ; timestamps(:,2)]; % First event start then event stop
                IntervalSeries.data = [ones(size(timestamps,1),1) ; ones(size(timestamps,1),1)* (-1)]; % First event start then event stop - 1 signifies events start, -1 event stop
                nwb.stimulus_presentation.set(name, IntervalSeries);
            else % Stimulation Onset only - The events need to be saved as AnnotationSeries
                AnnotationSeries = types.core.AnnotationSeries('data', repmat({name},length(timestamps),1),'timestamps',timestamps);
                nwb.stimulus_presentation.set(name, AnnotationSeries);
            end
            
            % Add the detectorinfo information
            
            detectorinfo = types.untyped.Set();
            
            
            
            FINISH THE DETECOTRING
            
            
            
            detector_fields = fields(events.(name).detectorinfo);
            for iField = 1:length(detector_fields)
                if strcmp(detector_fields{iField},'detectionintervals')
                    IntervalSeries = types.core.IntervalSeries();
                    IntervalSeries.timestamps = [events.(name).detectorinfo.(detector_fields{iField})(:,1) ; events.(name).detectorinfo.(detector_fields{iField})(:,2)]; % First event start then event stop
                    IntervalSeries.data = [ones(size(events.(name).detectorinfo.(detector_fields{iField}),1),1) ; ones(size(events.(name).detectorinfo.(detector_fields{iField}),1),1)* (-1)]; % First event start then event stop - 1 signifies events start, -1 event stop
                    detectorinfo.set(detector_fields{iField}, IntervalSeries);
                else
                	detectorinfo.set(detector_fields{iField}, events.(name).detectorinfo.(detector_fields{iField}));
                end
            end
            
            nwb.stimulus_presentation.set([name '_detectorinfo'], detectorinfo)
            
        end
        disp('Events added..')
    else
        disp('No *.events.mat files found')
    end
end




