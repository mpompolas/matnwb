classdef EventDetection < types.core.NWBDataInterface
% EventDetection Detected spike events from voltage trace(s).


% PROPERTIES
properties
    detection_method; % Description of how events were detected, such as voltage threshold, or dV/dT threshold, as well as relevant values.
    source_electricalseries; % HDF5 link to ElectricalSeries that this data was calculated from. Metadata about electrodes and their position can be read from that ElectricalSeries so it's not necessary to mandate that information be stored here
    source_idx; % Indices (zero-based) into source ElectricalSeries::data array corresponding to time of event. Module description should define what is meant by time of event (e.g., .25msec before action potential peak, zero-crossing time, etc). The index points to each event from the raw data
    times; % Timestamps of events, in Seconds
    times_unit; % The string ''Seconds''
end

methods
    function obj = EventDetection(varargin)
        % EVENTDETECTION Constructor for EventDetection
        %     obj = EVENTDETECTION(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % detection_method = char
        % source_electricalseries = ElectricalSeries
        % source_idx = int32
        % times = float64
        % times_unit = char
        varargin = [{'help' 'Detected spike events from voltage trace(s)' 'times_unit' 'Seconds'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'detection_method',[]);
        addParameter(p, 'source_electricalseries',[]);
        addParameter(p, 'source_idx',[]);
        addParameter(p, 'times',[]);
        addParameter(p, 'times_unit',[]);
        parse(p, varargin{:});
        obj.detection_method = p.Results.detection_method;
        obj.source_electricalseries = p.Results.source_electricalseries;
        obj.source_idx = p.Results.source_idx;
        obj.times = p.Results.times;
        obj.times_unit = p.Results.times_unit;
        if strcmp(class(obj), 'types.core.EventDetection')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.detection_method(obj, val)
        obj.detection_method = obj.validate_detection_method(val);
    end
    function obj = set.source_electricalseries(obj, val)
        obj.source_electricalseries = obj.validate_source_electricalseries(val);
    end
    function obj = set.source_idx(obj, val)
        obj.source_idx = obj.validate_source_idx(val);
    end
    function obj = set.times(obj, val)
        obj.times = obj.validate_times(val);
    end
    function obj = set.times_unit(obj, val)
        obj.times_unit = obj.validate_times_unit(val);
    end
    %% VALIDATORS
    
    function val = validate_detection_method(obj, val)
        val = types.util.checkDtype('detection_method', 'char', val);
    end
    function val = validate_source_electricalseries(obj, val)
        val = types.util.checkDtype('source_electricalseries', 'types.core.ElectricalSeries', val);
    end
    function val = validate_source_idx(obj, val)
        val = types.util.checkDtype('source_idx', 'int32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_times(obj, val)
        val = types.util.checkDtype('times', 'float64', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_times_unit(obj, val)
        val = types.util.checkDtype('times_unit', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.detection_method)
            if startsWith(class(obj.detection_method), 'types.untyped.')
                refs = obj.detection_method.export(fid, [fullpath '/detection_method'], refs);
            elseif ~isempty(obj.detection_method)
                io.writeDataset(fid, [fullpath '/detection_method'], class(obj.detection_method), obj.detection_method, false);
            end
        else
            error('Property `detection_method` is required.');
        end
        if ~isempty(obj.source_electricalseries)
            refs = obj.source_electricalseries.export(fid, [fullpath '/source_electricalseries'], refs);
        else
            error('Property `source_electricalseries` is required.');
        end
        if ~isempty(obj.source_idx)
            if startsWith(class(obj.source_idx), 'types.untyped.')
                refs = obj.source_idx.export(fid, [fullpath '/source_idx'], refs);
            elseif ~isempty(obj.source_idx)
                io.writeDataset(fid, [fullpath '/source_idx'], class(obj.source_idx), obj.source_idx, true);
            end
        else
            error('Property `source_idx` is required.');
        end
        if ~isempty(obj.times)
            if startsWith(class(obj.times), 'types.untyped.')
                refs = obj.times.export(fid, [fullpath '/times'], refs);
            elseif ~isempty(obj.times)
                io.writeDataset(fid, [fullpath '/times'], class(obj.times), obj.times, true);
            end
        else
            error('Property `times` is required.');
        end
        if ~isempty(obj.times_unit) && ~isempty(obj.times)
            io.writeAttribute(fid, [fullpath '/times/unit'], class(obj.times_unit), obj.times_unit, false);
        end
    end
end

end