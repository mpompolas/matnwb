classdef TimeSeries < types.core.NWBDataInterface
% TimeSeries General purpose time series.


% READONLY
properties(SetAccess=protected)
    starting_time_unit; % Value is 'Seconds'
    timestamps_interval; % Value is '1'
    timestamps_unit; % Value is 'Seconds'
end
% PROPERTIES
properties
    comments; % Human-readable comments about the TimeSeries. This second descriptive field can be used to store additional information, or descriptive information if the primary description field is populated with a computer-readable string.
    control; % Numerical labels that apply to each element in data[]. COMMENT: Optional field. If present, the control array should have the same number of elements as data[].
    control_description; % Description of each control value. COMMENT: Array length should be as long as the highest number in control minus one, generating an zero-based indexed array for control values.
    data; % Data values. Can also store binary data (eg, image frames) COMMENT: This field may be a link to data stored in an external file, especially in the case of raw data.
    data_conversion; % Scalar to multiply each element in data to convert it to the specified unit
    data_resolution; % Smallest meaningful difference between values in data, stored in the specified by unit. COMMENT: E.g., the change in value of the least significant bit, or a larger number if signal noise is known to be present. If unknown, use -1.0
    data_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    description; % Description of TimeSeries
    starting_time; % The timestamp of the first sample. COMMENT: When timestamps are uniformly spaced, the timestamp of the first sample can be specified and all subsequent ones calculated from the sampling rate.
    starting_time_rate; % Sampling rate, in Hz COMMENT: Rate information is stored in Hz
    timestamps; % Timestamps for samples stored in data.COMMENT: Timestamps here have all been corrected to the common experiment master-clock. Time is stored as seconds and all timestamps are relative to experiment start time.
end

methods
    function obj = TimeSeries(varargin)
        % TIMESERIES Constructor for TimeSeries
        %     obj = TIMESERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % data = any
        % data_unit = char
        % starting_time_rate = float32
        % starting_time_unit = char
        % timestamps_interval = int32
        % timestamps_unit = char
        % comments = char
        % control = uint8
        % control_description = char
        % data_conversion = float32
        % data_resolution = float32
        % description = char
        % starting_time = float64
        % timestamps = float64
        varargin = [{'comments' 'no comments' 'data_conversion' types.util.correctType(1.0, 'float32') 'data_resolution' types.util.correctType(0.0, 'float32') 'description' 'no description' 'help' 'General time series object' 'starting_time_unit' 'Seconds' 'timestamps_interval' types.util.correctType(1, 'int32') 'timestamps_unit' 'Seconds'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'data',[]);
        addParameter(p, 'data_unit',[]);
        addParameter(p, 'starting_time_rate',[]);
        addParameter(p, 'starting_time_unit',[]);
        addParameter(p, 'timestamps_interval',[]);
        addParameter(p, 'timestamps_unit',[]);
        addParameter(p, 'comments',[]);
        addParameter(p, 'control',[]);
        addParameter(p, 'control_description',[]);
        addParameter(p, 'data_conversion',[]);
        addParameter(p, 'data_resolution',[]);
        addParameter(p, 'description',[]);
        addParameter(p, 'starting_time',[]);
        addParameter(p, 'timestamps',[]);
        parse(p, varargin{:});
        obj.data = p.Results.data;
        obj.data_unit = p.Results.data_unit;
        obj.starting_time_rate = p.Results.starting_time_rate;
        obj.starting_time_unit = p.Results.starting_time_unit;
        obj.timestamps_interval = p.Results.timestamps_interval;
        obj.timestamps_unit = p.Results.timestamps_unit;
        obj.comments = p.Results.comments;
        obj.control = p.Results.control;
        obj.control_description = p.Results.control_description;
        obj.data_conversion = p.Results.data_conversion;
        obj.data_resolution = p.Results.data_resolution;
        obj.description = p.Results.description;
        obj.starting_time = p.Results.starting_time;
        obj.timestamps = p.Results.timestamps;
        if strcmp(class(obj), 'types.core.TimeSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.comments(obj, val)
        obj.comments = obj.validate_comments(val);
    end
    function obj = set.control(obj, val)
        obj.control = obj.validate_control(val);
    end
    function obj = set.control_description(obj, val)
        obj.control_description = obj.validate_control_description(val);
    end
    function obj = set.data(obj, val)
        obj.data = obj.validate_data(val);
    end
    function obj = set.data_conversion(obj, val)
        obj.data_conversion = obj.validate_data_conversion(val);
    end
    function obj = set.data_resolution(obj, val)
        obj.data_resolution = obj.validate_data_resolution(val);
    end
    function obj = set.data_unit(obj, val)
        obj.data_unit = obj.validate_data_unit(val);
    end
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.starting_time(obj, val)
        obj.starting_time = obj.validate_starting_time(val);
    end
    function obj = set.starting_time_rate(obj, val)
        obj.starting_time_rate = obj.validate_starting_time_rate(val);
    end
    function obj = set.timestamps(obj, val)
        obj.timestamps = obj.validate_timestamps(val);
    end
    %% VALIDATORS
    
    function val = validate_comments(obj, val)
        val = types.util.checkDtype('comments', 'char', val);
    end
    function val = validate_control(obj, val)
        val = types.util.checkDtype('control', 'uint8', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_control_description(obj, val)
        val = types.util.checkDtype('control_description', 'char', val);
    end
    function val = validate_data(obj, val)
    
    end
    function val = validate_data_conversion(obj, val)
        val = types.util.checkDtype('data_conversion', 'float32', val);
    end
    function val = validate_data_resolution(obj, val)
        val = types.util.checkDtype('data_resolution', 'float32', val);
    end
    function val = validate_data_unit(obj, val)
        val = types.util.checkDtype('data_unit', 'char', val);
    end
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_starting_time(obj, val)
        val = types.util.checkDtype('starting_time', 'float64', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_starting_time_rate(obj, val)
        val = types.util.checkDtype('starting_time_rate', 'float32', val);
    end
    function val = validate_timestamps(obj, val)
        val = types.util.checkDtype('timestamps', 'float64', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.comments)
            io.writeAttribute(fid, [fullpath '/comments'], class(obj.comments), obj.comments, false);
        end
        if ~isempty(obj.control)
            if startsWith(class(obj.control), 'types.untyped.')
                refs = obj.control.export(fid, [fullpath '/control'], refs);
            elseif ~isempty(obj.control)
                io.writeDataset(fid, [fullpath '/control'], class(obj.control), obj.control, true);
            end
        end
        if ~isempty(obj.control_description)
            if startsWith(class(obj.control_description), 'types.untyped.')
                refs = obj.control_description.export(fid, [fullpath '/control_description'], refs);
            elseif ~isempty(obj.control_description)
                io.writeDataset(fid, [fullpath '/control_description'], class(obj.control_description), obj.control_description, true);
            end
        end
        if ~isempty(obj.data)
            if startsWith(class(obj.data), 'types.untyped.')
                refs = obj.data.export(fid, [fullpath '/data'], refs);
            elseif ~isempty(obj.data)
                io.writeDataset(fid, [fullpath '/data'], class(obj.data), obj.data, true);
            end
        else
            error('Property `data` is required.');
        end
        if ~isempty(obj.data_conversion) && ~isempty(obj.data)
            io.writeAttribute(fid, [fullpath '/data/conversion'], class(obj.data_conversion), obj.data_conversion, false);
        end
        if ~isempty(obj.data_resolution) && ~isempty(obj.data)
            io.writeAttribute(fid, [fullpath '/data/resolution'], class(obj.data_resolution), obj.data_resolution, false);
        end
        if ~isempty(obj.data_unit) && ~isempty(obj.data)
            io.writeAttribute(fid, [fullpath '/data/unit'], class(obj.data_unit), obj.data_unit, false);
        elseif ~isempty(obj.data)
            error('Property `data_unit` is required.');
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], class(obj.description), obj.description, false);
        end
        if ~isempty(obj.starting_time)
            if startsWith(class(obj.starting_time), 'types.untyped.')
                refs = obj.starting_time.export(fid, [fullpath '/starting_time'], refs);
            elseif ~isempty(obj.starting_time)
                io.writeDataset(fid, [fullpath '/starting_time'], class(obj.starting_time), obj.starting_time, false);
            end
        end
        if ~isempty(obj.starting_time_rate) && ~isempty(obj.starting_time)
            io.writeAttribute(fid, [fullpath '/starting_time/rate'], class(obj.starting_time_rate), obj.starting_time_rate, false);
        elseif ~isempty(obj.starting_time)
            error('Property `starting_time_rate` is required.');
        end
        if ~isempty(obj.starting_time_unit) && ~isempty(obj.starting_time)
            io.writeAttribute(fid, [fullpath '/starting_time/unit'], class(obj.starting_time_unit), obj.starting_time_unit, false);
        elseif ~isempty(obj.starting_time)
            error('Property `starting_time_unit` is required.');
        end
        if ~isempty(obj.timestamps)
            if startsWith(class(obj.timestamps), 'types.untyped.')
                refs = obj.timestamps.export(fid, [fullpath '/timestamps'], refs);
            elseif ~isempty(obj.timestamps)
                io.writeDataset(fid, [fullpath '/timestamps'], class(obj.timestamps), obj.timestamps, true);
            end
        end
        if ~isempty(obj.timestamps_interval) && ~isempty(obj.timestamps)
            io.writeAttribute(fid, [fullpath '/timestamps/interval'], class(obj.timestamps_interval), obj.timestamps_interval, false);
        elseif ~isempty(obj.timestamps)
            error('Property `timestamps_interval` is required.');
        end
        if ~isempty(obj.timestamps_unit) && ~isempty(obj.timestamps)
            io.writeAttribute(fid, [fullpath '/timestamps/unit'], class(obj.timestamps_unit), obj.timestamps_unit, false);
        elseif ~isempty(obj.timestamps)
            error('Property `timestamps_unit` is required.');
        end
    end
end

end