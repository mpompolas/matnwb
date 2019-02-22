classdef DecompositionSeries < types.core.TimeSeries
% DecompositionSeries Holds spectral analysis of a timeseries. For instance of LFP or a speech signal


% PROPERTIES
properties
    bands; % A table for describing the bands that this series was generated from. There should be one row in this table for each band
    metric; % recommended: phase, amplitude, power
    source_timeseries; % HDF5 link to TimesSeries that this data was calculated from. Metadata about electrodes and their position can be read from that ElectricalSeries so it's not necessary to store that information here
end

methods
    function obj = DecompositionSeries(varargin)
        % DECOMPOSITIONSERIES Constructor for DecompositionSeries
        %     obj = DECOMPOSITIONSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % bands = DynamicTable
        % metric = char
        % source_timeseries = TimeSeries
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'bands',[]);
        addParameter(p, 'metric',[]);
        addParameter(p, 'source_timeseries',[]);
        parse(p, varargin{:});
        obj.bands = p.Results.bands;
        obj.metric = p.Results.metric;
        obj.source_timeseries = p.Results.source_timeseries;
        if strcmp(class(obj), 'types.core.DecompositionSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.bands(obj, val)
        obj.bands = obj.validate_bands(val);
    end
    function obj = set.metric(obj, val)
        obj.metric = obj.validate_metric(val);
    end
    function obj = set.source_timeseries(obj, val)
        obj.source_timeseries = obj.validate_source_timeseries(val);
    end
    %% VALIDATORS
    
    function val = validate_bands(obj, val)
        val = types.util.checkDtype('bands', 'types.core.DynamicTable', val);
    end
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf Inf Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_metric(obj, val)
        val = types.util.checkDtype('metric', 'char', val);
    end
    function val = validate_source_timeseries(obj, val)
        val = types.util.checkDtype('source_timeseries', 'types.core.TimeSeries', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.bands)
            refs = obj.bands.export(fid, [fullpath '/bands'], refs);
        else
            error('Property `bands` is required.');
        end
        if ~isempty(obj.metric)
            if startsWith(class(obj.metric), 'types.untyped.')
                refs = obj.metric.export(fid, [fullpath '/metric'], refs);
            elseif ~isempty(obj.metric)
                io.writeDataset(fid, [fullpath '/metric'], class(obj.metric), obj.metric, false);
            end
        else
            error('Property `metric` is required.');
        end
        if ~isempty(obj.source_timeseries)
            refs = obj.source_timeseries.export(fid, [fullpath '/source_timeseries'], refs);
        else
            error('Property `source_timeseries` is required.');
        end
    end
end

end