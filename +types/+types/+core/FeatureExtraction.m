classdef FeatureExtraction < types.core.NWBDataInterface
% FeatureExtraction Features, such as PC1 and PC2, that are extracted from signals stored in a SpikeEvent TimeSeries or other source.


% PROPERTIES
properties
    description; % Description of features (eg, ''PC1'') for each of the extracted features.
    electrodes; % the electrodes that this series was generated from
    features; % Multi-dimensional array of features extracted from each event.
    times; % Times of events that features correspond to (can be a link).
end

methods
    function obj = FeatureExtraction(varargin)
        % FEATUREEXTRACTION Constructor for FeatureExtraction
        %     obj = FEATUREEXTRACTION(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % electrodes = DynamicTableRegion
        % features = float32
        % times = float64
        varargin = [{'help' 'Container for salient features of detected events'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'electrodes',[]);
        addParameter(p, 'features',[]);
        addParameter(p, 'times',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        obj.electrodes = p.Results.electrodes;
        obj.features = p.Results.features;
        obj.times = p.Results.times;
        if strcmp(class(obj), 'types.core.FeatureExtraction')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.electrodes(obj, val)
        obj.electrodes = obj.validate_electrodes(val);
    end
    function obj = set.features(obj, val)
        obj.features = obj.validate_features(val);
    end
    function obj = set.times(obj, val)
        obj.times = obj.validate_times(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_electrodes(obj, val)
        val = types.util.checkDtype('electrodes', 'types.core.DynamicTableRegion', val);
    end
    function val = validate_features(obj, val)
        val = types.util.checkDtype('features', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf Inf Inf]};
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
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            if startsWith(class(obj.description), 'types.untyped.')
                refs = obj.description.export(fid, [fullpath '/description'], refs);
            elseif ~isempty(obj.description)
                io.writeDataset(fid, [fullpath '/description'], class(obj.description), obj.description, true);
            end
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.electrodes)
            refs = obj.electrodes.export(fid, [fullpath '/electrodes'], refs);
        else
            error('Property `electrodes` is required.');
        end
        if ~isempty(obj.features)
            if startsWith(class(obj.features), 'types.untyped.')
                refs = obj.features.export(fid, [fullpath '/features'], refs);
            elseif ~isempty(obj.features)
                io.writeDataset(fid, [fullpath '/features'], class(obj.features), obj.features, true);
            end
        else
            error('Property `features` is required.');
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
    end
end

end