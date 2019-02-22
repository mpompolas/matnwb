classdef SpatialSeries < types.core.TimeSeries
% SpatialSeries Direction, e.g., of gaze or travel, or position. The TimeSeries::data field is a 2D array storing position or direction relative to some reference frame. Array structure: [num measurements] [num dimensions]. Each SpatialSeries has a text dataset reference_frame that indicates the zero-position, or the zero-axes for direction. For example, if representing gaze direction, 'straight-ahead' might be a specific pixel on the monitor, or some other point in space. For position data, the 0,0 point might be the top-left corner of an enclosure, as viewed from the tracking camera. The unit of data will indicate how to interpret SpatialSeries values.


% PROPERTIES
properties
    reference_frame; % Description defining what exactly 'straight-ahead' means.
end

methods
    function obj = SpatialSeries(varargin)
        % SPATIALSERIES Constructor for SpatialSeries
        %     obj = SPATIALSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % reference_frame = char
        varargin = [{'data_unit' 'meter' 'help' 'Stores points in space over time. The data[] array structure is [num samples][num spatial dimensions]'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'reference_frame',[]);
        parse(p, varargin{:});
        obj.reference_frame = p.Results.reference_frame;
        if strcmp(class(obj), 'types.core.SpatialSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.reference_frame(obj, val)
        obj.reference_frame = obj.validate_reference_frame(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf], [Inf Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_data_unit(obj, val)
        val = types.util.checkDtype('data_unit', 'char', val);
    end
    function val = validate_reference_frame(obj, val)
        val = types.util.checkDtype('reference_frame', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.reference_frame)
            if startsWith(class(obj.reference_frame), 'types.untyped.')
                refs = obj.reference_frame.export(fid, [fullpath '/reference_frame'], refs);
            elseif ~isempty(obj.reference_frame)
                io.writeDataset(fid, [fullpath '/reference_frame'], class(obj.reference_frame), obj.reference_frame, false);
            end
        end
    end
end

end