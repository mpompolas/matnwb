classdef IndexSeries < types.core.TimeSeries
% IndexSeries Stores indices to image frames stored in an ImageSeries. The purpose of the ImageIndexSeries is to allow a static image stack to be stored somewhere, and the images in the stack to be referenced out-of-order. This can be for the display of individual images, or of movie segments (as a movie is simply a series of images). The data field stores the index of the frame in the referenced ImageSeries, and the timestamps array indicates when that image was displayed.


% PROPERTIES
properties
    indexed_timeseries; % HDF5 link to TimeSeries containing images that are indexed.
end

methods
    function obj = IndexSeries(varargin)
        % INDEXSERIES Constructor for IndexSeries
        %     obj = INDEXSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % indexed_timeseries = ImageSeries
        varargin = [{'help' 'A sequence that is generated from an existing image stack. Frames can be presented in an arbitrary order. The data[] field stores frame number in reference stack'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'indexed_timeseries',[]);
        parse(p, varargin{:});
        obj.indexed_timeseries = p.Results.indexed_timeseries;
        if strcmp(class(obj), 'types.core.IndexSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.indexed_timeseries(obj, val)
        obj.indexed_timeseries = obj.validate_indexed_timeseries(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'int', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_indexed_timeseries(obj, val)
        val = types.util.checkDtype('indexed_timeseries', 'types.core.ImageSeries', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.indexed_timeseries)
            refs = obj.indexed_timeseries.export(fid, [fullpath '/indexed_timeseries'], refs);
        else
            error('Property `indexed_timeseries` is required.');
        end
    end
end

end