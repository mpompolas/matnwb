classdef ImageSeries < types.core.TimeSeries
% ImageSeries General image data that is common between acquisition and stimulus time series. Sometimes the image data is stored in the HDF5 file in a raw format while other times it will be stored as an external image file in the host file system. The data field will either be binary data or empty. TimeSeries::data array structure: [frame] [y][x] or [frame][z][y][x].


% PROPERTIES
properties
    dimension; % Number of pixels on x, y, (and z) axes.
    external_file; % Path or URL to one or more external file(s). Field only present if format=external. NOTE: this is only relevant if the image is stored in the file system as one or more image file(s). This field should NOT be used if the image is stored in another HDF5 file and that file is HDF5 linked to this file.
    external_file_starting_frame; % Each entry is the frame number (within the full ImageSeries) of the first frame in the corresponding external_file entry. This serves as an index to what frames each file contains, allowing random access.Zero-based indexing is used.  (The first element will always be zero).
    format; % Format of image. If this is 'external' then the field external_file contains the path or URL information to that file. For tiff, png, jpg, etc, the binary representation of the image is stored in data. If the format is raw then the fields bit_per_pixel and dimension are used. For raw images, only a single channel is stored (eg, red).
end

methods
    function obj = ImageSeries(varargin)
        % IMAGESERIES Constructor for ImageSeries
        %     obj = IMAGESERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % external_file_starting_frame = int
        % dimension = int32
        % external_file = char
        % format = char
        varargin = [{'help' 'Storage object for time-series 2-D image data'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'external_file_starting_frame',[]);
        addParameter(p, 'dimension',[]);
        addParameter(p, 'external_file',[]);
        addParameter(p, 'format',[]);
        parse(p, varargin{:});
        obj.external_file_starting_frame = p.Results.external_file_starting_frame;
        obj.dimension = p.Results.dimension;
        obj.external_file = p.Results.external_file;
        obj.format = p.Results.format;
        if strcmp(class(obj), 'types.core.ImageSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.dimension(obj, val)
        obj.dimension = obj.validate_dimension(val);
    end
    function obj = set.external_file(obj, val)
        obj.external_file = obj.validate_external_file(val);
    end
    function obj = set.external_file_starting_frame(obj, val)
        obj.external_file_starting_frame = obj.validate_external_file_starting_frame(val);
    end
    function obj = set.format(obj, val)
        obj.format = obj.validate_format(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf Inf Inf], [Inf Inf Inf Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_dimension(obj, val)
        val = types.util.checkDtype('dimension', 'int32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_external_file(obj, val)
        val = types.util.checkDtype('external_file', 'char', val);
    end
    function val = validate_external_file_starting_frame(obj, val)
        val = types.util.checkDtype('external_file_starting_frame', 'int', val);
    end
    function val = validate_format(obj, val)
        val = types.util.checkDtype('format', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.dimension)
            if startsWith(class(obj.dimension), 'types.untyped.')
                refs = obj.dimension.export(fid, [fullpath '/dimension'], refs);
            elseif ~isempty(obj.dimension)
                io.writeDataset(fid, [fullpath '/dimension'], class(obj.dimension), obj.dimension, true);
            end
        end
        if ~isempty(obj.external_file)
            if startsWith(class(obj.external_file), 'types.untyped.')
                refs = obj.external_file.export(fid, [fullpath '/external_file'], refs);
            elseif ~isempty(obj.external_file)
                io.writeDataset(fid, [fullpath '/external_file'], class(obj.external_file), obj.external_file, true);
            end
        end
        if ~isempty(obj.external_file_starting_frame) && ~isempty(obj.external_file)
            io.writeAttribute(fid, [fullpath '/external_file/starting_frame'], class(obj.external_file_starting_frame), obj.external_file_starting_frame, true);
        elseif ~isempty(obj.external_file)
            error('Property `external_file_starting_frame` is required.');
        end
        if ~isempty(obj.format)
            if startsWith(class(obj.format), 'types.untyped.')
                refs = obj.format.export(fid, [fullpath '/format'], refs);
            elseif ~isempty(obj.format)
                io.writeDataset(fid, [fullpath '/format'], class(obj.format), obj.format, false);
            end
        end
    end
end

end