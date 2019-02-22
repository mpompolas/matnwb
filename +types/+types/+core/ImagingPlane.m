classdef ImagingPlane < types.core.NWBContainer
% ImagingPlane One of possibly many groups describing an imaging plane. COMMENT: Name is arbitrary but should be meaningful. It is referenced by TwoPhotonSeries and also ImageSegmentation and DfOverF interfaces


% PROPERTIES
properties
    description; % Description of image_plane_X
    device; % the device that was used to record from this electrode
    excitation_lambda; % Excitation wavelength in nm
    imaging_rate; % Rate images are acquired, in Hz.
    indicator; % Calcium indicator
    location; % Location of image plane
    manifold; % Physical position of each pixel. COMMENT: 'xyz' represents the position of the pixel relative to the defined coordinate space
    manifold_conversion; % Multiplier to get from stored values to specified unit (e.g., 1e-3 for millimeters)
    manifold_unit; % Base unit that coordinates are stored in (e.g., Meters)
    opticalchannel; % One of possibly many groups storing channel-specific data COMMENT: Name is arbitrary but should be meaningful
    reference_frame; % Describes position and reference frame of manifold based on position of first element in manifold. For example, text description of anotomical location or vectors needed to rotate to common anotomical axis (eg, AP/DV/ML). COMMENT: This field is necessary to interpret manifold. If manifold is not present then this field is not required
end

methods
    function obj = ImagingPlane(varargin)
        % IMAGINGPLANE Constructor for ImagingPlane
        %     obj = IMAGINGPLANE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % device = Device
        % excitation_lambda = float
        % imaging_rate = float
        % indicator = char
        % location = char
        % opticalchannel = OpticalChannel
        % description = char
        % manifold = float32
        % manifold_conversion = float
        % manifold_unit = char
        % reference_frame = char
        varargin = [{'help' 'Metadata about an imaging plane' 'manifold_conversion' types.util.correctType(1.0, 'float') 'manifold_unit' 'Meter'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        
        [obj.opticalchannel,ivarargin] = types.util.parseAnon('types.core.OpticalChannel', varargin{:});
        varargin(ivarargin) = [];
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'device',[]);
        addParameter(p, 'excitation_lambda',[]);
        addParameter(p, 'imaging_rate',[]);
        addParameter(p, 'indicator',[]);
        addParameter(p, 'location',[]);
        addParameter(p, 'description',[]);
        addParameter(p, 'manifold',[]);
        addParameter(p, 'manifold_conversion',[]);
        addParameter(p, 'manifold_unit',[]);
        addParameter(p, 'reference_frame',[]);
        parse(p, varargin{:});
        obj.device = p.Results.device;
        obj.excitation_lambda = p.Results.excitation_lambda;
        obj.imaging_rate = p.Results.imaging_rate;
        obj.indicator = p.Results.indicator;
        obj.location = p.Results.location;
        obj.description = p.Results.description;
        obj.manifold = p.Results.manifold;
        obj.manifold_conversion = p.Results.manifold_conversion;
        obj.manifold_unit = p.Results.manifold_unit;
        obj.reference_frame = p.Results.reference_frame;
        if strcmp(class(obj), 'types.core.ImagingPlane')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.device(obj, val)
        obj.device = obj.validate_device(val);
    end
    function obj = set.excitation_lambda(obj, val)
        obj.excitation_lambda = obj.validate_excitation_lambda(val);
    end
    function obj = set.imaging_rate(obj, val)
        obj.imaging_rate = obj.validate_imaging_rate(val);
    end
    function obj = set.indicator(obj, val)
        obj.indicator = obj.validate_indicator(val);
    end
    function obj = set.location(obj, val)
        obj.location = obj.validate_location(val);
    end
    function obj = set.manifold(obj, val)
        obj.manifold = obj.validate_manifold(val);
    end
    function obj = set.manifold_conversion(obj, val)
        obj.manifold_conversion = obj.validate_manifold_conversion(val);
    end
    function obj = set.manifold_unit(obj, val)
        obj.manifold_unit = obj.validate_manifold_unit(val);
    end
    function obj = set.opticalchannel(obj, val)
        obj.opticalchannel = obj.validate_opticalchannel(val);
    end
    function obj = set.reference_frame(obj, val)
        obj.reference_frame = obj.validate_reference_frame(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_device(obj, val)
        val = types.util.checkDtype('device', 'types.core.Device', val);
    end
    function val = validate_excitation_lambda(obj, val)
        val = types.util.checkDtype('excitation_lambda', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_imaging_rate(obj, val)
        val = types.util.checkDtype('imaging_rate', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_indicator(obj, val)
        val = types.util.checkDtype('indicator', 'char', val);
    end
    function val = validate_location(obj, val)
        val = types.util.checkDtype('location', 'char', val);
    end
    function val = validate_manifold(obj, val)
        val = types.util.checkDtype('manifold', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[3 Inf Inf], [3 Inf Inf Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_manifold_conversion(obj, val)
        val = types.util.checkDtype('manifold_conversion', 'float', val);
    end
    function val = validate_manifold_unit(obj, val)
        val = types.util.checkDtype('manifold_unit', 'char', val);
    end
    function val = validate_opticalchannel(obj, val)
        val = types.util.checkDtype('opticalchannel', 'types.core.OpticalChannel', val);
    end
    function val = validate_reference_frame(obj, val)
        val = types.util.checkDtype('reference_frame', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBContainer(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            if startsWith(class(obj.description), 'types.untyped.')
                refs = obj.description.export(fid, [fullpath '/description'], refs);
            elseif ~isempty(obj.description)
                io.writeDataset(fid, [fullpath '/description'], class(obj.description), obj.description, false);
            end
        end
        if ~isempty(obj.device)
            refs = obj.device.export(fid, [fullpath '/device'], refs);
        else
            error('Property `device` is required.');
        end
        if ~isempty(obj.excitation_lambda)
            if startsWith(class(obj.excitation_lambda), 'types.untyped.')
                refs = obj.excitation_lambda.export(fid, [fullpath '/excitation_lambda'], refs);
            elseif ~isempty(obj.excitation_lambda)
                io.writeDataset(fid, [fullpath '/excitation_lambda'], class(obj.excitation_lambda), obj.excitation_lambda, false);
            end
        else
            error('Property `excitation_lambda` is required.');
        end
        if ~isempty(obj.imaging_rate)
            if startsWith(class(obj.imaging_rate), 'types.untyped.')
                refs = obj.imaging_rate.export(fid, [fullpath '/imaging_rate'], refs);
            elseif ~isempty(obj.imaging_rate)
                io.writeDataset(fid, [fullpath '/imaging_rate'], class(obj.imaging_rate), obj.imaging_rate, false);
            end
        else
            error('Property `imaging_rate` is required.');
        end
        if ~isempty(obj.indicator)
            if startsWith(class(obj.indicator), 'types.untyped.')
                refs = obj.indicator.export(fid, [fullpath '/indicator'], refs);
            elseif ~isempty(obj.indicator)
                io.writeDataset(fid, [fullpath '/indicator'], class(obj.indicator), obj.indicator, false);
            end
        else
            error('Property `indicator` is required.');
        end
        if ~isempty(obj.location)
            if startsWith(class(obj.location), 'types.untyped.')
                refs = obj.location.export(fid, [fullpath '/location'], refs);
            elseif ~isempty(obj.location)
                io.writeDataset(fid, [fullpath '/location'], class(obj.location), obj.location, false);
            end
        else
            error('Property `location` is required.');
        end
        if ~isempty(obj.manifold)
            if startsWith(class(obj.manifold), 'types.untyped.')
                refs = obj.manifold.export(fid, [fullpath '/manifold'], refs);
            elseif ~isempty(obj.manifold)
                io.writeDataset(fid, [fullpath '/manifold'], class(obj.manifold), obj.manifold, true);
            end
        end
        if ~isempty(obj.manifold_conversion) && ~isempty(obj.manifold)
            io.writeAttribute(fid, [fullpath '/manifold/conversion'], class(obj.manifold_conversion), obj.manifold_conversion, false);
        end
        if ~isempty(obj.manifold_unit) && ~isempty(obj.manifold)
            io.writeAttribute(fid, [fullpath '/manifold/unit'], class(obj.manifold_unit), obj.manifold_unit, false);
        end
        if ~isempty(obj.opticalchannel)
            refs = obj.opticalchannel.export(fid, [fullpath '/'], refs);
        else
            error('Property `opticalchannel` is required.');
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