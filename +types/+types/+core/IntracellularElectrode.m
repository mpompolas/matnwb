classdef IntracellularElectrode < types.core.NWBContainer
% IntracellularElectrode One of possibly many. COMMENT: Name should be informative.


% PROPERTIES
properties
    description; % Recording description, description of electrode (e.g.,  whole-cell, sharp, etc)COMMENT: Free-form text (can be from Methods)
    device; % the device that was used to record from this electrode
    filtering; % Electrode specific filtering.
    initial_access_resistance; % Initial access resistance
    location; % Area, layer, comments on estimation, stereotaxis coordinates (if in vivo, etc)
    resistance; % Electrode resistance COMMENT: unit: Ohm
    seal; % Information about seal used for recording
    slice; % Information about slice used for recording
end

methods
    function obj = IntracellularElectrode(varargin)
        % INTRACELLULARELECTRODE Constructor for IntracellularElectrode
        %     obj = INTRACELLULARELECTRODE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % device = Device
        % filtering = char
        % initial_access_resistance = char
        % location = char
        % resistance = char
        % seal = char
        % slice = char
        varargin = [{'help' 'Metadata about an intracellular electrode'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'device',[]);
        addParameter(p, 'filtering',[]);
        addParameter(p, 'initial_access_resistance',[]);
        addParameter(p, 'location',[]);
        addParameter(p, 'resistance',[]);
        addParameter(p, 'seal',[]);
        addParameter(p, 'slice',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        obj.device = p.Results.device;
        obj.filtering = p.Results.filtering;
        obj.initial_access_resistance = p.Results.initial_access_resistance;
        obj.location = p.Results.location;
        obj.resistance = p.Results.resistance;
        obj.seal = p.Results.seal;
        obj.slice = p.Results.slice;
        if strcmp(class(obj), 'types.core.IntracellularElectrode')
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
    function obj = set.filtering(obj, val)
        obj.filtering = obj.validate_filtering(val);
    end
    function obj = set.initial_access_resistance(obj, val)
        obj.initial_access_resistance = obj.validate_initial_access_resistance(val);
    end
    function obj = set.location(obj, val)
        obj.location = obj.validate_location(val);
    end
    function obj = set.resistance(obj, val)
        obj.resistance = obj.validate_resistance(val);
    end
    function obj = set.seal(obj, val)
        obj.seal = obj.validate_seal(val);
    end
    function obj = set.slice(obj, val)
        obj.slice = obj.validate_slice(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_device(obj, val)
        val = types.util.checkDtype('device', 'types.core.Device', val);
    end
    function val = validate_filtering(obj, val)
        val = types.util.checkDtype('filtering', 'char', val);
    end
    function val = validate_initial_access_resistance(obj, val)
        val = types.util.checkDtype('initial_access_resistance', 'char', val);
    end
    function val = validate_location(obj, val)
        val = types.util.checkDtype('location', 'char', val);
    end
    function val = validate_resistance(obj, val)
        val = types.util.checkDtype('resistance', 'char', val);
    end
    function val = validate_seal(obj, val)
        val = types.util.checkDtype('seal', 'char', val);
    end
    function val = validate_slice(obj, val)
        val = types.util.checkDtype('slice', 'char', val);
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
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.device)
            refs = obj.device.export(fid, [fullpath '/device'], refs);
        else
            error('Property `device` is required.');
        end
        if ~isempty(obj.filtering)
            if startsWith(class(obj.filtering), 'types.untyped.')
                refs = obj.filtering.export(fid, [fullpath '/filtering'], refs);
            elseif ~isempty(obj.filtering)
                io.writeDataset(fid, [fullpath '/filtering'], class(obj.filtering), obj.filtering, false);
            end
        end
        if ~isempty(obj.initial_access_resistance)
            if startsWith(class(obj.initial_access_resistance), 'types.untyped.')
                refs = obj.initial_access_resistance.export(fid, [fullpath '/initial_access_resistance'], refs);
            elseif ~isempty(obj.initial_access_resistance)
                io.writeDataset(fid, [fullpath '/initial_access_resistance'], class(obj.initial_access_resistance), obj.initial_access_resistance, false);
            end
        end
        if ~isempty(obj.location)
            if startsWith(class(obj.location), 'types.untyped.')
                refs = obj.location.export(fid, [fullpath '/location'], refs);
            elseif ~isempty(obj.location)
                io.writeDataset(fid, [fullpath '/location'], class(obj.location), obj.location, false);
            end
        end
        if ~isempty(obj.resistance)
            if startsWith(class(obj.resistance), 'types.untyped.')
                refs = obj.resistance.export(fid, [fullpath '/resistance'], refs);
            elseif ~isempty(obj.resistance)
                io.writeDataset(fid, [fullpath '/resistance'], class(obj.resistance), obj.resistance, false);
            end
        end
        if ~isempty(obj.seal)
            if startsWith(class(obj.seal), 'types.untyped.')
                refs = obj.seal.export(fid, [fullpath '/seal'], refs);
            elseif ~isempty(obj.seal)
                io.writeDataset(fid, [fullpath '/seal'], class(obj.seal), obj.seal, false);
            end
        end
        if ~isempty(obj.slice)
            if startsWith(class(obj.slice), 'types.untyped.')
                refs = obj.slice.export(fid, [fullpath '/slice'], refs);
            elseif ~isempty(obj.slice)
                io.writeDataset(fid, [fullpath '/slice'], class(obj.slice), obj.slice, false);
            end
        end
    end
end

end