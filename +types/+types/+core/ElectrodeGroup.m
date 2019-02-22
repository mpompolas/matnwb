classdef ElectrodeGroup < types.core.NWBContainer
% ElectrodeGroup One of possibly many groups, one for each electrode group.


% PROPERTIES
properties
    description; % description of this electrode group
    device; % the device that was used to record from this electrode group
    location; % description of location of this electrode group
end

methods
    function obj = ElectrodeGroup(varargin)
        % ELECTRODEGROUP Constructor for ElectrodeGroup
        %     obj = ELECTRODEGROUP(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % location = char
        % device = Device
        varargin = [{'help' 'A physical grouping of channels'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'location',[]);
        addParameter(p, 'device',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        obj.location = p.Results.location;
        obj.device = p.Results.device;
        if strcmp(class(obj), 'types.core.ElectrodeGroup')
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
    function obj = set.location(obj, val)
        obj.location = obj.validate_location(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_device(obj, val)
        val = types.util.checkDtype('device', 'types.core.Device', val);
    end
    function val = validate_location(obj, val)
        val = types.util.checkDtype('location', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBContainer(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], class(obj.description), obj.description, false);
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.device)
            refs = obj.device.export(fid, [fullpath '/device'], refs);
        end
        if ~isempty(obj.location)
            io.writeAttribute(fid, [fullpath '/location'], class(obj.location), obj.location, false);
        else
            error('Property `location` is required.');
        end
    end
end

end