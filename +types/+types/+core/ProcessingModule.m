classdef ProcessingModule < types.core.NWBContainer
% ProcessingModule Module.  Name should be descriptive. Stores a collection of related data organized by contained interfaces.  Each interface is a contract specifying content related to a particular type of data.


% PROPERTIES
properties
    description; % Description of Module
    nwbdatainterface; % Interface objects containing data output from processing steps
end

methods
    function obj = ProcessingModule(varargin)
        % PROCESSINGMODULE Constructor for ProcessingModule
        %     obj = PROCESSINGMODULE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % nwbdatainterface = NWBDataInterface
        varargin = [{'help' 'A collection of analysis outputs from processing of data'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        [obj.nwbdatainterface, ivarargin] = types.util.parseConstrained(obj,'nwbdatainterface', 'types.core.NWBDataInterface', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        if strcmp(class(obj), 'types.core.ProcessingModule')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.nwbdatainterface(obj, val)
        obj.nwbdatainterface = obj.validate_nwbdatainterface(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_nwbdatainterface(obj, val)
        constrained = {'types.core.NWBDataInterface'};
        types.util.checkSet('nwbdatainterface', struct(), constrained, val);
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
        if ~isempty(obj.nwbdatainterface)
            refs = obj.nwbdatainterface.export(fid, fullpath, refs);
        end
    end
end

end