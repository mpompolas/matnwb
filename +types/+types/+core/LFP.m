classdef LFP < types.core.NWBDataInterface
% LFP LFP data from one or more channels. The electrode map in each published ElectricalSeries will identify which channels are providing LFP data. Filter properties should be noted in the ElectricalSeries description or comments field.


% PROPERTIES
properties
    electricalseries; % ElectricalSeries object containing LFP data for one or more channels
end

methods
    function obj = LFP(varargin)
        % LFP Constructor for LFP
        %     obj = LFP(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % electricalseries = ElectricalSeries
        varargin = [{'help' 'LFP data from one or more channels. Filter properties should be noted in the ElectricalSeries'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        [obj.electricalseries, ivarargin] = types.util.parseConstrained(obj,'electricalseries', 'types.core.ElectricalSeries', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        parse(p, varargin{:});
        if strcmp(class(obj), 'types.core.LFP')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.electricalseries(obj, val)
        obj.electricalseries = obj.validate_electricalseries(val);
    end
    %% VALIDATORS
    
    function val = validate_electricalseries(obj, val)
        constrained = {'types.core.ElectricalSeries'};
        types.util.checkSet('electricalseries', struct(), constrained, val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.electricalseries)
            refs = obj.electricalseries.export(fid, fullpath, refs);
        else
            error('Property `electricalseries` is required.');
        end
    end
end

end