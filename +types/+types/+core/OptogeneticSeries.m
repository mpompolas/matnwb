classdef OptogeneticSeries < types.core.TimeSeries
% OptogeneticSeries Optogenetic stimulus.  The data[] field is in unit of watts.


% PROPERTIES
properties
    site; % link to OptogeneticStimulusSite group that describes the site to which this stimulus was applied
end

methods
    function obj = OptogeneticSeries(varargin)
        % OPTOGENETICSERIES Constructor for OptogeneticSeries
        %     obj = OPTOGENETICSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % site = OptogeneticStimulusSite
        varargin = [{'data_unit' 'watt' 'help' 'Optogenetic stimulus'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'site',[]);
        parse(p, varargin{:});
        obj.site = p.Results.site;
        if strcmp(class(obj), 'types.core.OptogeneticSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.site(obj, val)
        obj.site = obj.validate_site(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_data_unit(obj, val)
        val = types.util.checkDtype('data_unit', 'char', val);
    end
    function val = validate_site(obj, val)
        val = types.util.checkDtype('site', 'types.core.OptogeneticStimulusSite', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.site)
            refs = obj.site.export(fid, [fullpath '/site'], refs);
        else
            error('Property `site` is required.');
        end
    end
end

end