classdef CurrentClampStimulusSeries < types.core.PatchClampSeries
% CurrentClampStimulusSeries Aliases to standard PatchClampSeries. Its functionality is to better tag PatchClampSeries for machine (and human) readability of the file.



methods
    function obj = CurrentClampStimulusSeries(varargin)
        % CURRENTCLAMPSTIMULUSSERIES Constructor for CurrentClampStimulusSeries
        %     obj = CURRENTCLAMPSTIMULUSSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        varargin = [{'help' 'Stimulus current applied during current clamp recording'} varargin];
        obj = obj@types.core.PatchClampSeries(varargin{:});
        if strcmp(class(obj), 'types.core.CurrentClampStimulusSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    
    %% VALIDATORS
    
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.PatchClampSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
    end
end

end