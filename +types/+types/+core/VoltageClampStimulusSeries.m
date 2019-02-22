classdef VoltageClampStimulusSeries < types.core.PatchClampSeries
% VoltageClampStimulusSeries Aliases to standard PatchClampSeries. Its functionality is to better tag PatchClampSeries for machine (and human) readability of the file.



methods
    function obj = VoltageClampStimulusSeries(varargin)
        % VOLTAGECLAMPSTIMULUSSERIES Constructor for VoltageClampStimulusSeries
        %     obj = VOLTAGECLAMPSTIMULUSSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        varargin = [{'help' 'Stimulus voltage applied during voltage clamp recording'} varargin];
        obj = obj@types.core.PatchClampSeries(varargin{:});
        if strcmp(class(obj), 'types.core.VoltageClampStimulusSeries')
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