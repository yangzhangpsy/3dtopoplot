function NewChanlocs = changeLocsForViews(chanlocs,viewPoint)

% argin:
%
% chanlocs:   [struct] A struct of the channel locs in EEGLAB format
% viewPoint:  [string] view point of the new locs, should be of {'back','top','front','left' or 'right'}
%
% argout:
%
% newChanlocs:  [struct] A struct of the channel locs in EEGLAB format
% Written by Yang Zhang at Soochow University
% 2018/5/14 14:03:05
if nargin < 1
    error('changeLocsForViews requires two input parameters, see help changeLocsForViews');
end
% in matlab XYZ Cartisan coordinates
% X is toward nose;
% Y is toward left ear;
% Z is toward vertex(the sky);

%

switch lower(viewPoint)
    case 'back'
        % new x:  <---  Z
        % new y:  <---  Y
        % new z:  <--- -X
        for iChan = 1:numel(chanlocs)
            NewChanlocs(iChan).labels = chanlocs(iChan).labels;
            NewChanlocs(iChan).X      = chanlocs(iChan).Z;
            NewChanlocs(iChan).Y      = chanlocs(iChan).Y;
            NewChanlocs(iChan).Z      = -chanlocs(iChan).X;
        end
        
    case 'top'
        % do nothing...
        NewChanlocs = chanlocs;
    case 'front'
        % new x:  <---  Z
        % new y:  <--- -Y
        % new z:  <---  X
        for iChan = 1:numel(chanlocs)
            NewChanlocs(iChan).labels = chanlocs(iChan).labels;
            NewChanlocs(iChan).X      = chanlocs(iChan).Z;
            NewChanlocs(iChan).Y      = -chanlocs(iChan).Y;
            NewChanlocs(iChan).Z      = chanlocs(iChan).X;
        end
        
    case 'left'
        % new x:  <---  Z
        % new y:  <---  X
        % new z:  <---  Y
        for iChan = 1:numel(chanlocs)
            NewChanlocs(iChan).labels = chanlocs(iChan).labels;
            NewChanlocs(iChan).X      = chanlocs(iChan).Z;
            NewChanlocs(iChan).Y      = chanlocs(iChan).X;
            NewChanlocs(iChan).Z      = chanlocs(iChan).Y;
        end
        
        
    case 'right'
        % new x:  <---  Z
        % new y:  <--- -X
        % new z:  <--- -Y
        for iChan = 1:numel(chanlocs)
            NewChanlocs(iChan).labels = chanlocs(iChan).labels; %#ok<*AGROW>
            NewChanlocs(iChan).X      = chanlocs(iChan).Z;
            NewChanlocs(iChan).Y      = -chanlocs(iChan).X;
            NewChanlocs(iChan).Z      = -chanlocs(iChan).Y;
        end
        
    otherwise
        error('"viewPoint" should be of {''back'',''top'',''front'',''left'', or ''right''}');
        
end



NewChanlocs = convertlocs(NewChanlocs,'cart2all');

end