function [handle,Zi,contourLineVals] = topoplot_bcl(bePlotedFileOrData,varargin)
% USAGE:
% [handle,Zi,contourLineVals] = topoplot_bcl(bePlotedFileOrData,'option',optionValue);
%
% e.g., [handle,Zi,contourLineVals] = topoplot_bcl('topoplotDemo.avg','negativePeriod',200,'plotAreas',[80 120],'viewpoint','back','colormap',2,'emarker2',{[60:63],'o',[1 1 1],10,1});
% argins:
%         bePlotedFileOrData:     [string]   filename of the to be plotted file including full file path.
%                              or [vector]   a vector of data to be plotted
% options:
%         'plotAreas':        [double vect] 1*2 vector of double, 1 and 2 col for start and end time (ms) respectively.
%         'EEGChans':         [cell string] (!!! only usefull if bePlotedFileOrData is a vector)the cell of strings comtains the chan names coresponding to each data point.
%         'negativePeriod'    [double]      the length of the inverval before the zero point in unit of ms (e.g., 200);
%         'maplimits'       - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0);
%                             'maxmin'   -> scale colors to the data range (makes green mid-range);
%                           - Or a double row vector:
%                             [lo,hi]    -> use user-definined lo/hi limits
%         'style'           - 'map'      -> plot colored map only
%                             'contour'  -> plot contour lines only
%                             'both'     -> plot both colored map and contour lines
%                             'fill'     -> plot constant color between contour lines
%                             'blank'    -> plot electrode locations only {default: 'fill'}
%         'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers'. To set the 'pts'
%                             marker,,see 'Plot detail options' below. {default: 'on' -> mark electrode
%                             locations with points ('.') unless more than 64 channels, then 'off'}.
%         'plotchans'       - [vector] channel numbers (indices) to use in making the head plot.
%                             {default: [] -> plot all chans}
%         'headrad'         - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head.
%                                  NOTE: Only headrad = 0.5 is anatomically correct! 0 -> don't draw head;
%                                 'rim' -> show cartoon head at outer edge of the plot {default: 0.71}
%         'numcontour'      - number of contour lines {default: 6}
%                             Or a vector of values for contour lines, each
%                             value corresponded to a contourline
%         'contourvals'     - values for contour {default: same as input values}
%        'contourLineWidth' - values for the width of the contourLines
%         'colormap'        -  (n,3) any size colormap {default: blue to red}
%                        or - a scaler : 1,2, and 3 for blue to red, hillyard lab's colormap, and Mishra (2008) respectively
%         'chanlocs'        - filename of the chanlocs, {default: 'elocsTopView.mat' }
%                             which is the average of elect locs measured
%                             from 88 chinese subject's brain.
%         'viewpoint'       - viewpoint of the topoplot map  {default: 'top'} ,'back','front','left' or 'right'
%         'emarker'         - Matlab marker char | {markerchar color size linewidth} char, else cell array
%                           specifying the electrode 'pts' marker. Ex: {'s','r',32,1} -> 32-point solid
%                           red square. {default: {'.','k',[],1} where marker size ([]) depends on the number
%                           of channels plotted}.
%        'emarker2'        - {markchans}|{markchans marker color size linewidth} cell array specifying
%                           an alternate marker for specified 'plotchans'. Ex: {[3 17],'s','g'}
%                           {default: none, or if {markchans} only are specified, then {markchans,'o','r',10,1}}
% argouts:
%
% handle          - figure handle of the topoplot
% Zi              - interpolated values
% contourLineVals - values of the contour lines
%
%     Written by Yang Zhang
%     Soochow University,china
%     zhangyang873@gmail.com
%     20150521
% 
% Rev by Yang. Fri May 22 22:10:00 2015
% Now will sort the chans in the locfile
% acorrding to the chan order contained in the avg file
%
% Rev by Yang at 2015/8/13 17:35:21
% added a pramemter contourLineWidth to control the width of the contour lines
% Rev. by Yang Zhang Thu Apr 19 22:39:52 2018
% 1) Added a prameter "emarker" to control the size of the electrodes
% 2) Added a prameter "emarker2" to control the highlighted elevtrodes (e.g., the slected electrodes for analsis)
% Rev. by Yang at 20180515
% Rev. by Yang Zhang Sun Nov 11 19:35:15 2018  Soochow University, China
% start to use a new methods to handle the input argins
%
% If you do think this function is usefull and have used it in your study, please cite our paper:
% Li A-S, Miao C-G, Han Y, He X and Zhang Y (2018) Electrophysiological Correlates of the Effect of Task Difficulty on Inhibition of Return. Front. Psychol. 9:2403.


helpStr = { 'USAGE:'
    ' [handle,Zi,contourLineVals] = topoplot_bcl(bePlotedFileOrData,''option'',optionValue);'
    ' '
    ' e.g., [handle,Zi,contourLineVals] = topoplot_bcl(''topoplotDemo.avg'',''negativePeriod'',200,''plotAreas'',[80 120],''viewpoint'',''back'',''colormap'',2,''emarker2'',{[55:60],''o'',[1 1 1],10,1});'
    ' argins:'
    '         bePlotedFileOrData:     [string]   filename of the be plotted file including full file path.'
    '                              or [vector]   a vector of data be plotted'
    ' options:'
    '         ''EEGChans'':         [cell string] (!!!only usefull if bePlotedFileOrData is a vector) the cell of strings comtains the chan names coresponding to each data point. '
    '         ''plotAreas'':        [double vect] 1*2 vector of double, 1 and 2 col for start and end time (ms) respectively.'
    '         ''negativePeriod''    [double]      the length of the inverval before the zero point in unit of ms (e.g., 200);'
    '         ''maplimits''       - ''absmax''   -> scale map colors to +/- the absolute-max (makes green 0); '
    '                               ''maxmin''   -> scale colors to the data range (makes green mid-range);'
    '                          or -  a double row vector:'
    '                                 [lo,hi]    -> use user-definined lo/hi limits'
    '         ''style''           - ''map''      -> plot colored map only'
    '                             ''contour''  -> plot contour lines only'
    '                             ''both''     -> plot both colored map and contour lines'
    '                             ''fill''     -> plot constant color between contour lines'
    '                             ''blank''    -> plot electrode locations only {default: ''fill''}'
    '         ''electrodes''      - ''on'',''off'',''labels'',''numbers'',''ptslabels'',''ptsnumbers''. To set the ''pts'' '
    '                             marker,,see ''Plot detail options'' below. {default: ''on'' -> mark electrode '
    '                             locations with points (''.'') unless more than 64 channels, then ''off''}.'
    '         ''plotchans''       - [vector] channel numbers (indices) to use in making the head plot. '
    '                             {default: [] -> plot all chans}'
    '         ''headrad''         - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head. '
    '                                  NOTE: Only headrad = 0.5 is anatomically correct! 0 -> don''t draw head; '
    '                                 ''rim'' -> show cartoon head at outer edge of the plot {default: 0.71}'
    '         ''numcontour''      - number of contour lines {default: 6} '
    '                             Or a vector of values for contour lines, each'
    '                             value corresponded to a contourline '
    '         ''contourvals''     - values for contour {default: same as input values} USE IT CAREFULLY!!'
    '        ''contourLineWidth'' - values for the width of the contourLines'
    '         ''colormap''        -  (n,3) any size colormap {default: existing colormap}'
    '         ''chanlocs''        - filename of the chanlocs, {default: ''elocsTopView.mat'' }'
    '                               which is the average of elect locs measured'
    '                               from 88 chinese subject''s brain. '
    '         ''viewpoint''        - viewpoint of the topoplot map  {default: ''top''} ,''back'',''front'',''left'' or ''right'''
    '         ''emarker''         - Matlab marker char | {markerchar color size linewidth} char, else cell array'
    '                            specifying the electrode ''pts'' marker. Ex: {''s'',''r'',32,1} -> 32-point solid'
    '                            red square. {default: {''.'',''k'',[],1} where marker size ([]) depends on the number'
    '                            of channels plotted}.'
    '        ''emarker2''        - {markchans}|{markchans marker color size linewidth} cell array specifying'
    '                            an alternate marker for specified ''plotchans''. Ex: {[3 17],''s'',''g''}'
    '                            {default: none, or if {markchans} only are specified, then {markchans,''o'',''r'',10,1}}'
    ' '
    ' argouts:'
    ' '
    ' handle          - figure handle of the topoplot'
    ' Zi              - interpolated values'
    ' contourLineVals - values of the contour lines'
    ' '
    '     Written by Yang Zhang'
    '     Soochow University,china'
    '     zhangyang873@gmail.com'
    '     20150521'
    ' '
    ' Rev by Yang. Fri May 22 22:10:00 2015'
    ' Now will sort the chans in the locfile'
    ' acorrding to the chan order contained in the avg file'
    ' '
    ' Rev by Yang at 2015/8/13 17:35:21'
    ' added a pramemter contourLineWidth to control the width of the contour lines'
    ' '
    ' Rev by Yang at 2016-09-25 1:48:31'
    'Added support the back view point topoplot map '
    'Rev. by Yang Zhang Sun Nov 11 19:35:15 2018  Soochow University, China'
    'start to use a new methods to handle the input argins'
    'If you do think this function is usefull and have used it in your study, please cite our paper:'
    'Li A-S, Miao C-G, Han Y, He X and Zhang Y (2018) Electrophysiological Correlates of the Effect of Task Difficulty on Inhibition of Return. Front. Psychol. 9:2403.'};

if nargin < 1
    for iRow = 1:numel(helpStr)
        disp(helpStr{iRow});
    end
    
    handle          = [];
    Zi              = [];
    contourLineVals = [];
    
    return;
end


mFileFolder = fileparts(mfilename('fullpath')); %#ok<*NASGU>

g = cell2struct(varargin(2:2:end),varargin(1:2:end),2);

isAvgFile = ischar(bePlotedFileOrData);

try g.viewpoint;        catch, g.viewpoint        = 'top';                end
try g.emarker2;         catch, g.emarker2         = {[],'o',[1 1 1],6,1}; end
try g.emarker;          catch, g.emarker          = {'.','k',20,1};       end
try g.chanlocs;         catch, g.chanlocs         = 'elocsTopView.mat';   end
try g.colormap;         catch, g.colormap         = 1;                    end % 1 for red gradually change to blue color map
try g.contourLineWidth; catch, g.contourLineWidth = 1;                    end
try g.contourvals;      catch, g.contourvals      = [];                   end
try g.numcontour;       catch, g.numcontour       = 6;                    end
try g.headrad;          catch, g.headrad          = 0.71;                 end
try g.plotchans;        catch, g.plotchans        = [];                   end
try g.electrodes;       catch, g.electrodes       = 'on';                 end
try g.style;            catch, g.style            = 'fill';               end
try g.maplimits;        catch, g.maplimits        = 'absmax';             end
try g.figureTitle;      catch, g.figureTitle      = '';                   end
try g.negativePeriod;   catch, g.negativePeriod   = 0;                    end
try g.EEGChans;         catch, g.EEGChans         = {};                   end
try g.plotAreas;        catch, g.plotAreas        = [];                   end




if numel(g.colormap) == 1
    switch g.colormap
        case 1
            %--- make default colormap ----/
            redMap         = repmat([1 0 0],8,1);
            blueMap        = repmat([0 0 1],8,1);
            alphaIdx       = repmat(linspace(0,1,8)',1,3);

            alphaedRed     = 1 - alphaIdx + redMap.*alphaIdx;
            alphaedBlue    = alphaIdx+ blueMap.*(1 - alphaIdx);

            g.colormap = [alphaedBlue;alphaedRed];

            %--free memory --/
            clear redMap blueMap alphaIdx alphaedRed alphaedBlue;
            %----------------\
            %------------------------------\
        case 2
            %  hillyard lab's colormap
            g.colormap = [...
            0.7216    0.2118    0.6314
            0.5490    0.2392    0.6392
            0.4275    0.2588    0.6471
            0.2471    0.3216    0.6745
            0.3725    0.7529    0.9373
            0.4745    0.8235    0.8588
            0.4431    0.8000    0.6235
            0.4196    0.7882    0.4157
            0.3725    0.7569    0.1961
            0.6980    0.8549    0.0941
            0.9294    0.9294    0.0588
            0.9725    0.8000    0.0275
            0.9569    0.4706    0.0824
            0.9255    0.1176    0.1412];
        case 3
            % Mishra, J., Martinez, A., & Hillyard, S. A. (2008). Cortical processes underlying sound-induced flash fusion. Brain Research, 1242, 102??15. 
            g.colormap = [...
             0 0 0
             255 0 0
             255 49 0
             255 116 0
             255 156 0
             255 192 0
             247 255 0
             158 255 0
             0 255 0
             154 248 255
             125 231 255
             72 217 255
             103 170 255
             78 126 229
             50 67 223
             70 0 255
             255 255 255]./255;
        otherwise
            error('for scaler input of ''colormap'', it should be of [1 2 3]!');
    end
end 

% g = finputcheck(varargin, {'plotAreas'          'real'    []                                [];
%     'EEGChans'           'cell'    []                                {};
%     'negativePeriod'     'real'    []                                 0;
%     'figureTitle'        'string'  ''                                '';
%     'maplimits'          'string'  {'absmax','maxmin'}          'absmax';
%     'maplimitsDouble'    'real'    []              [];
%     'style'              'string'  {'map','contour','both','fill','blank'} 'fill';
%     'electrodes'         'string'  { 'on','off','labels','numbers'}         'on';
%     'plotchans'          'real'    []              [];
%     'headrad'            'real'    []             0.71;
%     'numcontour'         'real'    []               6;
%     'contourvals'        'real'    []               [];
%     'contourLineWidth'   'real'    []               1;
%     'colormap'           'real'    []  colorMapMatrix;
%     'chanlocs'           'string'  []   'elocsTopView.mat';
%     'emarker'             'cell'   []   {'.','k',20,1};
%     'emarker2'            'cell'   []   {[],'o',[1 1 1],6,1};
%     'viewpoint'          'string'  {'front','back','top','left','right'}   'top'});
% 
% if ischar(g)
%     error(g);
% end

if isempty(g.headrad)
    g.headrad = 'rim';
end



if isAvgFile
    [filepath,Noused,EXT] = fileparts(bePlotedFileOrData); %#ok<*ASGLU>
    
    if strcmpi(EXT,'dat')
        error('We are not support the dat format of the data anymore :-( Please try to use the avg format data!');
    else
        
        if isempty(g.plotAreas)
            promptstr         = {'Define the topoplot areas (start,end)',...
                'Baseline period (before the zero point in ms)'};
            inistr            = {'300,500','200'};
            pop_title         = sprintf('Plot parameters');
            result            = inputdlg2(promptstr,pop_title,1,inistr,'');
            if size(result,1) == 0
                return;
            end;
            
            g.plotAreas      = eval(['[' result{1} ']']);
            g.negativePeriod = eval(['[' result{2} ']']);
        end
        
        if isempty(g.figureTitle)
            g.figureTitle = [num2str(g.plotAreas(1)),' - ',num2str(g.plotAreas(2)),' ms'];
        end
        
        
        [data,chan_names,Noused, pnts, fs]=loadavg_bcl(bePlotedFileOrData);
        
        EEGChans  = cell(size(chan_names,1),1);
        %----- transform the chan info into a cell of string ----/
        for iChan = 1:size(chan_names,1)
            EEGChans{iChan,1} = deblank(char(chan_names(iChan,:)));
        end
        %--------------------------------------------------------\
        
    end
    
else % for data vector inputs
    data     = bePlotedFileOrData;
    EEGChans = g.EEGChans;
    
    if numel(EEGChans) ~= numel(bePlotedFileOrData)
        errorStr = ['For vector input, the "EEGChans" paramter should be defined',char(10)];
        errorStr = [errorStr,'Numel of channames(EEGChans) should be equal to numel(data) '];

        error(errorStr);
    end
end % isAvgFile



EEGChans = lower(EEGChans);
isEog    = ismember(EEGChans,{'veo','heo','veog','heog','gfp'});

EEGChans(isEog) = [];

if isAvgFile
    % get the mean amplitudes
    timeinterval   = 1000/fs;
    
    plotAreas_down = (g.plotAreas(1)+g.negativePeriod)/timeinterval +1;
    plotAreas_up   = (g.plotAreas(2)+g.negativePeriod)/timeinterval +1;
    
    datavector     = mean(data(plotAreas_down:plotAreas_up,:));
else
    datavector     = data;
end

datavector(isEog) = [];

%----- get the default contour values ----/
if isempty(g.contourvals)
    g.contourvals = datavector;
end
%-----------------------------------------\

%-- load the chan locations info --/
chanInfo = load(g.chanlocs);
my_chan  = chanInfo.chanlocs;
my_chan  = changeLocsForViews(my_chan, g.viewpoint);
%----------------------------------\

%--- find out chans that has no location info ---/
locChanNames = lower({my_chan(:).labels}');
[isInChanDataBase,locsInChanDataBase] = ismember(EEGChans,locChanNames);
%-----------------------------------------------\

%--- remove chans that has no location info  ---/
my_chan    = my_chan(locsInChanDataBase);
datavector = datavector(isInChanDataBase);
%----------------------------------------------\

if isnumeric(g.maplimits)
    g.numcontour  =  min(g.maplimits):(max(g.maplimits) - min(g.maplimits))/(g.numcontour - 1):max(g.maplimits);
end






fprintf('---------- topoplot_bcl version: 2015/8/14 15:21:36 ------------\n');

if strcmpi(g.style,'blank')
    datavector = zeros(size(datavector));
    g.style = 'fill';
end

if strcmpi(g.viewpoint,'top')
    g.headrad = 0.71;
end


[handle,Zi,grid,Xi,Yi,contourLineVals,hp2b] = topoplotNew(... % added to topoplotNew by yang 2016-09-25 1:42:51
    datavector,         my_chan,...
    'maplimits',        g.maplimits, ...
    'style',            g.style,...
    'electrodes',       g.electrodes, ...
    'plotchans',        g.plotchans, ...
    'headrad',          g.headrad, ...
    'numcontour',       g.numcontour,...
    'contourvals',      g.contourvals,...
    'contourLineWidth', g.contourLineWidth,...
    'colormap',         g.colormap,...
    'emarker',          g.emarker,... % added by yang Thu Apr 19 21:50:28 2018
    'emarker2',         g.emarker2,...
    'viewpoint',        g.viewpoint); % added by yang 2016-09-25 1:42:25

if ~isempty(hp2b)
    set(hp2b,'Color',[0 0 0]);
end 

allValues                   = Zi;
allValues                   = allValues(:);
allValues(isnan(allValues)) = [];

title(g.figureTitle);

fprintf('interpolated Data: min: %-10.4f max: %-10.4f\n',min(allValues),max(allValues));


% if isempty(g.maplimitsDouble)
%     if strcmpi(g.maplimits,'absmax')

%         g.maplimits = [-max(abs(allValues)) max(abs(allValues))];

%     elseif strcmpi(g.maplimits,'maxmin')

%         g.maplimits = [min(allValues) max(allValues)];

%     else
%         error('maplimits should be of {''absmax'',''maxmin''}');
%     end
% end


% colorbarhandle = cbar(0,0,g.maplimits);

colorbarhandle = cbar;

colormap(g.colormap);

colorAxisRange = get(colorbarhandle,'YLim');

fprintf('Colormap range cf: min: %-10.4f max: %-10.4f\n',colorAxisRange([2 1]));

set(colorbarhandle,'GridLineStyle','-','YGrid','on','position',[0.92 0.11 0.013 0.48]);

fprintf('-------------- Yang: topoplot_bcl finished ---------------------\n');

fprintf('\n\n-------------------------------------------------------------------------\n');
fprintf('If you do think this function is usefull and have used it in your study, please cite our paper:\n');
fprintf('Li A-S, Miao C-G, Han Y, He X and Zhang Y (2018) Electrophysiological Correlates of the Effect of Task Difficulty on Inhibition of Return. Front. Psychol. 9:2403.\n');
fprintf('-------------------------------------------------------------------------\n');
end  % end of the main function






%%%%%%%%%%%%%%%%%%%%%%%
% sub function
%%%%%%%%%%%%%%%%%%%%%%%


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
            NewChanlocs(iChan).labels = chanlocs(iChan).labels; %#ok<*AGROW>
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
            NewChanlocs(iChan).labels = chanlocs(iChan).labels;
            NewChanlocs(iChan).X      = chanlocs(iChan).Z;
            NewChanlocs(iChan).Y      = -chanlocs(iChan).X;
            NewChanlocs(iChan).Z      = -chanlocs(iChan).Y;
        end
        
    otherwise
        error('"viewPoint" should be of {''back'',''top'',''front'',''left'', or ''right''}');
        
end

NewChanlocs = convertlocs(NewChanlocs,'cart2all');

end
