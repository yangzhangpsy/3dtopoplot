function [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps]=loadavg_bcl(FILENAME,chanNameList,noEEGChans,addGFP)
% written by Yang zhang based on EEGLAB's loadavg.m
% 2009-11-25
% Northeast Normal University
% Changchun, China
%
%
% Useage:
% [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps]=loadavg_bcl(FILENAME,chanNameList,noEEGChans,addGFP);
%
% input args:
%
% FILENAME               [string]             Filename with fullpath, e.g., 'C:\test.avg';
% chanNameList           [cell of string]     channel names, e.g, {'FP1','FP2'}; Default is {'all'}
% noEEGChans             [cell of string]     noEEG channel names that will be excluded from caculating GFP, e.g., {'heog','veog'};
% addGFP                 [double]             1 and [0] for product and not product GFP respectively when there is no gfp channel in the avg file
%
% update by Yang Zhang
% 2011-2-28
% Northeast Normal University
% Changchun, China

% revised by Yang Zhang
% 2011-4-6
% use an automatic routine to identifiy the type of the avg file to get the true nsweeps parameter
% rev. by YZ 2012-04-11 fix global field power problem
% bug fixation:
% now the output signal and chan_names will be sorted to correpond with the imput chanNameList.
% rev. by Yang Zhang
% 6-June-2012
%

% Rev. by Yang zhang
% Northeast Normal University
% Changchun, China
% 2012-11-19
% 1) Now, added a extra-input paramter to determine whether input the GFP or not if the avg file contains no gfp channel
% 2) completed the help infomation.
%
% Rev. by Yang Zhang Thu Jun  4 16:38:28 2015
% Soochow University, China
% added a output parameter "nsweeps" for NofSweep used to computer the avg file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             ARGs checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('chanNameList','var')
    chanNameList={'all'};
end

if ~exist('noEEGChans','var')||isempty(noEEGChans)
    noEEGChans={'gfp','veog','heog','heo','veo','ref'};
end

if ~exist('addGFP','var')||isempty(addGFP)
    addGFP=0;
end


if ischar(chanNameList)
    chanNameList = {chanNameList};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           BEGIN of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('==========================');
fprintf('\nDefined noEEG channels:');
%         fprintf('no eeg channel names = %s\n',noEEGChans{:});
fprintf(' %s',noEEGChans{:});
fprintf('\n');
%     disp('==========================');
try
    fid=fopen(FILENAME,'r','ieee-le');
    % read general part of the erp header and set variables
    % -----------------------------------------------------
    fseek(fid, 362, 'cof');% skip the firsts 362 from BOF (368 bytes in orginal comment?)
    % disp(ftell(fid));% 360 bytes
    % hdr.nsweeps = fread(fid, 1, 'ushort');
    % disp(ftell(fid));% 362
    % hdr.compsweeps = fread(fid, 1, 'ushort');% the exact sweep numbers for eeg and single subject avg file| in grand avg file it represents the number of subjects
    % hdr.acceptcnt = fread(fid, 1, 'ushort');% number of accepted sweeps also the exact sweep numbers for grand avg file
    % hdr.rejectcnt = fread(fid, 1, 'ushort');%number of rejected sweeps
    % disp(ftell(fid));% 368
    
    compsweeps = fread(fid, 1, 'ushort');% the exact sweeps numbers for eeg file| in Grand avg file it represented number of subjects
    
    acceptcnt  = fread(fid, 1, 'ushort');% number of accepted sweeps
    rejectcnt  = fread(fid, 1, 'ushort');%number of rejected sweeps
    
    
    % determine the type of avg file and choose the right value for nsweeps
    if (rejectcnt+acceptcnt)~=compsweeps
        fprintf('Grand avg : nSubs = %d; sweeps = %d\n',compsweeps,acceptcnt);
        %             disp(['Grand       avg: nSubs = ',num2str(compsweeps),'; sweeps = ',num2str(acceptcnt)]);
        nsweeps     = compsweeps;
        
        isSingleAvg = false;
    else
        % disp('It''s a single subject average file!!!');
        disp(['Single sub avg : nsweeps = ',num2str(compsweeps),'; Accepted sweeps = ',num2str(acceptcnt),'; Rejected sweeps = ',num2str(rejectcnt)]);
        nsweeps     = acceptcnt;
        
        isSingleAvg = true;
    end
    
    %         fprintf('nSweeps or nSubjects: %d\n',nsweeps);
    
    pnts          =fread(fid, 1, 'ushort');	% number of point per waveform
    chan          =fread(fid, 1, 'ushort');  % number of channels
    fseek(fid, 3, 'cof');	% skip 3 bytes
    variance_flag =fread(fid, 1, 'uchar');
    rate          =fread(fid, 1, 'ushort');  % sample rate (Hz)
    fseek(fid, 127, 'cof');	% skip 127 bytes
    xmin          =fread(fid, 1, 'float32'); % in s
    xmax          =fread(fid, 1, 'float32'); % in s
    fseek(fid, 387, 'cof');	% skip 387 bytes
    
    % read electrode configuration
    % ----------------------------
    for elec = 1:chan
        channel_label_tmp    = fread(fid, 10, 'uchar');
        electrodes(elec).tmp = channel_label_tmp;
        chan_names(elec,:)   = (channel_label_tmp)'; %#ok<*AGROW>
        for index = 2:9
            if chan_names(elec,index) == 0
                chan_names(elec,index)=' ';
            end;
        end;
        fseek(fid, 61, 'cof');%skip 61 bytes
        electrodes(elec).calib= fread(fid, 1, 'float32');
        
        %     	erp = fread(fid, 47-10, 'uchar');
        % 	baseline(elec) = fread(fid, 1, 'ushort');
        % 	erp = fread(fid, 10, 'uchar');
        % 	sensitivity(elec) = fread(fid, 1, 'float32');
        % 	erp = fread(fid, 8, 'uchar');
        % 	    electrodes(elec).calib= fread(fid, 1, 'float32');
        %
        % 	fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', ...
        % char(chan_names(elec,1:4)), baseline(elec), sensitivity(elec), electrodes(elec).calib);
    end;
    
    signal    = zeros(pnts, chan);
    variance  = zeros(pnts, chan);
    
    isEEGChan = true(size(chan_names,1),1);
    
    for elec = 1:chan
        
        if  ismember(lower(deblank(char(chan_names(elec,:)))),lower(noEEGChans))
            isEEGChan(elec)=false;
        end
        
        if  strcmpi(deblank(char(chan_names(elec,:))),'gfp')
            addGFP = elec;
        end
        % To scale a data point to
        % microvolts, multiply by the channel-specific calibration factor (i.e., for electrode j:
        % channel[j]->calib) and divide by the number of sweeps in the average (i.e.,
        % channel[j]->n);
        % skip sweeps header and read data
        % --------------------------------
        fseek(fid, 5, 'cof');
        
        if  strcmpi('gfp',deblank(char(chan_names(elec,:))))&&isSingleAvg
            signal(:, elec) =fread(fid, pnts, 'float32')*electrodes(elec).calib;
            disp(' gfp ');
        else
            signal(:, elec) =fread(fid, pnts, 'float32')*electrodes(elec).calib/nsweeps;
        end
        
    end;% FOR IELOC
    
    fprintf('NoEEG channels: ');
    fprintf('%d ',find(~isEEGChan)');
    fprintf('\n');
    
    if variance_flag
        for elec = 1:chan
            if strcmpi('gfp',deblank(char(chan_names(elec,:))))&&isSingleAvg
                variance(:, elec) = fread(fid, pnts, 'float32')*electrodes(elec).calib;
            else
                variance(:, elec) = fread(fid, pnts, 'float32')*electrodes(elec).calib/nsweeps;% not sure
            end
        end;
    else
        variance = 'novariance';
    end;
    
    
    %---- added the gfp chanel if not exist----/
    if ismember('gfp',lower(chanNameList))
        eegData=signal(:,isEEGChan);
        
        globalFieldPower=std(eegData,1,2);
        
        if addGFP
            disp('The input data alread have a GFP channel, we will use this data directly"');
            %                 signal(:,addGFP)=globalFieldPower;
        else
            disp('The input data have no GFP channel, so we will caculate the GFP as "gfp=std(eeg,1,2)"');
            signal=[signal,globalFieldPower];
            
            chan_names=[chan_names;chan_names(end,:)];
            chan_names(end,1:3)=double('GFP');
            
            if variance_flag
                variance=[variance,zeros(size(variance,1),1)];
            end
        end
        
    end
    %---------------------------------------\
    %   char(chan_names)
    %%
    if ~strcmpi(chanNameList{1},'all')
        inOrignalIdxs=[];
        
        for ichanList=1:numel(chanNameList)
            
            for elec=1:size(chan_names,1)
                if strcmpi(chanNameList{ichanList},char(chan_names(elec,1:numel(chanNameList{ichanList}))))
                    inOrignalIdxs=[inOrignalIdxs,elec];
                    break;
                end
            end
            
        end
        
        signal=signal(:,inOrignalIdxs);
        if variance_flag
            variance=variance(:,inOrignalIdxs);
        end
        chan_names=chan_names(inOrignalIdxs,:);
        
    end
    % signal = signal';
    % variance = variance';
    fclose(fid);
    
catch errorLOAD
    disp(FILENAME);
    rethrow(errorLOAD);
end
return;


