% Option 1: No GLM, just use SS coefficients = 1 to get individual trials
    % In this option, no CV is performed here. 
% Option 2: Use GLM only for SS coefficients since it's the same for all
% conditions. However, require building design matrix. Also not 100%
% independent of class conditions, still need to know trials labels.
    % In this option, CV is performed here.
% Option 3: Perform GLM on single trial using Gaussian basis and use HRF as
% input to classification.
% 4/5/2022: updated to work with HbT GLM fitting.

% TODO:
    % Specify beta weights
        % betaSS is nChn x 2
        % default all ones
    % Set beta weights to 0 if noisy SS channels.


% Note: data_raw is raw data. data_y is conc data.
function data_ynew = hmrR_ssBeta_CV(data_raw,data_y, probe, mlActAuto, tIncAuto,...
    betaSS)

    for iBlk=1:length(data_y)

        y      = data_y(iBlk).GetDataTimeSeries('reshape');
        y_raw  = data_raw(iBlk).GetDataTimeSeries('reshape');
        ml     = data_y(iBlk).GetMeasListSrcDetPairs();
        t      = data_y(iBlk).GetTime();
        data_ynew(iBlk)    = DataClass(data_y(iBlk));
        tInc = tIncAuto{iBlk};
        
        % default values
        SDrange = [0.0, 45.0];
        % dRange = [1e4 1e7];
        dRange = [5e3, 1e12];
        if ~exist('SNRthresh', 'var')
            SNRthresh = 1.6;
        end
        rhoSD_ssThresh = 15;
        
        SrcPos = probe.GetSrcPos();
        DetPos = probe.GetDetPos();
        mlAct = mlActAuto{iBlk};

        lst = 1:size(ml,1);
        rhoSD = zeros(length(lst),1);
        posM = zeros(length(lst),3);
        for iML = 1:length(lst)
            rhoSD(iML) = sum((SrcPos(ml(lst(iML),1),:) - DetPos(ml(lst(iML),2),:)).^2).^0.5;
            posM(iML,:) = (SrcPos(ml(lst(iML),1),:) + DetPos(ml(lst(iML),2),:)) / 2;
        end

        lstSS = lst(find(rhoSD<=rhoSD_ssThresh & mlAct(lst)==1));
        
        % clean lstSS of noisy channels
        % y is time x conc x chn
        % y_raw is time x (chn x wavelength)
        lstSS_raw = [lstSS lstSS*2];
        ySS = y_raw(:,lstSS_raw);
        ymean = squeeze(mean(ySS,1))';
        ystd = squeeze(std(ySS,[],1))';
        rhoSD_raw = [rhoSD(lstSS); rhoSD(lstSS)];
        lst2 = ymean>dRange(1) & ymean<dRange(2) & (ymean./ystd)>SNRthresh & rhoSD_raw>=SDrange(1) & rhoSD_raw<=SDrange(2);
        lst2Len = size(lst2,1);
        lst2Arr = [lst2(1:lst2Len/2) lst2(lst2Len/2+1:end)];
        lst2And = lst2Arr(:,1).*lst2Arr(:,2);
        
        lstSS(~lst2And) = [];

        dc = squeeze(y(:,1,:));
        dc = (dc-ones(length(dc),1)*mean(dc,1))./(ones(length(dc),1)*std(dc,[],1));
        cc(:,:,1) = dc'*dc / length(dc);

        % HbR
        dc = squeeze(y(:,2,:));
        dc = (dc-ones(length(dc),1)*mean(dc,1))./(ones(length(dc),1)*std(dc,[],1));
        cc(:,:,2) = dc'*dc / length(dc);
        
        dc = squeeze(y(:,3,:));
        dc = (dc-ones(length(dc),1)*mean(dc,1))./(ones(length(dc),1)*std(dc,[],1));
        cc(:,:,3) = dc'*dc / length(dc);

        clear dc
        % find short separation channel with highest correlation
        for iML = 1:size(cc,1)
            % HbO
            [foo,ii] = max(cc(iML,lstSS,1));
            iNearestSS(iML,1) = lstSS(ii);
            % HbR
            [foo,ii] = max(cc(iML,lstSS,2));
            iNearestSS(iML,2) = lstSS(ii);
            % iNearestSS is 36 x 2;
            % HbT
            [foo,ii] = max(cc(iML,lstSS,3));
            iNearestSS(iML,3) = lstSS(ii);
        end
        
        lstInc = find(tInc==1);
        
        nT = length(t);
        nCh = size(y,3);
        ynew    = zeros(nT,3,nCh);
        
        if ~exist('betaSS','var')
            betaSS = ones(nCh,3);
        end
        
        % Check noisy channels
    
        for conc = 1:3
            mlSSlst = unique(iNearestSS(:,conc));
            for iSS = 1:length(mlSSlst)
                lstMLtmp = 1:size(ml,1);
                lstML = find(iNearestSS(:,conc)==mlSSlst(iSS) & mlAct(lstMLtmp)==1);

                ytmp = y(lstInc,conc,lstML);
                
                ynew(lstInc,conc,lstML) = ytmp - permute(y(lstInc,conc,mlSSlst(iSS))*betaSS(lstML,conc)',[1 3 2]);

            end
        end
        
        %ynew(:,3,:) = ynew(:,1,:) + ynew(:,2,:);
        
        data_ynew(iBlk).SetDataTimeSeries(ynew);
        
    end
    
end