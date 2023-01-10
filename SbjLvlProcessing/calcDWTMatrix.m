%
% X is #time points x [HbO/HbR/HbT] x #channels
% dc.dataTimeSeries is time x channels
% When dealing with DWT, almost never deal directly with wavelet and
% support functions, only need high and low-pass filters and input
% signal
function WX = calcDWTMatrix(X)
    
    % X = ones(2^10,1);

    % Let matlab takes care of extension or manually extend the signal?
    dwtmode('sym');
    
    % M is support length of the (finest-level?) wavelet function

    wname = 'bior3.9';
    fb = dwtfilterbank('Wavelet',wname);
    spsi = waveletsupport(fb);
    % I still don't understand why the coarsest level is 20, not 4.
    %M = spsi(1);
    % Here it is 20, but in Jung's code, it's 9???
    M = spsi(2);
    % [psi, t] = wavelets(fb);
    % [phi, t2] = scalingfunctions(fb);
    % plot(t,psi');
    % grid on;
    % title('Time-Centered Wavelet');
    % xlabel('Time')
    % ylabel('Magnitude')
    
    Jp = floor(log2(size(X,1)/(M-1)));
    
    

    
    % reconstruction and decomposition filter, which are basically scaling
    % functions
    [RF,DF] = biorwavf(wname);
    % decomposition low/high-pass filter, reconstruction low/high-pass
    % filter. Low-pass has same shape as scaling func, high-pass => wavelet
    % func
    [LoD,HiD,LoR,HiR] = biorfilt(DF,RF);
    % sanity check
    % figure();stem(HiD);figure();stem(LoD);
    % perform single-level DWT, approximation and detail coefficients
    
    for chnIdx = 1:size(X,2)
    
        % really, the number of coefficients is determined by the
        % convolution process, which also depends on the filter length
        tempX = X(:,chnIdx);
        [c,l] = wavedec(tempX,Jp,LoD,HiD);

        cA = appcoef(c,l,wname,Jp);
        cD = cell(1,Jp);
        for i = 1:Jp
            cD{i} = detcoef(c,l,i);
        end

        tempW = cA;
        for i = Jp:-1:1
            tempW = [tempW; cD{i}];
        end
        
        if chnIdx == 1
            % W is time X #channels*#conc 
            WX = zeros(size(tempW,1),size(X,2));
        end
        
        WX(:,chnIdx) = tempW;
    end
    
%     for i = 1:size(X,2)*size(X,3)
%         for j = 1:J
%             [cA,cD] = dwt(X,LoD,HiD);
%             
%         end
%     end
    
    
    
end