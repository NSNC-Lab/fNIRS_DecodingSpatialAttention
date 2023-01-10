% learning wavelet
wname = 'bior3.9';
fb = dwtfilterbank('Wavelet',wname);
[Lo,Hi] = filters(fb);
[RF,DF] = biorwavf(wname);
[LoD,HiD,LoR,HiR] = biorfilt(DF,RF);
% I don't understand how different lvls of wavelets are obtained,
% interpolation?
[psi, t] = wavelets(fb);
[phi, t2] = scalingfunctions(fb);
spsi = waveletsupport(fb);
% wmaxlev give 6, seems like support length is 20 instead of 4. Why?
M2 = spsi(2);
M1 = spsi(1);
X = ones(1578,1);
Jp1 = floor(log2(size(X,1)/(M1-1)));
Jp2 = floor(log2(size(X,1)/(M2-1)));
 
Jmt = wmaxlev(size(X,1),wname);

% % perform convolution w upsampling
% % start with low-pass filter, get high-pass filter
% % then upsample by inserting 0s
% % then convolve repeatedly
% testHi = LoR(end:-1:1);
% for i = 1:length(testHi)
%     if mod(i,2) == 0
%         testHi(i) = -1*testHi(i);
%     end
% end
% 
% testHi = dyadup(testHi,0);
% testHi2 = conv(testHi,LoR);

% % construct finer wavelet one level up
psiTest1 = psi(end,:);
psiTest2 = 2^(1/2).*psiTest1(1:2:end);