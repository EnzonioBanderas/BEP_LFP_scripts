function [c_out,t_out,l_out] = corrgram2(varargin)
%   CORRGRAM2 Calculate windowed cross correlation between two signals.
%   C = CORRGRAM2(A,B,MAXLAG,WINDOW,NOVERLAP,FS) calculates the windowed cross   
%   correlation between the signals in vector A and vector B. CORRGRAM2 splits  
%   the signals into overlapping segments and forms the columns of C with 
%   their cross correlation values up to maximum lag specified by scalar 
%   MAXLAG. Each column of C contains the cross correlation coefficients between 
%   windows of signals A and B for different delays. Time increases linearly 
%   across the columns of C, from left to right.  Lag increases linearly down 
%   the rows, starting at -MAXLAG. If lengths of A and B differ, the shorter
%   signal is filled with zeros. If N is the length of the signals, C is
%   a matrix with 2*MAXLAG+1 rows and 
%         NWINDOW = floor((N-NOVERLAP)/(WINDOW-NOVERLAP)) 
%   columns. If the sampling rate FS is given as an input CORRGRAM2 uses 
%   seconds instead of sample numbers in its output.
%
%   [C,L,T] = CORRGRAM2(...) returns a column of lag L and one of time T
%   at which the correlation coefficients are computed. L has length equal 
%   to the number of rows of C, T has length NWINDOW.
%
%   C = CORRGRAM2(A,B) calculates windowed cross correlation using defeault
%   settings; the defeaults are MAXLAG = floor(0.1*N), WINDOW = floor(0.1*N), 
%   NOVERLAP = 0 and no FS (lags and time are given in sample numbers instead of seconds). 
%   You can tell CORRGRAM2 to use the default for any parameter by leaving 
%   it off or using [] for that parameter, e.g. 
%   CORRGRAM2(A,B,[],1000).
%
%   CORRGRAM2(A,B) with no output arguments plots the windowed cross 
%   correlation using the current figure.
%
%   EXAMPLE 1:
%       x = cos(0:.01:10*pi)';
%       y = sin(0:.01:10*pi)' + .5 * randn(length(x),1);
%       corrgram2(x,y)
%
%   EXAMPLE 2:
%       fs=1000; % sampling rate
%       t=0:1/fs:30; % 30 s signal
%       delay=.08; % 0.08 s delay
%       w(1)=(2*pi)/(1/3.5); % 3.5 Hz
%       w(2)=(2*pi)/(1/8); % 8 Hz
% 
%       % Input
%       Signal1=3.1.*sin(w(1)*t)+4.1*sin(w(2)*t); % sum of two sines with different frequencies and amplitudes
%       Signal2=1.2*sin(w(1)*(t-delay))+6.7*sin(w(2)*(t-delay)); % delayed and differently scaled version of Signal1
%       nMaxLag=1*fs; % 1 s maximum delay
%       nWindow=2*fs; % 2 s windows
%       nOverlap=0.5*nWindow; % 50% overlap
% 
%       corrgram2(Signal1,Signal2,nMaxLag,nWindow,nOverlap,fs)
%
%   See also CORRCOEF, CORR, XCORR, MIGRAM.
% 
% Copyright (c) 2007 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2007/06/14 08:52:19 $
% $Revision: 5.1 $
%
%   Revision by Enzo Nio 2018/06/26
%   The output and default input remains roughly the same as the previous corrgram,
%   but all the calculations in between are changed:
%
%   Changed calculation of cross correlation coefficients to be faster with 
%       C=ifft(conj(fft(A_WINDOW)).*fft(B_WINDOW))
%   Where A_WINDOW and B_WINDOW are possibly overlapping windows of A and B
%   defined by WINDOW and NOVERLAP. A_WINDOW is zeropadded with 2*MAXLAG
%   zeros, B_WINDOW is possibly zeropadded if no values are given by B 
%   (there are no values of B to compute cross correlation coefficients with 
%   A at the first and last windows for -MAXLAG:1 and 1:MAXLAG respectively)
%
%   This is equivalent to doing a circular convolution of a window of A 
%   (number of non-zero values is WINDOW) with a window of B 
%   (number of non-zero values is 2*MAXLAG++WINDOW) for
%   -MAXLAG:MAXLAG+WINDOW-1 where the last WINDOW-1 values can be discarded
%   after the convolution.
%
%   The cross correlation matrix is normalized by setting the
%   autocorrelation coefficients of the windows of A and B involved in 
%   computing the cross correlation coefficients for each delay to 1. 
%   The B autocorrelation coefficients change with delay and are in a matrix 
%   with dimensions corresponding to the dimensions of C, whereas the A 
%   autocorrelation coefficients remain constant with delay (because of the 
%   2*MAXLAG zero padding) and have 1 row and NWINDOW columns.
%
%   license.txt included at end of function
%
%% check input and inital setting of parameters
narginchk(2,6)
S1 = varargin{1}; S2 = varargin{2};
S1 = S1(:); S2 = S2(:);

nS1 = length(S1); nS2 = length(S2);
if nS1 < nS2    % zero-pad x if it has length less than y
    S1(nS2) = 0; nS1 = nS2;
end

if nS2 < nS1    % zero-pad y if it has length less than x
    S2(nS1) = 0;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    nMaxLag = varargin{3};
    if nMaxLag < 0, error('Requires positive integer value for maximum lag.'), end
    if length(nMaxLag) > 1, error('Requires MAXLAG to be a scalar.'), end
else
    nMaxLag = floor(nS1/10);
end

if length(varargin) > 3 && ~isempty(varargin{4})
    nWindowLength = varargin{4};
    if nWindowLength < 0, error('Requires positive integer value for window length.'), end
    if length(nWindowLength) > 1, error('Requires WINDOW to be a scalar.'), end
else
    nWindowLength = floor(nS1/10);
end

if length(varargin) > 4 && ~isempty(varargin{5})
    nOverlap = varargin{5};
    if nOverlap < 0, error('Requires positive integer value for NOVERLAP.'), end
    if length(nOverlap) > 1, error('Requires NOVERLAP to be a scalar.'), end
    if nOverlap >= nWindowLength, error('Requires NOVERLAP to be strictly less than the window length.'), end
else
    nOverlap = 0;
end

if length(varargin) > 5
    fs = varargin{6};
    if fs <= 0, error('Requires positive value for FS.'), end
else
    fs=[];
end

%% calculate cross correlation coefficient for each delay and each window

% lag and number of lag values
l=(-nMaxLag:nMaxLag)';
nLag=length(l); 

% overlap S1 segments
[S1,S1_trunc]=buffer(S1,nWindowLength,nOverlap,'nodelay');

% number of windows and center points of windows
nWindow=size(S1,2);
t=(floor(nWindowLength/2):nWindowLength-nOverlap:floor(nWindowLength/2)+(nWindow-1)*(nWindowLength-nOverlap))';

% add zeros to each S1 window for circular convolution
S1=[S1;zeros(2*nMaxLag,nWindow)];

% fill in missing values of S2 with zeros (e.g. nMaxLag on both sides)
S2=[zeros(nMaxLag,1);S2;zeros(nMaxLag-length(S1_trunc),1)];

% overlap S2 segments, zeros do not need to be added
[S2,~]=buffer(S2,nWindowLength+2*nMaxLag,nOverlap+2*nMaxLag,'nodelay');

% do circular convolution in frequency to compute non-normalized cross correlation coefficients
c=ifft(conj(fft(S1)).*(fft(S2)));

% truncate c to only include -nMaxLag:nMaxLag values
c=c(1:nLag,:);

% compute autocorrelation coefficients for S1
S1=S1.^2;
S1_auto=sum(S1);
S1_auto=repmat(S1_auto,[nLag,1]);
for i=1:nMaxLag
    S1_auto(nMaxLag+1-i,1)=S1_auto(nMaxLag+2-i,1)-S1(i,1);
    S1_auto(nMaxLag+1+i,end)=S1_auto(nMaxLag+i,end)-S1(nWindowLength+1-i,end);
end

% compute autocorrelation coefficients for nWindow windows of S2 for different delays
S2_auto=zeros(nLag,nWindow);
S2=S2.^2;
for i=1:nWindow
S2_auto(1,i)=sum(S2(1:nWindowLength,i)); 
    for ii=2:nLag
    S2_auto(ii,i)=S2_auto(ii-1,i)-S2(ii-1,i)+S2(ii-1+nWindowLength,i);  
    end
end

% normalize cross correlation coefficients
c=c./sqrt(S1_auto.*S2_auto);

%% display results if there are no output arguments and give output arguments otherwise
% if sampling rate is given convert sample number to seconds for lag and
% time vectors
if ~isempty(fs)
    l=l/fs;
    t=t/fs;
end

if nargout == 0
    newplot
    imagesc(t, l, c)
    xlabel('Time'), ylabel('Lag'), axis xy
    title('Windowed cross correlation', 'fontweight', 'bold')
    colorbar
    set(gca,'CLim',[-1,1])
elseif nargout == 1
    c_out = c;
elseif nargout == 2
    c_out = c;
    t_out = t;
elseif nargout == 3
    c_out = c;
    t_out = t;
    l_out = l;
end

end

% license.txt : 
% Copyright (c) 2009, Norbert Marwan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Potsdam Institute for Climate Impact Research Germany nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.