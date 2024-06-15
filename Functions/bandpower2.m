function BANDPOWERS = bandpower2(varargin)
% Approximate power in a frequency band by integrating the PSD using the 
% rectangle method assuming that the number of points used to compute the
% FFT is even.
% Input: [Power Spectral Density (PSD),frequency bins centers,frequency range]
% Output: Bandpower for each PSD in frequency range
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% error if not enough or too much input variables are given
narginchk(2,3)

% if no frequency range is given compute power over all frequencies
PSD=varargin{1};
fff=varargin{2};
if nargin<3
    freqrange=[min(fff),max(fff)];
else
    freqrange=varargin{3};
end

% error if frequency range is too low or high
if freqrange(1)<0
    error('Frequency range can not extend to negative values')
elseif freqrange(2)>max(fff)
    error('Frequency range can not extend to values above the Nyquist frequency')
end

% adjust frequency bin center vector if it is not in column form
fff=fff(:);

% adjust DC and Nyquist frequency components by multiplying by 2 so that
% the negatve part of the frequency bin is added
PSD([1,end],:)=PSD([1,end],:)*2;

% get frequency bin width
bin_width=fff(2)-fff(1);

% get widths of bins in selected range and frequency logical for PSD
fff(2:end)=fff(2:end)-bin_width/2;
if isempty(freqrange)
    freqlogical=true(size(fff,1),1);
else
    freqlogical=fff>freqrange(1)&fff<freqrange(2);    
end
index=find(fff<=freqrange(1),1,'last');
fff=fff(freqlogical);
freqlogical(index)=true;
fff=[freqrange(1);fff;freqrange(2)];
fff=diff(fff);

% get PSD values associated with bin widths
PSD=PSD(freqlogical,:);

% calculate area with rectangle method under PSD for selected frequency range
BANDPOWERS=sum(PSD.*fff,1);

end

