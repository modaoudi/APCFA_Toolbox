function [intesity_returned,tof] = H_fastsmooth(mz,intensity,param)

% Smoothing the signal using fastsmooth, 
% and changing the smoothing window 
% for each mz according to resolution
% The smoothing in done in 3 different windows separatly 
% then an average weighted smoothing is done

% So one needs:
% param.Range
% param.a
% param.b
% param.R

Range = param.Range;

a = param.a;
b = param.b;
R = param.R;
cases = 4;

tof = NaN(Range(2),1);
Window = round((Range(2) - Range(1))/(cases - 1));

h = 0;
intensity_tmp = NaN(length(intensity),cases);
weight = intensity_tmp;
for m = Range(1):Window:Range(2)
    h = h + 1;
    [~,Isel] = min(abs(mz - m));
    Width_estim = m/R; % Estimated peak width
    tof1 = sqrt(m - Width_estim)*a+b;
    tof2 = sqrt(m + Width_estim)*a+b;
    if m == 0
        m = m + 1;
    end
    tof(m) = max(round((tof2-tof1)/2),3);
    
    % Works better if window is odd
    if ~rem(tof(m),2)
        tof(m) = tof(m) + 1;
    end
    intensity_tmp(:,h) = tsmooth(intensity,tof(m));
    weight(Isel,:) = 0;
    weight(Isel,h) = 1;
end

weight(1,:) = 0;
weight(1,1) = 1;

weight(end,:) = 0;
weight(end,end) = 1;

weight = MDrepl1(weight);
if isnan(tof(1))
    tof(1) = tof(find(~isnan(tof),1));
end
intesity_returned = sum(intensity_tmp.*weight,2)./sum(weight,2);

