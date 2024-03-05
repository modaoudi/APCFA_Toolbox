function [enhanced_intensity,yd1,yd2,yd3,yd4] = tof_resolution_enhancement(intensity,enhanceFactor1,enhanceFactor2,remTails)
%
% Enhancing resolution of spectrum using derivative method
%
% [enhancedSpec,yd1,yd2,yd3,yd4] = tof_resolution_enhancement(spec',enhanceFactor1,enhanceFactor2,remTails)
%
% enhanceFactor1 - weighting factor for the second derivative addition.
% enhanceFactor2 - weighting factor for the forth derivative addition (when the second derivative is positive).
% remTails       - remove peak tails
% Bigger value increases the effect. Compromise between increased noise and
% enhancement have to be found experimentally
% enhancedSpec - resolution enhanced spectum
% yd1 - first derivative of the spectrum
% yd2 - second derivative of the spectrum
% yd3 - third derivative of the spectrum
% yd4 - forth derivative of the spectrum, only where yd2 is positive
%
% Resolution enhancement follows logic:
% enhancedSignal = signal - 2 nd Derivative*enhanceFactor1 + 4 th Derivative*enhanceFactor2
% 2 th Derivative = 2 th Derivative value where 2 nd Derivative is positive and
%                        signal value/enhanceFactor1 where negative (will remove peak tails)
% 4 th Derivative = 4 th Derivative value where 2 nd Derivative is positive
%
% based on ideas of T. C. O'Haver,2006
% Heikki Junninen
% March 2011

if nargin == 3
    remTails = 0;
end

% First derivative
yd0 = [0;intensity(2:end)-intensity(1:end-1)];
yd1 = tsmooth(yd0,3); % To reduce noise

% Second derivative
yd02 = [0;yd1(2:end)-yd1(1:end-1)];
yd2 = tsmooth(yd02,3); % To reduce noise

% Third derivative,
yd03 = [0;yd2(2:end)-yd2(1:end-1)];
yd3 = tsmooth(yd03,5); % To reduce noise

% Forth derivative,
yd04 = [0;yd3(2:end)-yd3(1:end-1)];
yd4 = tsmooth(yd04,5); % To reduce noise

% Chose only values where second derivative is positive
Ipos = yd2 > 0;
% Removing peak tails 
if remTails
    yd2(Ipos) = intensity(Ipos)/enhanceFactor1;
end

yd4(Ipos) = 0;

% Enhancing the signal
enhanced_intensity = intensity - yd2*enhanceFactor1 + yd4*enhanceFactor2;
