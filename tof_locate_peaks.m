function [PeakList, param, Isel, Coef_List_gauss, Soulders_count] = tof_locate_peaks(mz,intensity,param, text_size, posi, x0, tof)
% This function is used to locate peaks in a spectrum by detecting variation in signal derivative and amplitude
% It takes as:
% INPUTS:   
%           mz:     Mass to charge
%           spec:   Intensity values 
%           param:  Detection params defined in the main script
% 
% and provides as:
% OUTPUTS:  
% PeakList - [peak number, position, intensity, width, area, std, fitError]
%            columns:
%            1. peak number: The peak ID
%            2. position: Position of the peak in the mz axis
%            3. intensity of the peak
%            4. width: The full width at half maximum
%            5. area : The area of the peak 
%            6. std : The standard deviation
%            7. fitError: Norm of residuals
% 
% In case of shoulders width, area, std and fitError are atomatically 
% replaced by NaN
%
% param:  parameters used for fitting, if input parameters had 'auto'
% settings, the evaluated values are in output param structure.
%
% Resolution enhancement follows logic:
% enhancedSignal = signal - 1stDerivative*enhanceFactor1 + 
% 4thDerivative*enhanceFactor2
% 4thDerivative - 4thDerivative value where 1stDerivative is positive
%
% This will enhance peaks and shoulders, but will give artfacts in
% baseline. A compromise have to be found manually. The enhanced signal is
% only used for peak detection, not for any peak fitting.

% If the number of given params for this function is 2 :
% Consider the default params bellow
    if nargin == 6

        Range = [ceil(min(mz)),floor(max(mz))];
        Window = floor(Range(2)/10);
        SignalThreshold = 'auto';
        SlopeThreshold = 'auto'; 
        ScaleMedian = 10;
        R = 4000;
        doPlot = 0;
        enhanceFactor1 = 3;
        enhanceFactor2 = 50;
        DifLim = 0.01;
        remTails = 1;

    % Else, take the given params
    elseif nargin == 7

        Range = param.Range;
        Window = param.Window;
        SignalThreshold = param.SignalThreshold;
        SlopeThreshold = param.SlopeThreshold; 
        ScaleMedian = param.ScaleMedian;
        R = param.R;
        doPlot = param.doPlot;
        enhanceFactor1 = param.enhanceFactor1;
        enhanceFactor2 = param.enhanceFactor2;
        DifLim = param.DifLim;
        remTails = param.remTails;

    end
% Possible reading error message
    if max(mz) < Range(1)

        disp('Range exceeds the masses axes')

        if nargout == 1

            varargout{1} = [];

        elseif nargout == 2

            varargout{1} = [];
            varargout{2} = param;

        end

        return
    end
% Recalculating the calibration law coefficients: 
mz1 = mz(1);
mz2 = mz(end);
% tof1 = 0;
% tof2 = length(mz);
tof1 = 0;
tof2 = tof(end);
b = -0.5*(tof1^2/mz1 - tof2^2/mz2)/(tof1/mz1 - tof2/mz2);
a = sqrt(tof1^2/mz1 - 2*b*tof1/mz1);

clear mz1 mz2 tof1 tof2
param.a = a;
param.b = b;

mz_cons = round(mz);
Ind_mz = mz_cons >= Range(1) & mz_cons <= Range(2); % Indices of the considered mz in the specified range in param
mz_cons = mz(Ind_mz);
intensity_cons = intensity(Ind_mz);
clear Ind_mz intensity mz

% Smoothing the signal for initial peak finding
[intensity_filtered,swn] = H_fastsmooth(mz_cons,intensity_cons,param);

intensity_orig = intensity_cons;
intensity_cons = intensity_filtered;

[intensity_enhanced,~,~,~,~] = tof_resolution_enhancement(intensity_cons,enhanceFactor1,enhanceFactor2,remTails);
[~,yd1,yd2,~,~] = tof_resolution_enhancement(intensity_cons,0,0);
yd2 = -yd2;

aSteps = Range(1):Window:Range(2);
nrSteps = length(aSteps);

absDifs_1  = cell(nrSteps,1);

if ischar(SignalThreshold)
    % Auto set amplitude threshold: 5 times (ScaleMedian) non-zero values median difference of first
    % derivative
    ats = NaN(nrSteps,1);
    ats_1 = NaN(nrSteps,1);
    ats_2 = NaN(nrSteps,1);
    ats_3 = NaN(nrSteps,1);
    for i = 1:nrSteps - 1
        Isel = mz_cons > aSteps(i) & mz_cons <= aSteps(i+1);
        absDifs = abs(diff(intensity_orig(Isel)));
        absDifs(absDifs == 0) = [];
        absDifs_1{i} = absDifs;

        ats(i) = nanmedian(absDifs)*ScaleMedian;
        ats_1(i) = nanmedian(absDifs);
        ats_2(i) = nanmean(absDifs);    
        ats_3(i) = nanmean(absDifs)*ScaleMedian;    
    end
    try
        ats = MDrepl1(ats);
        ats_1 = MDrepl1(ats_1);
        ats_2 = MDrepl1(ats_2);
        ats_3 = MDrepl1(ats_3);
    catch
        ats = zeros(nrSteps,1);
        ats_1 = zeros(nrSteps,1);
        ats_2 = zeros(nrSteps,1);
        ats_3 = zeros(nrSteps,1);
    end
    ats(end) = ats(end-1);
    ats_1(end) = ats_1(end-1);
    ats_2(end) = ats_2(end-1);
else
    ats = SignalThreshold(ones(nrSteps,1),:);
    ats_1 = SignalThreshold(ones(nrSteps,1),:);
    ats_2 = SignalThreshold(ones(nrSteps,1),:);
    ats_3 = SignalThreshold(ones(nrSteps,1),:);
end

SignalThreshold = interp1(aSteps,ats,mz_cons);
clear absDifs

if ischar(SlopeThreshold)
    % Auto set slope treshold: 5 times median difference of seconds
    % derivative
    sts = NaN(nrSteps,1);
    for i = 1:nrSteps-1
        Isel = mz_cons > aSteps(i) & mz_cons <= aSteps(i+1);
        absDifs = abs(diff(yd2(Isel)));
        absDifs(absDifs == 0) = [];
        sts(i) = nanmedian(absDifs)*ScaleMedian;
    end
    try
        sts = MDrepl1(sts);
    catch
        sts = zeros(nrSteps,1);
    end
    sts(end) = sts(end-1);
else
    sts = SlopeThreshold(ones(nrSteps,1),:);
end
SlopeThreshold = interp1(aSteps,sts,mz_cons);

clear ats sts absDifs aSteps nrSteps

% Apply the different criteria
Iat = intensity_enhanced(1:end-1) > SignalThreshold(1:end-1);
clear intensity enhanced 
Izc = (yd1(1:end-1) >= 0 & yd1(2:end) <= 0);

% Use original data for slope detection
Ist = yd2(2:end) > SlopeThreshold(2:end);

Isel = find(Izc & Ist & Iat);
Idel = intensity_cons(Isel) < intensity_cons(Isel-1) & intensity_cons(Isel) < intensity_cons(Isel+1);
Isel(Idel)=[];
nrP = length(Isel);
clear Iat Izc Ist

%% Evaluate each peak

PeakList = NaN(nrP,7);
Coef_List = NaN(nrP,5);
Coef_List_gauss = NaN(nrP,3);
S_List = cell(nrP,1);

% Initializataion
Soulders_count = 0;

for p = 1:nrP
    
    m = mz_cons(Isel(p));
    pmPoints = max(round(swn(round(m))/2),1); % Plus minus points
    % Take points around the peak
    Ip = Isel(p) - pmPoints:Isel(p) + pmPoints;
    yy = intensity_orig(Ip);
    [~,Imx] = max(yy);
    Imx = Imx + Ip(1) - 1;
    Ip = Imx - pmPoints:Imx + pmPoints;
    yy = intensity_orig(Ip);
    
    xx = mz_cons(Ip);
    
    xx = xx(:);
    yy = yy(:);
    
    % Gaussian fit
    
    Tx_list = linspace(xx(1), xx(end), 100);
    gi = griddedInterpolant(xx, yy);   
    Ty_list = gi(Tx_list);

    [Coef_List_gauss(p,1) Coef_List_gauss(p,2) Coef_List_gauss(p,3)] = mygaussfit(Tx_list,Ty_list);
    
    [coef,S,MU] = polyfit(xx,log(abs(yy)),2);  % Fit parabola to log of sub-group with centering and scaling
    S_List{p} = S;
    if (isnan(Coef_List_gauss(p,1)) == 0) && (isnan(Coef_List_gauss(p,2)) == 0) && (isnan(Coef_List_gauss(p,3)) == 0) && (isreal(Coef_List_gauss(p,1)) == 1) && (isreal(Coef_List_gauss(p,1)) == 1) && (isreal(Coef_List_gauss(p,2)) == 1) && (isreal(Coef_List_gauss(p,3)) == 1)
        c1 = coef(3);
        c2 = coef(2);
        c3 = coef(1);
        PeakX = Coef_List_gauss(p,2);
        PeakY = exp(c1-c3*(c2/(2*c3))^2);
        width = Coef_List_gauss(p,3)*2*sqrt(2*log(2));
        std = Coef_List_gauss(p,3);
        Area = sum(Coef_List_gauss(p,1)*exp(-((mz_cons(Isel(p)) - PeakX).^2)/(2*std.^2)));
        PeakList(p,:) = [p PeakX PeakY width Area std Coef_List_gauss(p,1)]; % S.normr
        Coef_List(p,:) = [c1 c2 c3 MU(1) MU(2)]; 
        
        if abs(mz_cons(Isel(p))-PeakX) > DifLim
            Soulders_count = Soulders_count + 1;
            PeakList(p,:) = [p mz_cons(Isel(p)) intensity_orig(Isel(p)) NaN NaN NaN NaN];
            Coef_List(p,:) = [NaN NaN NaN NaN NaN];
            Coef_List_gauss(p,:) = [NaN NaN NaN];
        end
    else
        PeakList(p,:) = [p mz_cons(Isel(p)) intensity_orig(Isel(p)) NaN NaN NaN NaN];
        Coef_List(p,:) = [NaN NaN NaN NaN NaN];
        Coef_List_gauss(p,:) = [NaN NaN NaN];
    end  
end

clear p m pmPoints Ip yy Imx xx xx c1 c2 c3 S MU PeakX PeakY width std Area swn 

% G: Additional verifications
Ibad = isnan(PeakList(:,2))|isnan(PeakList(:,3))|isinf(PeakList(:,2))|isinf(PeakList(:,3))|PeakList(:,2)<0|PeakList(:,4)>=1|PeakList(:,3)<0;% |PeakList(:,3)>1e5
PeakList(Ibad,:) = [];
Coef_List(Ibad,:) = [];
Coef_List_gauss(Ibad,:) = [];
Isel(Ibad,:) = [];
clear Ibad
muh = size(PeakList);
i = 0;
while(i < muh(1)-1)
    i = i + 1;
    if PeakList(i,2) >= PeakList(i+1,2)
        PeakList(i,:) = [];
        Coef_List(i,:) = [];
        Coef_List_gauss(i,:) = [];
        Isel(i,:) = [];
        i = i - 1;
    end
    muh = size(PeakList);
end
clear muh

if doPlot
    
    mz_min = param.Range(1);
    mz_max = param.Range(2);
    
    tof_min = round((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_min))/(2*x0(3)));
    tof_max = round((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_max))/(2*x0(3)));
    
    %%%%%%%%%%%%%%%%%%
    number = 5;
    figure(),
    clf
    
    at = axes;
    
    xlabel(at(1), 'Time of flight [ns]', 'Interpreter', 'latex', 'FontSize', text_size - 2);
    ylabel(at(1), 'Ion count','Interpreter', 'latex', 'FontSize', text_size - 2);
    
    hold(at(1),'on')
    plot(at(1), ((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_cons))/(2*x0(3)))/10,intensity_orig,'-', 'color', 'b', 'linewidth', 2), hold on
    plot(at(1), ((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_cons))/(2*x0(3)))/10,SignalThreshold, 'color', 'k', 'linewidth', 2), hold on 
    plot(at(1), ((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_cons(Isel)))/(2*x0(3)))/10,intensity_orig(Isel),'o', 'color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
     
    legend(at(1), 'Original spectrum','Signal threshold','Detected preaks','Interpreter', 'latex', 'Location', 'NorthEast','FontSize',text_size-2) % ,'Filtered spectrum'
    
    at(1).XAxis.MinorTick = 'on';
    at(1).FontSize = text_size - 2;
    at(1).XLim = [(-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_min))/(2*x0(3)) (-x0(2)+ sqrt(x0(2)^2+4*x0(3)*mz_max))/(2*x0(3))]/10;
    grid(at(1), 'on')

    poso = get(at,'position');
    pos1 = poso(2);
    poso(2) = poso(2)*2; poso(4) = poso(4)-(pos1);
    set(at,'position',poso)
    poso(2) = pos1; poso(4) = .001;
    at(2) = axes('position',poso,'color','none'); 
    
    xlabel(at(2), 'm/z', 'Interpreter', 'latex', 'FontSize', text_size - 2);
    at(2).FontSize = text_size - 2;
    at(2).XAxis.MinorTick = 'on';
    at(2).XLim = [mz_min mz_max];
    grid(at(2), 'on')
    
    set(gca,'Fontname','Times','Fontsize',text_size - 2, 'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    set(gcf,'Position',posi)
    ay = gca;
    ay.GridAlpha = .25;
        
end
if doPlot
    disp(['Number of located peaks: ', num2str(nrP)])
end
clear intensity_cons intensity_filtered yd1 yd2 ydt2 mz_cons intensity_orig nrP


if nargout == 1
    varargout{1} = PeakList;
elseif nargout == 4
    
    param.Range = Range;
    param.Window = Window;
    param.SignalThreshold = SignalThreshold;
    param.SlopeThreshold = SlopeThreshold; 
    param.ScaleMedian = ScaleMedian;
    param.R = R;
    param.doPlot = doPlot;
    param.enhanceFactor1 = enhanceFactor1;
    param.enhanceFactor2 = enhanceFactor2;
    param.DifLim = DifLim;
    param.remTails = remTails;
    
    varargout{1} = PeakList;
    varargout{2} = param;
    varargout{3} = Isel;
    varargout{4} = Coef_List_gauss;
    clear intensity mz 
    
end

