function [Peak_List_each_spectra, Soulders_count] = Peak_locating_routine(tof_mz_intensity_aligned,param,i,text_size,posi,x0)

    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    % Before starting the calibration procedure, Peaks need to be located first.
    % This step uses peak locate feature provided by Heikki Junninen PhD
    % work. Locating peaks is based on detecting variation in signal derivative
    % and amplitude.
    % Peak_locate function takes as:
    % INPUTS:   
    %           mz          : Mass to charge
    %           Intensity   : Intensity (Total ion count) 
    %           Param       : Detection parameters % (Pre-defined at the begining of the script - see Experimental setup section)
    %
    % and provides as:
    % OUTPUTS:  
    % PeakList - [peak number, position, intensity, width, area, std, fitError]
    %            columns:
    %            1. Peak number: The peak ID
    %            2. Position: Position of the peak in the mz axis
    %            3. Intensity (Total ion count)
    %            4. Width: The full width at half maximum
    %            5. Area : The area of the peak 
    %            6. Std : The standard deviation
    %            7. FitError: Norm of residuals
    % 
    % This function applies a signal enhacement to try to get peaks and shoulders in tha mass spectra 
    % In case of shoulders width, area, std and fitError are atomatically replaced by NaN
    %
    % This will enhance peaks and shoulders, but will give artfacts in
    % baseline. A compromise have to be found manually. The enhanced signal is
    % only used for peak detection, not for any peak fitting.
    %
    % Resolution enhancement follows logic:
    % enhancedSignal = signal - 1stDerivative*enhanceFactor1 + 4thDerivative*enhanceFactor2
    % 4thDerivative - 4thDerivative value where 1stDerivative is positive
    if param.doPlot
        disp(['Mass spectrum number : ', num2str(i)])
    end
    [Peak_List_each_spectra,~, Isel_each_spectra, Coef_List_gauss, Soulders_count] = tof_locate_peaks(tof_mz_intensity_aligned(:,2),tof_mz_intensity_aligned(:,3),param, text_size, posi,x0, tof_mz_intensity_aligned(:,1));
        
    Peak_List_each_spectra = [Peak_List_each_spectra,...
                                tof_mz_intensity_aligned(Isel_each_spectra,1),...
                                round(Peak_List_each_spectra(:,2)),...
                                Peak_List_each_spectra(:,2) - round(Peak_List_each_spectra(:,2)),...
                                log10(Peak_List_each_spectra(:,3))*2];
                            
     % During the process of peak finding some rows have been deleted
     % creating a discontinuity in the Peaks ID in each table
     % This is corrected here
     if isequal(find(Peak_List_each_spectra(:,1)), Peak_List_each_spectra(:,1)) == 0
        if param.doPlot
            disp('ID rearranged this time')
        end
        Peak_List_each_spectra(:,1) = find(Peak_List_each_spectra(:,1));
        Coef_List_gauss(:,1) = find(Peak_List_each_spectra(:,1));
     end
     
end