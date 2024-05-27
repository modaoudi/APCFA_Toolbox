% This is the main routine of the APCFA toolbox 
%
% APCFA toolbox for MATLAB
% version 1.0 - March 2024
%
% LICENCE:
% APCFA toolbox © 2024 is licensed under CC BY-NC-SA 4.0
% 
% This APCFA toolbox is distributed with an Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) licence: https://creativecommons.org/licenses/by-nc-sa/4.0/
% 
% You are free to:
%
% Share — copy and redistribute the material in any medium or format
% Adapt — remix, transform, and build upon the material
% 
% The licensor cannot revoke these freedoms as long as you follow the license terms.
% Under the following terms:
%
% Attribution — You must give appropriate credit , provide a link to the license, and indicate if changes were made . You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
% NonCommercial — You may not use the material for commercial purposes .
% ShareAlike — If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.
% No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
%
% REFERENCE:
% The toolbox is freeware and may be used if proper reference is given to the authors, preferably refer to the following paper:
% A mass defect-based approach for the automatic construction of peak lists for databases of mass spectra with limited resolution. Application to ToF-SIMS data
%  https://doi.org/10.1002/rcm.9777
%
% By M. DAOUDI(1,2), N. NUNS(3), P. SCHIFFMANN(2), A.FROBERT(2), B.
% HANOUNE(1), P. DESGROUX(1), A. FACCINETTO(1)
%
% (1) Univ. Lille, CNRS, UMR 8522, PC2A, F59000 Lille, France
% (2) IFP Energies Nouvelles, Institut Carnot IFPEN TE, F-92852 Rueil-Malmaison, France
% (3) Univ. Lille, M. E. Chevreul Institute, F 59000 Lille, France
%
% HELP:
% For help, please adress your request to the corresponding author
% Corresponding author e-mail: alessandro.faccinetto@univ-lille.fr
%
% Credits
% Some features developed by Junninen et al., 2010 (doi:10.5194/amt-3-1039-2010) were used in this work

%% Clean previous work and reading/saving Directories paths
clearvars
close all force
clc

%% The post-processing starts from here  
% In this section paths must be filled in for :
% Scripts location for mass spectra post-processing
% Scripts location for principal component analysis
% Raw ".txt" data files location
% Storage location for mass spectra post-processing data

disp('..................................................................................')
disp('Reading paths...')

% Specify the directory for scripts
Script_path = '/APCFA_toolbox';
% Specify the directory for input data
Data_path = '/Data_Input';
% Specity the directory for output date (storage)
Storage_path = '/Data_Output';

% Saving paths
cd(Storage_path)                                                           % Going to storage location
disp('Saving paths...')
name_save = ('Paths');
save(name_save,'Script_path','Data_path','Storage_path')
disp('Saved')           
clc
diary('Report.txt')
disp('..................................................................................')
disp('APCFA post processing started')
disp('..................................................................................')

cd(Script_path)                                                            % Going to storage location
pause(1)
%% Setting general paramaters
% In this section the user may choose/select some general parameters 
% relative to plots, text font, saving data options and so on ...

disp('..................................................................................')
disp('Setting general parameters...')
% Specifying  numeric format
format long
% Defining general parameters for plots
text_size  = 18;                                                           % Font size for text in plots
left       = 1;                                                            % South east corner (x) 
bottom     = 1;                                                            % South west corner (y) 
width      = 700;                                                          % plot width (x)
height     = 650;                                                          % plot height (y)
% Position of plots in the screen
posi = [left bottom width height]; 
% Saving options
sv = 0;                                                                    % 0 : Not saving plots (.png and .fig format) and data (.mat format)
                                                                           % 1 : Saving plots (.png and .fig format) and data (.mat format)
       
sv_plus = 0;                                                               % 0 : Not saving data in the form of excell sheets
                                                                           % 1 : Saving data in the form of excell sheets

% Waiting time (in seconds), used during the aligmenet and calibration 
% routines under user assisted mode
time = 15;                                                                 % in seconds

% Saving general parameters
if sv 
    cd(Storage_path)
    disp('Saving general parameters...')                                   % Going to storage location
    name_save = ('General_parameters');
    save(name_save,'text_size','left','bottom','width','height','posi','sv','sv_plus','time')
    disp('Saved')
    cd(Script_path)                                                        % Going back to scripts location
end
pause(1)
%% Data parameters 
% In this section the user may choose parameters for peak picking/locating
% routine

disp('..................................................................................')
disp('Reading data characteristics...')

Binning = 4;                                                               % Data Binning

if sv 
    cd(Storage_path)
    disp('Saving data characteristics...')                                 % Going to storage location
    name_save = ('Data_characteristics');
    save(name_save,'Binning')
    disp('Saved')
    cd(Script_path)                                                        % Going back to scripts location
end
pause(1)
%% Peak picking parameters 
% In this section the user may choose parameters for peak picking/locating
% routine

disp('..................................................................................')
disp('Reading peak picking parameters...')
param = struct('Range',[0 600],...                                         % The range where to find peaks
               'Window',.2,...                                             % The window where peak and noise detection is calculated each time
               'SignalThreshold','auto',...                                % Signal threshold determined automatically in the window defined above. Threshold = noise(signal)*ScaleMedian
               'SlopeThreshold','auto',...                                 % Slope threshold ( = signal threshold of 1st derivative) determined automatically in the window defined above. Threshold = noise(signal 1st derivative)*ScaleMedian
               'ScaleMedian',10,...                                        % Multiplication of noise as threshold for peaks in the window
               'R',10000,...                                               % Used to evaluate smoothing window - it's a smoothing parameter and should be close as much as possible to the real resolution
               'doPlot',1,...                                              % Ploting or not the peak find raw result,
               'enhanceFactor1',3,...                                      % Factor for resolution enhancement that improves shoulder detection,
               'enhanceFactor2', 50,...                                    % Additional factor for resolution enhancement
               'DifLim',0.01,...                                           % If the polynome fpoluitted peak location is more than DifLim away from measured peak location: the peak is considerered a shoulder
               'remTails',0);                                              % remove peak tails in the signal enhance phase
% Saving peak locating parameters
if sv
    cd(Storage_path)                                                       % Going to storage location
    disp('Saving peak locating paramaters...')
    name_save = ('Peak_locating_parameters');
    save(name_save,'param')
    disp('Saved')
    cd(Script_path)                                                        % Going back to scripts location
end
pause(1)
%% Peak calibration parameters 
% In this section the user may choose parameters for peak calibration
% routine

disp('..................................................................................')
disp('Reading peak calibration parameters...')

% Calibration constants: Initialization

A = 0;                                                                     % Initial vale
B = 4.580203692436845e-07;                                                 % Initial value
C = 9.157148809116324e-10;                                                 % Initial value

% Intialization of the leastsq fit 
x0 = [A, B, C];

if sv 
    cd(Storage_path)                                                       % Going to storage location
    disp('Saving peak calibration paramaters...')
    name_save = ('Peak_calibration_initialization');
    save(name_save,'x0')
    disp('Saved')
    cd(Script_path)                                                        % Going back to scripts location
end
pause(1)
%% Peak finding parameters
% In this dsection the user may chose parameters for peak finding routine

disp('..................................................................................')
disp('Reading peak finding parameters...')

opt = 2;                                                                   % 1 : Interrogation window optimization
                                                                           % 2 : Suggest an attribution

max_iter_cluster_search = 100;                                             % Maximum iterations for local threshold variation

reg = 0;                                                                   % 1 : Regenerating the combination of elements 
                                                                           % 0 : Use the available (previously generated) combination of elements
% Saving peak locating parameters
if sv
    cd(Storage_path)                                                       % Going to storage location
    disp('Saving peak locating paramaters...')
    name_save = ('Peak_finding_parameters');
    save(name_save,'max_iter_cluster_search', 'opt', 'reg')
    disp('Saved')
    cd(Script_path)                                                        % Going back to scripts location
end
pause(1)
%% Reading data files names 
% In this section names of the data file are collected and ordered 

disp('..................................................................................')
% Select data directory
cd(Data_path)
disp('Reading names of available mass spectra files...')
% Getting the all files names
Files_list_temp = dir('*.TXT');
% Number of mass spectra to process 
n_spec = size(Files_list_temp,1); 
jetcustom = flip(jet(n_spec));                                             % Setting the color variation
disp(['Number of mass spectra detected : ', num2str(n_spec)])
% Sorting file names 
[~, reindex] = sort(str2double(regexp( {Files_list_temp.name}, '\d+', 'match', 'once' )));
Files_list = Files_list_temp(reindex);
clear Files_list_temp
% Saving some variable
if sv
    cd(Storage_path)
    disp('Saving files parameters...')
    name_save = ('Files_parameters');
    save(name_save,'Files_list','n_spec','jetcustom')
    disp('Saved')
    cd(Script_path)
end
pause(1)
%% Pre-conditioning routine
% In this section raw data is imported and pre-conditioned
% Raw data in 'tof_mz_intensity_raw' is structed as follows:
% Time of Flight = 1st column 
% m/z            = 2nd column (Pre-calibrated during the measurement phase)
% Intensity      = 3rd column

tic
disp('..................................................................................')
% Select data directory
cd(Data_path)
disp('Pre-conditioning routine launched...')
% Empty cell to store data
data = cell(size(Files_list));
pre_cal_params = zeros(size(Files_list,1),3);
% Empty cell to store raw data
tof_mz_intensity_raw = cell(n_spec,1);
data_reshape = cell(n_spec,1);
line_skip = 3;                                                             % Initialization : 3 since the first 3 lines are text lines
line_skip = inputdlg({'Please enter the number of text rows to skip in raw .txt files'},...
                                  'Number of rows to skip', [1 50]);
line_skip = str2num(cell2mat(line_skip));
% Raw data do not begin from 0 in time of flight scale
% zeros are added at the begining. The number of missing points might be 
% either even or odd, hence the if condition in the loop
disp('Filling missing values (at the beginning)...')
for i = 1 : n_spec 
    disp(['File: ' Files_list(i).name])
    data{i} = importdata( Files_list(i).name,'\t',line_skip);              
    data_reshape{i} = [data{i}.data(:,1), data{i}.data(:,2), data{i}.data(:,3)];
    pre_cal_params(i,:) = x0;
    clear temp
    if data{i}.data(1,1) == 0
        disp('Not needed this time')
    else
        if mod(data{i}.data(1,1),2) == 0
            tof = [floor(linspace(0,data{i}.data(1,1)-Binning,floor(data{i}.data(1,1)/Binning))');data{i}.data(:,1)];   % For time of flight column
            mz = [polyval(flip(pre_cal_params(i,:)), floor(linspace(0,data{i}.data(1,1)-Binning,floor(data{i}.data(1,1)/Binning))'));data{i}.data(:,2)]; 
            intensity = [zeros(floor(data{i}.data(1,1)/Binning),1);data{i}.data(:,3)];                                  % For intensity column
        else
            tof = [0;floor(linspace(1,data{i}.data(1,1)-Binning,floor(data{i}.data(1,1)/Binning))');data{i}.data(:,1)]; % For time of flight column
            mz = [0;polyval(flip(pre_cal_params(i,:)), floor(linspace(1,data{i}.data(1,1)-Binning,floor(data{i}.data(1,1)/Binning))'));data{i}.data(:,2)];
            intensity = [0;zeros(floor(data{i}.data(1,1)/Binning),1);data{i}.data(:,3)];                                % For intensity column
        end  
        tof_mz_intensity_raw{i} = [tof, mz, intensity];   
        disp('Needed this time')
    end
    clear tof mz intensity
end
% All mass spectra do not have the same data length 
% Zeros are added at the end as well to each spectrum
% to ease the post-processing later
data_length = zeros(n_spec,1);
for i = 1 : n_spec
    data_length(i) = size(tof_mz_intensity_raw{i},1);
end
max_length = max(data_length);
disp('Filling with zeros (at the end)...')
for i = 1 : n_spec
    disp(['File' Files_list(i).name])
    if size(tof_mz_intensity_raw{i},1) == max_length
        disp('No need this time')
    else
        m = size(tof_mz_intensity_raw{i},1);
        tof = linspace(tof_mz_intensity_raw{i}(end,1)+ Binning,tof_mz_intensity_raw{i}(end,1) + Binning*(max_length - m + 1),max_length - m + 1);
        mz = polyval(flip(pre_cal_params(i,:)), tof);
        intensity = zeros(max_length - m + 1,1);
        added = [tof', mz', intensity];
        clear tof mz intensity
        tof_mz_intensity_raw{i} = [tof_mz_intensity_raw{i}; added(1:end-1,:)];
        clear added 
        disp('Needed this time')
        clear m

    end
end
if sv_plus
    disp('Saving files...')
    disp('This might take few moments ...')
    cd(Storage_path)
    % Saving tof_mz_intensity_raw array in .xlsx format
    for i = 1 : n_spec
        disp(['Saving spectrum number: ', num2str(i)])
        xlswrite(['tof_mz_intensity_raw_' Files_list(i).name(1:end-4) '.xlsx'], cat(1, {'tof','mz','intensity'}, num2cell(tof_mz_intensity_raw{i})));
    end 
    disp('Saved')
    cd(Script_path)
end
clear max_length data_length Files_list_temp reindex
if sv
    cd(Storage_path)
    disp('Saving preconditioning routine parameters...')
    name_save = ('Pre_conditioning_parameters');
    save(name_save,'line_skip','pre_cal_params')
    disp('Saved')
    cd(Script_path)
end
disp('Pre-conditioning routine done')
disp('..................................................................................')
toc
pause(1)
%% Quick view of raw data
% In this section raw data is plotted to have a quick view 
tic
disp('..................................................................................')
disp('Quick view of raw data')
cd(Storage_path)
figure()
clf
for i = 1 : n_spec
    plot(tof_mz_intensity_raw{i}(:,1)/10, tof_mz_intensity_raw{i}(:,3), 'Color', jetcustom(i,:), 'Linewidth', 1.5);    % dividing by 10 since 10 points /ns
    hold all
end
grid on
ylabel('Ion count','Interpreter','latex');
xlabel('Time-of-flight [ns]','Interpreter','latex');
title('Raw spectra','Interpreter','latex')
set(gca,'Fontname','Times','Fontsize',text_size)
set(gcf,'Position',posi)
ax = gca;
ax.GridAlpha = .25;
hold off
disp('Inspecting data...')
% Communication message
answ = 1;
while answ == 1
    shg                                                                    % Calling back the plot
    pause(time)                                                            % Delay of 'time' in seconds
    % Question for the user
    answer = questdlg('Do you need more time ?', ...
        'Question', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            answ = 1;
        case 'No'
            answ = 0;
    end
end
if sv
    disp('Saving plot...')
    disp('This might take few moments ...')
    cd(Storage_path)
    % Saving the plot in .fig format
    savefig('Quick_view_raw_data.fig')
    % Saving the plot in .png format
    name_save = 'Quick_view_raw_data.png';
    print(name_save,'-dpng','-r300');
    % Saving tof_mz_intensity_raw array in .mat format
    name_save = ('tof_mz_intensity_raw');
    save(name_save,'tof_mz_intensity_raw')
    disp('Saved')
    cd(Script_path)
    close()
else
    close()
end
cd(Script_path)
disp('Quick view of raw data done')
disp('..................................................................................')
toc
pause(1)
%% Alignment routine
% In this section raw mass spectra are aligned taken into account
% parameters chosen by the user at the start
tof_mz_intensity_raw_temp = tof_mz_intensity_raw;
tic
disp('..................................................................................')
disp('Alignment routine launched...')
status = 0;
count = 1;
while status == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Alignment parameters
    % In this section the user is asked to choose mass spectra alignment
    % routine parameters 
    disp('..................................................................................')
    disp('Setting for alignment routine...')
    % Choices for spectra alignement step operating mode: user-assisted or
    % automated
    answer = questdlg('Operating mode:', ...
        'Alignment routine', ...
        'User-assisted', 'Automated', 'User-assisted');
    % Handle response
    switch answer
        case 'User-assisted'
            Option_Alignment = 0;
            disp('User-assisted option selected')
        case 'Automated'
            Option_Alignment = 1;
            disp('Automated option selected')
    end
    if Option_Alignment == 0
        % Choice how peak location is done for peaks used for alignment
        % routine: Raw maximum or maximum after a gaussian fit
        answer = questdlg('In this case, how would you like to locate peaks ?', ...
            'Alignment routine', ...
            'Raw maximum', 'Maximum after a gaussian fit', 'Maximum after a gaussian fit');
        % Handle response
        switch answer
            case 'Raw maximum'
                option_al = 0;
                disp('Raw maximum option selected')
            case 'Maximum after a gaussian fit'
                option_al = 1;
                disp('Maximum after a gaussian fit option selected')
        end
        % Choice for the number of peaks for alignment routine
        n_peak_align = 100;                                                % Initialization 
        while n_peak_align > 31  || n_peak_align < 5
             n_peak_align = inputdlg({'Please select the number of peaks (between 5 and 30) for alignment routine '},...
                                      'Number of peaks', [1 50]);
             n_peak_align = str2num(cell2mat(n_peak_align));
        end
        disp(['Alignment routine will consider ',num2str(n_peak_align),' peaks'])
        
        % Saving alignment parameters
        cd(Storage_path)
        disp('Saving alignment routine parameters...')
        name_save = ('Alignment_routine_parameters');
        save(name_save,'Option_Alignment','option_al','n_peak_align')
        disp('Saved')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd(Script_path)
        [tof_mz_intensity_aligned, Alignment_coef, values] = Alignment_routine(tof_mz_intensity_raw_temp, n_spec, n_peak_align, option_al, Option_Alignment, text_size, posi, time, jetcustom, pre_cal_params, sv_plus, Script_path, Storage_path, Binning);
        % Checking the quality of the alignment routine
        % In this subsection raw aligned data is plotted to check the quality of the alignment routine 
        disp('..................................................................................')
        disp('Checking the quality of alignment routine...')
        cd(Storage_path)
        figure()
        for i = 1 : n_spec
            plot(tof_mz_intensity_aligned{i}(:,1)/10, tof_mz_intensity_aligned{i}(:,3), 'Color', jetcustom(i,:), 'Linewidth', 1.5) % Dividing by 10 since 10 points/ ns
            hold on
        end
        grid on
        ylabel('Ion count','Interpreter','latex');
        xlabel('Time of flight [ns]','Interpreter','latex');
        title('Aligned spectra','Interpreter','latex')
        set(gca,'Fontname','Times','Fontsize',text_size)
        set(gcf,'Position',posi)
        ax = gca;
        ax.GridAlpha = .25;
        hold off
        % Communication message
        answ = 1;
        while answ == 1
            shg                                                            % Going back to the plot
            pause(time)                                                    % Delay of 'time' in seconds
            % Question for the user
            answer = questdlg('Do you need more time ?', ...
                'Question', ...
                'Yes','No','No');
            % Handle response
            switch answer
                case 'Yes'
                    answ = 1;
                case 'No'
                    answ = 0;
            end
        end
        % Choice for spectra calibration step operating mode: 1st Question
        answer = questdlg('Are you satisfied with the aligment result', ...
            'Alignment routine: evaluation', ...
            'Satisfied', 'Not satisfied: Re-aligne', 'Not satisfied: Re-aligne');
        % Handle response
        switch answer
            case 'Satisfied'
                status = 1;
                disp('Alignment approved')
                if sv == 1
                    cd(Storage_path)
                    disp('Saving files and plots...')
                    % Saving the plot in .fig format
                    savefig('Quick_view_aligned_data.fig')
                    % Saving the plot in .png format
                    name_save = 'Quick_view_aligned_data.png';
                    print(name_save,'-dpng','-r300');

                    % saving tof_mz_intensity_aligned array in .mat format
                    name_save = ('tof_mz_intensity_aligned');
                    save(name_save,'tof_mz_intensity_aligned')
                    % Saving linear regression values applyed to each raw mass spectra 
                    name_save = ('Alignment_coef');
                    save(name_save,'Alignment_coef', 'values')
                    if sv_plus
                        % Saving tof_mz_intensity_aligned array in .xlsx format
                        disp('This might take few moments ...')
                        for i = 1 : n_spec  
                            disp(['Saving aligned spectrum number: ', num2str(i)])
                            xlswrite(['tof_mz_intensity_aligned_' Files_list(i).name(1:end-4) '.xlsx'], cat(1, {'tof','mz','intensity'}, num2cell(tof_mz_intensity_aligned{i})));
                        end
                    end
                    disp('saved')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure()
                ax1 = subplot(2,1,1);
                plot([1:1:n_spec], Alignment_coef(:,1),'wo', 'Markersize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
                grid on
                ylabel('a','Interpreter','latex');
                set(gca,'Fontname','Times','Fontsize',text_size)
                set(gcf,'Position',posi)
                ax2 = subplot(2,1,2);
                plot([1:1:n_spec], Alignment_coef(:,2),'wo', 'Markersize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
                grid on
                ylabel('b','Interpreter','latex');
                xlabel('Number of mass spectra','Interpreter','latex');
                set(gca,'Fontname','Times','Fontsize',text_size)
                set(gcf,'Position',posi)
                ax = gca;
                ax.GridAlpha = .25;
                hold off

                if sv
                    cd(Storage_path)
                    % Saving the peak detection plot in .fig format
                    savefig(['Alignment_coef_', num2str(count),'.fig'])  
                    % Saving the peak detection plot in .png format
                    name_save = ['Alignment_coef_', num2str(count),'.png'];
                    print(name_save,'-dpng','-r300');
                    cd(Script_path)
                    close()
                else
                    close()
                end

            case 'Not satisfied: Re-aligne'
                status = 0;
                tof_mz_intensity_raw_temp = tof_mz_intensity_aligned;
                disp('Relaunching alignment routine...')
                if count == 1
                   Alignment_coef_start = Alignment_coef;
                   values_start = values;
                   
                   % Saving linear regression values applyed to each raw
                   % mass spectra for the first attempt
                    name_save = ('Alignment_coef_first_attempt');
                    save(name_save,'Alignment_coef_start', 'values_start')
                   
                end
                count = count + 1;
        end
        close()
    elseif Option_Alignment == 1                                           % If automated mode 
        disp('Automated mode not available for the moment. Please switch to user-assisted mode ...')
        status = 0; 
    end

end

clear status
cd(Script_path)
disp('Alignment routine done')
disp('..................................................................................')
toc
pause(1)

%% Peaks picking/locating routine
% In this section peaks are located in each spectrum
tic
disp('..................................................................................')
disp('Peak locating routine launched...')
cd(Script_path)
if Option_Alignment == 0 % In this case data have a given structure
    % Empty array
    Peak_List_each_spectra = cell(n_spec,1);
    Soulders_count = cell(n_spec,1);
    VariableNames = {'ID','Accurate_mass','Intensity(TIC)','Width','Area','Std','Gaussp1','tof','Nominal_mass','Mass_defect','log_Intensity_pow_2'};

   for i = 1 : n_spec

         [Peak_List_each_spectra{i}, Soulders_count{i}] = Peak_locating_routine(tof_mz_intensity_aligned{i}, param,i, text_size, posi,x0);         
         if param.doPlot == 1 && sv == 1
             disp('Saving data...')
             cd(Storage_path)
             % Saving the peak detection plot in .fig format 
             savefig(['Peak_locating_' Files_list(i).name(1:end-4) '.fig'])  
             % Saving the peak detection plot in .png format
             name_save = ['Peak_locating_' Files_list(i).name(1:end-4) '.png'];
             print(name_save,'-dpng','-r300');
             % Saving Peak_List_each_spectra array in .xlsx format after being
             % aligned and pre-calibrated
             if sv_plus
                xlswrite(['Peak_List_each_spectra_precalibrated_' Files_list(i).name(1:end-4) '.xlsx'], cat(1, VariableNames, num2cell(Peak_List_each_spectra{i})));
             end
             close()
             disp('Saved')
             cd(Script_path)
         end    
         close()
         clc
    end

    if sv 
        disp('Saving data (Pre-calibrated)...')
        cd(Storage_path)
        % Saving Peak_List_each_spectra in .mat format
        name_save = ('Peak_List_each_spectra_precalibrated');
        save(name_save,'Peak_List_each_spectra', 'Soulders_count')
        disp('Saved')
        cd(Script_path)
    end
    cd(Script_path)
    
else
     Peak_List_each_spectra = cell(n_spec,1);
     for i = 1 : 1
        disp(['Mass spectrum number : ', num2str(i)])
        Peak_List_each_spectra_temp = cell(size(tof_mz_intensity_aligned,2),1);
        parfor k  = 2 : size(tof_mz_intensity_aligned,2) - 1
            local_param = param;
            local_param.doPlot = 0;
            Local_spectrum = tof_mz_intensity_aligned{i,k};
            if isempty(Local_spectrum(:,3)) == 0 && (k <= local_param.Range(2))
                disp(['Window number : ', num2str(k)])
                [Peak_List_each_spectra_temp{k}] = Peak_locating_routine(Local_spectrum, local_param, i, text_size, posi,x0);    
            end
        end
        for k = 2 : size(tof_mz_intensity_aligned,2)
            if isempty(Peak_List_each_spectra_temp{k}) == 0
                Peak_List_each_spectra{i} = [Peak_List_each_spectra{i}; Peak_List_each_spectra_temp{k}];
            end
        end
     end
end

Number_of_detected_peaks = zeros(n_spec,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()

for i = 1 : n_spec
    Number_of_detected_peaks(i) = size(Peak_List_each_spectra{i},1);
    plot(i, size(Peak_List_each_spectra{i},1),'wo', 'Markersize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
    hold on
end
grid on

ylabel('Number of detected peaks','Interpreter','latex');
xlabel('Number of mass spectra','Interpreter','latex');
set(gca,'Fontname','Times','Fontsize',text_size)
set(gcf,'Position',posi)
ax = gca;
ax.GridAlpha = .25;
hold off

if sv
    cd(Storage_path)
    % Saving the peak detection plot in .fig format 
    savefig('Number_of_detected_peaks.fig')  
    % Saving the peak detection plot in .png format
    name_save = 'Number_of_detected_peaks.png';
    print(name_save,'-dpng','-r300');
    cd(Script_path)
    close()
else
    close()
end

disp('Peak picking/locating routine done')
disp('..................................................................................')
toc
pause(1)
%% Calibration routine
% In this section mass spectra are calibrated
tic
disp('..................................................................................')
disp('Calibration routine launched...')
status = 0;
count = 1;
while status == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration parameters
    % In this section the user is asked to choose mass spectra calibration
    % routine parameters 
    disp('..................................................................................')
    disp('Setting for calibration routine...')
    % Choice for spectra calibration step operating mode: 1st Question
    answer = questdlg('For the calibration process, would you like to ', ...
        'Calibration routine: Operating mode', ...
        'Use available constants', 'Revaluate constants', 'Revaluate constants');
    % Handle response
    switch answer
        case 'Use available constants'
            force = 0;
            disp('Use available constants option selected')
        case 'Revaluate constants'
            force = 1;
            disp('Revaluate constants option selected')
    end
    if force == 1                                                          % In case the user wants to revaluate the calibration constants
        % Choice for spectra calibration step operating mode: 2nd Question
        answer = questdlg('Please make a choice for mass spectra calibration operating mode', ...
            'Calibration routine: Operating mode', ...
            'User-assisted', 'Automated', 'Automated');
        % Handle response
        switch answer
            case 'User-assisted'
                Option_Calibration = 0;
                disp('User-assisted option selected')
            case 'Automated'
                Option_Calibration = 1;
                disp('Automated option selected')
        end
        
        if Option_Calibration == 1
          disp('Automated mode not available for the moment. Please switch to user-assisted mode...')
          status = 0;  
        else
            % Choice how peak location is considered for peaks used for alignment
            % routine
            answer = questdlg('For the calibration routine, would you like to consided a set of calibration constants for:', ...
                'Calibration routine: Operating mode', ...
                'All mass spectra', 'Each mass spectrum', 'Each mass spectrum');
            % Handle response
            switch answer
                case 'All mass spectra'
                    op_calib = 0;
                    disp('All mass spectra option selected')
                case 'Each mass spectrum'
                    op_calib = 1;
                    disp('Each mass spectrum option selected')
            end
                % Choice for the number of points selected on mass defect plots for 
                % calibration routine
                n_calib = 100;         % Initialization
                while n_calib >= 21
                     n_calib = inputdlg({'Please select a number below 20 !'},...
                              'Number of points used for mass spectra calibration', [1 50]);
                     n_calib = str2num(cell2mat(n_calib));
                end
                disp(['Calibration routine will consider ', num2str(n_calib), ' points clusters']) 
               
                % Calibration law: to convert the time of flight (tof) scale into a m/z scale of
                % the form: m/z = A + B*tof + C*tof^2
                Cal_Law = @(x,y) x(1)*y + x(2)*y.^2;
                % Saving calibration parameters
                cd(Storage_path)
                disp('Saving calibration routine parameters...')
                name_save = ('Calibration_routine_parameters');
                save(name_save,'Option_Calibration','op_calib', 'n_calib','Cal_Law','x0', 'force')
                disp('Saved')

                cd(Script_path)
                [tof_mz_intensity_calibrated, param_calib, Peak_List_each_spectra, tof_cal, init_cal_mz_list, mz_prev_each] = Calibration_routine(tof_mz_intensity_aligned, Peak_List_each_spectra, n_calib, Cal_Law, n_spec, x0, text_size, posi, time, op_calib, Option_Calibration, force);
                close()
                
                % Checking the quality of the calibration routine
                % In this subsection raw aligned data is plotted to check the quality of the alignment routine 
                disp('..................................................................................')
                disp('Checking the quality of calibration routine...')            
                cd(Storage_path)

                figure();
                for i = 1 : n_spec
                    mat = Peak_List_each_spectra{i};
                    L = find(mat(:,11) > 0);
                    scatter(mat(L,9),mat(L,10),mat(L,11),'filled','ok')
                    hold on
                    clear mat L
                end
                xlabel('m','Interpreter','latex');
                ylabel('$\Delta$m','Interpreter','latex');
                set(gca,'Fontname','Times','Fontsize',text_size)
                set(gcf,'Position',posi)
                ax = gca;
                ax.GridAlpha = .25;
                grid on
                hold off

                if sv == 1 
                    cd(Storage_path)
                    disp('Saving files and plots...')
                    % Saving the plot in .fig format
                    savefig('Mass_defect_plot_calibrated.fig')
                    % Saving the plot in .png format
                    name_save = 'Mass_defect_plot_calibrated.png';
                    print(name_save,'-dpng','-r300');

                    cd(Script_path)
                end
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fig = figure();
                plot(tof_mz_intensity_aligned{1}(1:2000:end,1)/10, tof_mz_intensity_aligned{1}(1:2000:end,2), 'bs', 'Markersize', 6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b')
                hold on
                plot(tof_cal/10, init_cal_mz_list,'hk', 'Markersize', 10, 'MarkerFaceColor', 'y')
                hold on
                plot(linspace(min(tof_cal), max(tof_cal),100)/10, Cal_Law(param_calib(2:end), linspace(min(tof_cal), max(tof_cal),100)),'-r', 'LineWidth', 2)
                hold on
                plot(linspace(max(tof_cal), 11.5*1E5,100)/10, Cal_Law(param_calib(2:end), linspace(max(tof_cal), 11.5*1E5,100)),'--r', 'LineWidth', 1.5)
                grid on
                xlim([0, 11.5*1E4])
                ylim([param.Range(1), param.Range(2)])
                pause(1)
                legend('Aligned raw data',...
                       'Points used for the calibration routine',...
                       ['Rendered calibration law ' newline '$m$/$z$ = 0.0 + ', num2str(10*param_calib(2),'%.4e') , ' $\times$ tof + ', num2str(100*param_calib(3),'%.4e') , ' $\times$ tof$^2$'],...
                       'Extended calibration law',...
                       'Location','NorthWest','Interpreter','latex','numColumns', 1)

                ylabel('$m$/$z$','Interpreter','latex');
                xlabel('Time of flight [ns]','Interpreter','latex');
                title('Calibration curve','Interpreter','latex')
                set(gca,'Fontname','Times','Fontsize',text_size)
                set(gcf,'Position',posi)
                ax = gca;
                ax.GridAlpha = .25;
                ylim([param.Range(1), param.Range(2)])

                hold off
                if sv 
                    cd(Storage_path)
                    disp('Saving files and plots...')
                    % Saving the plot in .fig format
                    savefig(['Calibration_curve_', num2str(count),'.fig'])
                    % Saving the plot in .png format
                    name_save = ['Calibration_curve_', num2str(count),'.png'];
                    print(name_save,'-dpng','-r300');
                    close()
                    cd(Script_path)
                end
                % Communication message
                answ = 1;
                while answ == 1
                    shg                                                    % Going back to the plot
                    pause(time)                                            % Delay of 'time' in seconds
                    % Question for the user
                    answer = questdlg('Do you need more time ?', ...
                        'Question', ...
                        'Yes','No','No');
                    % Handle response
                    switch answer
                        case 'Yes'
                            answ = 1;
                        case 'No'
                            answ = 0;
                    end
                end
                % Choice for spectra calibration step operating mode: 1st Question
                x0 = param_calib';
                close()
                answer = questdlg('Are you satisfied with the calibration result', ...
                    'Calibration routine: evaluation', ...
                    'Satisfied', 'Not satisfied: Re-calibrate', 'Not satisfied: Re-calibrate');
                % Handle response
                switch answer
                    case 'Satisfied'
                        status = 1;
                        disp('Calibration approved')
                        if sv == 1 
                            cd(Storage_path)
                            disp('Saving files ...')
                            % saving tof_mz_intensity_calibrated array in .mat format
                            name_save = ('tof_mz_intensity_calibrated');
                            save(name_save,'tof_mz_intensity_calibrated')
                            % Saving Peak_List_each_spectra in .mat format
                            name_save = ('Peak_List_each_spectra_Calibrated');
                            save(name_save,'Peak_List_each_spectra')
                            % Saving calibration parameters 
                            name_save = ('param_calib');
                            save(name_save,'param_calib', 'tof_cal', 'init_cal_mz_list', 'mz_prev_each')
                            % Saving tof_mz_intensity_calibrated array in .xlsx format
                            if sv_plus 
                                disp('This might take few moments...')
                                for i = 1 : n_spec  
                                    disp(['Saving calibrated spectrum number: ', num2str(i)])
                                    xlswrite(['tof_mz_intensity_calibrated_' Files_list(i).name(1:end-4) '.xlsx'], cat(1, {'tof','mz','intensity'}, num2cell(tof_mz_intensity_calibrated{i})));
                                end
                            end
                            disp('Saved')
                        end
                    case 'Not satisfied: Re-calibrate'
                        status = 0;
                        disp('Restarting calibration routine...')
                        close all
                        count = count + 1;
                end
               
        end
    else
        Option_Calibration = 2;                                            
        n_calib = 10;                                                      
        % In this case the 'All mass spectra' option is taken by default
        op_calib = 0;
        disp('The following calibration constants will be used : ')
        disp(['A = ', num2str(x0(1))])
        disp(['B = ', num2str(x0(2))])
        disp(['C = ', num2str(x0(3))])
   
        % Calibration law: to convert the time of flight(tof) scale into a m/z scale of
        % the form: m/z = A + B*tof + C*tof^2
        Cal_Law = @(x,y) x(1)*y + x(2)*y.^2;
        % Saving calibration parameters
        cd(Storage_path)
        disp('Saving calibration routine parameters...')
        name_save = ('Calibration_routine_parameters');
        save(name_save,'Option_Calibration','op_calib', 'n_calib','Cal_Law','x0', 'force')
        disp('Saved')

        cd(Script_path)
        [tof_mz_intensity_calibrated, param_calib, Peak_List_each_spectra, tof_cal, init_cal_mz_list] = Calibration_routine(tof_mz_intensity_aligned, Peak_List_each_spectra, n_calib, Cal_Law, n_spec, x0, text_size, posi, time, op_calib, Option_Calibration, force);
        
        % Checking the quality of the calibration routine
        % In this subsection raw aligned data is plotted to check the quality of the alignment routine 
        disp('..................................................................................')
        disp('Checking the quality of calibration routine...')            
        cd(Storage_path)

        figure();
        for i = 1 : n_spec
            mat = Peak_List_each_spectra{i};
            L = find(mat(:,11) > 0);
            scatter(mat(L,9),mat(L,10),mat(L,11),'filled','ok')
            hold on
            clear mat L
        end
        xlabel('m','Interpreter','latex');
        ylabel('$\Delta$m','Interpreter','latex');
        set(gca,'Fontname','Times','Fontsize',text_size)
        set(gcf,'Position',posi)
        ax = gca;
        ax.GridAlpha = .25;
        grid on
        hold off

        if sv == 1 
            cd(Storage_path)
            disp('Saving files and plots...')
            % Saving the plot in .fig format
            savefig('Mass_defect_plot_calibrated.fig')
            % Saving the plot in .png format
            name_save = 'Mass_defect_plot_calibrated.png';
            print(name_save,'-dpng','-r300');
            cd(Script_path)
        end

        % Communication message
        answ = 1;
        while answ == 1
            shg % Going back to the plot
            pause(time) % Delai of 'time' seconds
            % Question for the user
            answer = questdlg('Do you need more time ?', ...
                'Question', ...
                'Yes','No','No');
            % Handle response
            switch answer
                case 'Yes'
                    answ = 1;
                case 'No'
                    answ = 0;
            end
        end
        % Choice for spectra calibration step operating mode: 1st Question
        x0 = param_calib';
        close()
        answer = questdlg('Are you satisfied with the calibration result', ...
            'Calibration routine: evaluation', ...
            'Satisfied', 'Not satisfied: Re-calibrate', 'Not satisfied: Re-calibrate');
        % Handle response
        switch answer
            case 'Satisfied'
                status = 1;
                disp('Calibration approved')
                if sv == 1 
                    cd(Storage_path)
                    disp('Saving files ...')
                    % saving tof_mz_intensity_calibrated array in .mat format
                    name_save = ('tof_mz_intensity_calibrated');
                    save(name_save,'tof_mz_intensity_calibrated')
                    % Saving Peak_List_each_spectra in .mat format
                    name_save = ('Peak_List_each_spectra_Calibrated');
                    save(name_save,'Peak_List_each_spectra')
                    % Saving calibration parameters 
                    name_save = ('param_calib');
                    save(name_save,'param_calib', 'tof_cal', 'init_cal_mz_list')
                    % Saving tof_mz_intensity_calibrated array in .xlsx format
                    if sv_plus 
                        disp('This might take few moments...')
                        for i = 1 : n_spec  
                            disp(['Saving calibrated spectrum number: ', num2str(i)])
                            xlswrite(['tof_mz_intensity_calibrated_' Files_list(i).name(1:end-4) '.xlsx'], cat(1, {'tof','mz','intensity'}, num2cell(tof_mz_intensity_calibrated{i})));
                        end
                    end
                    disp('Saved')
                end
            case 'Not satisfied: Re-calibrate'
                status = 0;
                disp('Restarting calibration routine...')
                close all
                count = count + 1;
        end
    end
end
clear status
cd(Script_path)
disp('Calibration routine done')
disp('..................................................................................')
toc
pause(1)
%% Quick view of aligned and calibrated data
% In this section raw data is plotted to have a quick view 
tic
disp('..................................................................................')
disp('Quick view of aligned and calibrated data')
cd(Storage_path)
figure()
clf
for i = 1 : n_spec
    plot(tof_mz_intensity_calibrated{i}(:,2), tof_mz_intensity_calibrated{i}(:,3), 'Color', jetcustom(i,:), 'Linewidth', 1.5);    
    hold all
end
grid on
ylabel('Ion count','Interpreter','latex');
xlabel('$m$/$z$','Interpreter','latex');
title('Aligned $\&$ calibrated spectra','Interpreter','latex')
set(gca,'Fontname','Times','Fontsize',text_size)
set(gcf,'Position',posi)
ax = gca;
ax.GridAlpha = .25;
hold off
disp('Inspecting data...')
% Communication message
answ = 1;
while answ == 1
    shg                                                                    % Calling back the plot
    pause(time)                                                            % Delay of 'time' in seconds
    % Question for the user
    answer = questdlg('Do you need more time ?', ...
        'Question', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            answ = 1;
        case 'No'
            answ = 0;
    end
end
if sv
    disp('Saving plot...')
    disp('This might take few moments ...')
    cd(Storage_path)
    % Saving the plot in .fig format
    savefig('Quick_view_aligned_calibrated_data.fig')
    % Saving the plot in .png format
    name_save = 'Quick_view_aligned_calibrated_data.png';
    print(name_save,'-dpng','-r300');
    disp('Saved')
    cd(Script_path)
    close()
else
    close()
end
cd(Script_path)
disp('Quick view of aligned and calibrated data done')
disp('..................................................................................')
toc
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating some statistics
tic
disp('Getting data statistics...')
disp('This might take few minutes...')
x_table = [];
y_table = [];
y1_table = [];
z_table = [];
for i = 1 : n_spec
    mat = Peak_List_each_spectra{i};
    for k = 1 : size(mat,1)
        if isnan(mat(k,5)) == 0 
            [~,idx] = min(abs(tof_mz_intensity_calibrated{i}(:,2) - mat(k,2)));
            local_x = tof_mz_intensity_calibrated{i}(idx - 20:idx + 20,2);
            
            A = mat(k,7);
            mu = mat(k,2);
            sigma = mat(k,6);
            
            fwhm = 2*sqrt(2*log(2))*sigma;
            x_table = [x_table; mat(k,3)]; % Intensity
            y_table = [y_table; fwhm];     % fwhm
            y1_table = [y1_table; A];      % area
            z_table = [z_table; mu]; % m/z

            clear fwhm sigma mu A local_x
        end
    end
    clear mat
end
if sv
    % Saving data tables 
    cd(Storage_path)
    disp('Saving data statistics')
    name_save = ('Peaks_statistics');
    save(name_save,'x_table', 'y_table', 'y1_table', 'z_table')
    disp('Saved')
    cd(Script_path)
end
toc
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Peak finding routine
% In this section clusters of points are identified as being related (not
% specified yet) using the peak finding algorithm
tic
disp('..................................................................................')
disp('Peak finding routine launched...')
cd(Script_path)
vars_Elements  = [{'m_z'};{'H'};{'B'}; {'C'}; {'[13C]'}; {'N'}; {'[15N]'}; {'O'}; {'F'}; {'Na'}; {'Al'}; {'Si'}; {'[29Si]'}; {'[30Si]'}; {'P'}; {'S'}; {'Cl'}; {'[37Cl]'}; {'K'}; {'[41K]'}; {'Ca'}; {'[42Ca]'}; {'[44Ca]'}; {'[46Ti]'}; {'[47Ti]'}; {'Ti'}; {'[49Ti]'}; {'[50Ti]'}; {'Cr'}; {'Fe'}; {'[58Fe]'}; {'Cu'}; {'[65Cu]'}; {'Cs'}];
vars_Elements = vars_Elements';

if sv 
    % Need also to save vars_Elements
    disp('Saving data...')
    cd(Storage_path)
    name_save = ('Detailed_Elements_first_row');
    save(name_save,'vars_Elements')
    disp('Saved')
    cd(Script_path)

end
[Peak_List, Sum_Peak_list, Proposed_Species, Proposed_Species_latex, Cluster_stat_gstds, Cluster_stat_gmns] = Peak_finding_routine(Peak_List_each_spectra, param, n_spec, text_size, posi, Storage_path ,Script_path, sv, max_iter_cluster_search, opt, sv_plus, reg, vars_Elements(2:end));

Sum_Peak_list = [Sum_Peak_list, zeros(size(Sum_Peak_list,1),1)];

if sv
    cd(Storage_path)
    disp('Saving files and plots...')
    xlim auto
    ylim auto
    % Saving the plot in .fig format
    savefig('Peak_finding.fig')
    % Saving the plot in .png format
    name_save = 'Peak_finding.png';
    print(name_save,'-dpng','-r300');
    close()
    % Saving the Peak_List cell in .mat format
    name_save = ('Peak_List');
    save(name_save,'Peak_List', 'Sum_Peak_list', 'Proposed_Species', 'Proposed_Species_latex', 'Cluster_stat_gstds', 'Cluster_stat_gmns')
    All_names = [''];
    All_names_sum = [];
    All_Data = [];
    for i = 1 : n_spec
        disp(['Saving data for spectrum number: ', num2str(i)])
        All_names = [All_names; {['Accurate mass S', num2str(i)]}; {['Nominal mass S', num2str(i)]}; {['Mass Defect S', num2str(i)]}; {['Time of Flight S', num2str(i)]}; {['Intensity S', num2str(i)]}];
        All_Data  = [All_Data, Peak_List{1,3}(:,i), Peak_List{1,6}(:,i), Peak_List{1,5}(:,i), Peak_List{1,4}(:,i), Peak_List{1,2}(:,i)];
    end
    All_names_sum = [All_names_sum; {'m_z'}; {'m'}; {'Delta_m'}; {'Delta_m_std'}; {'Disperssion_ppm'}; {'flag'}; {'Accuracy'}];
    All_names = All_names';
    All_names_sum = All_names_sum';
    % Saving the Peak_List in .xlsx format
    disp('This might take few moments ...')
    xlswrite('Peak_List_not_attributed.xlsx', cat(1, All_names, num2cell(All_Data)));
    xlswrite('Peak_List_sum_not_attributed.xlsx', cat(1, All_names_sum, num2cell(Sum_Peak_list)));
    disp('Saved')
    clear All_Data
    cd(Script_path)
else
    close()
end
disp('Peak finding routine done')
disp('..................................................................................')
toc

%% Plotting some statistical indicators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()

plot(Sum_Peak_list(find(Sum_Peak_list(:,6)==0),1), Sum_Peak_list(find(Sum_Peak_list(:,6)==0),1)./Sum_Peak_list(find(Sum_Peak_list(:,6)==0),4),'o', 'Markersize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
hold on
plot(Sum_Peak_list(find(Sum_Peak_list(:,6)==1),1), Sum_Peak_list(find(Sum_Peak_list(:,6)==1),1)./Sum_Peak_list(find(Sum_Peak_list(:,6)==1),4),'o', 'Markersize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
hold on
plot(Sum_Peak_list(find(Sum_Peak_list(:,6)==2),1), Sum_Peak_list(find(Sum_Peak_list(:,6)==2),1)./Sum_Peak_list(find(Sum_Peak_list(:,6)==2),4),'o', 'Markersize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')

ylabel('Resolution (m/$\Delta$m)','Interpreter', 'latex');
xlabel('$m$/$z$','Interpreter', 'latex');
set(gca,'Fontname','Times','Fontsize',text_size)
set(gcf,'Position',posi)
grid on
ax = gca;
ax.GridAlpha = .25;
set(gca,'Yscale','log')

ylim([1E2, 1E8])

if sv
    cd(Storage_path)
    disp('Saving plot...')
    savefig('Resolution.fig')
    % Saving the plot in .png format
    name_save = 'Resolution.png';
    print(name_save,'-dpng','-r300');
    disp('Saved')
    cd(Script_path)
    close()
else
    close()
end
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
scatter(z_table,y_table,1,'filled','ok')  
hold on
i = 1 ;
while i < 600
    if size(Cluster_stat_gstds{i},1) > 0
        for j = 1 : size(Cluster_stat_gstds{i},1)
           plot(i+Cluster_stat_gmns{i}(j), 2*Cluster_stat_gstds{i}(j),'o', 'Markersize', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r')
           hold on
        end
    end
i = i+1;
end
ylabel('$\Delta$m','Interpreter', 'latex');
xlabel('$m$/$z$','Interpreter', 'latex');
set(gca,'Fontname','Times','Fontsize',text_size)
set(gcf,'Position',posi)
grid on
ax = gca;
ax.GridAlpha = .25;
ylim([0, 0.045])
legend('fwhm', '2 $\times \sigma$', ...
        'Interpreter', 'latex', 'Location', 'NorthEast','FontSize',text_size - 4)


if sv
    cd(Storage_path)
    disp('Saving plot...')
    savefig('std_vs_fwhm.fig')
    % Saving the plot in .png format
    name_save = 'std_vs_fwhm.png';
    print(name_save,'-dpng','-r300');
    disp('Saved')
    cd(Script_path)
    close()
else
    close()
end
pause(1)

%% Species attribution routine (Main peak list)
% In this section clusters may be attributed. A user-interface is used
% to ease this process
disp('..................................................................................')
disp('Species attribution routine launched...')
cd(Storage_path)    
clear t
global t
t = readtable('Peak_List_sum_not_attributed.xlsx');

cd(Script_path) 

ID = reshape(reshape([1:1:size(t,1)],size(t,1),1),1,[])';
if isempty(Proposed_Species) == 0
    Species = reshape(Proposed_Species,1,[])';
else
    Species = reshape({''},1,[])';
end
Accuracy = reshape(repmat(0,size(t,1),1),1,[])';

t = addvars(t, Species,'Before','m_z');
t = addvars(t, ID,'Before','Species');
t = addvars(t, Accuracy,'After','flag');
vars = [{'ID'};{'Species'}; {'m_z'}; {'m'}; {'Delta_m'}; {'Delta_m_std'}; {'Disperssion_ppm'}; {'flag'}; {'Accuracy'}];
vars = vars';

t = t(1:end,vars);
Peak_attribution_routine(t, 0, Peak_List, text_size, Sum_Peak_list, vars_Elements)
disp('Species attribution routine done (may be temporary)')
disp('..................................................................................')

diary('off');
pause(1)
%% Saving species attribution process need to be done seperatly
disp('..................................................................................')
if sv_plus
    cd(Storage_path)
    disp('Saving attribution table...')
    xlswrite('Peak_List_sum_attributed_temp.xlsx', cat(1, vars, table2cell(t)))
    disp('Saved')
    cd(Script_path)
end
disp('..................................................................................')
pause(1)
%% If the user want to continue the attribution routine later on run this section 
% The script will get the Peak_List_attributed.xlsx that was uncompleted (For this 'sv_plus' option should have being set to 1 from the beginning)
clc
disp('..................................................................................')
disp('Reloading parameters and data...')
cd(Storage_path)
load('Paths.mat')
load('General_parameters.mat')
load('Data_characteristics.mat')
load('Peak_locating_parameters.mat')
load('Peak_calibration_initialization.mat')
load('Peak_finding_parameters.mat')
load('Files_parameters.mat')
load('Pre_conditioning_parameters.mat')
load('tof_mz_intensity_raw.mat')
load('Alignment_routine_parameters.mat')
load('tof_mz_intensity_aligned.mat')
load('Alignment_coef.mat')
load('Peak_List_each_spectra_precalibrated.mat')
load('Calibration_routine_parameters.mat')
load('tof_mz_intensity_calibrated.mat')
load('Peak_List_each_spectra_Calibrated.mat')
load('param_calib.mat')
load('Peak_finding_parameters.mat')
load('Peak_List.mat')
load('Peaks_statistics.mat')
load('Detailed_Elements_first_row')

disp('Parameters and data loaded')
cd(Script_path)
%%
disp('Former attibuted species table reloaded...')
cd(Storage_path)
clear t
global t
t = readtable('Peak_List_sum_attributed.xlsx');
cd(Script_path)

vars = [{'ID'};{'Species'}; {'m_z'}; {'m'}; {'Delta_m'}; {'Delta_m_std'}; {'Disperssion_ppm'}; {'flag'}; {'Accuracy'}];
vars = vars';

t = t(1:end,vars);
Peak_attribution_routine(t, 0, Peak_List, text_size, Sum_Peak_list, vars_Elements)
disp('Species attribution routine done (may be temporary)')
disp('..................................................................................')
pause(1)
%% Saving species attribution process need to be done seperatly
disp('..................................................................................')
if sv_plus
    disp('Saving...')
    cd(Storage_path)
    xlswrite('Peak_List_sum_attributed.xlsx', cat(1, vars, table2cell(t)))
    disp('Saved')
    cd(Script_path)
end
disp('..................................................................................')
pause(1)
%% Plotting mass defect plot after attributing species
% In this section a mass defect plot is displayed after species attribution
% process
disp('..................................................................................')
disp('Display of mass defect plot after species attribution...')
cd(Storage_path)
% t = readtable('Peak_List_attributed.xlsx');
t = readtable('Peak_List_sum_attributed.xlsx');

%% Exporting the peak list with each formula 

cd(Script_path)
[mz, namEl, numEl] = tof_exact_mass(t.Species);
mat_comb_rev = name_transform(mz, namEl, numEl);

% Exporting attributed species detailing each element
if sv_plus
    disp('Attributed species detailing each element...')
    cd(Storage_path)
    
    disp('Saving detailed attribution table...')
    mat_comb_rev(:,1) = t.m_z;
    xlswrite('Peak_List_sum_attributed_detailed.xlsx', cat(1, vars_Elements, num2cell(mat_comb_rev)))
    % Need also to save vars_Elements namE1 numE1
    name_save = ('Detailed_Elements');
    save(name_save, 'mat_comb_rev', 'namEl', 'numEl')
    disp('Saved')
    cd(Script_path)
end

disp('..................................................................................')
disp('APCFA post processing ended')
disp('..................................................................................')
