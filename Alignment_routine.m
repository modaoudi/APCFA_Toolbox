function [tof_mz_intensity_aligned, Alignment_coef, values] = Alignment_routine(tof_mz_intensity_raw, n_spec, n_peak_align, option_al, Option_Alignment, text_size, posi,time, jetcustom, pre_calib_param, sv, Script_path, Storage_path, Binning)

    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    % When mass spectra are pre-conditioned, first they undergo an 
    % alignment procedure. This phase is either user-assisted/automated
    % and consists of selecting a number of group of peaks in the mass 
    % spactra (in time of flight scale) raw data.
    % For each selected window, the script selects the time of flight of
    % each peak of the each mass spectrum that lays within the selected
    % window. All mass spectra are aligned against the first mass spectra
    % (considered arbitrary as reference) using a linear regression fit of
    % the form tof(raw) = a.tof(aligned) + b

    
    % INPUTS: 
    % tof_mz_intensity_raw  :   Raw data structured as follows :    % Time-of-Flight = 1st column 
                                                                    % m/z            = 2nd column (Pre-calibrated during the measurement phase)
                                                                    % Intensity      = 3rd column 
    % n_spec                :   Number of mass spectra
    % n_peak_align          :   Number of peaks used for the alignment procedure    
    % option_al             :   Option for peak location (Raw maximum/Maximum after a gaussian fit)
    % Option_Alignment      :   Operating mode option (User-assisted/Automated)
    % text_size             :   Text parameters
    % posi                  :   Plot position
    % time                  :   Waiting period
    % Plot_al               :   Plot option (under automated mode only, since the peak location procedure is used for peak detection)
    % param                 :   Peak locating parameters (under automated mode only)
    % jetcustom             :   Color setting
    % x0                    :   Calibration paramaters (default values used under automated mode only)  
    
    % OUTPUTS:
    % tof_mz_intensity_aligned  :   Aligned data structured as follows :   % Time of Flight = 1st column 
                                                                           % m/z            = 2nd column (Pre-calibrated during the measurement phase)
                                                                           % Intensity      = 3rd column 
    % Alignment_coef            :   'a' and 'b' coefficients used for spectra alignment  
    
    % Creating an empty cell
    tof_mz_intensity_aligned = cell(n_spec,1);
    if Option_Alignment == 0 % User-assisted case
        disp('User-assisted mode launched...')
        % The aligment of the the mass spectra is done using a linear regression of
        % the form y = a*x + b, applyed on data in the time of flight scale
        % considering the first mass spectra as a reference.
        Aligment_points = cell(n_spec,2,n_peak_align);
        values = zeros(n_spec,n_peak_align);
        % Aligment coefficient 'a' and 'b' for each spectrum
        Alignment_coef = zeros(n_spec,2);
        % Plot all mass spectra all at once
        fig = figure();
        for j = 1 : n_spec
            plot(tof_mz_intensity_raw{j}(:,1), tof_mz_intensity_raw{j}(:,3), 'Color', jetcustom(j,:), 'Linewidth', 1.5);
            hold on
        end
        hold off
        grid on
        ylabel('Ion count','Interpreter','latex');
        xlabel('Channel','Interpreter','latex');
        title('Raw spectra','Interpreter','latex')
        set(gca,'Fontname','Times','Fontsize',text_size)
        set(gcf,'Position',posi)
        ax = gca;
        ax.GridAlpha = .25;
        % For each alignement peak
        for i = 1: n_peak_align
            shg                                                            % Calling back the figure
            % Communication message
            answ = 1;
            while answ == 1
                wb = waitbar(0.67,['Please zoom in, peak number ', num2str(i) ,' out of ', num2str(n_peak_align),' : you have ', num2str(time),' secondes !!']); % Waiting bar
                pause(3)                                                                                                                                         % Closing waiting bar after 3 seconds
                close(wb)                                                                                                                                        % Closing waiting bar
                pause(time)                                                                                                                                      % Delai of 'time' seconds
                % Question for the user
                answer = questdlg('Have you zoomed in on a proper peak ?', ...
                    'Question', ...
                    'Yes','More Time', 'More Time');
                % Handle response
                switch answer
                    case 'Yes'
                        answ = 0;
                    case 'More Time'
                        answ = 1;
                end
            end
                % Waiting bar
                wb = waitbar(0.67,'Select a rectangle encompassing a pack of peaks');
                pause(3)                                                   % Closing waiting bar after 3 seconds
                close(wb)                                                  % Closing waiting bar
                % Recovering the rectangle coordinates
                selec_rect = getrect(fig.CurrentAxes);
                x1 = selec_rect(1);                                        % Time of flight lower value 
                x2 = selec_rect(1) + selec_rect(3);                        % Time of flight highest value
                % Communication message
                % Waiting bar
                wb = waitbar(0.9,['Peak between ', num2str(x1),' (ns) and ', num2str(x2),' (ns), selected']);  
                pause(3)                                                                                       % Closing waiting bar after 3 seconds
                close(wb)                                                                                      % Closing waiting bar
            hold off
            % Identifying time of flight and intensity of each peak within the selected rectangle 
            count = 0;
            for j = 1 : n_spec
                Aligment_points{j,1,i} = tof_mz_intensity_raw{j}(near_value(tof_mz_intensity_raw{j}(:,1),x1):near_value(tof_mz_intensity_raw{j}(:,1),x2),1); % All points between x1 and x2 (storing time of flight data)
                Aligment_points{j,2,i} = tof_mz_intensity_raw{j}(near_value(tof_mz_intensity_raw{j}(:,1),x1):near_value(tof_mz_intensity_raw{j}(:,1),x2),3); % All points between x1 and x2 (storing intensity data)
                if option_al == 0                                          % Raw maximum case
                    [~,ind] = max(Aligment_points{j,2,i});                 % Getting the index of the maximum intensity
                    count = count + 1;
                    values(j,i) = tof_mz_intensity_raw{j}(ind + near_value(tof_mz_intensity_raw{j}(:,1),x1),1);
                    clear ind
                elseif option_al == 1                                      % Maximum after a gaussian fit case
                    count = count + 1;
                    Tx_list = linspace(Aligment_points{j,1,i}(1), Aligment_points{j,1,i}(end), 100);
                    gi = griddedInterpolant(Aligment_points{j,1,i}, Aligment_points{j,2,i});   
                    Ty_list = gi(Tx_list);
                    [p1 p2 p3] = mygaussfit(Tx_list,Ty_list);              % Need a specific toolbox
                    [~,ind] = max(Aligment_points{j,2,i});                 % Getting the index of the maximum intensity

                    if (isnan(p1) == 1) || (isnan(p2) == 1) || (isnan(p3) == 1) || (isreal(p1) == 1) || (isreal(p2) == 1) || (isreal(p3) == 1) || (abs(p2 - tof_mz_intensity_raw{j}(ind + near_value(tof_mz_intensity_raw{j}(:,1),x1),1)) > Binning)
                        disp('Local switch to raw maximum option')
                        values(j,i) = tof_mz_intensity_raw{j}(ind + near_value(tof_mz_intensity_raw{j}(:,1),x1),1);
                    else
                        
                        if sv 

                            figure()
                            plot(Aligment_points{j,1,i}, Aligment_points{j,2,i}, 'sk')
                            hold on
                            shade(Tx_list, p1*exp(-((Tx_list - p2).^2)/(2*p3.^2)),...
                                          Tx_list, zeros(size(Tx_list,2),1),...
                                          'FillType',[1 2;2 1], 'FillColor', 'b','FillAlpha',.4);

                            ylabel('Ion count','Interpreter', 'latex');
                            xlabel('Time of flight [ns]','Interpreter', 'latex');
                            set(gca,'Fontname','Times','Fontsize',10)
                            set(gcf,'Position',posi)
                            ax = gca;
                            ax.GridAlpha = .25;
                            legend('Raw mass spectra', 'Gaussian fit', 'Location', 'best')
                            grid on

                            cd(Storage_path)
                            % Saving the plot in .fig format
                            savefig(['Peak_alignment_manual_fit_',num2str(i),'_', num2str(j),'.fig'])
                            % Saving the plot in .png format
                            name_save = ['Peak_alignment_manual_fit_',num2str(i),'_', num2str(j),'.png'];
                            print(name_save,'-dpng','-r300');
                            close()
                        end
                        values(j,i) = p2;
                        
                    end
                    clear p1 p2 p3 ind Tx_list Ty_list
                end
            end
            % Getting the number of detected peaks (indicator)
            disp(['Interrogation window number: ', num2str(i),' (out of ', num2str(n_peak_align),') '])
            disp(['Number of selected points: ', num2str(count),' '])
            clear ind f x1 x2
         
        end
        close(fig)
        
        disp('Regression ...')
        for i = 1 : n_spec
            % Deducing the alignement coefficient for each spectrum
            % Considering the first one as reference 
            Alignment_coef(i,:) = polyfit(values(i,:),values(1,:),1);
        end
        % Applying the alignement to raw data
        for i = 1 : n_spec
            tof_mz_intensity_aligned{i}(:,1) = polyval(Alignment_coef(i,:), tof_mz_intensity_raw{i}(:,1));                        % Time of flight
            tof_mz_intensity_aligned{i}(:,2) = pre_calib_param(i,1) + pre_calib_param(i,2)*tof_mz_intensity_aligned{i}(:,1) + pre_calib_param(i,3)*tof_mz_intensity_aligned{i}(:,1).^2;      % m/z
            tof_mz_intensity_aligned{i}(:,3) = tof_mz_intensity_raw{i}(:,3);                                                      % Intensity

        end
        disp('Mass spectra aligned')
    elseif Option_Alignment == 1                                            % Automated case
        
        disp('Automated mode launched...')

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Applying cross correlation with a local window
        half_Window_size = 120;                                            % Initializing half size of tof interrogation window               
        max_mz = round(tof_mz_intensity_raw{1}(end,2)) - 1;                % Getting the upper limite m/z range
        values = zeros(n_spec, max_mz);                                    % Initializing lag matrix

        store_window_center = [];                                          % Initialization
        x0 = pre_calib_param(1,:);                                         % Initialization
        
        disp('This might take few minutes...')
        for j  = 2 : 10 : max_mz                                                 % Loop over interrogation windows
                   
            Tof_window_center = round((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*j))/(2*x0(3)));    % Getting the center of tof interrogation window
            [~, Tof_window_center_pos] = min(abs(tof_mz_intensity_raw{1}(:,1) - Tof_window_center)); % Getting the index of center of tof interrogation window
            Tof_window_pos = [ Tof_window_center_pos - half_Window_size : 1 :  Tof_window_center_pos + half_Window_size];% Setting indices of tof interrogation window
            store_window_center = [store_window_center, Tof_window_center]; % Window indices
            
            if isequal(tof_mz_intensity_raw{1}(Tof_window_pos,3), zeros(size(tof_mz_intensity_raw{1}(Tof_window_pos,3)))) == 0
                
                [~,ind] = max(tof_mz_intensity_raw{1}(Tof_window_pos,3));  % Getting the index of the maximum intensity
                values(1,j) = tof_mz_intensity_raw{1}(ind + near_value(tof_mz_intensity_raw{1}(:,1),Tof_window_pos(1)),1);
                    
                for i = 2 : n_spec
                
                    [~,ind] = max(tof_mz_intensity_raw{i}(Tof_window_pos,3)); % Getting the index of the maximum intensity
                    values(i,j) = tof_mz_intensity_raw{i}(ind + near_value(tof_mz_intensity_raw{i}(:,1),Tof_window_pos(1)),1);
                    if abs(values(i,j) - values(1,j)) > 20*Binning
                        values(i,j) = NaN;
                    end
                
                end
                
            else
                values(:,j) = NaN;
            end
            
            clear Tof_window_center Tof_window_center_pos Tof_window_pos

        end 
        
        Alignment_coef = zeros(n_spec,2);
        L_ref    = find(values(1,:));
        L_temp   = find(isnan(values(1,:)) == 0);
        L_ref    = intersect(L_ref, L_temp);
        clear L_temp
        
        for i = 2: n_spec
            
            L_signal = find(values(i,:));
            L_temp   = find(isnan(values(i,:)) == 0);
            L_signal = intersect(L_signal, L_temp);
            clear L_temp
            
            L        = intersect(L_ref, L_signal);
            
            
            Alignment_coef(i,:) = polyfit(values(i,L), values(1,L),1);
            
            clear L_signal
        end
        clear L_ref
        
        % Applying the alignement to raw data
        for i = 1 : n_spec
            tof_mz_intensity_aligned{i}(:,1) = polyval(Alignment_coef(i,:), tof_mz_intensity_raw{i}(:,1));                         % Time of flight
            tof_mz_intensity_aligned{i}(:,2) = pre_calib_param(i,1) + pre_calib_param(i,2)*tof_mz_intensity_aligned{i}(:,1) + pre_calib_param(i,3)*tof_mz_intensity_aligned{i}(:,1).^2;      % m/z
            tof_mz_intensity_aligned{i}(:,3) = tof_mz_intensity_raw{i}(:,3);                                                       % Intensity

        end
        disp('Mass spectra aligned')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    end
    
end
