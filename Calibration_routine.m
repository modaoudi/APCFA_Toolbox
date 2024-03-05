function [tof_mz_intensity_calibrated, param_calib, Peak_List_each_spectra, tof_cal_each, init_cal_mz_list_each, mz_prev_each] = Calibration_routine(tof_mz_intensity_aligned, Peak_List_each_spectra, n_calib, Cal_Law, n_spec, x0, text_size, posi, time, op_calib, Option_Calibration, force)

    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    % Once all mass spectra are aligned, as the variable of interest is m/z
    % instead of time of flight (tof), all mass spactra are calibrated by
    % selecting a number of easyly identified chemical species.
    % The calibration consist of transforming the tof scale into m/z scale 
    % using the following law : m/z = A + B*tof +C*tof^2 (careful it is channele instead of tof)
    % To determine accurate calibration coefficients, A, B and C a mass 
    % defect plot is genetared considering a intial values of these 
    % calibration coefficients (Pre-defined at the begining of the script -
    % see Experimental setup section) this post processing step is 
    % user-assisted, to select a number of calibration points (Pre-defined 
    % at the begining of the script - see Experimental setup section)
    % During this process the function tof_exacte_mass, provided by Heikki 
    % Junninen PhD work, is used.

    % Before procceding with the next processing steps, all relevant
    % columns resulting from peak locating procedure are corrected with the
    % new calibration coefficients.

    if force == 0 % if the user wants to apply the calibration constants fed to the function
        disp('Applying calibration constants...')
        param_calib = x0;
        tof_cal_each = [];
        mz_prev_each = [];
        init_cal_mz_list_each = [];
        for i = 1 : n_spec
            % Params for all spectra (Applying calibration constants for all mass spectra is the only available option)
            Peak_List_each_spectra{i}(:,2) = param_calib(1) + param_calib(2)*Peak_List_each_spectra{i}(:,8) + param_calib(3)*Peak_List_each_spectra{i}(:,8).^2;
            tof_mz_intensity_aligned{i}(:,2) = param_calib(1) + param_calib(2)*tof_mz_intensity_aligned{i}(:,1) + param_calib(3)*tof_mz_intensity_aligned{i}(:,1).^2;
        end
        tof_mz_intensity_calibrated = tof_mz_intensity_aligned;
        disp('Mass spectra calibrated')
    else % if the user wants to revaluate the calibration constants
        if Option_Calibration == 0
            disp('Please wait for a few moments...')
            fig = figure();
            % Calculating CnHm pattern in a mass defect plot to help identifying 
            % some species fssssor the calibration process

            for i = 1 : 25                                                 %  For each C
                for j = 0 : 30                                             % For each H 
                    mz = tof_exact_mass(['C', num2str(i), 'H', num2str(j)]);
                    scatter(round(mz),mz-round(mz),'or')
                    hold on
                    if j == 0
                        if i == 1
                        text(round(mz)+.25,(mz-round(mz)),'C$^+$', 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                        else
                            text(round(mz)+.25,(mz-round(mz)),['C$_{', num2str(i),'}^+$'], 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                        end
                    elseif j == 1 && i == 1
                        text(round(mz)+.25,(mz-round(mz)),'CH$^+$', 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                    elseif i == 1 && j > 1 
                        text(round(mz)+.25,(mz-round(mz)),['CH$_{', num2str(j),'}^+$'], 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                    elseif i > 1 && j == 1
                        text(round(mz)+.25,(mz-round(mz)),['C$_{', num2str(i),'}$H$^+$'], 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                    else
                        text(round(mz)+.25,(mz-round(mz)),['C$_{', num2str(i),'}$H$_{', num2str(j),'}^+$'], 'Interpreter', 'latex', 'Fontsize',text_size,'Color','k');
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1 : n_spec
                mat = Peak_List_each_spectra{i};  
                L = find(mat(:,11) > 0);
                scatter(mat(L,9),mat(L,10),mat(L,11),'filled','ok')
                ylabel('$\Delta$m','Interpreter', 'latex');
                xlabel('m','Interpreter', 'latex');
                set(gca,'Fontname','Times','Fontsize',text_size)
                set(gcf,'Position',posi)
                grid on
                ay = gca;
                ay.GridAlpha = .25;
                hold on
                clear mat L
            end
            if op_calib == 0 
                init_cal_mz_list_each = [];
                tof_cal_each = [];
                mz_prev_each = [];
            elseif op_calib == 1
                init_cal_mz_list_each = cell(n_spec,1);
                tof_cal_each = cell(n_spec,1);
                mz_prev_each = cell(n_spec,1);
            end
            for i = 1 : n_calib
                shg
                answ = 1;
                while answ == 1
                    wb = waitbar(0.5,['Please, zoom in to select point number ', num2str(i) ,': you have ', num2str(time),' secondes !!']);
                    pause(2)
                    close(wb)
                    pause(time)
                    answer = questdlg('Have you zoomed in on a proper cluster of points ?', ...
                        'So !', ...
                        'Yes','More Time', 'More Time');
                    % Handle response
                    switch answer
                        case 'Yes'
                            answ = 0;
                        case 'More Time'
                            answ = 1;
                    end
                end
                wb = waitbar(0.7,'Select a rectangle encompassing the points of interest');
                pause(1)
                close(wb)
                shg
                selec_rect = getrect(fig.CurrentAxes);
                x = inputdlg({'Specie'},...
                          'Choice', [1 50]);
                disp(['Selected specie: ' x])
                if op_calib == 0 
                    temp_size = size(init_cal_mz_list_each,1);
                end

                for l = 1 : n_spec
                    L = find(selec_rect(1) < Peak_List_each_spectra{l}(:,9) & Peak_List_each_spectra{l}(:,9) < selec_rect(1) + selec_rect(3) & selec_rect(2) < Peak_List_each_spectra{l}(:,10) & Peak_List_each_spectra{l}(:,10) < selec_rect(2) + selec_rect(4));
                    if L~=0
                        if op_calib == 0
                            init_cal_mz_list_each = [init_cal_mz_list_each; tof_exact_mass(x)];
                            tof_cal_each = [tof_cal_each; Peak_List_each_spectra{l}(L(1),8)];
                            mz_prev_each = [mz_prev_each; Peak_List_each_spectra{l}(L(1),2)];
                        elseif op_calib == 1
                            init_cal_mz_list_each{l} = [init_cal_mz_list_each{l}; tof_exact_mass(x)];
                            tof_cal_each{l} = [tof_cal_each{l}; Peak_List_each_spectra{l}(L(1),8)];
                            mz_prev_each = [mz_prev_each; Peak_List_each_spectra{l}(L(1),2)];
                        end
                    end
                end
                disp(['Calibration point number: ', num2str(i),' (out of ', num2str(n_calib),') '])
                if op_calib == 0
                    disp(['Number of selected points: ', num2str(size(init_cal_mz_list_each,1) - temp_size), ' ']) 
                end
                clear x selec_rect temp_size answ
            end
            if op_calib == 0
                % Calculating calibration parameters
                temp = lsqcurvefit(Cal_Law, x0(2:3), tof_cal_each, init_cal_mz_list_each);
                param_calib = [0;reshape(temp,2,1)];
                clear temp
            elseif op_calib == 1
                param_calib = zeros(n_spec, 3);
                for i = 1 : n_spec
                      temp = flip(polyfit(tof_cal_each{i}, init_cal_mz_list_each{i},2));
                      param_calib(i,:) = [0, temp(2:end)];
                      clear temp
                end
            end
            % Data Update
            for i = 1 : n_spec
                if op_calib == 0
                    % Params for all spectra
                    Peak_List_each_spectra{i}(:,2) = param_calib(1) + param_calib(2)*Peak_List_each_spectra{i}(:,8) + param_calib(3)*Peak_List_each_spectra{i}(:,8).^2;
                    tof_mz_intensity_aligned{i}(:,2) = param_calib(1) + param_calib(2)*tof_mz_intensity_aligned{i}(:,1) + param_calib(3)*tof_mz_intensity_aligned{i}(:,1).^2;
                    Peak_List_each_spectra{i}(:,4) = param_calib(2)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,4)))/(2*x0(3))) + param_calib(3)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,4)))/(2*x0(3))).^2;
                    Peak_List_each_spectra{i}(:,5) = param_calib(2)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,5)))/(2*x0(3))) + param_calib(3)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,5)))/(2*x0(3))).^2;
                    Peak_List_each_spectra{i}(:,6) = param_calib(2)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,6)))/(2*x0(3))) + param_calib(3)*((-x0(2)+ sqrt(x0(2)^2+4*x0(3)*Peak_List_each_spectra{i}(:,6)))/(2*x0(3))).^2;

                elseif op_calib == 1 
                    % Params for each spectra
                    Peak_List_each_spectra{i}(:,2) = param_calib(i,1) + param_calib(i,2)*Peak_List_each_spectra{i}(:,8) + param_calib(i,3)*Peak_List_each_spectra{i}(:,8).^2;
                    tof_mz_intensity_aligned{i}(:,2) = param_calib(i,1) + param_calib(i,2)*tof_mz_intensity_aligned{i}(:,1) + param_calib(i,3)*tof_mz_intensity_aligned{i}(:,1).^2;
                end
                % Peak_List_each_spectra updated
                Peak_List_each_spectra{i}(:,9) = round(Peak_List_each_spectra{i}(:,2));
                Peak_List_each_spectra{i}(:,10) = Peak_List_each_spectra{i}(:,2) - Peak_List_each_spectra{i}(:,9);
            end
            tof_mz_intensity_calibrated = tof_mz_intensity_aligned;
            disp('Mass spectra calibrated')
        elseif Option_Calibration == 1
            disp('Automated mode not available for the moment. Please switch to user-assisted mode ...')
            param_calib = x0;
            tof_cal_each = [];
            mz_prev_each = [];
            init_cal_mz_list_each = [];
            
            for i = 1 : n_spec
                % Params for all spectra (Applying calibration constants for all mass spectra is the only available option)
                Peak_List_each_spectra{i}(:,2) = param_calib(1) + param_calib(2)*Peak_List_each_spectra{i}(:,8) + param_calib(3)*Peak_List_each_spectra{i}(:,8).^2;
                tof_mz_intensity_aligned{i}(:,2) = param_calib(1) + param_calib(2)*tof_mz_intensity_aligned{i}(:,1) + param_calib(3)*tof_mz_intensity_aligned{i}(:,1).^2;
            end
            tof_mz_intensity_calibrated = tof_mz_intensity_aligned;
            
        end
    end
 
end