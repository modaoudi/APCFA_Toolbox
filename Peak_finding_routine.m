function [Peak_List, Sum_Peak_list, Proposed_Species, Proposed_Species_latex, Cluster_stat_gstds, Cluster_stat_gmns] = Peak_finding_routine(Peak_List_each_spectra, param, n_spec, text_size, posi, Storage_path, Script_path, sv, max_iter_cluster_search, opt, sv_plus, reg, vars_Elements)
            
    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    % At this point of the post processing each peak list is independent in
    % which each row contains a certain information for a given peak, still
    % not attributed to a chemical specie. The main issue here is to be 
    % able to get a global Peak list for all mass spectra in witch 'i' th 
    % row, of one of the peak lists (i.e, one of the mass spectra) is 
    % associated to the 'j' th row of another peak list (i.e, another mass 
    % spectra).

    % Untill now, this was done manually using a mass defect plot. 
    % A mass defect plot allows to represent mass spectra data by plotting
    % located peaks intensities, for each a mass spectra, in a 2D 
    % coordinates. The first one (x-axis) is the nominal mass (i.e, 
    % round(mz)) and the second one (y-axis) is the mass defect (i.e., mz - 
    % round(mz)).

    % The difficulty here is that a given nominal mass more than one specie
    % exhibiting different mass defect for each, can be present in each 
    % mass spectra. The main issu is to be able to distinguish and 
    % associate, at a given nominal mass, peaks (i.e, rows in peak list 
    % table) that correspond to the same chemical specie. Knowing they do 
    % not have exactly the same mass defect and that more than one chemical
    % specie can be present at the same nominal mass. The other issu is 
    % that even with a well alignment and calibration procedure, for a 
    % given specie, present in all mass spactra for instance, it would 
    % still exhibit a mass defect dispersion and this mass defect 
    % dispersion is more obvious for higher nominal mass.
    % The idea here is to exploite mass defect graphic representation and 
    % determine for each nominal mass the number of cluseters that would 
    % probably describe the same chemical specie. Based on that clusters 
    % detection an interogation window is set to search for the position 
    % (ID) of the peak, if it exists, for each mass spactra (i.e, Peak list
    % of each spectra) and associated them. The interogation window heigth 
    % (mass defect dispersion) is dynamic as it varie considering on the 
    % one hand the dispersion of the point within each cluster and on the 
    % other hand the fact that the dispesrion is more significant for hight
    % nominal mass.

    % Having this feature permits a quick rows association (particularly 
    % beacause a dynamic interogation window is depolyed) and limits 
    % significatly the number of times peak attribution procedure (the next
    % step) the user will have to undergo. Additionnaly, is gives a quick 
    % acces to evaluat which mass spectra contain/or not a given chemical 
    % specie.

    % For the mass defect axis a dynamic range is set based on cluster 
    % point detection and standard deviation elongation along the nominal
    % mass axis as calibration was based on low mass range at this point
    
    % generate combination of elements and save to amu based
    
    Proposed_Species = [];
    Proposed_Species_latex = [];
    
    if reg
        % Need to call the fucntion here
        Elements_data_base_generation
    else
       
        % Loading 
        load('mat_comb.mat')
        
    end

    
    nm_L = [];
    nm_M = [];
    
    dx = 1;
    
    nm_min = 0.5;
    nm_max = param.Range(2) + 0.5;

    Np_nm = round((nm_max-nm_min)/dx);

    nm_grid = linspace(nm_min, nm_max, Np_nm+1);

    Temp_Peak_List = cell(1,6);
    Temp_Peak_List_eval = cell(1,6);
    reshaped_data = cell(param.Range(2),1);
    reshaped_data_1 = cell(param.Range(2),1);

    % Plotting the mass defect plot for all mass spectra
    figure();
    for k = 1 : n_spec
        mat = Peak_List_each_spectra{k};
        TR = setdiff(find(mat(:,11)), find(mat(:,11)<0));
        scatter(mat(TR,9),mat(TR,10),mat(TR,11),'filled','ok')
        hold on
    end
    ylabel('$\Delta$m','Interpreter', 'latex');
    xlabel('m','Interpreter', 'latex');
    set(gca,'Fontname','Times','Fontsize',text_size)
    set(gcf,'Position',posi)
    ay = gca;
    ay.GridAlpha = .25;
    
    grid on
    
    if sv_plus 
        cd(Storage_path)
        vidObj = VideoWriter('Peak_finding'); 
        vidObj.FrameRate = 10;
        open(vidObj);
        cd(Script_path)
    end
    
    for j = 1 : Np_nm
        reshaped_data{j} = zeros(4,1);
        reshaped_data{j} = reshaped_data{j}';
        reshaped_data_1{j} = zeros(2,1);
        reshaped_data_1{j} = reshaped_data_1{j}';
    end
    
    
    % Variable to save clusters data
    nom_mass            = zeros(Np_nm,1);
    Cluster_stat_nGrps  = zeros(Np_nm,1);
    Cluster_stat_gmns   = cell(Np_nm,1);
    Cluster_stat_gstds  = cell(Np_nm,1);
    Cluster_stat_gmin   = cell(Np_nm,1);
    Cluster_stat_gmax   = cell(Np_nm,1);
    
    thrsh_table = zeros(Np_nm, 1);
    
    for j = 1 : Np_nm
        
        for i = 1 : n_spec
            L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1));
            if L1~=0
                ID = Peak_List_each_spectra{i}(L1,1);                      % ID
                width = Peak_List_each_spectra{i}(L1,4);                   % gaussian fit width
                reshaped_data_1{j} = [reshaped_data_1{j}; [i*ones(size(ID,1),1), width]];
            end
            clear L1
        end

        if (nanmean(reshaped_data_1{j}(:,2)) ~= 0) && (isequal(reshaped_data_1{j}(:,2), real(reshaped_data_1{j}(:,2))))
            thrsh_table(j) = nanmean(reshaped_data_1{j}(:,2));   
        else
            thrsh_table(j) = NaN;
        end
        
    end
    
    nm_table = [1:1:Np_nm]';
    fit_param = polyfit(nm_table(~isnan(thrsh_table)), thrsh_table(~isnan(thrsh_table)),2);
    
    Threshold = polyval(fit_param, nm_table);
    
    figure()
    plot([1:1:Np_nm], thrsh_table, 'wo', 'Markersize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    hold on
    grid on
    plot([1:1:Np_nm], Threshold, 'Color', 'r', 'LineWidth', 3)
    
    ylabel('fwhm','Interpreter', 'latex');
    xlabel('$m$/$z$','Interpreter', 'latex');
    set(gca,'Fontname','Times','Fontsize',text_size)
    set(gcf,'Position',posi)
    ay = gca;
    ay.GridAlpha = .25;
    
    legend('Peaks fwhm - gaussian fit', 'Threshold law', 'Location', 'NorthWest', 'FontSize', text_size, 'Interpreter', 'latex')
    
    if sv 
        cd(Storage_path)
        % Saving the plot in .fig format
        savefig('Peak_finding_treshold_law.fig')
        % Saving the plot in .png format
        name_save = 'Peak_finding_treshold_law.png';
        print(name_save,'-dpng','-r300');
        cd(Script_path)
        close()
    else 
        close()
    end
    
    win_thrsh = flip(logspace(0,1,max_iter_cluster_search)/10);
    
    figure()
    plot([1:1:Np_nm], thrsh_table, 'wo', 'Markersize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    hold on
    grid on
    plot([1:1:Np_nm], Threshold, 'Color', 'r', 'LineWidth', 3)
    hold on
    
    for i = 1 : max_iter_cluster_search
        
        plot([1:1:Np_nm], Threshold*win_thrsh(i), 'Color', [17 17 17]/255, 'Linewidth', .5)
        hold on
        
    end
    
    ylabel('fwhm','Interpreter', 'latex');
    xlabel('$m$/$z$','Interpreter', 'latex');
    set(gca,'Fontname','Times','Fontsize',text_size)
    set(gcf,'Position',posi)
    ay = gca;
    ay.GridAlpha = .25;
    
    legend('Peaks fwhm - gaussian fit', 'Threshold law', 'Location', 'NorthWest', 'FontSize', text_size, 'Interpreter', 'latex')
    
    if sv 
        cd(Storage_path)
        % Saving the plot in .fig format
        savefig('Peak_finding_treshold_variation.fig')
        % Saving the plot in .png format
        name_save = 'Peak_finding_treshold_variation.png';
        print(name_save,'-dpng','-r300');
        cd(Script_path)
        close()
    else 
        close()
    end
    
    count = 0;
    flag_counter = [];
    for j = 1 : Np_nm 

        % For each nominal mass point clusters are detected to be able limit
        % the number of interogation windows along the mass defect axis
        % First data is reshaped
        for i = 1 : n_spec
            L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1));
            if L1~=0
                ID = Peak_List_each_spectra{i}(L1,1);  % ID
                md = Peak_List_each_spectra{i}(L1,10); % Mass defect
                int = Peak_List_each_spectra{i}(L1,3); % Intensity
                reshaped_data{j} = [reshaped_data{j}; [i*ones(size(ID,1),1), ID, md, int]];
            end
            clear ID L1 md int
        end

        thrsh = Threshold(j);
        [nGrps, gmns, gstds, gmin, gmax] = get_clusters(reshaped_data{j}(2:end,3),thrsh, 0);
        
        % continue only if a cluster of point are detected at a given nominal mass
        if nGrps > 0
            iter = nGrps;
            k = 1;
            while iter ~=0
                % storing current row in the Temp_Peak_List
                count = count + 1;
                for i = 1 : 6
                    % Filling with zeros for each interogation window
                    % 1 st column for ID
                    % 2nd column for Intensity
                    % 3rd column for m/z
                    % 4th column for time of flight
                    % 5th olumn for mass defect
                    % 6th column for nominal mass
                    Temp_Peak_List{1,i} = [Temp_Peak_List{1,i}; zeros(1,n_spec)];
                    Temp_Peak_List_eval{1,i} = [Temp_Peak_List_eval{1,i}; zeros(1,n_spec)];
                end
                
                overlap = 0;                                               % Initialization
                flag = 0 ;                                                 % Initialization
                
                L0 = cell(n_spec,1);
                for i = 1 : n_spec
                    % Here the interogation window is set
                    % L0{i} = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmns(k) - 3*gstds(k) < Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) < gmns(k) + 3*gstds(k));
                    L0{i} = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));
                    if size(L0{i},1) > 1
                        flag = flag + 1;
                    end
                end
                clear L0
                
                if flag > 0
                    overlap = 1;
                end 
                    
                overlap_flag_for_rect = 0;                                 % Initialization
                while overlap == 1
                    temp_flag = flag - 1 ;
                    temp_flag_iter = 1;
                    while (temp_flag > 0 ) && (temp_flag_iter < max_iter_cluster_search/10) % flag ~=0
                        old_nGrps = nGrps;
                        flag_2 = 1;                                        % Initialization
                        % If the overlap exists at the end
                        if k == nGrps
                            disp('Overlap exists at the end !')  
                            ii = 1;                                        % initalization
                            while flag_2 == 1 && ii < max_iter_cluster_search
                                thrsh_var = win_thrsh(ii)*thrsh;
                                % Here the data that should be processed should
                                % include only the remaining cluster of
                                % points
                                temporary_data = reshaped_data{j}(2:end,3);
                                L0_1 = find( gmin(k) <= temporary_data & temporary_data <= gmax(k));
                                [nGrps_add, gmns_add, gstds_add, gmin_add, gmax_add] = get_clusters(temporary_data(L0_1),thrsh_var, 1);
                                if nGrps_add > 1
                                    nGrps = nGrps + nGrps_add - 1;
                                    gmns  = [gmns(1:end-1);gmns_add]; 
                                    gstds = [gstds(1:end-1);gstds_add];
                                    gmin  = [gmin(1:end-1);gmin_add];
                                    gmax  = [gmax(1:end-1);gmax_add]; 
                                    flag_2 = 0;
                                    clear L0_1 temporary_data gmns_add gmax_add gmin_add gstds_add gmns_add thrsh_var
                                else
                                    flag_2 = 1; 
                                end
                                ii = ii + 1;
                            end
                            
                            if ii < max_iter_cluster_search
                                disp(['Number of clusters revaluating to ', num2str(nGrps),' after ', num2str(ii), ' attempts'])
                            end
                            
                            ii = 1;
                            iter = iter + (nGrps - old_nGrps);
                            clear old_nGrps

                        % If the overlap exists at the beginning
                        elseif k == 1 
                            disp('Overlap exists at the beginning !')                            
                            ii = 1;                                        % initalization
                            while flag_2 == 1 && ii < max_iter_cluster_search

                                thrsh_var = win_thrsh(ii)*thrsh;
                                % Here the data that should be processed should
                                % include only the remaining cluster of
                                % points
                                temporary_data = reshaped_data{j}(2:end,3);
                                L0_1 = find( gmin(k) <= temporary_data & temporary_data <= gmax(k));
                                [nGrps_add, gmns_add, gstds_add, gmin_add, gmax_add] = get_clusters(temporary_data(L0_1),thrsh_var, 1);
                                if nGrps_add > 1
                                    nGrps = nGrps + nGrps_add - 1;
                                    gmns  = [gmns_add;gmns(2:end)]; 
                                    gstds = [gstds_add;gstds(2:end)];
                                    gmin  = [gmin_add;gmin(2:end)];
                                    gmax  = [gmax_add;gmax(2:end)]; 
                                    flag_2 = 0;
                                    clear L0_1 temporary_data gmns_add gmax_add gmin_add gstds_add gmns_add thrsh_var
                                else
                                    flag_2 = 1;
                                end
                                ii = ii + 1;
                            end
                            
                            if ii < max_iter_cluster_search
                                disp(['Number of clusters revaluating to ', num2str(nGrps),' after ', num2str(ii), ' attempts'])
                            end
                            ii = 1;
                            iter = iter + (nGrps - old_nGrps);
                            clear old_nGrps
                           
                        % If the overlap exists at the middle
                        else
                            disp('Overlap exists at the middle !')                            
                            ii = 1;      % initalization
                            while flag_2 == 1 && ii < max_iter_cluster_search
                                thrsh_var = win_thrsh(ii+1)*thrsh;                                
                                % Here the data that should be processed and should
                                % include only the remaining cluster of
                                % points
                                temporary_data = reshaped_data{j}(2:end,3);
                                L0_1 = find( gmin(k) <= temporary_data & temporary_data <= gmax(k));
                                [nGrps_add, gmns_add, gstds_add, gmin_add, gmax_add] = get_clusters(temporary_data(L0_1),thrsh_var, 1);
                                if nGrps_add > 1
                                    nGrps = nGrps + nGrps_add - 1;
                                    gmns  = [gmns(1:k-1);gmns_add;gmns(k+1:end)]; 
                                    gstds = [gstds(1:k-1);gstds_add;gstds(k+1:end)];
                                    gmin  = [gmin(1:k-1);gmin_add;gmin(k+1:end)];
                                    gmax  = [gmax(1:k-1);gmax_add;gmax(k+1:end)]; 
                                    flag_2 = 0;
                                    clear L0_1 temporary_data gmns_add gmax_add gmin_add gstds_add gmns_add thrsh_var
                                else
                                    flag_2 = 1;
                                end
                                ii = ii + 1;
                            end
                            
                            if ii < max_iter_cluster_search
                                disp(['Number of clusters revaluating to ', num2str(nGrps),' after ', num2str(ii), ' attempts'])
                            end
                            ii = 1;
                            iter = iter + (nGrps - old_nGrps);
                            clear old_nGrps

                        end 
                        % Revaluate the flag
                        disp('Revaluating the overlap...')
                        flag = 0;
                        L0 = cell(n_spec,1);
                        for i = 1 : n_spec
                            % Here the interogation window is set
                            L0{i} = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));
                            if size(L0{i},1) > 1
                                flag = flag + 1;
                            end
                        end
                        temp_flag = flag - 1 ;
                        clear L0
                        temp_flag_iter = temp_flag_iter + 1 ;

                    end
                    if (temp_flag_iter >= max_iter_cluster_search/10)
                        disp('....................... Overlap still exists ........................')
                        % Should discard the cluster
                        overlap_flag_for_rect = 1;
                    else
                        disp('.......................   Overlap solved     ........................')
                    end
                    overlap = 0;
                    temp_flag_iter = 1;
                end
                % if overlap == 0
                % Plotting the interrogation windows positions
                Counter = 0; % Intialization
                for i = 1 : n_spec
                    L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));

                    if L1~=0
                        Counter = Counter + 1; 
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Before setting the interrogation window it needs to be
                % optimized in order to look for a candidate specie
                if (opt == 2)
                    % Switching to local variables
                    local_gmin  = gmin(k);
                    local_gmax  = gmax(k);
                    local_gstds = gstds(k);
                    local_gmns  = gmns(k);
                    
                    if (Counter >= 3)                        
                        % Geting cluster data 
                        for i = 1 : n_spec
                            % Here the interogation window is set
                            L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & local_gmin <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= local_gmax);

                            if L1~=0 
                                ID = Peak_List_each_spectra{i}(L1,1);
                                [~,m] = max(Peak_List_each_spectra{i}(ID,3)); % if more than one ID is detected consider the one with maximum intensity
                                Temp_Peak_List_eval{1,1}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),1); % ID
                                Temp_Peak_List_eval{1,2}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),3); % Intensity
                                Temp_Peak_List_eval{1,3}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),2); % m/z
                                Temp_Peak_List_eval{1,4}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),8); % tof
                                Temp_Peak_List_eval{1,5}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),10);% mass defect
                                Temp_Peak_List_eval{1,6}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),9); % nominal mass

                            end
                            clear ID L1
                        end 
                        if j == 62 
                           display('I ama here') 
                        end
                        figure();
                        H = histogram(Temp_Peak_List_eval{1,3}(count,find(Temp_Peak_List_eval{1,3}(count,:))),... 
                                               'Normalization','count',...
                                               'DisplayStyle','bar',...
                                               'BinMethod','fd', 'Facecolor','k');
                                           
                        ylabel('Count','Interpreter', 'latex');
                        xlabel('$m$/$z$','Interpreter', 'latex');
                        set(gca,'Fontname','Times','Fontsize',text_size)
                        set(gcf,'Position',posi)
                        grid on
                        ay = gca;
                        ay.GridAlpha = .25;
    
                        % Extracting histogram data
                        xx = H.BinEdges;
                        yy = [H.BinCounts,0];
                        if sv_plus
                            pause(2)
                            hold on
                            plot([j + local_gmin, j + local_gmin],[0, max(yy)], 'Color', 'k', 'LineWidth', 2)
                            hold on
                            plot([j + local_gmax, j + local_gmax],[0, max(yy)], 'Color', 'k', 'LineWidth', 2)
                        else
                            close()
                        end
                        Tx_list = linspace(xx(1), xx(end), 100);
                        gi = griddedInterpolant(xx, yy);   
                        Ty_list = gi(Tx_list);
                    
                        [o1, o2, o3] = mygaussfit(Tx_list, Ty_list);       % Perform a gaussian fit 
                        
                        if (isnan(o1) == 0) && (isnan(o2) == 0) && (isnan(o3) == 0) && (isreal(o1) == 1) && (isreal(o2) == 1) && (isreal(o3) == 1)

                            if  (max([o2 - round(o2) - 3*o3, local_gmin]) > min([o2 - round(o2) + 3*o3, local_gmax]))  ||  min([o2 - round(o2) + 3*o3, local_gmax]) > local_gmax || (max([o2 - round(o2) - 3*o3, local_gmin]) < local_gmin)

                                disp('Bad fitting')

                                local_gmin = local_gmns - 3*local_gstds;
                                local_gmax = local_gmns + 3*local_gstds;

                            else
                                if local_gmin >0
                                    local_gmin = 0.999*max([o2 - round(o2) - 3*o3, local_gmin]);
                                else
                                    local_gmin = 1.001*max([o2 - round(o2) - 3*o3, local_gmin]);
                                end
                                if local_gmax > 0
                                    local_gmax = 1.001*min([o2 - round(o2) + 3*o3, local_gmax]);
                                else
                                    local_gmax = 0.999*min([o2 - round(o2) + 3*o3, local_gmax]);
                                end
                                
                                gmns(k)  = o2 - round(o2);
                                gstds(k) = o3;
                                
                                if sv_plus
                                    hold on
                                    shade(Tx_list, o1*exp(-((Tx_list - o2).^2)/(2*o3.^2)),...
                                          Tx_list, zeros(size(Tx_list,2),1),...
                                          'FillType',[1 2;2 1], 'FillColor', 'b','FillAlpha',.4);
                                end
                                clear xx H o1 o2 o3 

                            end
                            
                        end
                        
                        if sv_plus
                            hold on 
                            plot([j + local_gmin, j + local_gmin],[0 max(yy)], '--', 'Color', 'r', 'LineWidth', 2)
                            hold on 
                            plot([j + local_gmax, j + local_gmax],[0 max(yy)], '--', 'Color', 'r', 'LineWidth', 2)
                        end
                    else
                        
                        if local_gmax < local_gmin                         % check if min < max 
                                   
                            temp = local_gmin;
                            local_gmin = local_gmax;
                            local_gmax = temp;
                            clear temp

                        end
                        
                        if local_gmin > 0
                            local_gmin = 0.999*local_gmin;
                        else
                            local_gmin = 1.001*local_gmin;
                        end
                        if local_gmax > 0
                            local_gmax = 1.001*local_gmax;
                        else
                            local_gmax = 0.999*local_gmax;
                        end
                    
                    end  
                    
                    % At this point the interrogation window have been
                    % optimized 
                    
                    Searched_nm = j;                                        % nominal mass to look for 
                    Searched_dm = (local_gmax + local_gmin)/2;              % Mass defect to look for 
                    
                    Searched_mz = Searched_nm + Searched_dm;                % Accurate mass to look for 
                    Tol = Threshold(j);                                     % Tolerence
                    [elemComp, Isel] = tof_find_elemental_composition(Searched_mz,Tol, mat_comb);
                    
                    % Considering the first best guess
                                
                    if isempty(elemComp) == 0 && (Counter >= 3)
                 
                        [~, choice] = min(abs(elemComp(:,1) - Searched_mz)); % To add : the one with minimun  number in atoms example: take CH not CHN (nombre de type d'atome)

                        [name, stored_shown_name] = Element_Composition_Reader_recent(elemComp(:,2:end), choice, vars_Elements);
                           
                        if sv_plus && (Counter >= 3)
                            
                            shg
                            hold on
                            plot([j + gmns(k), j + gmns(k)],[0 max(yy)], '-.', 'Color', 'w', 'LineWidth', 2)
                            hold on
                            plot([round(elemComp(choice,1)) + elemComp(choice,1) - round(elemComp(choice,1)), round(elemComp(choice,1)) + elemComp(choice,1) - round(elemComp(choice,1))], [0 max(yy)- 2], ':', 'Color', 'r', 'LineWidth', 2)
                            hold on
                            
                            text(round(elemComp(choice,1)) + elemComp(choice,1) - round(elemComp(choice,1)), max(yy) - 2, stored_shown_name, 'Interpreter', 'latex', 'FontSize', text_size,'Color','r');
                            
                            clear yy
                            
                            cd(Storage_path)
                            % Saving the plot in .fig format
                            savefig(['Peak_finding_histogram_',num2str(j),'_', num2str(k),'.fig'])
                            % Saving the plot in .png format
                            name_save = ['Peak_finding_histogram_',num2str(j),'_', num2str(k),'.png'];
                            print(name_save,'-dpng','-r300');
                            cd(Script_path)
                            close()

                        end
                        
                        hold on
                        scatter(round(elemComp(choice,1)), elemComp(choice,1) - round(elemComp(choice,1)),'^r', 'filled')
                        hold on
                        text(round(elemComp(choice,1)) + .25, elemComp(choice,1) - round(elemComp(choice,1)), stored_shown_name, 'Interpreter', 'latex', 'Fontsize',text_size-2,'Color','r');
                        
                        Proposed_Species = [Proposed_Species; cellstr(name)];
                        Proposed_Species_latex = [Proposed_Species_latex; cellstr(stored_shown_name)];
                        
                        clear name stored_shown_name 
                        
                    else
                        
                        if sv_plus && (Counter >= 3)
                            
                            cd(Storage_path)
                            % Saving the plot in .fig format
                            savefig(['Peak_finding_histogram_',num2str(j),'_', num2str(k),'.fig'])
                            % Saving the plot in .png format
                            name_save = ['Peak_finding_histogram_',num2str(j),'_', num2str(k),'.png'];
                            print(name_save,'-dpng','-r300');
                            cd(Script_path)
                            close()
                            
                        end
                        
                        Proposed_Species = [Proposed_Species; cellstr([''])];
                        Proposed_Species_latex = [Proposed_Species_latex; cellstr([''])];
                        
                    end 
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                shg
                rect = [nm_grid(j) + 0.25, gmin(k), dx/2, gmax(k) - gmin(k)];
                tf = 0;
                if (Counter >= 3) && (overlap_flag_for_rect == 0)
                    flag_counter = [flag_counter, 0];
                    rectangle('Position',rect,'EdgeColor','b','Linewidth',1)
                    for i = 1 : n_spec
                        % Here the interogation window is set
                        L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));
                        
                        if L1~=0 
                            ID = Peak_List_each_spectra{i}(L1,1);
                            [~,m] = max(Peak_List_each_spectra{i}(ID,3)); % if more than one ID is detected consider the one with maximum intensity
                            Temp_Peak_List{1,1}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),1); % ID
                            Temp_Peak_List{1,2}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),3); % Intensity
                            Temp_Peak_List{1,3}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),2); % m/z
                            Temp_Peak_List{1,4}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),8); % tof
                            Temp_Peak_List{1,5}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),10);% mass defect
                            Temp_Peak_List{1,6}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),9); % nominal mass

                        end
                        clear ID L1
                    end
                elseif (overlap_flag_for_rect == 0)
                    flag_counter = [flag_counter, 2];
                    rectangle('Position',rect,'EdgeColor','r','Linewidth',1)
                    for i = 1 : n_spec
                        % Here the interrogation window is set
                        L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));

                        if L1~=0 
                            ID = Peak_List_each_spectra{i}(L1,1);
                            [~,m] = max(Peak_List_each_spectra{i}(ID,3));  % if more than one ID is detected consider the one with maximum intensity
                            Temp_Peak_List{1,1}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),1); % ID
                            Temp_Peak_List{1,2}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),3); % Intensity
                            Temp_Peak_List{1,3}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),2); % m/z
                            Temp_Peak_List{1,4}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),8); % tof
                            Temp_Peak_List{1,5}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),10);% mass defect
                            Temp_Peak_List{1,6}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),9); % nominal mass

                        end
                        clear ID L1
                    end
                end
                if overlap_flag_for_rect == 1
                    rectangle('Position',rect,'EdgeColor','g','Linewidth',1)
                    flag_counter = [flag_counter, 1];
                    L0 = cell(n_spec,1);
                    for i = 1 : n_spec
                        % Here the interrogation window is set
                        L0{i} = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));
                    end
                    for i = 1 : n_spec
                        % Here the interrogation window is set
                        L1 = find( nm_grid(j) < Peak_List_each_spectra{i}(:,9) & Peak_List_each_spectra{i}(:,9) < nm_grid(j+1) & gmin(k) <= Peak_List_each_spectra{i}(:,10) & Peak_List_each_spectra{i}(:,10) <= gmax(k));
                        
                        if size(L1,1) > 1
                           % in this case Select the most intense peak
                           [~, d] = max(Peak_List_each_spectra{i}(L1,3));
                           L2 = L1(~d); % The one discarded
                           L1 = L1(d); % The one kept
                           hold on
                           scatter(Peak_List_each_spectra{i}(L2,9),Peak_List_each_spectra{i}(L2,10),Peak_List_each_spectra{i}(L2,11),'filled','om')
                           clear L2
                           hold on
                        end

                        if L1~=0 
                            ID = Peak_List_each_spectra{i}(L1,1);
                            [~,m] = max(Peak_List_each_spectra{i}(ID,3)); % if more than one ID is detected consider the one with maximum intensity
                            Temp_Peak_List{1,1}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),1); % ID
                            Temp_Peak_List{1,2}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),3); % Intensity
                            Temp_Peak_List{1,3}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),2); % m/z
                            Temp_Peak_List{1,4}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),8); % tof
                            Temp_Peak_List{1,5}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),10);% mass defect
                            Temp_Peak_List{1,6}(count,i) = Peak_List_each_spectra{i}(Peak_List_each_spectra{i}(:,1) == ID(m),9); % nominal mass

                        end
                        clear ID L1
                    end
                end
                
                % Making the video 
                    shg
                    if (j >= 12.5) &&(j <= param.Range(2) - 12.5)
                        xlim([j - 12.5  j + 12.5])
                    elseif j < 12.5 
                        xlim([0 25])
                    elseif j > param.Range(2) - 12.5
                        xlim([param.Range(2) - 25 param.Range(2)])
                    end
                    ylim([ ((-0.4 + 0.05))*j/(param.Range(2) - (param.Range(1) + 1)) - 0.05 - ((-0.4 + 0.05))/(param.Range(2) - (param.Range(1) + 1)) ,...
                           ((0.4 - 0.05))*j/(param.Range(2) - (param.Range(1) + 1)) + 0.05 - ((0.4 - 0.05))/(param.Range(2) - (param.Range(1) + 1))])
                    
                if sv_plus
                    
                    cd(Storage_path)
                    % Write each frame to the file.
                    currFrame = getframe(gcf);
                    writeVideo(vidObj,currFrame);
                    cd(Script_path)
                end
                %%%%%%%%%%%%%%%%%%
                
                % end
                iter = iter - 1;
                k = k + 1;
            end
            nom_mass(j)             = j;
            Cluster_stat_nGrps(j)   = nGrps;
            Cluster_stat_gmns{j}    = gmns;
            Cluster_stat_gstds{j}   = gstds;
            Cluster_stat_gmin{j}    = gmin;
            Cluster_stat_gmax{j}    = gmax;
        end
    end
    
    if sv_plus
        
        cd(Storage_path)
        close(vidObj);
        
        xlim([param.Range(1) param.Range(2)])
        ylim auto      
        
        cd(Script_path)
    end
    
    Peak_List = cell(1,6);
    % Loop over all samples
    for i = 1 : n_spec
        % Recognizing empty rows
        A = ismember(Temp_Peak_List{1},zeros(1, n_spec),'rows');
        for j = 1 : 6
            % 1 st column for ID
            % 2nd column for Intensity
            % 3rd column for m/z
            % 4th column for time of flight
            % 5th olumn for mass defect
            % 6th column for nominal mass

            % Keeping only non emty rows
            Peak_List{1,j} = Temp_Peak_List{1,j}(A==0,:);
            
        end
    end
    % Exporting a main peak list
    % 1st column : accurate mass (mz)
    % 2nd column : nominal mass (m/u)
    % 3rd column : mass defect
    % 4rd column : mass defect variation within a cluster
    % 5th column : Accuracy in ppm if the mean of the cluster corresponds
    % to a specie
    % 6th column : flag : 0 if peak exists in more than 3 mass spectra 
    %                     1 if peak exists in more than 3 mass spectra for
    %                     witch overlap was detected and the highest peak
    %                     was considered within the final interogation window 
    %                     2 if peak exists in less than 3 mass spectra
    disp('Preparing a summarized peak list...')
    Sum_Peak_list = zeros(size(Peak_List{1,1},1),6); % Initialization
   
    latest_count = 0;
    for j = 1 : Np_nm
       if Cluster_stat_nGrps(j) ~=0
           Table = [latest_count + 1 : 1 : Cluster_stat_nGrps(j) + latest_count];
           
           Sum_Peak_list(Table,1) = Cluster_stat_gmns{j}' + repmat(j, 1, Cluster_stat_nGrps(j));    % Cluster mean accurate mass
           Sum_Peak_list(Table,2) = repmat(j, 1, Cluster_stat_nGrps(j));                            % Cluster nominal mass 
           Sum_Peak_list(Table,3) = Cluster_stat_gmns{j};                                           % Cluster mean mass defect
           Sum_Peak_list(Table,4) = (Cluster_stat_gmax{j} - Cluster_stat_gmin{j})/2;                % Cluster dispression in m/z

           latest_count = Cluster_stat_nGrps(j) + latest_count;
           clear Table 
       end
    end
    Sum_Peak_list(:,5) = ((Sum_Peak_list(:,4)./Sum_Peak_list(:,1)))*1E6;   % Cluster dispression in ppm 
    Sum_Peak_list(:,6) = flag_counter;                                     % Cluster flag
    disp('Done')
    
end