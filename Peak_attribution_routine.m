function Peak_attribution_routine(t, stop, Peak_List, text_size, Sum_Peak_list, vars_Elements)

    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    % This step is user assisted as it requires the user to decide what
    % chemical specie to be attributed to a given point in the mass defet 
    % plot. To ease this process an interface is proposed with a global
    % peak list, with the first column editable to put the appropriate 
    % chemical specie along with a live update of the mass defect plot, to
    % indicate for the user the accuracy of his choice. Additionnaly, 
    % aligned & calibrated mass spectra is plotted as well so, the user can
    % inspect the data if needed. At the end of this routine a Peak list 
    % is generated for each mass spectra for further statistical post 
    % processing. Data is saved in . dat. xlsx and .txt files for storage 
    % During this process the function tof_exacte_mass, provided by Heikki 
    % Junninen PhD work, is used.

    if stop == 0
        % Create table array
        global t
        % Create UI figure
        fig = uifigure('Name','Species attribution - GUI ');
        fig.Position(1:4) = get(0,'ScreenSize')/2;

        btn = uibutton('Parent',fig,...
            'Position',[20 950 100 22],... % Finish button position
            'Text','Finished',...
            'ButtonPushedFcn', @(btn,event) Stopbutton(btn));

        % Create table UI component
        uit = uitable(fig);
                
        uit.Data = t;
        uit.ColumnSortable = false;
        uit.ColumnEditable = [false true false false false false false false];
        uit.Position(1:4) = [100 100 400 800]; % Table position
        uit.DisplayDataChangedFcn = @updatePlot;
        
        
        % Create scatter plot
        ax1 = uiaxes(fig);
        % Creating 'ax' twin of 'ax1' for other permanent plots, 'ax1' is
        % the dynamics one
        i = 1;
        while  i == 1
            ax = ax1;            
            i = 0;
        end
        ax.Position(1:4) = [550 100 400 800]; % mass defect plot Figure Position
        ax.XLabel.String = 'Nominal mass';
        ax.YLabel.String = 'Mass defect';

        x1 = round(tof_exact_mass(t.Species));
        y1 = tof_exact_mass(t.Species) - round(tof_exact_mass(t.Species));
        sz1 = 100;
        
        mat1 = Peak_List{6};
        mat2 = Peak_List{5};
        mat3 = Peak_List{2};
        
        scatter(ax,mat1(mat3~=0),mat2(mat3~=0),log10(mat3(mat3~=0))*2,'filled','ob')
        hold(ax, 'on')
        xlim(ax,[0 25])
        ylim(ax, 'auto')
        
        % Getting those that have been attributed
        % 0: if attributted
        % 1: if still not attributted
        loc = cellfun('isempty', t{:,'Species'});
        loc_1 = repmat(loc,1,size(mat3,2));
        loc_1 = ~loc_1;
        loc_2 = mat3~=0;
        loc_3 = loc_1.*loc_2;
        
        grid(ax, 'on')
        h = scatter(ax1,x1,y1,sz1,'filled','^m');
        hold(ax, 'on')
        h_2 = scatter(ax1,mat1(loc_3~=0),mat2(loc_3~=0),log10(mat3(loc_3~=0))*2,'filled','or'); 
        
        
        [mz, namEl, numEl] = tof_exact_mass(t.Species);
        mat_comb_rev = name_transform(mz, namEl, numEl);
        
        stored_shown_name   = [''];

        for i = 1 : size(mat_comb_rev,1)

            [~, temp_2]    = Element_Composition_Reader_recent(mat_comb_rev(i,2:end),1, vars_Elements(2:end));
            stored_shown_name   = [stored_shown_name;{temp_2}];

            clear temp_1 temp_2
        end

%         a = text(ax1,x1+0.25,y1,t.Species, 'Interpreter', 'latex','Fontsize',14,'Color','r');
        a = text(ax1,x1+0.25,y1,stored_shown_name,'Interpreter', 'latex', 'FontSize', 14,'Color','r');

    end
    % Create the function for the ButtonPushedFcn callback
    function Stopbutton(btn)
        Peak_attribution_routine(t, 1)
        close all force
        disp('Species attribution process finished')
        return
    end

    % Create an updating function
    % Update the scatter plot when table data changes
    function updatePlot(src,event)
        delete(h)
        delete(h_2)
        delete(a)
        
        t = uit.DisplayData;
        
        t.Accuracy(find(~cellfun(@isempty,t.Species))) = ((tof_exact_mass(t.Species(find(~cellfun(@isempty,t.Species)))) - t.m_z(find(~cellfun(@isempty,t.Species))))./t.m_z(find(~cellfun(@isempty,t.Species))))*1E6;
        uit.Data = t;
        
        x1 = round(tof_exact_mass(t.Species));
        y1 = tof_exact_mass(t.Species) - round(tof_exact_mass(t.Species));
        sz1 = 100;
        
        % Getting those that have been attributed
        % 0: if attributted
        % 1: if still not attributted
        loc = cellfun('isempty', t{:,'Species'});
        loc_1 = repmat(loc,1,size(mat3,2));
        loc_1 = ~loc_1;
        loc_2 = mat3~=0;
        loc_3 = loc_1.*loc_2;
        
        h = scatter(ax1,x1,y1,sz1, 'filled', '^m');
        hold(ax1, 'on')
        h_2 = scatter(ax1,mat1(loc_3~=0),mat2(loc_3~=0),log10(mat3(loc_3~=0))*2,'filled','or');         

        [mz, namEl, numEl] = tof_exact_mass(t.Species);
        mat_comb_rev = name_transform(mz, namEl, numEl);
        
        stored_shown_name   = [''];

        for i = 1 : size(mat_comb_rev,1)

            [~, temp_2]    = Element_Composition_Reader_recent(mat_comb_rev(i,2:end),1, vars_Elements(2:end));
            stored_shown_name   = [stored_shown_name;{temp_2}];

            clear temp_1 temp_2
        end

%         a = text(ax1,x1+0.25,y1,t.Species, 'Interpreter', 'latex','Fontsize',14,'Color','r');
        a = text(ax1,x1+0.25,y1,stored_shown_name,'Interpreter', 'latex', 'FontSize', 14,'Color','r');

    end

    return

end