function [name, stored_shown_name] = Element_Composition_Reader_recent(elemComp, choice, vars_Elements)

        % This is an internal routine for the APCFA toolbox.
        % The main routine to execute the entire toolbox is APCFA_toolbox

        name = [''];                                                       % Initialization
        stored_shown_name = [''];                                          % Initialization                
        
        Ind = find(elemComp(choice,:));      
        
        if (isempty(Ind) == 0)
            
            stored_shown_name = [stored_shown_name, ''];
            
            % Getting a given order
            Order = [33, 32, 31, 30, 29, 28, 27, 26, 24, 23, 25, 22, 21, 20, 19, 18, 9, 17, 16, 13, 12, 11, 14, 15, 10, 8, 2, 4, 3, 1, 6, 5, 7];
            [~, Ordered_Ind] = ismember(Ind, Order);
            
            for i = sort(Ordered_Ind)

                if (elemComp(choice,Order(i)) == 1)
                    name  = strcat(name ,vars_Elements(Order(i)));
                    stored_shown_name = strcat(stored_shown_name, vars_Elements(Order(i)));
                else 
                    stored_shown_name = strcat(stored_shown_name , vars_Elements(Order(i)));
                    stored_shown_name = strcat(stored_shown_name , '$$_{', num2str(elemComp(choice,Order(i))),'}$$');

                    name  = strcat(name ,vars_Elements(Order(i)));
                    name  = strcat(name , num2str(elemComp(choice,Order(i))));
                end

            end
            
            stored_shown_name = strcat(stored_shown_name, '$$^+$$');
            
        end
        
        clear Ind
        
        
end