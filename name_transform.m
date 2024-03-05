function [mat_comb_rev] = name_transform(mz, namEl, numEl)

    mat_comb_rev = [];

    for i = 1: size(namEl,1)
        
        % Initialization
        h    = 0;
        b    = 0;
        c    = 0;
        c13  = 0;
        n    = 0;
        n15  = 0;
        o    = 0;
        f    = 0;
        na   = 0;
        al   = 0;
        si   = 0;
        si29 = 0;
        si30 = 0;
        p    = 0;
        s    = 0;
        cl   = 0;
        cl37 = 0;
        k    = 0;
        k41  = 0;
        ca   = 0;
        ca42 = 0;
        ca44 = 0;
        ti46 = 0;
        ti47 = 0;
        ti   = 0;
        ti49 = 0;
        ti50 = 0;
        cr   = 0;
        fe   = 0;
        fe58 = 0;
        cu   = 0;
        cu65 = 0;
        cs   = 0;
        
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[13C]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                c13    = numEl{i}(find(loc));
            end
        else
            c13    = 0;
        end
        
        [atome,loc] = ismember(namEl{i}, 'C');
        if (isempty(atome) == 0) 
            if (find(loc) ~=0)
                c    = numEl{i}(find(loc));
            end
        else
            c    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'H');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                h    = numEl{i}(find(loc));
            end
        else
            h    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'O');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                o    = numEl{i}(find(loc));
            end
        else
            o    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[15N]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                n15    = numEl{i}(find(loc));
            end
        else
            n15    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'N');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                n    = numEl{i}(find(loc));
            end
        else
            n    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'P');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                p    = numEl{i}(find(loc));
            end
        else
            p    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Ti');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti    = numEl{i}(find(loc));
            end
        else
            ti    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[46Ti]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti46     = numEl{i}(find(loc));
            end
        else
            ti46     = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[47Ti]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti47     = numEl{i}(find(loc));
            end
        else
            ti47     = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Ti');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti    = numEl{i}(find(loc));
            end
        else
            ti    = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[49Ti]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti49     = numEl{i}(find(loc));
            end
        else
            ti49     = 0;
        end

        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[50Ti]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ti50     = numEl{i}(find(loc));
            end
        else
            ti50     = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'B');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                b    = numEl{i}(find(loc));
            end
        else
            b    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'S');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                s    = numEl{i}(find(loc));
            end
        else
            s    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Na');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                na    = numEl{i}(find(loc));
            end
        else
            na    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Cl');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cl    = numEl{i}(find(loc));
            end
        else
            cl    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[37Cl]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cl37    = numEl{i}(find(loc));
            end
        else
            cl37    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[41K]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                k41    = numEl{i}(find(loc));
            end
        else
            k41    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'K');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                k    = numEl{i}(find(loc));
            end
        else
            k    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[42Ca]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ca42    = numEl{i}(find(loc));
            end
        else
            ca42    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[44Ca]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ca44    = numEl{i}(find(loc));
            end
        else
            ca44    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Ca');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                ca    = numEl{i}(find(loc));
            end
        else
            ca    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'F');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                f    = numEl{i}(find(loc));
            end
        else
            f    = 0;
        end
        
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Al');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                al    = numEl{i}(find(loc));
            end
        else
            al    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Si');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                si    = numEl{i}(find(loc));
            end
        else
            si    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[29Si]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                si29    = numEl{i}(find(loc));
            end
        else
            si29    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[30Si]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                si30    = numEl{i}(find(loc));
            end
        else
            si30    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Cs');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cs    = numEl{i}(find(loc));
            end
        else
            cs    = 0;
        end
        
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Cr');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cr    = numEl{i}(find(loc));
            end
        else
            cr    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[58Fe]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                fe58    = numEl{i}(find(loc));
            end
        else
            fe58    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Fe');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                fe    = numEl{i}(find(loc));
            end
        else
            fe    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, '[65Cu]');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cu65    = numEl{i}(find(loc));
            end
        else
            cu65    = 0;
        end
        
        clear atome loc
        [atome,loc] = ismember(namEl{i}, 'Cu');
        if (isempty(atome) == 0)
            if (find(loc) ~=0)
                cu    = numEl{i}(find(loc));
            end
        else
            cu    = 0;
        end
        
        mat_comb_rev =  [mat_comb_rev;[mz(i) h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

    end
end