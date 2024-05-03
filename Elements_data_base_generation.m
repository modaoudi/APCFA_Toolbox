% This is an internal routine for the APCFA toolbox.
% The main routine to execute the entire toolbox is APCFA_toolbox 

disp('Generating combination of elements ...')
disp('This might take few minutes...')
tic

maxN  = 1;
maxO  = 1;
maxC  = 60;
maxH  = 60;
maxTi = 7;
maxP  = 1;

nrCombination=(maxN+1)*(maxO+1)*(maxC+1)*(maxH+1);

H    = tof_exact_mass('H');
B    = tof_exact_mass('B');
C    = tof_exact_mass('C');
C13  = tof_exact_mass('[13C]');
N    = tof_exact_mass('N');
N15  = tof_exact_mass('[15N]');
O    = tof_exact_mass('O');
F    = tof_exact_mass('F');
Na   = tof_exact_mass('Na');
Al   = tof_exact_mass('Al');
Si   = tof_exact_mass('Si');
Si29 = tof_exact_mass('[29Si]');
Si30 = tof_exact_mass('[30Si]');
P    = tof_exact_mass('P');
S    = tof_exact_mass('S');
Cl   = tof_exact_mass('Cl');
Cl37 = tof_exact_mass('[37Cl]');
K    = tof_exact_mass('K');
K41  = tof_exact_mass('[41K]');
Ca   = tof_exact_mass('Ca');
Ca42 = tof_exact_mass('[42Ca]');
Ca44 = tof_exact_mass('[44Ca]');
Ti46 = tof_exact_mass('[46Ti]');
Ti47 = tof_exact_mass('[47Ti]');
Ti   = tof_exact_mass('Ti');
Ti49 = tof_exact_mass('[49Ti]');
Ti50 = tof_exact_mass('[50Ti]');
Cr   = tof_exact_mass('Cr');
Fe  = tof_exact_mass('Fe');
Fe58  = tof_exact_mass('[58Fe]');
Cu   = tof_exact_mass('Cu');
Cu65 = tof_exact_mass('[65Cu]');
Cs   = tof_exact_mass('Cs');

el   = 0.00054858;

maxAmu = param.Range(2);
mat_comb = [];                                                   % Initialization

% Adding CnH2n+2 groupe

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

for c = 0 : maxC 
    if c == 0
        for h = 1 : 2*c + 2
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

            end

        end
    else
        for h = 0 : 2*c + 2
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

            end

        end        
        
    end
end

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

% Adding CnH2n+2 groupe with 13C isotope
for c13 = 2 : maxC 
  
    for h = 4 : 2*c13 + 2
        M = H*h + B*b + C*(c13-1) + C13*1 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
        iamu = round(M);
        if iamu < maxAmu & iamu > 0

            mat_comb =  [mat_comb;[M h b c13-1 1 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

        end

    end        

end

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

% Adding CnH2nO groupe 

for c = 0 : maxC 
    if c == 0
        for h = 1 : 2*c
            for o = 1 : maxO
                M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
                iamu = round(M);
                if iamu < maxAmu & iamu > 0

                    mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

                end
            end
        end
    else
        for h = 0 : 2*c
            for o = 1 : maxO
                M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
                iamu = round(M);
                if iamu < maxAmu & iamu > 0

                    mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

                end
            end

        end        
        
    end
end

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

% Adding and TinPO2n+2 groupe with 46 Ti iosotope

for ti46 = 1:maxTi 
    for p = 0 : maxP
        for o = 0: 2*ti46 + 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*1 + Ti47*ti47 + Ti*(ti46-1) + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 1 ti47 ti46-1 ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2n+2 and TinPO2n+2 groupe with 47 Ti iosotope

for ti47 = 1:maxTi 
    for p = 0 : maxP
        for o = 0: 2*ti47 + 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*1 + Ti*(ti47-1) + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 1 (ti47-1) ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2n+2 and TinPO2n+2 groupe 

for ti = 1:maxTi 
    for p = 0 : maxP
        for o = 0: 2*ti + 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2n+2 and TinPO2n+2 groupe with 49 Ti iosotope

for ti49 = 1:maxTi 
    for p = 0 : maxP
        for o = 0: 2*ti49 + 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*(ti49-1) + Ti49*1 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti49-1 1 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2n+2 and TinPO2n+2 groupe with 50 Ti iosotope

for ti50 = 1:maxTi 
    for p = 0 : maxP
        for o = 0: 2*ti50 + 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*(ti50-1) + Ti49*ti49 + Ti50*1 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti50-1 ti49 1 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2nHn-2 group with 46 Ti isotope

for ti46 = 1:maxTi 
    for o = 1 : 2*ti46 + 2 
        for h = 0 : 2*ti46 - 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*1 + Ti47*ti47 + Ti*(ti46-1) + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 1 ti47 ti46-1 ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2nHn-2 group with 47 Ti isotope

for ti47 = 1:maxTi 
    for o = 1 : 2*ti47 + 2 
        for h = 0 : 2*ti47 - 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*1 + Ti*(ti47-1) + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 1 ti47-1 ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2nHn-2 group 

for ti = 1:maxTi 
    for o = 1 : 2*ti + 2 
        for h = 0 : 2*ti - 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*ti + Ti49*ti49 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti ti49 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2nHn-2 group with 49 Ti isotope

for ti49 = 1:maxTi 
    for o = 1 : 2*ti49 + 2 
        for h = 0 : 2*ti49 - 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*(ti49-1) + Ti49*1 + Ti50*ti50 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti49-1 1 ti50 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

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

% Adding TinO2nHn-2 group with 50 Ti isotope

for ti50 = 1:maxTi 
    for o = 1 : 2*ti50 + 2 
        for h = 0 : 2*ti50 - 2 
            M = H*h + B*b + C*c + C13*c13 + N*n + N15*n15 + O*o + F*f + Na*na + Al*al + Si*si + Si29*si29 + Si30*si30 + P*p + S*s + Cl*cl + Cl37*cl37 +...
                K*k + k41*k41 + Ca*ca + Ca42*ca42 + Ca44*ca44 + Ti46*ti46 + Ti47*ti47 + Ti*(ti50-1) + Ti49*ti49 + Ti50*1 + Cr*cr + Fe*fe + Fe58*fe58 + Cu*cu + Cu65*cu65 + Cs*cs - el;
            iamu = round(M);
            if iamu < maxAmu & iamu > 0

                mat_comb =  [mat_comb;[M h b c c13 n n15 o f na al si si29 si30 p s cl cl37 k k41 ca ca42 ca44 ti46 ti47 ti50-1 ti49 1 cr fe fe58 cu cu65 cs]];

            end
        end
    end
end

disp('Done')
toc
% Saving 
% Saving the Peak_List cell in .mat format
name_save = ('mat_comb');
save(name_save,'mat_comb')

