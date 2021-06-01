function Input_menu()
    clc; clear;
    disp('--------------------------------------------------------')
    disp('PropSim is a simulator for solid propellant rocket motors');
    disp('Written by Daniel Cantos')
    disp('Last edition: 3 January 2021')
    disp('Version 1.1.1')
    disp('--------------------------------------------------------')
    out = false;
    manual_data = {}; % Void manual data
    f_change = ones(1,10); 
    % Vectorial flag which marks which preset fields to edit
    f_csv = false;
    
    disp('In which format do you want to save the results? (either .csv or .mat)');
    file_type = input('.mat is a native Matlab format while .csv is more compatible: ','s');
    if strcmpi('mat',file_type) || strcmpi('.mat',file_type) || strcmpi('m',file_type) || strcmpi('ma',file_type)
        f_csv = false;
    else % As default csv will be used
        f_csv = true;
    end
    
    while out == false % Repeat until valid choice
        disp('Choose data for problem definition:')
        input_type = input('1 - Manual data introduction / 2 - Introduce case preset: ');
        switch input_type
            case 1 % Manual loading
                disp('--------------------------------------------------------')
                disp('MANUAL LOADING SELECTED')
                manual_data = manual_loading(f_change,manual_data,f_csv);
                % Calls manual data to create presets
                data = manual_data;              
                out = true;
            case 2 % Preset loading
                disp('--------------------------------------------------------')
                disp('PRESET LOADING SELECTED')
                if f_csv
                    [file,path] = uigetfile('.\Libraries\Case\*.csv');
                    fullpath = cat(2,path,file);  
                    preset_data = readtable(fullpath);
                    % The preset is read
                    data.ambient = preset_data{1,1:2};
                    data.chamber.geometry = preset_data{2,1:4};
                    data.chamber.material.properties = preset_data{3,1:5};
                    data.chamber.material.yield = preset_data{5:6,1:preset_data{4,1}};
                    data.chamber.material.uts = preset_data{8:9,1:preset_data{7,1}};
                    data.propellant.geometry = preset_data{10,1:6};
                    data.propellant.properties = preset_data{11,1:17};
                    data.nozzle.geometry = preset_data{12,1:6};
                    data.nozzle.material.properties = preset_data{13,1:5};
                    data.nozzle.material.yield = preset_data{15:16,1:preset_data{14,1}};
                    data.nozzle.material.uts = preset_data{18:19,1:preset_data{17,1}};
                    data.liner = preset_data{20,1:3};
                    data.insulator = preset_data{21,1:3};
                    data.simulation = preset_data{22,1:4};
                    h = 1;
                else
                    [file,path] = uigetfile('.\Libraries\Case\');
                    fullpath = cat(2,path,file);
                    preset_data = load(fullpath);
                    data = preset_data.manual_data;               
                end     
                change_preset_edit = input('Do you want to edit the case? (Y/N): ','s');
                if change_preset_edit == 'y' || change_preset_edit == 'Y' 
                    f_change = zeros(1,10);
                    [f_change,out] = final_changes(change_preset_edit,data,f_change,f_csv);
                    manual_data = manual_loading(f_change,data,f_csv);
                    % Reuses manual_loading to edit specific fields
                    data = manual_data;
                else
                    out = true;
                end     
            otherwise
                disp('Wrong selection')
                clc; clear;
                disp('--------------------------------------------------------')
                out = false;
        end
    end
    
    out = false;
    
      
    
    function manual_data = manual_loading(f_change,manual_data,f_csv)
        out2 = false;  
        while out2 == false % Out when you confirm you don't want to change anything 
            if f_change(1) == 1 % Ambient properties
                disp('--------------------------------------------------------')
                disp('Ambient properties')
                disp('--------------------------------------------------------')
                manual_data.ambient(1) = input('Ambient temperature (K): ');
                manual_data.ambient(2) = input('Ambient pressure (bar): ');
            end
            
            if f_change(2) == 1 % Chamber properties
                disp('--------------------------------------------------------')
                disp('Chamber properties')
                disp('--------------------------------------------------------')
                out3 = false;
                while out3 == false % Dimension sense
                    manual_data.chamber.geometry(1) = input('External diameter (mm): ');
                    manual_data.chamber.geometry(2) = input('Internal diameter (mm): ');
                    if manual_data.chamber.geometry(1) <= manual_data.chamber.geometry(2)
                        % Internal diameter must be smaller than external
                        disp('Invalid dimensions!'); out3 = false;
                    else
                        out3 = true;
                    end
                end
                manual_data.chamber.geometry(3) = input('Chamber length (mm): ');
                manual_data.chamber.geometry(4) = input('Bulkhead union resistance (N): '); 
                
                out3 = false;
                while out3 == false % Material selection loop
                    disp('Select the material.');
                    material_load = input('1 - Load material preset / 2 - Create new material: ');
                    switch material_load
                        case 1 
                            if f_csv
                                [file2,path2] = uigetfile('.\Libraries\Mat\*.csv'); % Load preset
                                fullpath2 = cat(2,path2,file2);
                                mat_load = readtable(fullpath2);
                                manual_data.chamber.material.properties = mat_load{2,1:mat_load{1,1}};
                                manual_data.chamber.material.yield = mat_load{4:5,1:mat_load{3,1}};
                                manual_data.chamber.material.uts = mat_load{7:8,1:mat_load{6,1}};
                            else
                                [file2,path2] = uigetfile('.\Libraries\Mat\'); % Load preset
                                fullpath2 = cat(2,path2,file2); 
                                preset_mat = load(fullpath2);
                                manual_data.chamber.material = preset_mat.material;                               
                            end     
                            out3 = true;
                        case 2
                            while out3 == false % Create material
                                [a,b,c] = material_creation; % Calls function, importing and exporting the data
                                manual_data.chamber.material.properties = a;
                                manual_data.chamber.material.yield = b;
                                manual_data.chamber.material.uts = c;
 
                                mat_name = input('Provide a name for the material: ','s');
                                change = input('Do you want to change anything? (Y/N): ','s');
                                if change == 'y' || change == 'Y'
                                    out3 = false;
                                else % We'll assume n, N or any other thing
                                    out3 = true;
                                    path2 = '.\Libraries\Mat\';
                                    if f_csv
                                        full_mat_name = cat(2,path2,mat_name,'.csv');
                                        mat_save = zeros(8,max([length(a), length(b(1,:)), length(c(1,:))]));
                                        mat_save(1,1) = length(a(1,:));
                                        mat_save(2,1:mat_save(1,1)) = a;
                                        mat_save(3,1) = length(b(1,:));
                                        mat_save(4:5,1:mat_save(3,1)) = b(1:2,:);
                                        mat_save(6,1) = length(c(1,:));
                                        mat_save(7:8,1:mat_save(6,1)) = c(1:2,:);
                                        writematrix(mat_save,full_mat_name);
                                    else
                                        full_mat_name = cat(2,path2,mat_name,'.mat');
                                        material = manual_data.chamber.material;
                                        save(full_mat_name, 'material');  
                                    end   
                                end
                                
                            end
                            out3 = true;
                        otherwise
                            disp('Invalid selection')
                            out3 = false;
                    end
                end
                out3 = false;
            end
            if f_change(3) == 1 % Propellant properties
                disp('--------------------------------------------------------')
                disp('Propellant properties')
                disp('--------------------------------------------------------')
                out3 = false;
                while out3 == false
                    manual_data.propellant.geometry(1) = input('Dext (mm): ');
                    manual_data.propellant.geometry(2) = input('Dint (mm): ');
                    manual_data.propellant.geometry(3) = input('L0 (one grain) (mm): ');
                    manual_data.propellant.geometry(4) = input('Number of grains - N: ');
                    if manual_data.chamber.geometry(3) < manual_data.propellant...
                            .geometry(3)*manual_data.propellant.geometry(4) ||...
                            manual_data.propellant.geometry(1) > manual_data.chamber.geometry(2)
                        disp('Invalid grain geometry')
                        out3 = false;
                    else
                        out3 = true;
                    end
                end
                out3 = false;
                disp('Inhibited sides must be determined.')
                manual_data.propellant.geometry(5) = input...
                    ('Is the outside inhibited? (1 - Yes / 2 - No): ');
                manual_data.propellant.geometry(6) = input...
                    ('Are the ends inhibited? (1 - Yes / 2 - No): ');
                
                while out3 == false %Propellant selection loop
                    disp('Select the propellant.')
                    propellant_load = input('1 - Load propellant preset / 2 - Create new propellant: ');
                    switch propellant_load
                        case 1
                            if f_csv
                                [file3,path3] = uigetfile('.\Libraries\Prop\*.csv');
                                fullpath3 = cat(2,path3,file3);
                                prop_load = readtable(fullpath3);                                
                                manual_data.propellant.properties = ...
                                    prop_load{2,1:prop_load{1,1}}; % Just copies the full vector
                            else
                                [file3,path3] = uigetfile('.\Libraries\Prop\');
                                fullpath3 = cat(2,path3,file3);
                                preset_prop = load(fullpath3);
                                manual_data.propellant.properties = preset_prop.prop;
                            end      
                            out3 = true;
                        case 2 
                            while out3==false % Loop for propellant properties
                                prop_name = input('Provide a name for the propellant: ','s');
                                manual_data.propellant.properties(1) = input('Temperature coefficient a (mm/s·MPa^n): ');
                                manual_data.propellant.properties(2) = input('Combustion index n: ');
                                manual_data.propellant.properties(3) = input('Ideal density (kg/m3): ');
                                manual_data.propellant.properties(4) = input('Density ratio (0-1): ');
                                manual_data.propellant.properties(5) = input('Thermal conductivity propellant (W/mK): ');
                                manual_data.propellant.properties(6) = input('Heat capacity propellant (J/kgK): '); % c -> dQ = c·dT
                                manual_data.propellant.properties(7) = input('Mean molar mass gas (g/mol): ');
                                manual_data.propellant.properties(8) = input('Specific heat at constant volume gas (J/molK): ');
                                manual_data.propellant.properties(9) = input('Ratio of specific heats gas: ');
                                manual_data.propellant.properties(10) = input('Dynamic viscosity gas (Pa·s): ');
                                manual_data.propellant.properties(11) = input('Conduction coefficient gas (W/mK): ');
                                manual_data.propellant.properties(12) = input('Adiabatic flame temperature (K): ');
                                manual_data.propellant.properties(13) = input('Autoignition temperature (K): ');
                                manual_data.propellant.properties(14) = input('Erosive burn rate coefficient r (m/s, Pa): ');
                                manual_data.propellant.properties(15) = 53; % Always true
                                manual_data.propellant.properties(16) = input('Combustion efficiency (0-1): ');
                                manual_data.propellant.properties(17) = input('Mole fraction water in exhaust: ');
                                change = input('Do you want to change anything? (Y/N): ','s');
                                if change == 'y' || change == 'Y'
                                    out3 = false;
                                else % We'll assume n, N or any other thing
                                    out3 = true;
                                    path3 = '.\Libraries\Prop\';
                                    if f_csv
                                        full_prop_name = cat(2,path3,prop_name,'.csv');
                                        prop_save = zeros(2,17);
                                        prop_save(1,1) = 17;
                                        prop_save(2,1:length(prop_save(1,:))) = ...
                                            manual_data.propellant.properties;
                                        writematrix(prop_save,full_prop_name);
                                    else
                                        full_prop_name = cat(2,path3,prop_name,'.mat');
                                        prop = manual_data.propellant.properties;
                                        save(full_prop_name, 'prop');                                        
                                    end
                                end                              
                            end
                        otherwise
                            disp('Invalid selection')
                            out3 = false;
                    end
                end
            end
            
            if f_change(4) == 1 % Nozzle properties
                disp('--------------------------------------------------------')
                disp('Nozzle properties')
                disp('--------------------------------------------------------')  
                manual_data.nozzle.geometry(1) = input('Inlet diameter (mm): ');
                manual_data.nozzle.geometry(2) = input('Throat diameter (mm): ');
                manual_data.nozzle.geometry(3) = input('Exit diameter (mm): ');
                manual_data.nozzle.geometry(4) = input('Angle convergent section (º): ');
                manual_data.nozzle.geometry(5) = input('Angle divergent section (º): ');
                manual_data.nozzle.geometry(6) = input('Nozzle union resistance (N): ');
                
                out3 = false;
                while out3 == false % Material selection loop
                    disp('Select the material.');
                    material_load = input('1 - Load material preset / 2 - Create new material: ');
                    switch material_load
                        case 1 
                            if f_csv
                                [file2,path2] = uigetfile('.\Libraries\Mat\*.csv'); % Load preset
                                fullpath2 = cat(2,path2,file2);
                                mat_load = readtable(fullpath2);
                                manual_data.nozzle.material.properties = mat_load{2,1:mat_load{1,1}};
                                manual_data.nozzle.material.yield = mat_load{4:5,1:mat_load{3,1}};
                                manual_data.nozzle.material.uts = mat_load{7:8,1:mat_load{6,1}};                             
                            else
                                [file2,path2] = uigetfile('.\Libraries\Mat\'); % Load preset
                                fullpath2 = cat(2,path2,file2); 
                                preset_mat = load(fullpath2);
                                manual_data.nozzle.material = preset_mat.material;
                            end 
                            out3 = true;
                        case 2
                            while out3 == false % Create material
                                [a,b,c] = material_creation; 
                                % Calls function, to determine material
                                % properties
                                manual_data.nozzle.material.properties = a;
                                manual_data.nozzle.material.yield = b;
                                manual_data.nozzle.material.uts = c;
                                  
                                mat_name = input('Provide a name for the material: ','s');
                                change = input('Do you want to change anything? (Y/N): ','s');
                                if change == 'y' || change == 'Y'
                                    out3 = false;
                                else % We'll assume n, N or any other thing
                                    out3 = true;
                                    path2 = '.\Libraries\Mat\'; 
                                    if f_csv
                                        full_mat_name = cat(2,path2,mat_name,'.csv');
                                        mat_save = zeros(8,max([length(a), length(b(1,:)), length(c(1,:))]));
                                        mat_save(1,1) = length(a(1,:));
                                        mat_save(2,1:mat_save(1,1)) = a;
                                        mat_save(3,1) = length(b(1,:));
                                        mat_save(4:5,1:mat_save(3,1)) = b(1:2,:);
                                        mat_save(6,1) = length(c(1,:));
                                        mat_save(7:8,1:mat_save(6,1)) = c(1:2,:);
                                        writematrix(mat_save,full_mat_name);
                                    else
                                        full_mat_name = cat(2,path2,mat_name,'.mat');
                                        material = manual_data.chamber.material;
                                        save(full_mat_name, 'material');  
                                    end                                     
                                end  
                            end
                            out3 = true;
                        otherwise
                            disp('Invalid selection')
                            out3 = false;
                    end
                end
                out3 = false;
            end
            
            
            if f_change(5) == 1 % Liner & insulator properties
                out3 = false;
                disp('--------------------------------------------------------')
                disp('Liner & insulator properties')
                disp('--------------------------------------------------------')
                while out3 == false        
                    manual_data.liner(1) = input('Thickness liner (mm): ');
                    %manual_data.liner(2) = input('Heat capacity liner (J/kgK): ');
                    %manual_data.liner(3) = input('Thermal conductivity liner (W/mK): ');  
                    manual_data.liner(2) = -1; % UNUSED
                    manual_data.liner(3) = -1; % UNUSED
                    
                    manual_data.insulator(1) = input('Thickness insulator (mm): ');
                    %manual_data.insulator(2) = input('Heat capacity insulator (J/kgK): ');
                    %manual_data.insulator(3) = input('Thermal conductivity insulator (W/mK): ');
                    manual_data.insulator(2) = -1; % UNUSED
                    manual_data.insulator(3) = -1; % UNUSED
                    
                    if manual_data.propellant.geometry(1)+2*manual_data.liner(1)+...
                            2*manual_data.insulator(1) > manual_data.chamber.geometry(2)
                        disp('Invalid dimensions.')
                        out3 = false;
                    else
                        change = input('Do you want to change anything? (Y/N): ','s');
                        if change == 'y' || change == 'Y'
                            out3 = false;
                        else % We'll assume n, N or any other thing 
                            % These values are only saved as part of case
                            out3 = true;                            
                        end
                    end
                    
                end
            end
            
            if f_change(6) == 1 % Simulation control
                disp('--------------------------------------------------------')
                disp('Simulation parameters')
                disp('--------------------------------------------------------')
                out3 = false;
                while out3 == false
                    manual_data.simulation(1) = input('Time step (s): ');
                    manual_data.simulation(2) = input('Number of FVM elements: ');
                    %manual_data.simulation(3) = input('Number of propellant elements: ');
                    manual_data.simulation(3) = -1; % UNUSED
                    manual_data.simulation(4) = input('Simulated time (s) (should be as close as possible to the expected burn time): ');                 
                    change = input('Do you want to change anything? (Y/N): ','s');
                    if change == 'y' || change == 'Y'
                        out3 = false;
                    else % We'll assume n, N or any other thing
                        out3 = true;
                    end
                end
                out3 = false;
            end
            disp('--------------------------------------------------------')
            disp('Final changes');
            disp('--------------------------------------------------------')
            change = input('Do you want to change anything? (Y/N): ','s');
            f_change = zeros(1,10);
            [f_change,out2] = final_changes(change,manual_data,f_change,f_csv);
        end
    end

    function [a,b,c] = material_creation % Only adds the material properties
        a(1) = input('Young Modulus (GPa): ');
        a(2) = input('Density (kg/m3): ');
        a(3) = input('Poisson ratio: ');
        a(4) = input('Melt temp (ºC): ');
        a(5) = input('Thermal conductivity (W/mK): ');
        disp('Yield and UTS as function of T must be inputted now.')
                                
        n_data_yield = input('Number of yield data points: ');
        b = zeros(2,n_data_yield);
        for i=1:n_data_yield
            b(1,i) = input('Temperature (ºC): ');
            b(2,i) = input('Yield strength (MPa): ');
        end
                                
        n_data_uts = input('Number of UTS data points: ');
        c = zeros(2,n_data_uts);
        for i=1:n_data_uts
            c(1,i) = input('Temperature (ºC): ');
            c(2,i) = input('UTS (MPa): ');
        end        
    end

    function [f_change,out2] = final_changes(change,manual_data,f_change,f_csv) 
        % Allows user to revisit selected fields and rewrite parameters
        if change == 'y' || change == 'Y'
            out3 = false; out2 = false;
            while out3 == false
                f_change_sel = input('Select data to change: 1-Ambient, 2-Chamber, 3-Propellant, 4-Nozzle, 5-Insulation, 6-Simulation, 7-Done: ');
                % Selects which fields to change. Once done input 7
                switch f_change_sel
                    case 1
                        f_change(1) = 1;
                    case 2
                        f_change(2) = 1;
                    case 3
                        f_change(3) = 1;
                    case 4
                        f_change(4) = 1;
                    case 5
                        f_change(5) = 1;
                    case 6
                        f_change(6) = 1;
                    case 7
                        out3 = true;
                    otherwise
                        disp('Invalid selection.')
                        out3 = false;
                end
            end
        else % We'll assume n, N or any other thing
            disp('No change selected.')
            out2 = true;
            save_case = input('Do you want to save the new case? (Y/N): ','s');
            if save_case == 'y' || save_case == 'Y'
                case_name = input('Provide a name for the case: ','s');
                path4 = '.\Libraries\Case\'; 
                if f_csv
                    full_case_name = cat(2,path4,case_name,'.csv');
                    case_save = zeros(22,max([length(manual_data.propellant.properties),...
                        length(manual_data.chamber.material.yield(1,:)),...
                        length(manual_data.nozzle.material.uts(1,:))]));
                    % Ambient
                    case_save(1,1:2) = manual_data.ambient;
                    % Chamber
                    case_save(2,1:4) = manual_data.chamber.geometry;
                    case_save(3,1:5) = manual_data.chamber.material.properties;
                    case_save(4,1) = length(manual_data.chamber.material.yield(1,:));
                    case_save(5:6,1:case_save(4,1)) = manual_data.chamber.material.yield(1:2,:);
                    case_save(7,1) = length(manual_data.chamber.material.uts(1,:));
                    case_save(8:9,1:case_save(7,1)) = manual_data.chamber.material.uts(1:2,:);
                    % Propellant
                    case_save(10,1:6) = manual_data.propellant.geometry;
                    case_save(11,1:17) = manual_data.propellant.properties;
                    % Nozzle
                    case_save(12,1:6) = manual_data.nozzle.geometry;
                    case_save(13,1:5) = manual_data.nozzle.material.properties;
                    case_save(14,1) = length(manual_data.nozzle.material.yield(1,:));
                    case_save(15:16,1:case_save(14,1)) = manual_data.nozzle.material.yield(1:2,:);
                    case_save(17,1) = length(manual_data.nozzle.material.uts(1,:));
                    case_save(18:19,1:case_save(17,1)) = manual_data.nozzle.material.uts(1:2,:);
                    % Liner & Insulator
                    case_save(20,1:3) = manual_data.liner;
                    case_save(21,1:3) = manual_data.insulator;
                    case_save(22,1:4) = manual_data.simulation;
                    
                    writematrix(case_save,full_case_name);
                else
                    full_case_name = cat(2,path4,case_name,'.mat');
                    save(full_case_name, 'manual_data');                    
                end
            end
        end  
    end
    
    save('_data.mat', 'data')
end