function [list_BC, list_Gravity, list_Pressure, list_Surftraction, list_Concentrated]...
    = sub_BuildList_Load(Data, IMO, line_range_BC, line_range_Loads)
% 2022/01/21
% initialize list----------------------------------------------------------
list_BC           = cell(IMO, 3);
list_Gravity      = cell(2, 5);
list_Pressure     = cell(IMO, 5);
list_Surftraction = cell(IMO, 4);
list_Concentrated = cell(IMO, 4);

list_BC{1, 1} = 'object';
list_BC{1, 2} = 'first degree of freedom';
list_BC{1, 3} = 'last degree of freedom';
list_Gravity{1, 1} = 'name';
list_Gravity{1, 2} = 'object';
list_Gravity{1, 3} = 'g-const';
list_Gravity{1, 4} = 'comp1';
list_Gravity{1, 5} = 'comp2';
list_Gravity{1, 6} = 'comp3';
list_Pressure{1, 1} = 'name';
list_Pressure{1, 2} = 'option';
list_Pressure{1, 3} = 'object';
list_Pressure{1, 4} = 'type';
list_Pressure{1, 5} = 'magnitude';
list_Surftraction{1, 1} = 'name';
list_Surftraction{1, 2} = 'option';
list_Surftraction{1, 3} = 'object';
list_Surftraction{1, 4} = 'type';
list_Surftraction{1, 5} = 'magnitude';
list_Surftraction{1, 6} = 'comp1';
list_Surftraction{1, 7} = 'comp2';
list_Surftraction{1, 8} = 'comp3';
list_Concentrated{1, 1} = 'name';
list_Concentrated{1, 2} = 'object';
list_Concentrated{1, 3} = 'degree of freedom';
list_Concentrated{1, 4} = 'magnitude';

% Partition data-----------------------------------------------------------
Data_BC = Data(line_range_BC(1): line_range_BC(2));
Data_Loads = Data(line_range_Loads(1): line_range_Loads(2));

% Complete list (from data)------------------------------------------------
count_bc = 1;
count_pressure = 1;
count_surftraction = 1;
count_concentrated = 1;

for k = 1: length(Data_BC)
    if contains(Data_BC(k), '** Name') == 1
        if contains(Data_BC(k), 'Displacement/Rotation') == 1
            for p = k+2: length(Data_BC)
                if contains(Data_BC(p), '**') == 1
                    break
                end
                count_bc = count_bc + 1;
                data_bc_temp = split(Data_BC(p), ',');
                objec_temp = data_bc_temp(1);
                list_BC{count_bc, 1} = objec_temp(end);
                list_BC{count_bc, 2} = str2double(data_bc_temp(2));
                list_BC{count_bc, 3} = str2double(data_bc_temp(3));
            end
        elseif contains(Data_BC(k), 'Symmetry/Antisymmetry/Encastre') == 1
            for p = k+2: length(Data_BC)
                if contains(Data_BC(p), '**') == 1
                    break
                end
                count_bc = count_bc + 1;
                data_bc_temp = split(Data_BC(p), ',');
                objec_temp = data_bc_temp(1);
                list_BC{count_bc, 1} = objec_temp(end);
                list_BC{count_bc, 2} = 1;
                list_BC{count_bc, 3} = 2;
            end
        end
    end
end
for k = 1: length(Data_Loads)
    if contains(Data_Loads(k), '** Name') == 1
        if contains(Data_Loads(k), 'Gravity') == 1
            load_name_temp = split(Data_Loads(k), ' ');
            data_loads_temp = split(Data_Loads(k+2), ',');
            objec_temp = data_loads_temp(1);
            gravity_const_temp = str2double(data_loads_temp(3));
            list_Gravity{2, 1} = load_name_temp(3);
            list_Gravity{2, 2} = objec_temp(end);
            list_Gravity{2, 3} = gravity_const_temp;
            list_Gravity{2, 4} = str2double(data_loads_temp(4));
            list_Gravity{2, 5} = str2double(data_loads_temp(5));
            list_Gravity{2, 6} = str2double(data_loads_temp(6));
        elseif contains(Data_Loads(k), 'Pressure') == 1
            load_name_temp = split(Data_Loads(k), ' ');
            if contains(Data_Loads(k+1), 'Dload') ==1
                option_temp = 'Dload';
            elseif contains(Data_Loads(k+1), 'Dsload') ==1
                option_temp = 'Dsload';
            end

            for p = k+2: length(Data_Loads)
                if contains(Data_Loads(p), '**')
                    break
                end
                data_loads_temp = split(Data_Loads(p), ',');
                objec_temp = data_loads_temp(1);

                count_pressure = count_pressure + 1;
                list_Pressure{count_pressure, 1} = load_name_temp(3);
                list_Pressure{count_pressure, 2} = option_temp;
                list_Pressure{count_pressure, 3} = objec_temp(end);
                list_Pressure{count_pressure, 4} = data_loads_temp(2);
                list_Pressure{count_pressure, 5} = str2double(data_loads_temp(3));

            end

        elseif contains(Data_Loads(k), 'Surface traction') == 1
            load_name_temp = split(Data_Loads(k), ' ');
            if contains(Data_Loads(k+1), 'Dload') == 1
                option_temp = 'Dload';
            elseif contains(Data_Loads(k+1), 'Dsload') == 1
                option_temp = 'Dsload';
            end
            
            for p = k+2: length(Data_Loads)
                if contains(Data_Loads(p), '**')
                    break
                end

                data_loads_temp = split(Data_Loads(p), ',');
                objec_temp = data_loads_temp(1);

                count_surftraction = count_surftraction + 1;
                list_Surftraction{count_surftraction, 1} = load_name_temp(3);
                list_Surftraction{count_surftraction, 2} = option_temp;
                list_Surftraction{count_surftraction, 3} = objec_temp(end);
                list_Surftraction{count_surftraction, 4} = data_loads_temp(2);
                list_Surftraction{count_surftraction, 5} = str2double(data_loads_temp(3));
                list_Surftraction{count_surftraction, 6} = str2double(data_loads_temp(4));
                list_Surftraction{count_surftraction, 7} = str2double(data_loads_temp(5));
                list_Surftraction{count_surftraction, 8} = str2double(data_loads_temp(6));
            end
        elseif contains(Data_Loads(k), 'Concentrated force') == 1
            count_concentrated = count_concentrated + 1;
            load_name_temp = split(Data_Loads(k), ' ');
            data_loads_temp = split(Data_Loads(k+2), ',');
            objec_temp = data_loads_temp(1);
            list_Concentrated{count_concentrated, 1} = load_name_temp(3);
            list_Concentrated{count_concentrated, 2} = objec_temp(end);
            list_Concentrated{count_concentrated, 3} = str2double(data_loads_temp(2));
            list_Concentrated{count_concentrated, 4} = str2double(data_loads_temp(3));
        end
    end
end

% delete space of list
list_BC(count_bc+1:end, :) = [];
list_Pressure(count_pressure+1:end, :) = [];
list_Surftraction(count_surftraction+1:end, :) = [];
list_Concentrated(count_concentrated+1:end, :) = [];
end