function [list_Material, Instance_name] = sub_BuildList_Assembly...
    (Data, line_range_Material, line_Instance)
% 2022/01/21
% initialize list----------------------------------------------------------
list_Material = cell(2, 4);

list_Material{1, 1} = 'name';
list_Material{1, 2} = 'density';
list_Material{1, 3} = 'E';
list_Material{1, 4} = 'Poisson ratio';

% Partition data-----------------------------------------------------------
Data_Material = Data(line_range_Material(1):line_range_Material(2));
Data_Instance = Data(line_Instance);

% Complete list (from data)------------------------------------------------
load_name_temp = split(Data_Material{1}, '=');
list_Material{2, 1} = load_name_temp(end);
for k = 2: length(Data_Material)
    if contains(Data_Material{k}, '*Density') ==1
        density_temp = split(Data_Material{k+1}, ',');
        list_Material{2, 2} = str2double(density_temp(1));
    end
    if contains(Data_Material{k}, '*Elastic') ==1
        elastic_temp = split(Data_Material{k+1}, ',');
        list_Material{2, 3} = str2double(elastic_temp(1));
        list_Material{2, 4} = str2double(elastic_temp(2));
    end
end

% find instance name
Instance_name = split(Data_Instance, ',');
Instance_name = split(Instance_name(2), '=');
Instance_name = Instance_name(end);

end