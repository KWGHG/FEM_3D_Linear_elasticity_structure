function [NEN, num_Node, num_Element, num_Nset, num_Elset, num_Set, num_Surface, ...
    num_Surftraction, num_Pressure, num_Concentrated, num_BC, ...
    list_Node, list_Element, list_Material, Instance_name, list_Surface, list_Set, ...
    list_BC, list_F_Gravity, list_F_Pressure, list_F_Surftraction, list_F_Concentrated] = ...
    sub_Data_Prepropose(inp_name)

Data = readlines(inp_name);

% Initial Memory Option----------------------------------------------------
IMO = 0.25;   % defult value of initial memory for variable
% If you wanna modify value of IMO, doing it at upper value
IMO = round(length(Data) * IMO);   % Don't modify this
% Key words definition-----------------------------------------------------
str_Node     = '*Node';
str_Element  = '*Element';
str_Nset     = '*Nset';
str_Elset    = '*Elset';
str_Surface  = '*Surface';
str_Instance = '*Instance';
str_Material = '*Material';
str_BC       = '** BOUNDARY CONDITIONS';
str_Loads    = '** LOADS';
str_End      = '** OUTPUT REQUESTS';

% Scan Data for lines that contains Key words------------------------------
% "line_range" record the range of line contains related information
% "line" record the # of line that contains related information
line_range_Node    = zeros(2, 1);
line_range_Element = zeros(2, 1);
line_range_BC      = zeros(2, 1);
line_range_Loads   = zeros(2, 1);
line_range_Material  = zeros(2, 1);
line_Instance  = zeros(1, 1);
line_Nset      = zeros(IMO, 1);
line_Elset     = zeros(IMO, 1);
line_Surface   = zeros(IMO, 1);

% "count" is used for record times that key word appeard
count_nset     = 0;
count_elset    = 0;
count_surface  = 0;

% scaning all data
for k = 1: length(Data)
    if contains(Data(k), str_Node) == 1
        line_range_Node(1) = k+1;
    elseif contains(Data(k), str_Element) ==1
        Element_type = split(Data(k), '=');
        Element_type = Element_type(end);
        line_range_Node(2) = k-1;
        line_range_Element(1) = k+1;
        for p = k+1: length(Data)
            if contains(Data(p), '*') ==1
                line_range_Element(2) = p-1;
                break
            end
        end
    elseif contains(Data(k), str_Nset) ==1
        count_nset = count_nset + 1;
        line_Nset(count_nset) = k;
    elseif contains(Data(k), str_Elset) ==1
        count_elset = count_elset + 1;
        line_Elset(count_elset) = k;
    elseif contains(Data(k), str_Surface) ==1
        count_surface = count_surface + 1;
        line_Surface(count_surface) = k;
    elseif contains(Data(k), str_Instance) ==1
        line_Instance(1) = k;
    elseif contains(Data(k), str_Material) ==1
        line_range_Material(1) = k;
        for p = k+2: length(Data)
            if contains(Data(p), '** ') ==1
                line_range_Material(2) = p-1;
                break
            end
        end
    elseif contains(Data(k), str_BC) == 1
        line_range_BC(1) = k+1;
    elseif contains(Data(k), str_Loads) ==1
        line_range_BC(2) = k-1;
        line_range_Loads(1) = k+1;
    elseif contains(Data(k), str_End) ==1
        line_range_Loads(2) = k-1;
    end
end

% determine element type
if strcmp(Element_type, 'C3D8') || strcmp(Element_type, 'C3D8R')
    NEN = 8;
else
    error('**Error! Unknown element type. Check your input file!')
end

% Bulid list---------------------------------------------------------------
% "num_xxx": the number of xxx
num_Node     = line_range_Node(2) - line_range_Node(1) + 1;
num_Element  = line_range_Element(2) - line_range_Element(1) + 1;
num_Nset     = count_nset;
num_Elset    = count_elset;
num_Set      = num_Nset + num_Elset;
num_Surface  = count_surface;
% if you need to modify the data structure of list, modify the sub-function
[list_Node, list_Element] = ...
    sub_BuildList_Mesh(Data, NEN, line_range_Node, line_range_Element, num_Node, num_Element);
[list_Material, Instance_name] = ...
    sub_BuildList_Assembly(Data, line_range_Material, line_Instance);
[list_Surface] = ...
    sub_BuildList_Surface(Data, num_Surface, line_Surface, line_Instance, Instance_name);
[list_Set] = ...
    sub_BuildList_Set(Data, IMO, num_Set, num_Nset, num_Elset, line_Nset, line_Elset, line_Instance, Instance_name);
[list_BC, list_F_Gravity, list_F_Pressure, list_F_Surftraction, list_F_Concentrated] = ...
    sub_BuildList_Load(Data, IMO, line_range_BC, line_range_Loads);

% "num_xxx": the number of xxx
[num_Surftraction, ~] = size(list_F_Surftraction);
num_Surftraction = num_Surftraction -1;
[num_Pressure, ~] = size(list_F_Pressure);
num_Pressure = num_Pressure -1;
[num_Concentrated, ~] = size(list_F_Concentrated);
num_Concentrated = num_Concentrated -1;
[num_BC, ~] = size(list_BC);
num_BC = num_BC -1;

disp('The model has been imported from an input file.')
end