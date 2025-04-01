function [list_Surface] = sub_BuildList_Surface...
    (Data, num_Surface, line_Surface, line_Instance, Instance_name)
% 2022/01/21
% initialize list----------------------------------------------------------
Data_Surface  = cell(num_Surface, 1);
surface_content = cell(num_Surface, 1);

list_Surface  = cell(1+num_Surface, 4);
list_Surface{1, 1} = 'name';
list_Surface{1, 2} = '# of content';
list_Surface{1, 3} = 'elset';
list_Surface{1, 4} = 'face#';

if num_Surface ~= 0
% Partition data-----------------------------------------------------------
endPart = line_Surface(num_Surface) +1;

for k = 1: num_Surface
    if line_Surface(k) > line_Instance
        endPart = k;
    end
    Data_Surface{k, 1} = Data(line_Surface(k));
    tt = 0;
    for p = line_Surface(k)+1: length(Data)
        if contains(Data(p), '*') ==1
            break
        end
        tt = tt + 1;
        surface_content{k, 1}{tt} = Data(p);
    end
end

% Complete list (from data)------------------------------------------------
tt = 1;
for k = 1: num_Surface
    data_Surface_temp = split(Data_Surface{k}, ',');
    surf_name_temp = split(data_Surface_temp(3), '=');
    surf_name_temp = surf_name_temp(end);
    tt = tt + 1;
    if k < endPart
        list_Surface{tt, 1} = strcat(Instance_name,'.',surf_name_temp(end));
    else
        list_Surface{tt, 1} = surf_name_temp(end);
    end
    num_content = length(surface_content{k});
    list_Surface{tt, 2} = num_content;

    for p = 1: num_content
        cc = tt + p - 1;
        surface_temp = surface_content{k}{p};
        surface_temp = split(surface_temp, ',');
        if k < endPart
            list_Surface{cc, 3} = strcat(Instance_name,'.',surface_temp(1));
        else
            list_Surface{cc, 3} = surface_temp(1);
        end
        list_Surface{cc, 4} = surface_temp(2);
    end
end
end

end