function [list_Set] = sub_BuildList_Set(Data, IMO, num_Set, num_Nset, ...
    num_Elset, line_Nset, line_Elset, line_Instance, Instance_name)
% 2022/01/21
% initialize list----------------------------------------------------------
Data_Nset     = cell(num_Nset, 1);
Data_Elset    = cell(num_Elset, 1);
Nset_content  = cell(num_Nset, 1);
Elset_content = cell(num_Elset, 1);

list_Set      = cell(1+num_Set, IMO*10);
list_Set{1, 1} = 'type';
list_Set{1, 2} = 'name';
list_Set{1, 3} = '# of content';
list_Set{1, 4} = 'content~';

if num_Set ~=0
% Partition data-----------------------------------------------------------
endPart_n = line_Nset(num_Nset) +1;
endPart_e = line_Elset(num_Elset) +1;

for k = 1: num_Nset
    Data_Nset{k, 1} = Data(line_Nset(k));
    tt = 0;
    for p = line_Nset(k)+1: length(Data)
        if contains(Data(p), '*') ==1
            break
        end
        tt = tt + 1;
        Nset_content{k, 1}{tt} = Data(p);
    end
end
for k = 1: num_Elset
    Data_Elset{k, 1} = Data(line_Elset(k));
    tt = 0;
    for p = line_Elset(k)+1: length(Data)
        if contains(Data(p), '*') ==1
            break
        end
        tt = tt + 1;
        Elset_content{k, 1}{tt} = Data(p);
    end
end

% endPart is used to determine if the set's name need adding part's name
for k = 1: num_Nset
    if line_Nset(k) > line_Instance
        endPart_n = k;
        break
    end
end
for k = 1: num_Elset
    if line_Elset(k) > line_Instance
        endPart_e = k;
        break
    end
end

% Complete list (from data)------------------------------------------------
set_col_max = 0;
for k = 1: num_Nset
    data_Nset_temp = split(Data_Nset{k}, ',');
    set_name_temp = split(data_Nset_temp(2), '=');
    generate_key = any(contains(data_Nset_temp, 'generate'));
    list_Set{k+1, 1} = 'Nset';
    if k < endPart_n
        list_Set{k+1, 2} = strcat(Instance_name,'.',set_name_temp(end));
    else
        list_Set{k+1, 2} = set_name_temp(end);
    end
    
    if generate_key ==1
        nset_temp = str2double(split(Nset_content{k}{1}, ','));
        count = 3;
        for p = nset_temp(1):nset_temp(3):nset_temp(2)
            count = count + 1;
            list_Set{k+1, count} = p;
            list_Set{k+1, 3} = count - 3;
            if count > set_col_max
                set_col_max = count;
            end
        end
    else
        num_content = length(Nset_content{k});
        tt = 3;
        for m = 1: num_content
            nset_temp = str2double(split(Nset_content{k}{m}, ','));
            for p = 1:length(nset_temp)
                if isnan(nset_temp(p)) ~=1
                    tt = tt + 1;
                    list_Set{k+1, tt} = nset_temp(p);
                    list_Set{k+1, 3} = tt - 3;
                    if p+3 > set_col_max
                        set_col_max = p+3;
                    end
                end
            end
        end
    end
end
for k = 1: num_Elset
    data_Elset_temp = split(Data_Elset{k}, ',');
    set_name_temp = split(data_Elset_temp(2), '=');
    generate_key = any(contains(data_Elset_temp, 'generate'));
    list_Set{k+1+num_Nset, 1} = 'Elset';
    if k < endPart_e
        list_Set{k+1+num_Nset, 2} = strcat(Instance_name,'.',set_name_temp(end));
    else
        list_Set{k+1+num_Nset, 2} = set_name_temp(end);
    end
    
    if generate_key ==1
        elset_temp = str2double(split(Elset_content{k}{1}, ','));
        count = 3;
        for p = elset_temp(1):elset_temp(3):elset_temp(2)
            count = count + 1;
            list_Set{k+1+num_Nset, count} = p;
            list_Set{k+1+num_Nset, 3} = count - 3;
            if count > set_col_max
                set_col_max = count;
            end
        end
    else
        num_content = length(Elset_content{k});
        tt = 3;
        for m = 1: num_content
            elset_temp = str2double(split(Elset_content{k}{m}, ','));
            for p = 1:length(elset_temp)
                if isnan(elset_temp(p)) ~=1
                    tt = tt + 1;
                    list_Set{k+1+num_Nset, tt} = elset_temp(p);
                    list_Set{k+1+num_Nset, 3} = tt - 3;
                    if p+3 > set_col_max
                        set_col_max = p+3;
                    end
                end
            end
        end
    end
end

% delete space of list
list_Set(:, set_col_max+1:end) = [];

end

end