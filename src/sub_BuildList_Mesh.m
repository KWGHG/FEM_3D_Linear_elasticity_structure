function [list_Node, list_Element]= sub_BuildList_Mesh(Data, NEN, line_range_Node, line_range_Element, num_Node, num_Element)
% 2022/01/21
% initialize list----------------------------------------------------------
list_Node     = cell(1+num_Node, 4);
list_Element  = cell(1+num_Element, NEN+1);

list_Node{1, 1} = 'node#';
list_Node{1, 2} = 'x-coord.';
list_Node{1, 3} = 'y-coord.';
list_Node{1, 4} = 'z-coord.';
list_Element{1, 1} = 'element#';
list_Element{1, 2} = 'node#_1';
list_Element{1, 3} = 'node#_2';
list_Element{1, 4} = 'node#_3';
list_Element{1, 5} = 'node#_4';
list_Element{1, 6} = 'node#_5';
list_Element{1, 7} = 'node#_6';
list_Element{1, 8} = 'node#_7';
list_Element{1, 9} = 'node#_8';

% Partition data-----------------------------------------------------------
Data_Node     = split(Data(line_range_Node(1):line_range_Node(2)), ',');
Data_Element  = split(Data(line_range_Element(1):line_range_Element(2)), ',');

% Complete list (from data)------------------------------------------------
for k = 1: num_Node
    for p = 1:4
        list_Node{k+1, p} = str2double(Data_Node(k, p));
    end
end
for k = 1: num_Element
    for p = 1:NEN+1
        list_Element{k+1, p} = str2double(Data_Element(k, p));
    end
end

end