clear all; close all; clc;
%% Read inp file data
% File Name
inp_name = './input/Mesh_005.inp';
disp("Read .inp file name : " + inp_name);

% Data Pre-process function
[NEN,num_Node, num_Element, num_Nset, num_Elset, num_Set, num_Surface, ...
    num_Surftraction, num_Pressure, num_Concentrated, num_BC, ...
    list_Node, list_Element, list_Material, Instance_name, list_Surface, list_Set, ...
    list_BC, list_F_Gravity, list_F_Pressure, list_F_Surftraction, list_F_Concentrated] = ...
    sub_Data_Prepropose(inp_name);

%% Basic parameters for FEM calculation
% material property
elastic_modulus = list_Material{2, 3};
poisson_ratio   = list_Material{2, 4};
density         = list_Material{2, 2};

% node & element & integration point
Node = cell2mat(list_Node(2:end, :));
Element = cell2mat(list_Element(2:end, 1:end));
% integration_points' position in physical domain
num_INT_physical = num_Element * NEN;
% INT = [INT#, element#, x, y, z]
INT_physical = zeros(num_INT_physical, 5);   

% Element information (C3D8)
line_queue = [1,2;2,3;3,4;4,1;1,5;2,6;3,7;4,8;5,6;6,7;7,8;8,5];
patch_queue = [1,2,3,4;5,8,7,6;1,5,6,2;2,6,7,3;3,7,8,4;4,8,5,1];

%% Loading
% pressure
Pressure_cell = {};
for i = 2:size(list_F_Pressure, 1)
    P_keyword = string(list_F_Pressure(i, 1));
    for j = 2:size(list_Set, 1)
        if contains(string(list_Set{j, 2}), P_keyword)
            Pressure_cell{end+1, 1} = list_Set(j, :);
        end
    end
end

% define pressure at which face, element, value
Pressure_element_info = {};
num_normal = 0;
for i = 1:size(Pressure_cell,1)
    for j = 2:size(list_F_Pressure,1)
        keyword = string(list_F_Pressure(j, 1));
        cell_temp = Pressure_cell{i,1};
        num_normal = num_normal+cell_temp{3};
        P_info_temp = zeros(cell_temp{3},2);
        P_info_temp(:,1) = cell2mat(cell_temp(4:end));
        if contains(string(cell_temp(2)),keyword)
            P_info_temp(:,2) = cell2mat(list_F_Pressure(j, 5));
            if contains(string(cell_temp(2)),'S1')
                P_Surface = 1;
            elseif contains(string(cell_temp(2)),'S2')
                P_Surface = 2;
            elseif contains(string(cell_temp(2)),'S3')
                P_Surface = 3;
            elseif contains(string(cell_temp(2)),'S4')
                P_Surface = 4;
            elseif contains(string(cell_temp(2)),'S5')
                P_Surface = 5;
            elseif contains(string(cell_temp(2)),'S6')
                P_Surface = 6;
            end
        end
    end
    
    Pressure_element_info{i,1} = P_Surface; 
    Pressure_element_info{i,2} = P_info_temp;
end

% gravity
gravity = -list_F_Gravity{2,3};

%% boundary condition of displacement
for k = 2: num_BC+1
    object = list_BC{k, 1};
    for p = 2: num_Set+1
        if strcmp(list_Set{p, 2}, object) && strcmp(list_Set{p, 1}, 'Nset')
            nset_num = list_Set{p, 3};
            nset = cell2mat(list_Set(p, 4:4+nset_num-1));
        end
    end
end

% bc_dof = [node_number, Dof];
bc_dof = zeros(size(nset,2)*3,2);
ic = 0;
for i = 1:length(nset)
    for j = 1:3
        ic = ic+1;
        bc_dof(ic,:) = [nset(i),j];
    end
end

%% Define matrix

% declarae U, E, S---------------------------------------------------------
% U = [node#, U1, U2, U3, U_magnitude]
U_deformed = zeros(num_Node, 5);    

% strain and stress on integration point
% E_INT = [INT#, E11, E22, E33, E12, E13, E23]
% S_INT = [INT#, S11, S22, S33, S12, S13, S23]
E_INT = zeros(NEN*num_Element, 7);
S_INT = zeros(NEN*num_Element, 7);

% strain and stress on element
% E_e = [element#, E11, E22, E33, E12, E13, E23]
% S_e = [element#, S11, S22, S33, S12, S13, S23]
E_element = zeros(num_Element, 7);    
S_element = zeros(num_Element, 7);    

% von-Mises stress on integration point 
% Sv_INT = [INT#, S_Mises]
Sv_INT = zeros(NEN*num_Element, 2);  
% von-Mises stress on element
% Sv_element = [element#, S_Mises]
Sv_element = zeros(num_Element, 2);    

% principle stress on integration point
% Sp_INT = [INT#, MaxP, MinP]
Sp_INT = zeros(NEN*num_Element, 3);  
% principle stress on element 
% Sp_element = [element#, sigma1, sigma2, sigma3]
Sp_element = zeros(num_Element, 4);   

% Shear stress on integration point
% Ss_INT = [INT#, MaxSs, MinSs]
Ss_INT = zeros(NEN*num_Element, 3);  
% Shear stress on element
% Ss_element = [element#, MaxSsX, MinSsX, MaxSsY MinSsY MaxSsZ MinSsz]
Ss_element = zeros(num_Element, 7);   

disp('The job has been created and submitted for analysis.')

%%Main code----------------------------------------------------------------
lamda = poisson_ratio * elastic_modulus / (1+poisson_ratio) / (1-2*poisson_ratio);
mu = 0.5 * elastic_modulus / (1+poisson_ratio);
D_matrix = zeros(6, 6);
D_matrix(1, 1) = lamda + 2*mu; D_matrix(1, 2) = lamda; D_matrix(1, 3) = lamda;
D_matrix(2, 1) = lamda; D_matrix(2, 2) = lamda + 2*mu; D_matrix(2, 3) = lamda;
D_matrix(3, 1) = lamda; D_matrix(3, 2) = lamda; D_matrix(3, 3) = lamda + 2*mu;
D_matrix(4, 4) = mu; D_matrix(5, 5) = mu; D_matrix(6, 6) = mu;
X = zeros(8, 4); %%XYZ positions
K_localMatrix = zeros(24, 24);%%K element matrix
K_gobalMatrix = zeros(3 * size(Node, 1), 3 * size(Node, 1));
Fb_localMatrix = zeros(24, 1);%%degree of freedom * 8nodes
Fh_localMatrix = zeros(24, 1);%%degree of freedom * 8nodes
F_gobalMatrix = zeros(3 * size(Node, 1), 1);

%%K_matrix && Fq matrix----------------------------------------------------

for i = 1:size(Element, 1)
    bodyForce = [0; 0; gravity * density];%%bodyForce
    for j = 1:8
        X(j, 1) = Node(Element(i, j+1), 1);%%Node
        X(j, 2) = Node(Element(i, j+1), 2);%%x coordinate
        X(j, 3) = Node(Element(i, j+1), 3);%%y coordinate
        X(j, 4) = Node(Element(i, j+1), 4);%%z coordinate
    end
    % X = sortX(X);
    K_localMatrix = K_e(D_matrix, X);
    Fb_localMatrix = Fb_e(bodyForce, X);
    %%index for local to gobal matrix
    index = [3*X(1, 1)-2; 3*X(1, 1)-1; 3*X(1, 1);
             3*X(2, 1)-2; 3*X(2, 1)-1; 3*X(2, 1);
             3*X(3, 1)-2; 3*X(3, 1)-1; 3*X(3, 1);
             3*X(4, 1)-2; 3*X(4, 1)-1; 3*X(4, 1);
             3*X(5, 1)-2; 3*X(5, 1)-1; 3*X(5, 1);
             3*X(6, 1)-2; 3*X(6, 1)-1; 3*X(6, 1);
             3*X(7, 1)-2; 3*X(7, 1)-1; 3*X(7, 1);
             3*X(8, 1)-2; 3*X(8, 1)-1; 3*X(8, 1)];
             
    %%Assemble
    for j = 1:24
        for k = 1:24
            K_gobalMatrix(index(j), index(k)) = K_gobalMatrix(index(j), index(k)) + K_localMatrix(j, k);  
        end
        F_gobalMatrix(index(j), 1) = F_gobalMatrix(index(j), 1) + Fb_localMatrix(j, 1);
    end
end

%% Fh matrix----------------------------------------------------------------

for i = 1:size(Pressure_element_info, 1)
    indexFace = Pressure_element_info{i, 1};
    %disp('indexFace');
    for j = 1:size(Pressure_element_info{i, 2}, 1)
        indexElement = Pressure_element_info{i, 2}(j, 1);
        index= faceIndex(indexFace);
        pressure = Pressure_element_info{i, 2}(j, 2);
        elementNode = Element(indexElement, 2:9);
        X = zeros(8, 4);
        for k = 1:8
            X(k, 1) = Node(elementNode(k), 1);%%Node
            X(k, 2) = Node(elementNode(k), 2);%%x coordinate
            X(k, 3) = Node(elementNode(k), 3);%%y coordinate
            X(k, 4) = Node(elementNode(k), 4);%%z coordinate
        end
        NodeA = X(index(1), 2:4);
        NodeB = X(index(2), 2:4);
        NodeC = X(index(3), 2:4);
        NodeD = X(index(4), 2:4);
        normalVector1 = -cross(NodeA - NodeB, NodeC - NodeB);
        normalVector2 = -cross(NodeB - NodeC, NodeD - NodeC);
        normalVector3 = -cross(NodeC - NodeD, NodeA - NodeD);
        normalVector4 = -cross(NodeD - NodeA, NodeB - NodeA);
        Length_ = [norm(NodeA - NodeB) norm(NodeB - NodeC) norm(NodeC - NodeD) norm(NodeD - NodeA)];
        n1 = norm(normalVector1);
        n2 = norm(normalVector2);
        n3 = norm(normalVector3);
        n4 = norm(normalVector4);
        Area = 0.25 * (n1 + n2 + n3 + n4);
        normalVector1 = [(normalVector1(1) / n1) (normalVector1(2) / n1) (normalVector1(3) / n1)];
        normalVector2 = [(normalVector2(1) / n2) (normalVector2(2) / n2) (normalVector2(3) / n2)];
        normalVector3 = [(normalVector3(1) / n3) (normalVector3(2) / n3) (normalVector3(3) / n3)];
        normalVector4 = [(normalVector4(1) / n4) (normalVector4(2) / n4) (normalVector4(3) / n4)];
        normalVector = 0.25 * (normalVector1 + normalVector2 + normalVector3 + normalVector4);

        tractionForce = Pressure_element_info{i, 2}(j, 2) * normalVector.';
        % XX = sortX(X);
        XX = X;
        % newFace = sortFace(index, X, XX);
        Fh_localMatrix = Fh_e(tractionForce, XX, indexFace);

        %%index for local to gobal matrix
        index2 = [3*XX(1, 1)-2; 3*XX(1, 1)-1; 3*XX(1, 1);
                  3*XX(2, 1)-2; 3*XX(2, 1)-1; 3*XX(2, 1);
                  3*XX(3, 1)-2; 3*XX(3, 1)-1; 3*XX(3, 1);
                  3*XX(4, 1)-2; 3*XX(4, 1)-1; 3*XX(4, 1);
                  3*XX(5, 1)-2; 3*XX(5, 1)-1; 3*XX(5, 1);
                  3*XX(6, 1)-2; 3*XX(6, 1)-1; 3*XX(6, 1);
                  3*XX(7, 1)-2; 3*XX(7, 1)-1; 3*XX(7, 1);
                  3*XX(8, 1)-2; 3*XX(8, 1)-1; 3*XX(8, 1)];

    %%Assemble
        for k = 1:24
            F_gobalMatrix(index2(k), 1) = F_gobalMatrix(index2(k), 1) + Fh_localMatrix(k, 1);
        end
    end
end

%%penalty method to restrict degree of freedoms----------------------------
maxK_gobal = max(K_gobalMatrix);
index = zeros(size(bc_dof, 1), 1);
for i = 1:size(bc_dof, 1)
    if bc_dof(i, 2) == 1
        index (i)= 3 * bc_dof(i, 1) - 2;
    elseif bc_dof(i, 2) == 2
        index(i) = 3 * bc_dof(i, 1) - 1;
    else
        index(i) = 3 * bc_dof(i, 1);
    end
    % K_gobalMatrix(index, index) = 1*10^12;
    % F_gobalMatrix(index, 1) = 0.0;
end

for i = 1:size(index, 1)
    K_gobalMatrix(index(i), index(i)) = 1e10;
    F_gobalMatrix(index(i), 1) = 0.0;
end
displacement = K_gobalMatrix\F_gobalMatrix;
index = 1;
for i = 1:size(Node, 1)
    %%Node#
    U_deformed(i, 1) = i;

    %%Ux
    U_deformed(i, 2) = displacement(index, 1);
    U_deformed(i, 3) = displacement(index+1, 1);
    U_deformed(i, 4) = displacement(index+2, 1);
    U_deformed(i, 5) = sqrt(U_deformed(i, 2)^2 + U_deformed(i, 3)^2 + U_deformed(i, 4)^2);
    index = index + 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Max deformed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_U_deformed = sqrt(U_deformed.^2);
disp('U_deformed x-dir MAX');
disp(max(abs_U_deformed(:, 2)));
disp('U_deformed y-dir MAX');
disp(max(abs_U_deformed(:, 3)));
disp('U_deformed z-dir MAX');
disp(max(abs_U_deformed(:, 4)));
disp('U_deformed Magnitude MAX');
disp(max(abs_U_deformed(:, 5)));

%%Calculate strain and stress----------------------------------------------

for i = 1:size(Element, 1)
    X = zeros(8, 3); %%XYZ positions
    localDisplacement = zeros(24, 1);

    indexNatural_Coor = [-0.5773502692 -0.5773502692 0.5773502692;
                         0.5773502692 -0.5773502692 0.5773502692;
                         0.5773502692 0.5773502692 0.5773502692;
                         -0.5773502692 0.5773502692 0.5773502692;
                         -0.5773502692 -0.5773502692 -0.5773502692;
                         0.5773502692 -0.5773502692 -0.5773502692;
                         0.5773502692 0.5773502692 -0.5773502692;
                         -0.5773502692 +0.5773502692 -0.5773502692;];
    strain = zeros(6, 8);
    stress = zeros(6, 8);

    for j = 1:8
        X(j, 1) = Node(Element(i, j+1), 1);%%Node
        X(j, 2) = Node(Element(i, j+1), 2);%%x coordinate
        X(j, 3) = Node(Element(i, j+1), 3);%%y coordinate
        X(j, 4) = Node(Element(i, j+1), 4);%%z coordinate

        localDisplacement(3*j-2, 1) = U_deformed(X(j, 1), 2);
        localDisplacement(3*j-1, 1) = U_deformed(X(j, 1), 3);
        localDisplacement(3*j, 1) = U_deformed(X(j, 1), 4);
    end

    for j = 1:8
        indexNode = Node(Element(i, j+1), 1);%%Node
        B_ = B_matrix(indexNatural_Coor(j, 1), indexNatural_Coor(j, 2), indexNatural_Coor(j, 3), X);
        strain(:, j) = B_ * localDisplacement;
        stress(:, j) = D_matrix * strain(:, j);

        E_INT(i+j-1, 1) = indexNode;
        E_INT(i+j-1, 2) = strain(1, j);
        E_INT(i+j-1, 3) = strain(2, j);
        E_INT(i+j-1, 4) = strain(3, j);
        E_INT(i+j-1, 5) = strain(4, j);
        E_INT(i+j-1, 6) = strain(5, j);
        E_INT(i+j-1, 7) = strain(6, j);

        S_INT(i+j-1, 1) = indexNode;
        S_INT(i+j-1, 2) = stress(1, j);
        S_INT(i+j-1, 3) = stress(2, j);
        S_INT(i+j-1, 4) = stress(3, j);
        S_INT(i+j-1, 5) = stress(4, j);
        S_INT(i+j-1, 6) = stress(5, j);
        S_INT(i+j-1, 7) = stress(6, j);
    end
    strain_ = mean(strain, 2);
    stress_ = mean(stress, 2);

    E_element(i, 1) = i;
    E_element(i, 2) = strain_(1);
    E_element(i, 3) = strain_(2);
    E_element(i, 4) = strain_(3);
    E_element(i, 5) = strain_(4);
    E_element(i, 6) = strain_(5);
    E_element(i, 7) = strain_(6);

    

    S_element(i, 1) = i;
    S_element(i, 2) = stress_(1);
    S_element(i, 3) = stress_(2);
    S_element(i, 4) = stress_(3);
    S_element(i, 5) = stress_(4);
    S_element(i, 6) = stress_(5);
    S_element(i, 7) = stress_(6);
    
    %%principle stress-----------------------------------------------------
    stress_matrix = [stress_(1) stress_(4) stress_(5);
                     stress_(4) stress_(2) stress_(6);
                     stress_(5) stress_(6) stress_(3)];
    eStress_matrix = eig(stress_matrix);
    sigma1 = max(eStress_matrix);
    sigma3 = min(eStress_matrix);
    sigma2 = stress_matrix(1, 1) + stress_matrix(2, 2) + stress_matrix(3, 3) - sigma1 - sigma3;
    Sp_element(i, 1) = i;
    Sp_element(i, 2) = sigma1;
    Sp_element(i, 3) = sigma2;
    Sp_element(i, 4) = sigma3;

    %%von-Mises stress-----------------------------------------------------
    J2 = ((sigma1 - sigma2)^2 + (sigma2 - sigma3)^2 + (sigma3 - sigma1)^2)/6;
    sigmaV = sqrt(3 * J2);
    Sv_element(i, 1) = i;
    Sv_element(i, 2) = sigmaV;

    %%Max shear stress-----------------------------------------------------
    Max_TauX = 0.5 * abs(sigma2 - sigma3); Min_TauX = -0.5 * abs(sigma2 - sigma3);
    Max_TauY = 0.5 * abs(sigma1 - sigma3); Min_TauY = -0.5 * abs(sigma1 - sigma3);
    Max_TauZ = 0.5 * abs(sigma1 - sigma2); Min_TauZ = -0.5 * abs(sigma1 - sigma2);

    Ss_element(i, 1) = i;
    Ss_element(i, 2) = Max_TauX; Ss_element(i, 3) = Min_TauX;
    Ss_element(i, 4) = Max_TauY; Ss_element(i, 5) = Min_TauY;
    Ss_element(i, 6) = Max_TauZ; Ss_element(i, 7) = Min_TauZ;

end

%% result visualization----------------------------------------------------
scalingFactor = 500;

% figure('Name','Deformation contour');
% movegui('southeast');
% hold on; grid on; box on;
% axis equal;
% view(-45,30);xlim([-0.1,2.1]);ylim([-0.8,0.8]);zlim([-0.5,0.5]);
% xlabel('x (m)');
% ylabel('y (m)');
% zlabel('z (m)');
% for i = 1:size(Element, 1)
%     x_e = zeros(1,NEN);
%     y_e = zeros(1,NEN);
%     z_e = zeros(1,NEN);
%     x_eU = zeros(1,NEN);
%     y_eU = zeros(1,NEN);
%     z_eU = zeros(1,NEN);
%     for j = 1:NEN
%         x_e(j) = Node(Element(i,j+1),2);
%         y_e(j) = Node(Element(i,j+1),3);
%         z_e(j) = Node(Element(i,j+1),4);
%         x_eU(j) = Node(Element(i,j+1),2) + scalingFactor * U_deformed(Element(i,j+1), 2);
%         y_eU(j) = Node(Element(i,j+1),3) + scalingFactor * U_deformed(Element(i,j+1), 3);
%         z_eU(j) = Node(Element(i,j+1),4) + scalingFactor * U_deformed(Element(i,j+1), 4);
%     end
%     for j = 1:size(patch_queue,1)
%         fill3(x_e(patch_queue(j,:)),y_e(patch_queue(j,:)),z_e(patch_queue(j,:)),...
%             'c','FaceAlpha',0.5);
%         fill3(x_eU(patch_queue(j,:)),y_eU(patch_queue(j,:)),z_eU(patch_queue(j,:)),...
%             'r','FaceAlpha',0.5);
%     end
% end

% % plot principle stress on specify path------------------------------------
figure('Name', 'Principle stress in path');
movegui('northeast');
hold on;
x_path = linspace(0,2.00579,100);
y_path = 0.0505*x_path.^2+0.0085*x_path-0.04;
z_path = 0.0018*x_path.^2-0.0104*x_path+0.015;
X_location = [0 -0.04 0.015;
              2.00579 0.180220487967050 0.00138153234338000];
distance = zeros(size(Element, 1), 10);
store_ElementIndex = zeros(1, size(x_path, 2));
store_ElementIndex_Converge = zeros(1, 2);
S_gX = zeros(1, size(x_path, 2));
S_gX_Converge = zeros(1, 2);
%%Path
for i = 1:size(x_path, 2)
    for j = 1:size(Element, 1)
        x_e = zeros(1,NEN);
        y_e = zeros(1,NEN);
        z_e = zeros(1,NEN);
        for k = 1:NEN
            x_e(k) = Node(Element(j,k+1),2);
            y_e(k) = Node(Element(j,k+1),3);
            z_e(k) = Node(Element(j,k+1),4);
        end

        x_eSort = sort(x_e);
        y_eSort = sort(y_e);
        z_eSort = sort(z_e);

        if x_path(i) >= x_eSort(1) && x_path(i) <= x_eSort(8) && y_path(i) >= y_eSort(1) && y_path(i) <= y_eSort(8) && z_path(i) >= z_eSort(1) && z_path(i) <= z_eSort(8)
            store_ElementIndex(i) = j;
        end

    end
end

for i = 1:size(store_ElementIndex, 2)
    X = zeros(8, 4);
    targetX = [x_path(i) y_path(i) z_path(i)];
    for j = 1:NEN
        X(j, 1) = Node(Element(store_ElementIndex(i),j+1),1);
        X(j, 2) = Node(Element(store_ElementIndex(i),j+1),2);
        X(j, 3) = Node(Element(store_ElementIndex(i),j+1),3);
        X(j, 4) = Node(Element(store_ElementIndex(i),j+1),4);
    end

    [natural_coordinate, err] = NewtonMethod(targetX, X);

    for j = 1:8
        localDisplacement(3*j-2, 1) = U_deformed(X(j, 1), 2);
        localDisplacement(3*j-1, 1) = U_deformed(X(j, 1), 3);
        localDisplacement(3*j, 1) = U_deformed(X(j, 1), 4);
    end
    B_ = B_matrix(natural_coordinate(1), natural_coordinate(2), natural_coordinate(3), X);
    strain = B_ * localDisplacement;
    stress = D_matrix * strain;
    %%principle stress-----------------------------------------------------
    stress_matrix = [stress(1) stress(4) stress(5);
                     stress(4) stress(2) stress(6);
                     stress(5) stress(6) stress(3)];
    eStress_matrix = eig(stress_matrix);
    sigma1 = max(eStress_matrix);
    sigma3 = min(eStress_matrix);
    sigma2 = stress_matrix(1, 1) + stress_matrix(2, 2) + stress_matrix(3, 3) - sigma1 - sigma3;
    S_gX(i) = sigma3;

end
plot(x_path, S_gX);
% 
% %%Converge-----------------------------------------------------------------
figure('Name', 'Principle stress in select points');
movegui('east');
for i = 1:size(X_location, 1)
    for j = 1:size(Element, 1)
        x_e = zeros(1,NEN);
        y_e = zeros(1,NEN);
        z_e = zeros(1,NEN);
        for k = 1:NEN
            x_e(k) = Node(Element(j,k+1),2);
            y_e(k) = Node(Element(j,k+1),3);
            z_e(k) = Node(Element(j,k+1),4);
        end
        x_eSort = sort(x_e);
        y_eSort = sort(y_e);
        z_eSort = sort(z_e);

        if X_location(i, 1) >= x_eSort(1) && X_location(i, 1) <= x_eSort(8) && X_location(i, 2) >= y_eSort(1) && X_location(i, 2) <= y_eSort(8) && X_location(i, 3) >= z_eSort(1) && X_location(i, 3) <= z_eSort(8)
            store_ElementIndex_Converge(i) = j;
        end

    end
end

for i = 1:size(X_location, 1)
    X = zeros(8, 4);
    targetX = [X_location(i, 1) X_location(i, 2) X_location(i, 3)];
    for j = 1:NEN
        X(j, 1) = Node(Element(store_ElementIndex_Converge(i),j+1),1);
        X(j, 2) = Node(Element(store_ElementIndex_Converge(i),j+1),2);
        X(j, 3) = Node(Element(store_ElementIndex_Converge(i),j+1),3);
        X(j, 4) = Node(Element(store_ElementIndex_Converge(i),j+1),4);
    end
    [natural_coordinate,err] = NewtonMethod(targetX, X);

    for j = 1:8
        localDisplacement(3*j-2, 1) = U_deformed(X(j, 1), 2);
        localDisplacement(3*j-1, 1) = U_deformed(X(j, 1), 3);
        localDisplacement(3*j, 1) = U_deformed(X(j, 1), 4);
    end
    B_ = B_matrix(natural_coordinate(1), natural_coordinate(2), natural_coordinate(3), X);
    strain = B_ * localDisplacement;
    stress = D_matrix * strain;
    %%principle stress-----------------------------------------------------
    stress_matrix = [stress(1) stress(4) stress(5);
                     stress(4) stress(2) stress(6);
                     stress(5) stress(6) stress(3)];
    eStress_matrix = eig(stress_matrix);
    sigma1 = max(eStress_matrix);
    sigma3 = min(eStress_matrix);
    sigma2 = stress_matrix(1, 1) + stress_matrix(2, 2) + stress_matrix(3, 3) - sigma1 - sigma3;
    scatter(X_location(i, 1), sigma1, 100, 'filled');
    hold on;grid on;
end
% 
% 
% %%plot element that pass the line------------------------------------------
figure('Name','element that pass the line ');
hold on; grid on; box on;
axis equal;
view(-45,30);xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
plot3(x_path,y_path,z_path,'LineWidth',3,'Color','k');
plot3(X_location(1, 1), X_location(1, 2), X_location(1, 3),'o');
plot3(X_location(2, 1), X_location(2, 2), X_location(2, 3),'o');
for i = 1:size(store_ElementIndex, 2)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);

    for j = 1:NEN
        x_e(j) = Node(Element(store_ElementIndex(i),j+1),2);
        y_e(j) = Node(Element(store_ElementIndex(i),j+1),3);
        z_e(j) = Node(Element(store_ElementIndex(i),j+1),4);
    end
    for j = 1:size(patch_queue,1)
        fill3(x_e(patch_queue(j,:)),y_e(patch_queue(j,:)),z_e(patch_queue(j,:)),...
            'c','FaceAlpha',0.5);
    end
end
for i = 1:size(store_ElementIndex_Converge, 2)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);

    for j = 1:NEN
        x_e(j) = Node(Element(store_ElementIndex_Converge(i),j+1),2);
        y_e(j) = Node(Element(store_ElementIndex_Converge(i),j+1),3);
        z_e(j) = Node(Element(store_ElementIndex_Converge(i),j+1),4);
    end
    for j = 1:size(patch_queue,1)
        fill3(x_e(patch_queue(j,:)),y_e(patch_queue(j,:)),z_e(patch_queue(j,:)),...
            'r','FaceAlpha',0.5);
    end
end

%Plot 3 DOF deformation contour plot--------------------------------------
figure('Name','ux');
movegui('northwest');
ax1 = nexttile;title(ax1,'ux');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','uy');
movegui('west');
ax2 = nexttile;title(ax2,'uy');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','uz');
movegui('southwest');
ax3 = nexttile;title(ax3,'uz');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','u magnitude');
movegui('southwest');
ax4 = nexttile;title(ax4,'u magnitude');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    x_eU = zeros(1,NEN);
    y_eU = zeros(1,NEN);
    z_eU = zeros(1,NEN);
    U_magnitude = zeros(1, NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
        x_eU(j) = U_deformed(Element(i,j+1), 2);
        y_eU(j) = U_deformed(Element(i,j+1), 3);
        z_eU(j) = U_deformed(Element(i,j+1), 4);
        U_magnitude(j) = U_deformed(Element(i,j+1), 5);
    end
    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [x_eU(1),x_eU(2),x_eU(3),x_eU(4),x_eU(5),x_eU(6),x_eU(7),x_eU(8)];
    c2 = [y_eU(1),y_eU(2),y_eU(3),y_eU(4),y_eU(5),y_eU(6),y_eU(7),y_eU(8)];
    c3 = [z_eU(1),z_eU(2),z_eU(3),z_eU(4),z_eU(5),z_eU(6),z_eU(7),z_eU(8)];
    c4 = [U_magnitude(1),U_magnitude(2),U_magnitude(3),U_magnitude(4),U_magnitude(5),U_magnitude(6),U_magnitude(7),U_magnitude(8)];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));

        fill3(ax2, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c2(patch_queue(j,:)));

        fill3(ax3, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c3(patch_queue(j,:)));

        fill3(ax4, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c4(patch_queue(j,:)));
    end


end

%%Plot 6 directions strian contour plot------------------------------------
figure('Name','strainXX');
movegui([200 670]);
ax1 = nexttile;title(ax1,'strainXX');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','strainYY');
movegui([200 350]);
ax2 = nexttile;title(ax2,'strainYY');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressZZ');
movegui([200 10]);
ax3 = nexttile;title(ax3,'strainZZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','strainXY');
movegui([300 670]);
ax4 = nexttile;title(ax4,'strainXY');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','strainXZ');
movegui([300 350]);
ax5 = nexttile;title(ax5,'strainXZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','strainYZ');
movegui([300 10]);
ax6 = nexttile;title(ax6,'strainYZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    x_eU = zeros(1,NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
    end
    strainXX = E_element(i, 2);
    strainYY = E_element(i, 3);
    strainZZ = E_element(i, 4);
    strainXY = E_element(i, 5);
    strainXZ = E_element(i, 6);
    strainYZ = E_element(i, 7);

    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [strainXX,strainXX,strainXX,strainXX,strainXX,strainXX,strainXX,strainXX];
    c2 = [strainYY,strainYY,strainYY,strainYY,strainYY,strainYY,strainYY,strainYY];
    c3 = [strainZZ,strainZZ,strainZZ,strainZZ,strainZZ,strainZZ,strainZZ,strainZZ];
    c4 = [strainXY,strainXY,strainXY,strainXY,strainXY,strainXY,strainXY,strainXY];
    c5 = [strainXZ,strainXZ,strainXZ,strainXZ,strainXZ,strainXZ,strainXZ,strainXZ];
    c6 = [strainYZ,strainYZ,strainYZ,strainYZ,strainYZ,strainYZ,strainYZ,strainYZ];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));

        fill3(ax2, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c2(patch_queue(j,:)));

        fill3(ax3, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c3(patch_queue(j,:)));

        fill3(ax4, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c4(patch_queue(j,:)));

        fill3(ax5, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c5(patch_queue(j,:)));

        fill3(ax6, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c6(patch_queue(j,:)));
    end
end

%%Plot 6 directions stress contour plot--------------------------------
figure('Name','stressXX');
movegui([400 670]);
ax1 = nexttile;title(ax1,'stressXX');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressYY');
movegui([400 350]);
ax2 = nexttile;title(ax2,'stressYY');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressZZ');
movegui([400 10]);
ax3 = nexttile;title(ax3,'stressZZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressXY');
ax4 = nexttile;title(ax4,'stressXY');
movegui([500 670]);
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressXZ');
movegui([500 350]);
ax5 = nexttile;title(ax5,'stressXZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;
figure('Name','stressYZ');
movegui([500 10]);
ax6 = nexttile;title(ax6,'stressYZ');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
    end
    stressXX = S_element(i, 2);
    stressYY = S_element(i, 3);
    stressZZ = S_element(i, 4);
    stressXY = S_element(i, 5);
    stressXZ = S_element(i, 6);
    stressYZ = S_element(i, 7);

    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [stressXX,stressXX,stressXX,stressXX,stressXX,stressXX,stressXX,stressXX];
    c2 = [stressYY,stressYY,stressYY,stressYY,stressYY,stressYY,stressYY,stressYY];
    c3 = [stressZZ,stressZZ,stressZZ,stressZZ,stressZZ,stressZZ,stressZZ,stressZZ];
    c4 = [stressXY,stressXY,stressXY,stressXY,stressXY,stressXY,stressXY,stressXY];
    c5 = [stressXZ,stressXZ,stressXZ,stressXZ,stressXZ,stressXZ,stressXZ,stressXZ];
    c6 = [stressYZ,stressYZ,stressYZ,stressYZ,stressYZ,stressYZ,stressYZ,stressYZ];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));

        fill3(ax2, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c2(patch_queue(j,:)));

        fill3(ax3, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c3(patch_queue(j,:)));

        fill3(ax4, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c4(patch_queue(j,:)));

        fill3(ax5, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c5(patch_queue(j,:)));

        fill3(ax6, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c6(patch_queue(j,:)));
    end

end

% %%Principle stress contour-------------------------------------------------
figure('Name','principle stress1');
movegui([600 670]);
ax1 = nexttile;title(ax1,'principle stress1');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','principle stress2');
movegui([600 350]);
ax2 = nexttile;title(ax2,'principle stress2');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','principle stress3');
movegui([600 10]);
ax3 = nexttile;title(ax3,'principle stress3');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
    end
    principleStress1 = Sp_element(i, 2);
    principleStress2 = Sp_element(i, 3);
    principleStress3 = Sp_element(i, 4);
    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [principleStress1,principleStress1,principleStress1,principleStress1,principleStress1,principleStress1,principleStress1,principleStress1];
    c2 = [principleStress2,principleStress2,principleStress2,principleStress2,principleStress2,principleStress2,principleStress2,principleStress2];
    c3 = [principleStress3,principleStress3,principleStress3,principleStress3,principleStress3,principleStress3,principleStress3,principleStress3];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));

        fill3(ax2, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c2(patch_queue(j,:)));

        fill3(ax3, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c3(patch_queue(j,:)));
    end
end

% %%Von-Mises stress contour-------------------------------------------------
figure('Name','Von-Mises stress');
movegui([700 670]);
ax1 = nexttile;title(ax1,'Von-Mises stress');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
    end
    Von_Mises = Sv_element(i, 2);
    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [Von_Mises,Von_Mises,Von_Mises,Von_Mises,Von_Mises,Von_Mises,Von_Mises,Von_Mises];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));
    end
end

% %%Shear stress contour-----------------------------------------------------
figure('Name','Shear stressX_Max');
movegui([800 670]);
ax1 = nexttile;title(ax1,'Shear stressX_M_a_x');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','Shear stressX_Min');
movegui([800 350]);
ax2 = nexttile;title(ax2,'Shear stressX_M_i_n');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','Shear stressY_Max');
movegui([800 10]);
ax3 = nexttile;title(ax3,'Shear stressY_M_a_x');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','Shear stressY_Min');
movegui([900 670]);
ax4 = nexttile;title(ax4,'Shear stressY_M_i_n');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','Shear stressZ_Max');
movegui([900 350]);
ax5 = nexttile;title(ax5,'Shear stressZ_M_a_x');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

figure('Name','Shear stressZ_Min');
movegui([900 10]);
ax6 = nexttile;title(ax6,'Shear stressZ_M_i_n');
hold on; box on; grid on;axis equal;
view(-45,30);xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');xlim([-0.1,2.1]);ylim([-0.25,0.25]);zlim([-0.25,0.25]);
colorbar;colormap jet;

for i = 1:size(Element, 1)
    x_e = zeros(1,NEN);
    y_e = zeros(1,NEN);
    z_e = zeros(1,NEN);
    for j = 1:NEN
        x_e(j) = Node(Element(i,j+1),2);
        y_e(j) = Node(Element(i,j+1),3);
        z_e(j) = Node(Element(i,j+1),4);
    end
    MaxSsX = Ss_element(i, 2);
    MinSsx = Ss_element(i, 3);
    MaxSsY = Ss_element(i, 4);
    MinSsY = Ss_element(i, 5);
    MaxSsZ = Ss_element(i, 6);
    MinSsZ = Ss_element(i, 7);

    x = [x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),x_e(6),x_e(7),x_e(8)];
    y = [y_e(1),y_e(2),y_e(3),y_e(4),y_e(5),y_e(6),y_e(7),y_e(8)];
    z = [z_e(1),z_e(2),z_e(3),z_e(4),z_e(5),z_e(6),z_e(7),z_e(8)];
    c1 = [MaxSsX,MaxSsX,MaxSsX,MaxSsX,MaxSsX,MaxSsX,MaxSsX,MaxSsX];
    c2 = [MinSsx,MinSsx,MinSsx,MinSsx,MinSsx,MinSsx,MinSsx,MinSsx];
    c3 = [MaxSsY,MaxSsY,MaxSsY,MaxSsY,MaxSsY,MaxSsY,MaxSsY,MaxSsY];
    c4 = [MinSsY,MinSsY,MinSsY,MinSsY,MinSsY,MinSsY,MinSsY,MinSsY];
    c5 = [MaxSsZ,MaxSsZ,MaxSsZ,MaxSsZ,MaxSsZ,MaxSsZ,MaxSsZ,MaxSsZ];
    c6 = [MinSsZ,MinSsZ,MinSsZ,MinSsZ,MinSsZ,MinSsZ,MinSsZ,MinSsZ];

    for j = 1:size(patch_queue,1)
        fill3(ax1, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c1(patch_queue(j,:)));

        fill3(ax2, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c2(patch_queue(j,:)));

        fill3(ax3, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c3(patch_queue(j,:)));

        fill3(ax4, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c4(patch_queue(j,:)));

        fill3(ax5, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c5(patch_queue(j,:)));

        fill3(ax6, x(patch_queue(j,:)),y(patch_queue(j,:)),z(patch_queue(j,:)),...
        c6(patch_queue(j,:)));
    end
end