clc; clear all;
syms ksi eta zeta

dN1_matrix = [0.125 * (1-eta) * (1-zeta);
             -0.125 * (1+ksi) * (1-zeta);
             -0.125 * (1+ksi) * (1-eta)];

dN2_matrix = [0.125 * (1+eta) * (1-zeta);
             0.125 * (1+ksi) * (1-zeta);
             -0.125 * (1+ksi) * (1+eta)];

dN3_matrix = [-0.125 * (1+eta) * (1-zeta);
             0.125 * (1-ksi) * (1-zeta);
             -0.125 * (1-ksi) * (1+eta)];

dN4_matrix = [-0.125 * (1-eta) * (1-zeta);
             -0.125 * (1-ksi) * (1-zeta);
             -0.125 * (1-ksi) * (1-eta)];

dN5_matrix = [0.125 * (1-eta) * (1+zeta);
             -0.125 * (1+ksi) * (1+zeta);
             0.125 * (1+ksi) * (1-eta);];

dN6_matrix = [0.125 * (1+eta) * (1+zeta);
             0.125 * (1+ksi) * (1+zeta);
             0.125 * (1+ksi) * (1+eta)];

dN7_matrix = [-0.125 * (1+eta) * (1+zeta);
             0.125 * (1-ksi) * (1+zeta);
             0.125 * (1-ksi) * (1+eta)];

dN8_matrix = [-0.125 * (1-eta) * (1+zeta);
             -0.125 * (1-ksi) * (1+zeta);
             0.125 * (1-ksi) * (1-eta)];

Jacobi_matrix = [dx_dksi() dy_dksi() dz_dksi();
                  dx_deta() dy_deta() dz_deta();
                  dx_dzeta() dy_dzeta() dz_dzeta()];


Jv_exact = det(Jacobi_matrix);

[cof, det] = cofactor_Det(Jacobi_matrix);
Jv_ = det;
% B1 = inv(Inverse_Matrix) * dN1_matrix;
% B2 = Inverse_Matrix \ dN2_matrix;
% B3 = Inverse_Matrix \ dN3_matrix;
% B4 = Inverse_Matrix \ dN4_matrix;
% B5 = Inverse_Matrix \ dN5_matrix;
% B6 = Inverse_Matrix \ dN6_matrix;
% B7 = Inverse_Matrix \ dN7_matrix;
% B8 = Inverse_Matrix \ dN8_matrix;
% 
% A1 = [dx_dksi() dy_dksi() dz_dksi()]; 
% B1 = [dx_deta() dy_deta() dz_deta()]; 
% temp = cross(A1, B1);
% Js1 = sqrt(dot(temp, temp));
% 
% A3 = [dx_dzeta() dy_dzeta() dz_dzeta()]; 
% B3 = [dx_deta() dy_deta() dz_deta()]; 
% temp = cross(A3, B3);
% Js3 = sqrt(dot(temp, temp));
% 
% A4 = [dx_dksi() dy_dksi() dz_dksi()]; 
% B4 = [dx_dzeta() dy_dzeta() dz_dzeta()]; 
% temp = cross(A4, B4);
% Js4 = sqrt(dot(temp, temp));

function Ans = dx_dksi()
syms eta zeta x11 x12 x13 x0
    Ans = 0.125 * x11..., 
    + 0.125 * eta * x12..., 
    + 0.125 * zeta * x13...,
    + 0.125 * eta * zeta * x0;
end

function Ans = dy_dksi()
syms eta zeta y11 y12 y13 y0
    Ans = 0.125 * y11...,
    + 0.125 * eta * y12...,
    + 0.125 * zeta * y13...,
    + 0.125 * eta * zeta * y0;
end

function Ans = dz_dksi()
syms eta zeta z11 z12 z13 z0
    Ans = 0.125 * z11...,
    + 0.125 * eta * z12...,
    + 0.125 * zeta * z13...,
    + 0.125 * eta * zeta * z0;
end

function Ans = dx_deta()
syms ksi  zeta x21 x22 x23 x0
    Ans = 0.125 * x21...,
    + 0.125 * ksi * x22...,
    + 0.125 * zeta * x23...,
    + 0.125 * ksi * zeta * x0;
end

function Ans = dy_deta()
syms ksi zeta y21 y22 y23 y0
    Ans = 0.125 * y21...,
    + 0.125 * ksi * y22...,
    + 0.125 * zeta * y23...,
    + 0.125 * ksi * zeta * y0;
end

function Ans = dz_deta()
syms ksi zeta z21 z22 z23 z0
    Ans = 0.125 * z21...,
    + 0.125 * ksi * z22...,
    + 0.125 * zeta * z23...,
    + 0.125 * ksi * zeta * z0;
end

function Ans = dx_dzeta()
syms ksi eta x31 x32 x33 x0
    Ans = 0.125 * x31...,
    + 0.125 * ksi * x32...,
    + 0.125 * eta * x33...,
    + 0.125 * ksi * eta * x0;
end

function Ans = dy_dzeta()
syms ksi eta y31 y32 y33 y0
   Ans = 0.125 * y31...,
   + 0.125 * ksi * y32...,
   + 0.125 * eta * y33...,
   + 0.125 * ksi * eta * y0;
end

function Ans = dz_dzeta()
syms ksi eta z31 z32 z33 z0
    Ans = 0.125 * z31...,
    + 0.125 * ksi * z32...,
    + 0.125 * eta * z33...,
    + 0.125 * ksi * eta * z0;
end

function detM = detminor(M)
   
    n = size(M, 1);
    if n == 1
        detM = M(1,1);
    else
        detM = 0;
        for i = 1:n
            subM = M([1:i-1,i+1:n],2:n);
            detM = detM + ((-1)^(i+1)) * M(i,1) * detminor(subM);
        end
    end
end

function [cof, det] = cofactor_Det(M)
    % Function to calculate cofactor matrix and determinant of matrix M
    
    n = size(M, 1); % Assuming M is a square matrix
    cof = zeros(n);
    det = 0;
    
    for i = 1:n
        for j = 1:n
            % Calculate cofactor of M(i,j)
            subM = M([1:i-1,i+1:n],[1:j-1,j+1:n]);
            cof(i,j) = ((-1)^(i+j)) * detminor(subM);
        end
    end
    det = sum(M(1,:) .* cof(1,:));
end
