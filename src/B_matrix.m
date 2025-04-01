function B = B_matrix(ksi, eta, zeta, X)
    x = X(:, 2);
    y = X(:, 3);
    z = X(:, 4);
    %%X
    x11 = x(1) + x(2) - x(3) - x(4) + x(5) + x(6) - x(7) - x(8);
    x12 = -x(1) + x(2) - x(3) + x(4) - x(5) + x(6) - x(7) + x(8);
    x13 = -x(1) - x(2) + x(3) + x(4) + x(5) + x(6) - x(7) - x(8);
    x21 = -x(1) + x(2) + x(3) - x(4) - x(5) + x(6) + x(7) - x(8);
    x22 = -x(1) + x(2) - x(3) + x(4) - x(5) + x(6) - x(7) + x(8);
    x23 = x(1) - x(2) - x(3) + x(4) - x(5) + x(6) + x(7) - x(8);
    x31 = -x(1) - x(2) - x(3) - x(4) + x(5) + x(6) + x(7) + x(8);
    x32 = -x(1) - x(2) + x(3) + x(4) + x(5) + x(6) - x(7) - x(8);
    x33 = x(1) - x(2) - x(3) + x(4) - x(5) + x(6) + x(7) - x(8);
    x0 = x(1) - x(2) + x(3) - x(4) - x(5) + x(6) - x(7) + x(8);

    %%Y
    y11 = y(1) + y(2) - y(3) - y(4) + y(5) + y(6) - y(7) - y(8);
    y12 = -y(1) + y(2) - y(3) + y(4) - y(5) + y(6) - y(7) + y(8);
    y13 = -y(1) - y(2) + y(3) + y(4) + y(5) + y(6) - y(7) - y(8);
    y21 = -y(1) + y(2) + y(3) - y(4) - y(5) + y(6) + y(7) - y(8);
    y22 = -y(1) + y(2) - y(3) + y(4) - y(5) + y(6) - y(7) + y(8);
    y23 = y(1) - y(2) - y(3) + y(4) - y(5) + y(6) + y(7) - y(8);
    y31 = -y(1) - y(2) - y(3) - y(4) + y(5) + y(6) + y(7) + y(8);
    y32 = -y(1) - y(2) + y(3) + y(4) + y(5) + y(6) - y(7) - y(8);
    y33 = y(1) - y(2) - y(3) + y(4) - y(5) + y(6) + y(7) - y(8);
    y0 = y(1) - y(2) + y(3) - y(4) - y(5) + y(6) - y(7) + y(8);

    %Z
    z11 = z(1) + z(2) - z(3) - z(4) + z(5) + z(6) - z(7) - z(8);
    z12 = -z(1) + z(2) - z(3) + z(4) - z(5) + z(6) - z(7) + z(8);
    z13 = -z(1) - z(2) + z(3) + z(4) + z(5) + z(6) - z(7) - z(8);
    z21 = -z(1) + z(2) + z(3) - z(4) - z(5) + z(6) + z(7) - z(8);
    z22 = -z(1) + z(2) - z(3) + z(4) - z(5) + z(6) - z(7) + z(8);
    z23 = z(1) - z(2) - z(3) + z(4) - z(5) + z(6) + z(7) - z(8);
    z31 = -z(1) - z(2) - z(3) - z(4) + z(5) + z(6) + z(7) + z(8);
    z32 = -z(1) - z(2) + z(3) + z(4) + z(5) + z(6) - z(7) - z(8);
    z33 = z(1) - z(2) - z(3) + z(4) - z(5) + z(6) + z(7) - z(8);
    z0 = z(1) - z(2) + z(3) - z(4) - z(5) + z(6) - z(7) + z(8);

    %%[dN_dksi; dN_deta; dN_dzeta] matrix
    dN1_matrix = [0.125 * (1-eta) * (1-zeta);-0.125 * (1+ksi) * (1-zeta);-0.125 * (1+ksi) * (1-eta)];
    dN2_matrix = [0.125 * (1+eta) * (1-zeta);0.125 * (1+ksi) * (1-zeta);-0.125 * (1+ksi) * (1+eta)];
    dN3_matrix = [-0.125 * (1+eta) * (1-zeta);0.125 * (1-ksi) * (1-zeta);-0.125 * (1-ksi) * (1+eta)];
    dN4_matrix = [-0.125 * (1-eta) * (1-zeta);-0.125 * (1-ksi) * (1-zeta);-0.125 * (1-ksi) * (1-eta)];
    dN5_matrix = [0.125 * (1-eta) * (1+zeta);-0.125 * (1+ksi) * (1+zeta);0.125 * (1+ksi) * (1-eta)];
    dN6_matrix = [0.125 * (1+eta) * (1+zeta);0.125 * (1+ksi) * (1+zeta);0.125 * (1+ksi) * (1+eta)];            
    dN7_matrix = [-0.125 * (1+eta) * (1+zeta);0.125 * (1-ksi) * (1+zeta);0.125 * (1-ksi) * (1+eta)];
    dN8_matrix = [-0.125 * (1-eta) * (1+zeta);-0.125 * (1-ksi) * (1+zeta);0.125 * (1-ksi) * (1-eta)];
    
    dx_dksi = 0.125 * x11 + 0.125 * eta * x12 + 0.125 * zeta * x13 + 0.125 * eta * zeta * x0; 
    dy_dksi = 0.125 * y11 + 0.125 * eta * y12 + 0.125 * zeta * y13 + 0.125 * eta * zeta * y0;
    dz_dksi = 0.125 * z11 + 0.125 * eta * z12 + 0.125 * zeta * z13 + 0.125 * eta * zeta * z0; 
    
    dx_deta = 0.125 * x21 + 0.125 * ksi * x22 + 0.125 * zeta * x23 + 0.125 * ksi * zeta * x0;
    dy_deta = 0.125 * y21 + 0.125 * ksi * y22 + 0.125 * zeta * y23 + 0.125 * ksi * zeta * y0;
    dz_deta = 0.125 * z21 + 0.125 * ksi * z22 + 0.125 * zeta * z23 + 0.125 * ksi * zeta * z0;
    
    dx_dzeta = 0.125 * x31 + 0.125 * ksi * x32 + 0.125 * eta * x33 + 0.125 * ksi * eta * x0;
    dy_dzeta = 0.125 * y31 + 0.125 * ksi * y32 + 0.125 * eta * y33 + 0.125 * ksi * eta * y0;
    dz_dzeta = 0.125 * z31 + 0.125 * ksi * z32 + 0.125 * eta * z33 + 0.125 * ksi * eta * z0;
    
    Jacobi_matrix = [dx_dksi dy_dksi dz_dksi;
                     dx_deta dy_deta dz_deta;
                     dx_dzeta dy_dzeta dz_dzeta];

    
    %invJacobi_matrix = inv(Jacobi_matrix);

    B1 = Jacobi_matrix \ dN1_matrix;
    B2 = Jacobi_matrix \ dN2_matrix;
    B3 = Jacobi_matrix \ dN3_matrix;
    B4 = Jacobi_matrix \ dN4_matrix;
    B5 = Jacobi_matrix \ dN5_matrix;
    B6 = Jacobi_matrix \ dN6_matrix;
    B7 = Jacobi_matrix \ dN7_matrix;
    B8 = Jacobi_matrix \ dN8_matrix;
    
    %%B1
    B(1, 1) = B1(1, 1);%%dN/dx
    B(2, 2) = B1(2, 1);%%dN/dy
    B(3, 3) = B1(3, 1);%%dN/dz
    B(4, 1) = B1(2, 1);%%dN/dy
    B(4, 2) = B1(1, 1);%%dN/dx
    B(5, 1) = B1(3, 1);%%dN/dz
    B(5, 3) = B1(1, 1);%%dN/dx
    B(6, 2) = B1(3, 1);%%dN/dz
    B(6, 3) = B1(2, 1);%%dN/dy

    %%B2
    B(1, 4) = B2(1, 1);%%dN/dx
    B(2, 5) = B2(2, 1);%%dN/dy
    B(3, 6) = B2(3, 1);%%dN/dz
    B(4, 4) = B2(2, 1);%%dN/dy
    B(4, 5) = B2(1, 1);%%dN/dx
    B(5, 4) = B2(3, 1);%%dN/dz
    B(5, 6) = B2(1, 1);%%dN/dx
    B(6, 5) = B2(3, 1);%%dN/dz
    B(6, 6) = B2(2, 1);%%dN/dy

    %%B3
    B(1, 7) = B3(1, 1);%%dN/dx
    B(2, 8) = B3(2, 1);%%dN/dy
    B(3, 9) = B3(3, 1);%%dN/dz
    B(4, 7) = B3(2, 1);%%dN/dy
    B(4, 8) = B3(1, 1);%%dN/dx
    B(5, 7) = B3(3, 1);%%dN/dz
    B(5, 9) = B3(1, 1);%%dN/dx
    B(6, 8) = B3(3, 1);%%dN/dz
    B(6, 9) = B3(2, 1);%%dN/dy

    %%B4
    B(1, 10) = B4(1, 1);%%dN/dx
    B(2, 11) = B4(2, 1);%%dN/dy
    B(3, 12) = B4(3, 1);%%dN/dz
    B(4, 10) = B4(2, 1);%%dN/dy
    B(4, 11) = B4(1, 1);%%dN/dx
    B(5, 10) = B4(3, 1);%%dN/dz
    B(5, 12) = B4(1, 1);%%dN/dx
    B(6, 11) = B4(3, 1);%%dN/dz
    B(6, 12) = B4(2, 1);%%dN/dy
    
    %%B5
    B(1, 13) = B5(1, 1);%%dN/dx
    B(2, 14) = B5(2, 1);%%dN/dy
    B(3, 15) = B5(3, 1);%%dN/dz
    B(4, 13) = B5(2, 1);%%dN/dy
    B(4, 14) = B5(1, 1);%%dN/dx
    B(5, 13) = B5(3, 1);%%dN/dz
    B(5, 15) = B5(1, 1);%%dN/dx
    B(6, 14) = B5(3, 1);%%dN/dz
    B(6, 15) = B5(2, 1);%%dN/dy

    %%B6
    B(1, 16) = B6(1, 1);%%dN/dx
    B(2, 17) = B6(2, 1);%%dN/dy
    B(3, 18) = B6(3, 1);%%dN/dz
    B(4, 16) = B6(2, 1);%%dN/dy
    B(4, 17) = B6(1, 1);%%dN/dx
    B(5, 16) = B6(3, 1);%%dN/dz
    B(5, 18) = B6(1, 1);%%dN/dx
    B(6, 17) = B6(3, 1);%%dN/dz
    B(6, 18) = B6(2, 1);%%dN/dy

    %%B7
    B(1, 19) = B7(1, 1);%%dN/dx
    B(2, 20) = B7(2, 1);%%dN/dy
    B(3, 21) = B7(3, 1);%%dN/dz
    B(4, 19) = B7(2, 1);%%dN/dy
    B(4, 20) = B7(1, 1);%%dN/dx
    B(5, 19) = B7(3, 1);%%dN/dz
    B(5, 21) = B7(1, 1);%%dN/dx
    B(6, 20) = B7(3, 1);%%dN/dz
    B(6, 21) = B7(2, 1);%%dN/dy

    %%B8
    B(1, 22) = B8(1, 1);%%dN/dx
    B(2, 23) = B8(2, 1);%%dN/dy
    B(3, 24) = B8(3, 1);%%dN/dz
    B(4, 22) = B8(2, 1);%%dN/dy
    B(4, 23) = B8(1, 1);%%dN/dx
    B(5, 22) = B8(3, 1);%%dN/dz
    B(5, 24) = B8(1, 1);%%dN/dx
    B(6, 23) = B8(3, 1);%%dN/dz
    B(6, 24) = B8(2, 1);%%dN/dy

end