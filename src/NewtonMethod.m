function [Ans, ERR] = NewtonMethod(targerX, X)

    %%initial guass value and setting
    targetErr = 1e-10;
    At = zeros(3, 1);
    %%
    xt = targerX(1);
    yt = targerX(2);
    zt = targerX(3);
    %%
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
    
    err = 10;
    while err > targetErr
        ksi = At(1);
        eta = At(2);
        zeta = At(3);

        %%N1
        N1 = 0.125 * (1 + ksi) * (1 - eta) * (1 - zeta);
        %%N2
        N2 = 0.125 * (1 + ksi) * (1 + eta) * (1 - zeta);
        %%N3
        N3 = 0.125 * (1 - ksi) * (1 + eta) * (1 - zeta);
        %%N4
        N4 = 0.125 * (1 - ksi) * (1 - eta) * (1 - zeta);  
        %%N5
        N5 = 0.125 * (1 + ksi) * (1 - eta) * (1 + zeta);
        %%N6
        N6 = 0.125 * (1 + ksi) * (1 + eta) * (1 + zeta);
        %%N7
        N7 = 0.125 * (1 - ksi) * (1 + eta) * (1 + zeta);
        %%N8
        N8 = 0.125 * (1 - ksi) * (1 - eta) * (1 + zeta);

        %%Jacobi matrix
        dx_dksi = 0.125 * x11 + 0.125 * eta * x12 + 0.125 * zeta * x13 + 0.125 * eta * zeta * x0; 
        dy_dksi = 0.125 * y11 + 0.125 * eta * y12 + 0.125 * zeta * y13 + 0.125 * eta * zeta * y0;
        dz_dksi = 0.125 * z11 + 0.125 * eta * z12 + 0.125 * zeta * z13 + 0.125 * eta * zeta * z0; 
        
        dx_deta = 0.125 * x21 + 0.125 * ksi * x22 + 0.125 * zeta * x23 + 0.125 * ksi * zeta * x0;
        dy_deta = 0.125 * y21 + 0.125 * ksi * y22 + 0.125 * zeta * y23 + 0.125 * ksi * zeta * y0;
        dz_deta = 0.125 * z21 + 0.125 * ksi * z22 + 0.125 * zeta * z23 + 0.125 * ksi * zeta * z0;
        
        dx_dzeta = 0.125 * x31 + 0.125 * ksi * x32 + 0.125 * eta * x33 + 0.125 * ksi * eta * x0;
        dy_dzeta = 0.125 * y31 + 0.125 * ksi * y32 + 0.125 * eta * y33 + 0.125 * ksi * eta * y0;
        dz_dzeta = 0.125 * z31 + 0.125 * ksi * z32 + 0.125 * eta * z33 + 0.125 * ksi * eta * z0;

        Jacobi_matrix = [-dx_dksi -dx_deta -dx_dzeta;
                         -dy_dksi -dy_deta -dy_dzeta;
                         -dz_dksi -dz_deta -dz_dzeta];

        f1 = xt - (N1 * x(1) + N2 * x(2) + N3 * x(3) + N4 * x(4) + N5 * x(5) + N6 * x(6) + N7 * x(7) + N8 * x(8));
        f2 = yt - (N1 * y(1) + N2 * y(2) + N3 * y(3) + N4 * y(4) + N5 * y(5) + N6 * y(6) + N7 * y(7) + N8 * y(8));
        f3 = zt - (N1 * z(1) + N2 * z(2) + N3 * z(3) + N4 * z(4) + N5 * z(5) + N6 * z(6) + N7 * z(7) + N8 * z(8));
        F = [f1; f2; f3];

        tmp = Jacobi_matrix\F;

        At_1 = At - tmp;
        err = mean(abs(At_1 - At));
        

        At = At_1;

    end
    
    if At(1) >1
        At(1) = 1;
    end

    if At(2) >1
        At(2) = 1;
    end

    if At(3) >1
        At(3) = 1;
    end
    
    if At(1) < -1
        At(1) = -1;
    end

    if At(2) < -1
        At(2) = -1;
    end

    if At(3) < -1
        At(3) = -1;
    end
    
    Ans = At;
    ERR = err;
end