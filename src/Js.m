function Ans = Js(ksi, eta, zeta, X, faceIndex)
    x = X(:, 2);
    y = X(:, 3);
    z = X(:, 4);
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

    if faceIndex == 1
        dx_dksi = 0.25 * (1-eta)*x(1) + 0.25 * (1+eta)*x(2) - 0.25*(1+eta)*x(3) - 0.25*(1-eta)*x(4);
        dy_dksi = 0.25 * (1-eta)*y(1) + 0.25 * (1+eta)*y(2) - 0.25*(1+eta)*y(3) - 0.25*(1-eta)*y(4);
        dz_dksi = 0.25 * (1-eta)*z(1) + 0.25 * (1+eta)*z(2) - 0.25*(1+eta)*z(3) - 0.25*(1-eta)*z(4);
        
        dx_deta = -0.25 * (1+ksi)*x(1) + 0.25 * (1+ksi)*x(2) + 0.25*(1-ksi)*x(3) - 0.25*(1-ksi)*x(4);
        dy_deta = -0.25 * (1+ksi)*y(1) + 0.25 * (1+ksi)*y(2) + 0.25*(1-ksi)*y(3) - 0.25*(1-ksi)*y(4);
        dz_deta = -0.25 * (1+ksi)*z(1) + 0.25 * (1+ksi)*z(2) + 0.25*(1-ksi)*z(3) - 0.25*(1-ksi)*z(4);
        A_vector = [dx_dksi dy_dksi dz_dksi];
        B_vector = [dx_deta dy_deta dz_deta];
        Ans = norm(cross(A_vector, B_vector));
    end

    if faceIndex == 2
        dx_dksi = 0.25 * (1-eta)*x(5) + 0.25*(1+eta)*x(6) - 0.25*(1+eta)*x(7) - 0.25*(1-eta)*x(8);
        dy_dksi = 0.25 * (1-eta)*y(5) + 0.25*(1+eta)*y(6) - 0.25*(1+eta)*y(7) - 0.25*(1-eta)*y(8);
        dz_dksi = 0.25 * (1-eta)*z(5) + 0.25*(1+eta)*z(6) - 0.25*(1+eta)*z(7) - 0.25*(1-eta)*z(8);
        
        dx_deta = -0.25*(1+ksi)*x(5) + 0.25*(1+ksi)*x(6) + 0.25*(1-ksi)*x(7) - 0.25*(1-ksi)*x(8);
        dy_deta = -0.25*(1+ksi)*y(5) + 0.25*(1+ksi)*y(6) + 0.25*(1-ksi)*y(7) - 0.25*(1-ksi)*y(8);
        dz_deta = -0.25*(1+ksi)*z(5) + 0.25*(1+ksi)*z(6) + 0.25*(1-ksi)*z(7) - 0.25*(1-ksi)*z(8);
        A_vector = [dx_dksi dy_dksi dz_dksi];
        B_vector = [dx_deta dy_deta dz_deta];
        Ans = norm(cross(A_vector, B_vector));
    end

    if faceIndex == 3
        dx_deta = -0.25*(1-zeta)*x(1) + 0.25*(1-zeta)*x(2) -0.25*(1+zeta)*x(5) + 0.25*(1+zeta)*x(6);
        dy_deta = -0.25*(1-zeta)*y(1) + 0.25*(1-zeta)*y(2) -0.25*(1+zeta)*y(5) + 0.25*(1+zeta)*y(6);
        dz_deta = -0.25*(1-zeta)*z(1) + 0.25*(1-zeta)*z(2) -0.25*(1+zeta)*z(5) + 0.25*(1+zeta)*z(6);
        
        dx_dzeta = -0.25*(1-eta)*x(1) - 0.25*(1+eta)*x(2) + 0.25*(1-eta)*x(5) + 0.25*(1+eta)*x(6);
        dy_dzeta = -0.25*(1-eta)*y(1) - 0.25*(1+eta)*y(2) + 0.25*(1-eta)*y(5) + 0.25*(1+eta)*y(6);
        dz_dzeta = -0.25*(1-eta)*z(1) - 0.25*(1+eta)*z(2) + 0.25*(1-eta)*z(5) + 0.25*(1+eta)*z(6);
        A_vector = [dx_deta dy_deta dz_deta];
        B_vector = [dx_dzeta dy_dzeta dz_dzeta];
        Ans = norm(cross(A_vector, B_vector));
    end

    if faceIndex == 4
        dx_dksi = 0.25 * (1-zeta)*x(2) - 0.25*(1-zeta)*x(3) + 0.25*(1+zeta)*x(6) - 0.25 * (1+zeta)*x(7);
        dy_dksi = 0.25 * (1-zeta)*y(2) - 0.25*(1-zeta)*y(3) + 0.25*(1+zeta)*y(6) - 0.25 * (1+zeta)*y(7);
        dz_dksi = 0.25 * (1-zeta)*z(2) - 0.25*(1-zeta)*z(3) + 0.25*(1+zeta)*z(6) - 0.25 * (1+zeta)*z(7);
    
        dx_dzeta = -0.25*(1+ksi)*x(2) - 0.25*(1-ksi)*x(3) + 0.25*(1+ksi)*x(6) + 0.25*(1-ksi)*x(7);
        dy_dzeta = -0.25*(1+ksi)*y(2) - 0.25*(1-ksi)*y(3) + 0.25*(1+ksi)*y(6) + 0.25*(1-ksi)*y(7);
        dz_dzeta = -0.25*(1+ksi)*z(2) - 0.25*(1-ksi)*z(3) + 0.25*(1+ksi)*z(6) + 0.25*(1-ksi)*z(7);
        A_vector = [dx_dksi dy_dksi dz_dksi];
        B_vector = [dx_dzeta dy_dzeta dz_dzeta];
        Ans = norm(cross(A_vector, B_vector));
    end

    if faceIndex == 5
        dx_deta = 0.25*(1-zeta)*x(3) - 0.25*(1-zeta)*x(4) + 0.25 * (1+zeta)*x(7) - 0.25*(1+zeta)*x(8);
        dy_deta = 0.25*(1-zeta)*y(3) - 0.25*(1-zeta)*y(4) + 0.25 * (1+zeta)*y(7) - 0.25*(1+zeta)*y(8);
        dz_deta = 0.25*(1-zeta)*z(3) - 0.25*(1-zeta)*z(4) + 0.25 * (1+zeta)*z(7) - 0.25*(1+zeta)*z(8);
        
        dx_dzeta = -0.25*(1+eta)*x(3) - 0.25*(1-eta)*x(4) + 0.25*(1+eta)*x(7) + 0.25*(1-eta)*x(8);
        dy_dzeta = -0.25*(1+eta)*y(3) - 0.25*(1-eta)*y(4) + 0.25*(1+eta)*y(7) + 0.25*(1-eta)*y(8);
        dz_dzeta = -0.25*(1+eta)*z(3) - 0.25*(1-eta)*z(4) + 0.25*(1+eta)*z(7) + 0.25*(1-eta)*z(8);
        A_vector = [dx_deta dy_deta dz_deta];
        B_vector = [dx_dzeta dy_dzeta dz_dzeta];
        Ans = norm(cross(A_vector, B_vector));
    end

    if faceIndex == 6
        dx_dksi = 0.25 * (1-zeta)*x(1) - 0.25*(1-zeta)*x(4) + 0.25 * (1+zeta)*x(5) - 0.25*(1+zeta)*x(8);
        dy_dksi = 0.25 * (1-zeta)*y(1) - 0.25*(1-zeta)*y(4) + 0.25 * (1+zeta)*y(5) - 0.25*(1+zeta)*y(8);
        dz_dksi = 0.25 * (1-zeta)*z(1) - 0.25*(1-zeta)*z(4) + 0.25 * (1+zeta)*z(5) - 0.25*(1+zeta)*z(8);
    
        dx_dzeta = -0.25*(1+ksi)*x(1) - 0.25*(1-ksi)*x(4) + 0.25*(1+ksi)*x(5) + 0.25*(1-ksi)*x(8);
        dy_dzeta = -0.25*(1+ksi)*y(1) - 0.25*(1-ksi)*y(4) + 0.25*(1+ksi)*y(5) + 0.25*(1-ksi)*y(8);
        dz_dzeta = -0.25*(1+ksi)*z(1) - 0.25*(1-ksi)*z(4) + 0.25*(1+ksi)*z(5) + 0.25*(1-ksi)*z(8);
        A_vector = [dx_dksi dy_dksi dz_dksi];
        B_vector = [dx_dzeta dy_dzeta dz_dzeta];
        Ans = norm(cross(A_vector, B_vector));
    end
end