function Ans = N_matrix(ksi, eta, zeta)
    Ans = zeros(3, 24);
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

    Ans(1, 1) = N1;
    Ans(2, 2) = N1;
    Ans(3, 3) = N1;

    Ans(1, 4) = N2;
    Ans(2, 5) = N2;
    Ans(3, 6) = N2;

    Ans(1, 7) = N3;
    Ans(2, 8) = N3;
    Ans(3, 9) = N3;

    Ans(1, 10) = N4;
    Ans(2, 11) = N4;
    Ans(3, 12) = N4;

    Ans(1, 13) = N5;
    Ans(2, 14) = N5;
    Ans(3, 15) = N5;

    Ans(1, 16) = N6;
    Ans(2, 17) = N6;
    Ans(3, 18) = N6;

    Ans(1, 19) = N7;
    Ans(2, 20) = N7;
    Ans(3, 21) = N7;

    Ans(1, 22) = N8;
    Ans(2, 23) = N8;
    Ans(3, 24) = N8;

end