function Ans = K_e(D_matrix, X)
    [XG, WG] = sub_Gauss_info(2);
    K1 = zeros(24, 24);%%K element matrix
    Jv_ = zeros(1, 8);
    index = 1;
       for i = 1:2%%K
            for j = 1:2
                for k = 1:2
                    B = B_matrix(XG(i), XG(j), XG(k), X);
                    B_T = B.';
                    Jv_(index) = Jv(XG(i), XG(j), XG(k), X);
                    K1 = K1 + WG(i) * WG(j) * WG(k) * B_T * D_matrix * B;
                    index = index + 1;
                end
            end 
       end
       Jv_Mean = mean(Jv_);
       Ans = K1 * Jv_Mean;
end