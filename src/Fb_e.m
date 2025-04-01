function Ans = Fb_e(bodyForce, X)
    [XG, WG] = sub_Gauss_info(2);
    Fbe = zeros(24, 1);%%Fb_e element matrix
    Jv_ = zeros(1, 8);
    index = 1;
       for i = 1:2
            for j = 1:2
                for k = 1:2
                    N_ = N_matrix(XG(i), XG(j), XG(k));
                    Nt = N_.';
                    Jv_(index) = Jv(XG(i), XG(j), XG(k), X);
                    Fbe = Fbe + WG(i) * WG(j) * WG(k) * Nt * bodyForce;
                    index = index + 1;
                end
            end 
       end
       Jv_Mean = max(Jv_);
       Ans = Fbe * Jv_Mean;
end