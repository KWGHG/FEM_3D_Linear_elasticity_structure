function Ans = Fh_e(tractionForce, X, faceIndex)
    [XG, WG] = sub_Gauss_info(2);
    Fhe = zeros(24, 1);%%Fhe element matrix
    Js_ = zeros(1, 4);
    index = 1;
    if faceIndex == 1
        for j = 1:2
            for k = 1:2
                N_ = N_matrix(XG(j), XG(k), -1);
                Nt = N_.';
                Js_(index)= Js(XG(j), XG(k), -1, X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForce;
                index = index + 1;
            end
        end 
    end

    if faceIndex == 2
        for j = 1:2
            for k = 1:2
                N_ = N_matrix(XG(j), XG(k), 1);
                Nt = N_.';
                Js_(index) = Js(XG(j), XG(k), 1, X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForce;
                index = index + 1;
            end
        end 
    end

    if faceIndex == 3
        for j = 1:2
            for k = 1:2
                N_ = N_matrix(1, XG(j), XG(k));
                Nt = N_.';
                Js_(index) = Js(1, XG(j), XG(k), X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForce;
                index = index + 1;
            end
        end 
    end

    if faceIndex == 4
        for j = 1:2
            for k = 1:2
                N_ = N_matrix( XG(j), 1, XG(k));
                Nt = N_.';
                Js_(index) = Js(XG(j), 1, XG(k), X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForc_;
                index = index + 1;
            end
        end 
    end

    if faceIndex == 5
        for j = 1:2
            for k = 1:2
                N_ = N_matrix(-1, XG(j), XG(k));
                Nt = N_.';
                Js_(index) = Js(-1, XG(j), XG(k), X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForce;
                index = index + 1;
            end
        end 
    end

    if faceIndex == 6
        for j = 1:2
            for k = 1:2
                N_ = N_matrix( XG(j), -1, XG(k));
                Nt = N_.';
                Js_(index) = Js(XG(j), -1, XG(k), X, faceIndex);
                Fhe = Fhe + WG(j) * WG(k) * Nt * tractionForce;
                index = index + 1;
            end
        end 
    end

    Js_Mean = max(Js_);
    Ans = Fhe * Js_Mean;
end