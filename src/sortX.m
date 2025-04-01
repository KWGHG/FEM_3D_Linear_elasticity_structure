function Ans = sortX(X)
    Ans = zeros(8, 4);
    Z1 = zeros(4, 4);
    Z2 = zeros(4, 4);

    Z = sort(X(:, 4));
    index1 = 1; index2 = 1;
    for i = 1:8
       if X(i, 4) < Z(5)
           Z1(index1,:) = X(i, :);
           index1 = index1 + 1;
       else
           Z2(index2, :) = X(i, :);
           index2 = index2 + 1;
       end
    end
    
    Y1 = sort(Z1(:, 3));
    Y2 = sort(Z2(:, 3));
    Z11_1 = zeros(2, 4);
    Z11_2 = zeros(2, 4);
    Z22_1 = zeros(2, 4);
    Z22_2 = zeros(2, 4);

    for i = 1:4
       if Z1(i, 3) == Y1(1)
           Z11_1(1,:) = Z1(i, :);
       end
       if Z1(i, 3) == Y1(2)
           Z11_1(2,:) = Z1(i, :);
       end
       if Z1(i, 3) == Y1(3)
           Z11_2(1,:) = Z1(i, :);
       end
       if Z1(i, 3) == Y1(4)
           Z11_2(2,:) = Z1(i, :);
       end

       if Z2(i, 3) == Y2(1)
           Z22_1(1,:) = Z2(i, :);
       end
       if Z2(i, 3) == Y2(2)
           Z22_1(2, :) = Z2(i, :);
       end
       if Z2(i, 3) == Y2(3)
           Z22_2(1, :) = Z2(i, :);
       end
       if Z2(i, 3) == Y2(4)
           Z22_2(2, :) = Z2(i, :);
       end
    end

    X1_1 = sort(Z11_1(:, 2));
    X1_2 = sort(Z11_2(:, 2));
    X2_1 = sort(Z22_1(:, 2));
    X2_2 = sort(Z22_2(:, 2));

    for i = 1:2
        if Z11_1(i, 2) == X1_1(1)
            Ans(4, :) = Z11_1(i, :);
        else
            Ans(1, :) = Z11_1(i, :);
        end

        if Z11_2(i, 2) == X1_2(1)
            Ans(3, :) = Z11_2(i, :);
        else
            Ans(2, :) = Z11_2(i, :);
        end
        %%2
        if Z22_1(i, 2) == X2_1(1)
            Ans(8, :) = Z22_1(i, :);
        else
            Ans(5, :) = Z22_1(i, :);
        end

        if Z22_2(i, 2) == X2_2(1)
            Ans(7, :) = Z22_2(i, :);
        else
            Ans(6, :) = Z22_2(i, :);
        end

    end



end