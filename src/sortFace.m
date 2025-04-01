function Ans  = sortFace(beforeIndex, beforeX, X)
    face1 = [1 2 3 4];
    face2 = [5 6 7 8];
    face3 = [1 2 5 6];
    face4 = [2 3 6 7];
    face5 = [3 4 7 8];
    face6 = [1 4 5 8];
    newIndex = zeros(1, 4);
    oldX = beforeX(beforeIndex, :);

    oldX_Nodes = oldX(:, 1);
    for i = 1:8
        if X(i, 1) == oldX_Nodes(1)
            newIndex(1) = i;
        end
        if X(i, 1) == oldX_Nodes(2)
            newIndex(2) = i;
        end
        if X(i, 1) == oldX_Nodes(3)
            newIndex(3) = i;
        end
        if X(i, 1) == oldX_Nodes(4)
            newIndex(4) = i;
        end
    end
    sortX_ = sort(newIndex);
    if sortX_(:) == face1(:)
        Ans = 1;
    end

    if sortX_(:) == face2(:)
        Ans = 2;
    end

    if sortX_(:) == face3(:)
        Ans = 3;
    end

    if sortX_(:) == face4(:)
        Ans = 4;
    end

    if sortX_(:) == face5(:)
        Ans = 5;
    end

    if sortX_(:) == face6(:)
        Ans = 6;
    end
end