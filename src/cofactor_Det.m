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