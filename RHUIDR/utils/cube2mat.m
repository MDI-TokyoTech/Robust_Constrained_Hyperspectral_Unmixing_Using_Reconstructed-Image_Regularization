function Mat = cube2mat(Cube, n1, n2, n3)
    Mat = reshape(Cube, n1*n2, n3)';
end
    