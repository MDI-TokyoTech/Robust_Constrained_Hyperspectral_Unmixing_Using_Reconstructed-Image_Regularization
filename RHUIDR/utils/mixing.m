function HSI_rec = mixing(Mat_endmember_library, abundance)
    [n1, n2, m] = size(abundance);
    l = size(Mat_endmember_library, 1);
    Mat_abundance = cube2mat(abundance, n1, n2, m);
    HSI_rec = mat2cube(Mat_endmember_library*Mat_abundance, n1, n2, l);
end

