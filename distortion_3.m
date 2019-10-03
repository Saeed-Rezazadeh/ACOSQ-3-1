function D = distortion_3 (f , codebook , delta , Pr , T )
summation = 0;
parfor x = 1 : 8
    for y = 1 : 8
        u_index = find (T (: , 2) == x) ;
        u = T(u_index , 1) ;
        summation = summation + Pr (x , y) * delta * sum (f(u_index) .* (u - codebook (y)) .^ 2) ;
    end
end
D = summation ;
end