function [SDR , Distortion , T , codebook] = COSQ_3 (Pr , f , T ,  codebook , delta)

FileID = fopen ('REsults.txt' , 'a') ;
D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros (1 , 8) ;
        u = T(u_index , 1) ;
        for x = 1 : 8
            for y = 1 : 8
                summation = summation + Pr(x , y) * (u - codebook (y)) ^ 2 ;
            end
            temp (x) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min (temp) ;
        T_u (u_index) = partition_index ;
    end
    T (: , 2) = T_u ;
    
    %% Optimal Centroids
    parfor y = 1 : 8
        numerator = 0 ;
        denominator = 0 ;
        for x = 1 : 8
            u_index = find (T (: , 2) == x) ;
            u = T(u_index , 1) ;
            numerator = numerator + Pr (x , y) * sum (u .* f(u_index)) ;
            denominator = denominator + Pr (x , y) * sum (f(u_index)) ;
        end
        codebook(y) = numerator / denominator ;
    end
    %% Distortion 
    D(2) = distortion_3 (f , codebook , delta , Pr , T ) ;
    fprintf (FileID , 'Overall D_3 = %f\n' ,D(2)) ; 
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_3 = %f\n' , SDR) ; 
fclose (FileID) ; 
end