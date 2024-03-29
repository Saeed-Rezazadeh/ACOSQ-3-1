function codebook = init_codebook(f , delta , T , alpha)
F = zeros (length(f) , 1) ; 

for u_index = 1 : length(T)
    u = T(u_index , 1);
    index = find (T(: , 1) <= u) ;
    F(u_index) = sum(f(index)) * delta ;
end
Training_Set = zeros (alpha , 1) ; 
for u_index = 1 : alpha 
    u = rand ;
    index = find(u < F(: , 1) ,1,'first') ; 
    Training_Set(u_index , 1) = T(index , 1) ; 
end 
[~  , codebook] = kmeans (Training_Set , 2 , 'MaxIter',1000 , 'OnlinePhase','on') ;
end