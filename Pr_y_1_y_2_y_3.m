function Probability_y_1_2_3 = Pr_y_1_y_2_y_3 (y_1 , y_2 , y_3 , f , T , delta , Pr)
y = (y_1 - 1) * 4 + (y_2 - 1) * 2 + y_3 ;

summation = 0 ;
for x = 1 : 8
    u_index = find (T(: , 2) == x) ;
    summation = summation + Pr(x , y) * delta * sum(f(u_index)) ;
end
Probability_y_1_2_3 = summation ;
end