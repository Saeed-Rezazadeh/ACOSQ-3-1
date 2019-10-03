function [f_u_given_y_1_y_2_y_3] = generate_pdf_rate_4_I(y_1 , y_2 , y_3 , T , Pr , f , delta )
y = (y_1 - 1) * 4 + (y_2 - 1) * 2 + y_3 ;
numerator = zeros(length(T) , 1) ;
for u_index = 1 : length(T)
    x = T(u_index , 2) ;
    numerator (u_index) = Pr(x , y) * f(u_index) ;
end

summation = 0 ;

for x = 1 : 8
    u_index_x = find(T(: , 2) == x) ;
    summation = summation + Pr(x , y) * sum(f(u_index_x)) * delta ;
end
denominator = summation ;
f_u_given_y_1_y_2_y_3 = numerator ./ denominator ;
f_u_given_y_1_y_2_y_3 = f_u_given_y_1_y_2_y_3 ./ (sum(f_u_given_y_1_y_2_y_3) * delta ) ;
end