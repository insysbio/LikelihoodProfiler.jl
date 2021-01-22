
f_1p(x) = 5.0 + (x[1]-3.0)^2 # [100., missing, missing, missing]

f_2p_1im(x) = 5.0 + (x[1]-3.0)^2 + 0.0*x[2]  # [3., missing, missing, missing]

f_2p(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2  # [3., 4., missing, missing]

f_3p_1im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 # [3., missing, missing, missing]

f_3p_1im_dep(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2

f_4p_2im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4] # [3., 4., missing, missing]

f_4p_3im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 + 0.0*x[4] # [3., missing, missing, missing]

f_1p_ex(x) = 5.0 + (x[1]-1e-8)^2 # [1e-8, missing, missing, missing]

f_5p_3im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + (x[3]/x[4]-4.0)^2 + 0.0*x[5]

f_3p_im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3]
