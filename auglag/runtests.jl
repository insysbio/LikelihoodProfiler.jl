using Test

include("auglag2.jl")
include("cases_func.jl")

# test models
model_1p = plmodel(f_1p, [3.])
model_2p_1im = plmodel(f_2p_1im, [3.,1.])
model_2p = plmodel(f_2p, [3.,4.])
model_3p_1im = plmodel(f_3p_1im, [3.,4.,1.]) 
model_3p_1im_dep = plmodel(f_3p_1im_dep, [3.,2.,2.])
model_4p_2im = plmodel(f_4p_2im, [3.,4.,1.,1.])
model_4p_3im = plmodel(f_4p_3im,[3.,4.,1.,1.])
model_1p_ex = plmodel(f_1p_ex, [1e-8])
model_5p_3im = plmodel(f_5p_3im, [3., 0.1, 4, 1.1, 8.])
model_3p_im = plmodel(f_3p_im, [3.,1.,1,])

# identification results
res_1p = identify_auglag(model_1p)
res_2p_1im = identify_auglag(model_2p_1im)
res_2p = identify_auglag(model_2p)
res_3p_1im = identify_auglag(model_3p_1im)
res_3p_1im_dep = identify_auglag(model_3p_1im_dep)
res_4p_2im = identify_auglag(model_4p_2im)
res_4p_3im = identify_auglag(model_4p_3im)
res_1p_ex = identify_auglag(model_1p_ex)
res_5p_3im = identify_auglag(model_5p_3im)
res_3p_im = identify_auglag(model_3p_im)

# check results
#1
  @test isapprox(res_1p[1][1][1], 1.0, atol=1e-2)
  @test isapprox(res_1p[1][2][1], 5.0, atol=1e-2)
#2
  @test isapprox(res_2p_1im[1][1][1], 1.0, atol=1e-2)
  @test isapprox(res_2p_1im[1][2][1], 5.0, atol=1e-2)
  @test isequal(res_2p_1im[2][1][3], :unidentifiable)
  @test isequal(res_2p_1im[2][2][3], :unidentifiable)
#3
  @test isapprox(res_2p[1][1][1], 1.0, atol=1e-2)
  @test isapprox(res_2p[1][2][1], 5.0, atol=1e-2)
  @test isapprox(res_2p[2][1][1], 2.0, atol=1e-2)
  @test isapprox(res_2p[2][2][1], 6.0, atol=1e-2)
#4
  @test isapprox(res_3p_1im[1][1][1], 1.0, atol=1e-2) # wrong
  @test isapprox(res_3p_1im[1][2][1], 5.0, atol=1e-2) # wrong
  @test isequal(res_3p_1im[2][1][3], :unidentifiable)
  @test isequal(res_3p_1im[2][2][3], :unidentifiable)
  @test isequal(res_3p_1im[3][1][3], :unidentifiable)
  @test isequal(res_3p_1im[3][2][3], :unidentifiable)
#5
  @test isapprox(res_3p_1im_dep[1][1][1], 1.0, atol=1e-2)
  @test isapprox(res_3p_1im_dep[1][2][1], 5.0, atol=1e-2)
  @test isapprox(res_3p_1im_dep[2][1][1], 2.0-2.0*sqrt(2.), atol=1e-2)
  @test isapprox(res_3p_1im_dep[2][2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
  @test isequal(res_3p_1im_dep[3][1][3], :unidentifiable)
  @test isequal(res_3p_1im_dep[3][2][3], :unidentifiable)
#6
  @test isapprox(res_4p_2im[1][1][1], 1.0, atol=1e-2)
  @test isapprox(res_4p_2im[1][2][1], 5.0, atol=1e-2)
  @test isapprox(res_4p_2im[2][1][1], 2.0, atol=1e-2)
  @test isapprox(res_4p_2im[2][2][1], 6.0, atol=1e-2)
  @test isequal(res_4p_2im[3][1][3], :unidentifiable)
  @test isequal(res_4p_2im[3][2][3], :unidentifiable)
  @test isequal(res_4p_2im[4][1][3], :unidentifiable)
  @test isequal(res_4p_2im[4][2][3], :unidentifiable)
#7
  @test isapprox(res_4p_3im[1][1][1], 1.0, atol=1e-2) # wrong
  @test isapprox(res_4p_3im[1][2][1], 5.0, atol=1e-2) # wrong
  @test isequal(res_4p_3im[2][1][3], :unidentifiable)
  @test isequal(res_4p_3im[2][2][3], :unidentifiable)
  @test isequal(res_4p_3im[3][1][3], :unidentifiable)
  @test isequal(res_4p_3im[3][2][3], :unidentifiable) # failure
  @test isequal(res_4p_3im[4][1][3], :unidentifiable)
  @test isequal(res_4p_3im[4][2][3], :unidentifiable)
#8
  @test isapprox(res_1p_ex[1][1][1], -2.0+1e-8, atol=1e-2)
  @test isapprox(res_1p_ex[1][2][1], 2.0-1e-8, atol=1e-2)
#9
  @test isapprox(res_5p_3im[1][1][1], 1.0, atol=1e-2) # wrong
  @test isapprox(res_5p_3im[1][2][1], 5.0, atol=1e-2) # wrong
  @test isequal(res_5p_3im[2][1][3], :unidentifiable)
  @test isapprox(res_5p_3im[2][2][1], log(3), atol=1e-2)
  @test isequal(res_5p_3im[3][1][3], :unidentifiable)
  @test isequal(res_5p_3im[3][2][3], :unidentifiable) # wrong
  @test isequal(res_5p_3im[4][1][3], :unidentifiable)
  @test isequal(res_5p_3im[4][2][3], :unidentifiable) # failure
  @test isequal(res_5p_3im[5][1][3], :unidentifiable)
  @test isequal(res_5p_3im[5][2][3], :unidentifiable)
#10
  @test isapprox(res_3p_im[1][1][1], 1.0, atol=1e-2) # wrong
  @test isapprox(res_3p_im[1][2][1], 5.0, atol=1e-2) # wrong
  @test isequal(res_3p_im[2][1][3], :unidentifiable)
  @test isapprox(res_3p_im[2][2][1], log(3), atol=1e-2)
  @test isequal(res_3p_im[2][1][3], :unidentifiable)
  @test isequal(res_3p_im[3][2][3], :unidentifiable)
