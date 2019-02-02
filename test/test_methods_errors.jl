
for method = methods_list
  str = "Using " * String(method)
  @testset "$str" begin
      @testset "bad local_alg" begin
          @test_logs (:warn, r"Using local_alg .* may result in wrong output\.") get_endpoint(
              [3., 2., 2.1],
              1,
              f_3p_1im_dep,
              method;
              loss_crit = 9.,
              local_alg = :LN_BOBYQA
              )

          @test_throws ArgumentError get_endpoint(
              [3., 2., 2.1],
              1,
              f_3p_1im_dep,
              method;
              loss_crit = 9.,
              local_alg = :xxx
              )
      end

      @testset "warn for LN_NELDERMEAD" begin
          @test_logs (:warn, r"Close-to-zero parameters found when using :LN_NELDERMEAD\.") get_endpoint(
              [3., 0., 2.1],
              1,
              f_3p_1im_dep,
              method;
              loss_crit = 90.,
              local_alg = :LN_NELDERMEAD
              )
      end
  end
end
