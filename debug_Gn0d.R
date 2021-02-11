devtools::load_all()

dA = Beta$new(shape1 = 2, shape2 = 3)
dB = Binomial$new(size = 2, prob = 0.2)
dC = Binomial$new(size = 20, prob = 0.2)

pm_continuous_case(-0.9, dA, dA, n = 7) # -0.914
pm_continuous_case(-0.6, dA, dA, n = 7) # -0.611
pm_continuous_case(-0.3, dA, dA, n = 7) # -0.306
pm_continuous_case( 0.3, dA, dA, n = 7) #  0.304
pm_continuous_case( 0.6, dA, dA, n = 7) #  0.606
pm_continuous_case( 0.9, dA, dA, n = 7) #  0.904

pm_discrete_case(-0.5, dB, dB, n = 17)  # -0.937
pm_discrete_case(-0.3, dB, dB, n = 7)  # -0.501
pm_discrete_case(-0.2, dB, dB, n = 7)  # -0.322
pm_discrete_case( 0.3, dB, dB, n = 7)  #  0.418
pm_discrete_case( 0.6, dB, dB, n = 7)  #  0.769
pm_discrete_case( 0.8, dB, dB, n = 17)  #  0.944

pm_discrete_case(-0.9, dC, dC, n = 7) # -0.939
pm_discrete_case(-0.6, dC, dC, n = 7) # -0.624
pm_discrete_case(-0.3, dC, dC, n = 7) # -0.311
pm_discrete_case( 0.3, dC, dC, n = 7) #  0.310
pm_discrete_case( 0.6, dC, dC, n = 7) #  0.618
pm_discrete_case( 0.9, dC, dC, n = 7) #  0.925
