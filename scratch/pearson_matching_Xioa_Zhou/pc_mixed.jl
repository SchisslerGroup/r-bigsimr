include("common.jl")

margins = [
    Beta(2, 3),
    Normal(12, 2),
    Binomial(20, 0.2)
]

typeof(margins[1]) <: ContinuousUnivariateDistribution
typeof(margins[2]) <: ContinuousUnivariateDistribution
typeof(margins[3]) <: ContinuousUnivariateDistribution

typeof(margins[1]) <: DiscreteUnivariateDistribution
typeof(margins[2]) <: DiscreteUnivariateDistribution
typeof(margins[3]) <: DiscreteUnivariateDistribution

[typeof(m) <: ContinuousUnivariateDistribution for m in margins]
