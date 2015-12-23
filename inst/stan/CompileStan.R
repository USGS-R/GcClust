library( rstan )

tr <- stanc( file = "inst\\stan\\MixtureModel.stan",
             model_name = "MixtureModel" )
sm <- stan_model( stanc_ret = tr, verbose = FALSE )
save( sm, file = "inst\\stan\\MixtureModel.bin" )
