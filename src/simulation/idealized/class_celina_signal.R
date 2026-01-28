celina_like_simulator_signal <- setRefClass(
  "celina_simulator_signal",
  contains = "SimDataset",
  
  fields = list(
    scenario = "character",
    pattern = "character",
    pi = "numeric",
    fold_change = "numeric",
    center = "data.frame",
    radius = "numeric"
  ),
  
  methods = list(
    initialize = function(seed = 1, dispersion = 0.95,
                          scenario = "I",
                          pattern = "hotspot",
                          pi = 0.4,
                          fold_change = c(2, 0.5)) {
      .self$seed = seed
      .self$dispersion = dispersion
      .self$scenario = scenario
      .self$pattern = pattern
      .self$pi = pi
      .self$fold_change = fold_change
      .self$simulate()
    },
    
    simulate = function() {
      set.seed(.self$seed)
      stop("CELINA-like signal simulation not implemented yet")
    }
  )
)
