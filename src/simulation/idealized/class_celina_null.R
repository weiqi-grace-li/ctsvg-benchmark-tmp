celina_simulator_null <- setRefClass(
  "celina_simulator_null",
  contains = "SimDataset",
  
  fields = list(
    scenario = "character"
  ),
  
  methods = list(
    initialize = function(seed = 1, dispersion = 0.95, scenario = "I") {
      .self$seed = seed
      .self$dispersion = dispersion
      .self$scenario = scenario
      .self$simulate()
    },
    
    simulate = function() {
      set.seed(.self$seed)
      stop("CELINA-like null simulation not implemented yet")
    }
  )
)
