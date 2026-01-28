stance_simulator_signal <- setRefClass(
  "stance_simulator_signal",
  contains = "SimDataset",
  
  fields = list(
    center = "data.frame",
    radius = "numeric"
  ),
  
  methods = list(
    initialize = function(seed = 1, dispersion = 0.7,
                          cell_type_proportion = c(Type1=0.1, Type2=0.3, Type3=0.6)) {
      .self$seed = seed
      .self$dispersion = dispersion
      .self$cell_type_proportion = cell_type_proportion
      .self$simulate()
    },
    
    simulate = function() {
      set.seed(.self$seed)
      stop("STANCE-like signal simulation not implemented yet")
    }
  )
)
