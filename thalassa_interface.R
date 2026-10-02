library(ggplot2)

source('model.R')

get.overview <- function() {
  # Markdown shown in the Overview tab (see get.overview.markdown() in shiny-cea).
  # The text lives in overview.md and describes what simulate() in model.R and
  # this interface do, so it must be kept in sync with them.
  return(paste(readLines('overview.md'), collapse='\n'))
}

# Grouped by screening modality: the app reads a named list of entries as the
# headings the strategies hang from (see normalize.strategies() in shiny-cea).
# The name is still what run.simulation() is given, so grouping changes nothing
# but how the strategies are shown.
get.strategies <- function() {
  return(list(
    `No screening`=list(
      list(name='no_screening', display.name='No screening')
    ),
    Test=list(
      list(name='test_biennial', display.name='Biennial test, 50-74')
    ),
    Procedure=list(
      list(name='procedure_10y', display.name='Procedure every 10 years, 50-79'),
      list(name='procedure_45', display.name='Procedure every 10 years, 45-79')
    )
  ))
}

# Base values of the three age-dependent (and calibrated) transition
# probabilities, one per stratum of MODEL.STRATA.
BASE.LGL.ONSET <- c(0.0050, 0.0100, 0.0180, 0.0270, 0.0320, 0.0340, 0.0370)
BASE.FPL.ONSET <- c(0.00018, 0.00035, 0.00056, 0.00077, 0.00088, 0.00095, 0.00102)
BASE.HGL.PROGRESS <- c(0.0053, 0.0064, 0.0094, 0.0139, 0.0188, 0.0225, 0.0251)

# Descriptors of every model parameter, in the shape the app expects from
# get.parameters(). model.R only consumes parameter values, so the base case is
# declared here.
stratified.parameter <- function(name, display.name, values, class) {
  lapply(seq_along(MODEL.STRATA), function(i) {
    list(name = name,
         display.name = display.name,
         base.value = values[i],
         stratum = MODEL.STRATA[i],
         class = class)
  })
}

MODEL.PARAMETERS <- c(
  stratified.parameter('p.lgl.onset',
                       'Annual probability of developing a low-grade lesion',
                       BASE.LGL.ONSET, 'Natural history (slow pathway)'),
  stratified.parameter('p.fpl.onset',
                       'Annual probability of developing a fast-pathway lesion',
                       BASE.FPL.ONSET, 'Natural history (fast pathway)'),
  stratified.parameter('p.hgl.progress',
                       'Annual probability that a high-grade lesion becomes cancer',
                       BASE.HGL.PROGRESS, 'Natural history (slow pathway)'),
  list(
    list(name = 'p.lgl.progress',
         display.name = 'Annual probability that a low-grade lesion becomes high-grade',
         base.value = 0.015, class = 'Natural history (slow pathway)'),
    list(name = 'p.lgl.regress',
         display.name = 'Annual probability that a low-grade lesion regresses',
         base.value = 0.015, class = 'Natural history (slow pathway)'),
    list(name = 'p.hgl.regress',
         display.name = 'Annual probability that a high-grade lesion regresses',
         base.value = 0.005, class = 'Natural history (slow pathway)'),
    list(name = 'p.fpl.progress',
         display.name = 'Annual probability that a fast-pathway lesion becomes cancer',
         base.value = 0.040, class = 'Natural history (fast pathway)'),
    list(name = 'p.fpl.regress',
         display.name = 'Annual probability that a fast-pathway lesion regresses',
         base.value = 0.020, class = 'Natural history (fast pathway)'),
    list(name = 'p.stage.progress.slow',
         display.name = 'Annual probability of localized to advanced stage, slow-pathway cancer',
         base.value = 0.300, class = 'Cancer progression'),
    list(name = 'p.stage.progress.fast',
         display.name = 'Annual probability of localized to advanced stage, fast-pathway cancer',
         base.value = 0.450, class = 'Cancer progression'),
    list(name = 'p.symptomatic.early.slow',
         display.name = 'Annual probability of symptomatic diagnosis, localized slow-pathway cancer',
         base.value = 0.270, class = 'Cancer progression'),
    list(name = 'p.symptomatic.late.slow',
         display.name = 'Annual probability of symptomatic diagnosis, advanced slow-pathway cancer',
         base.value = 0.550, class = 'Cancer progression'),
    list(name = 'p.symptomatic.early.fast',
         display.name = 'Annual probability of symptomatic diagnosis, localized fast-pathway cancer',
         base.value = 0.120, class = 'Cancer progression'),
    list(name = 'p.symptomatic.late.fast',
         display.name = 'Annual probability of symptomatic diagnosis, advanced fast-pathway cancer',
         base.value = 0.450, class = 'Cancer progression'),
    list(name = 'p.cure.early',
         display.name = 'Annual probability of cure after a localized cancer diagnosis',
         base.value = 0.550, class = 'Cancer survival'),
    list(name = 'p.cancer.death.early',
         display.name = 'Annual probability of cancer death, localized cancer',
         base.value = 0.030, class = 'Cancer survival'),
    list(name = 'p.cure.late',
         display.name = 'Annual probability of cure after an advanced cancer diagnosis',
         base.value = 0.120, class = 'Cancer survival'),
    list(name = 'p.cancer.death.late',
         display.name = 'Annual probability of cancer death, advanced cancer',
         base.value = 0.250, class = 'Cancer survival'),
    list(name = 'p.survivor.death',
         display.name = 'Annual probability of cancer death after recovery (late recurrence)',
         base.value = 0.010, class = 'Cancer survival'),
    list(name = 'mortality.other.base',
         display.name = 'Other-cause mortality at age 20',
         base.value = 0.0005, class = 'General'),
    list(name = 'mortality.other.rate',
         display.name = 'Gompertz rate of other-cause mortality',
         base.value = 0.085, class = 'General'),
    list(name = 'adherence.test',
         display.name = 'Proportion attending an invitation to the screening test',
         base.value = 0.650, class = 'Screening (test)'),
    list(name = 'sens.test.lgl',
         display.name = 'Test sensitivity for low-grade lesion',
         base.value = 0.040, class = 'Screening (test)'),
    list(name = 'sens.test.hgl',
         display.name = 'Test sensitivity for high-grade lesion',
         base.value = 0.240, class = 'Screening (test)'),
    list(name = 'sens.test.fpl',
         display.name = 'Test sensitivity for fast-pathway lesion',
         base.value = 0.080, class = 'Screening (test)'),
    list(name = 'sens.test.cancer',
         display.name = 'Test sensitivity for localized cancer',
         base.value = 0.550, class = 'Screening (test)'),
    list(name = 'spec.test',
         display.name = 'Test specificity',
         base.value = 0.950, class = 'Screening (test)'),
    list(name = 'adherence.procedure',
         display.name = 'Proportion attending an invitation to the screening procedure',
         base.value = 0.550, class = 'Screening (procedure)'),
    list(name = 'sens.procedure.lgl',
         display.name = 'Procedure sensitivity for low-grade lesion',
         base.value = 0.750, class = 'Screening (procedure)'),
    list(name = 'sens.procedure.hgl',
         display.name = 'Procedure sensitivity for high-grade lesion',
         base.value = 0.920, class = 'Screening (procedure)'),
    list(name = 'sens.procedure.fpl',
         display.name = 'Procedure sensitivity for fast-pathway lesion',
         base.value = 0.550, class = 'Screening (procedure)'),
    list(name = 'sens.procedure.cancer',
         display.name = 'Procedure sensitivity for localized cancer',
         base.value = 0.950, class = 'Screening (procedure)'),
    list(name = 'rr.detection.fast',
         display.name = 'Relative detection of fast-pathway versus slow-pathway cancer',
         base.value = 0.850, class = 'Screening (procedure)'),
    list(name = 'p.procedure.complication',
         display.name = 'Probability of a serious complication per procedure',
         base.value = 0.002, class = 'Screening (procedure)'),
    list(name = 'cost.test',
         display.name = 'Cost per screening test performed',
         base.value = 30, class = 'Costs'),
    list(name = 'cost.procedure',
         display.name = 'Cost per procedure performed',
         base.value = 800, class = 'Costs'),
    list(name = 'cost.removal',
         display.name = 'Additional cost per lesion removed',
         base.value = 250, class = 'Costs'),
    list(name = 'cost.complication',
         display.name = 'Cost per serious procedure complication',
         base.value = 3000, class = 'Costs'),
    list(name = 'cost.treatment.early',
         display.name = 'Annual cost of treating a localized cancer',
         base.value = 18000, class = 'Costs'),
    list(name = 'cost.treatment.late',
         display.name = 'Annual cost of treating an advanced cancer',
         base.value = 40000, class = 'Costs'),
    list(name = 'cost.followup',
         display.name = 'Annual follow-up cost of a cancer survivor',
         base.value = 1200, class = 'Costs'),
    list(name = 'utility.cancer.early',
         display.name = 'Utility of a year with treated localized cancer',
         base.value = 0.720, class = 'Utilities'),
    list(name = 'utility.cancer.late',
         display.name = 'Utility of a year with treated advanced cancer',
         base.value = 0.500, class = 'Utilities'),
    list(name = 'utility.survivor',
         display.name = 'Utility of a year as a cancer survivor',
         base.value = 0.920, class = 'Utilities'),
    list(name = 'disutility.complication',
         display.name = 'Disutility of a serious procedure complication',
         base.value = 0.020, class = 'Utilities'),
    list(name = 'discount',
         display.name = 'Discount rate for costs and effects',
         base.value = 0.030, class = 'General')
  )
)

# Base case parameter values in the shape simulate() expects: a named list where
# stratified parameters are themselves named lists over MODEL.STRATA.
default.parameters <- function() {
  pars <- list()
  for (descriptor in MODEL.PARAMETERS) {
    if (is.null(descriptor$stratum)) {
      pars[[descriptor$name]] <- descriptor$base.value
    } else {
      if (is.null(pars[[descriptor$name]])) pars[[descriptor$name]] <- list()
      pars[[descriptor$name]][[descriptor$stratum]] <- descriptor$base.value
    }
  }
  return(pars)
}

get.parameters <- function() {
  return(MODEL.PARAMETERS)
}

get.strata <- function() {
  return(MODEL.STRATA)
}

get.model.settings <- function() {
  return(list(
    cost.unit='€',
    effectiveness.unit='QALY'
  ))
}

# State diagrams shown in the Overview tab, one tab each. They describe
# transition.matrix() and apply.screening() in model.R and must be kept in sync
# with them. Other-cause death is left out of both: every living state leads to
# it, and an arrow from each into one node would hide the rest.
MODEL.STATE.NODES <- data.frame(
  id=c('normal', 'lgl', 'hgl', 'fpl',
       'pre.early.slow', 'pre.late.slow', 'pre.early.fast', 'pre.late.fast',
       'clin.early', 'clin.late', 'survivor', 'dead.cancer'),
  label=c('Normal', 'Low-grade\nlesion', 'High-grade\nlesion', 'Fast-pathway\nlesion',
          'Preclinical\nlocalized (slow)', 'Preclinical\nadvanced (slow)',
          'Preclinical\nlocalized (fast)', 'Preclinical\nadvanced (fast)',
          'Diagnosed\nlocalized', 'Diagnosed\nadvanced', 'Survivor', 'Cancer death'),
  group=c('No lesion', 'Slow pathway', 'Slow pathway', 'Fast pathway',
          'Slow pathway', 'Slow pathway', 'Fast pathway', 'Fast pathway',
          'Diagnosed cancer', 'Diagnosed cancer', 'Survivor', 'Death'),
  description=c(
    'No lesion. The whole cohort starts here.',
    'Low-grade lesion, the first step of the slow pathway.',
    'High-grade lesion, the precursor of slow-pathway cancer.',
    'Fast-pathway lesion, harder to detect than the slow-pathway ones.',
    'Undiagnosed localized cancer, slow pathway.',
    'Undiagnosed advanced cancer, slow pathway.',
    'Undiagnosed localized cancer, fast pathway.',
    'Undiagnosed advanced cancer, fast pathway.',
    'Diagnosed localized cancer, under treatment.',
    'Diagnosed advanced cancer, under treatment.',
    'Recovered after treatment, with a small annual risk of late recurrence.',
    'Death from cancer.'
  ),
  # Laid out by hand: one column per step of the disease, the slow pathway on the
  # upper row and the fast one on the lower.
  x=290 * c(0, 1, 2, 2, 3, 4, 3, 4, 5, 5, 6, 7),
  y=75 * c(0, -1, -1, 1, -1, -1, 1, 1, -1, 1, 0, 0)
)

get.model.states <- function() {
  natural.history <- data.frame(
    from=c('normal', 'normal', 'lgl', 'lgl', 'hgl', 'hgl', 'fpl', 'fpl',
           'pre.early.slow', 'pre.early.slow', 'pre.late.slow',
           'pre.early.fast', 'pre.early.fast', 'pre.late.fast',
           'clin.early', 'clin.early', 'clin.late', 'clin.late', 'survivor'),
    to=c('lgl', 'fpl', 'hgl', 'normal', 'pre.early.slow', 'lgl', 'pre.early.fast', 'normal',
         'clin.early', 'pre.late.slow', 'clin.late',
         'clin.early', 'pre.late.fast', 'clin.late',
         'survivor', 'dead.cancer', 'survivor', 'dead.cancer', 'dead.cancer'),
    label=c('Onset', 'Onset', 'Progression', 'Regression', 'Progression', 'Regression',
            'Progression', 'Regression',
            'Symptoms', 'Stage progression', 'Symptoms',
            'Symptoms', 'Stage progression', 'Symptoms',
            'Cure', 'Death', 'Cure', 'Death', 'Late recurrence'),
    description=c('p.lgl.onset', 'p.fpl.onset', 'p.lgl.progress', 'p.lgl.regress',
                  'p.hgl.progress', 'p.hgl.regress', 'p.fpl.progress', 'p.fpl.regress',
                  'p.symptomatic.early.slow', 'p.stage.progress.slow', 'p.symptomatic.late.slow',
                  'p.symptomatic.early.fast', 'p.stage.progress.fast', 'p.symptomatic.late.fast',
                  'p.cure.early', 'p.cancer.death.early', 'p.cure.late', 'p.cancer.death.late',
                  'p.survivor.death')
  )

  screening.states <- c('normal', 'lgl', 'hgl', 'fpl',
                        'pre.early.slow', 'pre.late.slow', 'pre.early.fast', 'pre.late.fast',
                        'clin.early', 'clin.late')
  screening <- data.frame(
    from=c('lgl', 'hgl', 'fpl',
           'pre.early.slow', 'pre.early.fast', 'pre.late.slow', 'pre.late.fast'),
    to=c('normal', 'normal', 'normal',
         'clin.early', 'clin.early', 'clin.late', 'clin.late'),
    label=c('Removal', 'Removal', 'Removal',
            'Detection', 'Detection', 'Detection', 'Detection'),
    description=c(
      'Detected and removed: sens.procedure.lgl, after a positive test (sens.test.lgl) in a test round.',
      'Detected and removed: sens.procedure.hgl, after a positive test (sens.test.hgl) in a test round.',
      'Detected and removed: sens.procedure.fpl, after a positive test (sens.test.fpl) in a test round.',
      'Diagnosed at the localized stage: sens.procedure.cancer, after a positive test (sens.test.cancer) in a test round.',
      'Diagnosed at the localized stage, as the slow-pathway one but scaled by rr.detection.fast.',
      'Diagnosed at the advanced stage: sens.procedure.cancer, after a positive test (sens.test.cancer) in a test round.',
      'Diagnosed at the advanced stage, as the slow-pathway one but scaled by rr.detection.fast.'
    )
  )

  return(list(
    `Natural history`=list(
      description='The annual transitions of the two pathways, from a first lesion to the outcome of a diagnosed cancer. Every living state is also exposed to other-cause death, which is not drawn. Hover a transition for the parameter behind it.',
      nodes=MODEL.STATE.NODES,
      edges=natural.history
    ),
    `Screening round`=list(
      description='What a screening round moves, before the natural history of the cycle is applied: detected lesions are removed and detected preclinical cancers are diagnosed at the stage they had reached.',
      nodes=MODEL.STATE.NODES[MODEL.STATE.NODES$id %in% screening.states, ],
      edges=screening
    )
  ))
}

run.simulation <- function(strategies, pars) {
  # simulate() takes the parameters as the named list the app already builds,
  # and copes with any of them arriving as one value per stratum.
  return(simulate(strategies, pars))
}

# Strata used for calibration: all of them.
#
# The cohort starts at age 20 with no lesions, so the first stratum is burn-in
# and its cancer incidence is essentially zero whatever the parameters. It is
# excluded from the targets, which is the documented way of leaving a stratum
# out of the error, rather than from the scheme's strata: the app builds the
# initial guess from the base values over *every* stratum
# (`unlist(calib.pars[scheme$parameters])` in shiny_calibration.R), so a scheme
# that calibrates a subset of the strata gets a longer vector than it consumes
# and every value ends up on the wrong stratum.
CALIB.STRATA <- MODEL.STRATA

# Targets, in the range the literature reports for a population without
# screening: cancer incidence and the stage distribution at diagnosis are
# registry-like figures, lesion prevalence is what a study examining everyone
# for lesions would find at the same ages. All three refer to the same cohort, which is
# what ties the two pathways together.
#
# They were generated by running the model on a known parameter set and rounding
# (see tests/reference_solution.R), so the exercise has a global optimum with an
# error of about 1e-5 and an optimizer that stops well above it has genuinely
# been trapped rather than run out of model.
CALIB.TARGET <- list(
  `Cancer incidence`=list(
    `30-39`=0.000114,
    `40-49`=0.000365,
    `50-59`=0.000867,
    `60-69`=0.001720,
    `70-79`=0.002820,
    `80-89`=0.003980
  ),
  `Lesion prevalence`=list(
    `30-39`=0.0742,
    `40-49`=0.1490,
    `50-59`=0.2510,
    `60-69`=0.3560,
    `70-79`=0.4430,
    `80-89`=0.5090
  ),
  `Late-stage share`=list(
    `30-39`=0.649,
    `40-49`=0.629,
    `50-59`=0.593,
    `60-69`=0.567,
    `70-79`=0.555,
    `80-89`=0.549
  )
)

# The calibration vector is laid out parameter by parameter, and within a
# parameter stratum by stratum. That is the order the app builds its initial
# guess in and the order the calibrated values are reported in, so the initial
# guess below and the constraints on the vector follow it too.
initial.guess.for <- function(params) {
  base <- default.parameters()
  unname(unlist(lapply(params, function(param) {
    sapply(CALIB.STRATA, function(stratum) base[[param]][[stratum]])
  })))
}

CALIB.PARAMS.CORE <- c('p.lgl.onset', 'p.fpl.onset')
CALIB.PARAMS.FULL <- c('p.lgl.onset', 'p.fpl.onset', 'p.hgl.progress')

get.calibration.schemes <- function() {
  return(list(
    natural_history=list(
      description=CALIB.DESCRIPTION.FULL,
      parameters=CALIB.PARAMS.FULL,
      target=CALIB.TARGET,
      strata=CALIB.STRATA,
      initial_guess=initial.guess.for(CALIB.PARAMS.FULL),
      error_function=calibration.error,
      latent_space_training_set=generate.training.dataset,
      other.plots=list(plot.pathway.split)
    ),
    constrained=list(
      description=CALIB.DESCRIPTION.CONSTRAINED,
      parameters=CALIB.PARAMS.FULL,
      target=CALIB.TARGET,
      strata=CALIB.STRATA,
      initial_guess=initial.guess.for(CALIB.PARAMS.FULL),
      error_function=calibration.error,
      constraints=calibration.constraints,
      constraints_description=calibration.constraints.llm,
      latent_space_training_set=generate.constrained.training.dataset,
      other.plots=list(plot.pathway.split)
    )
  ))
}

# Turn a calibration vector into the full parameter list the model runs on. The
# app does this itself before it calls the error function; this copy is what
# tests/calibration_landscape.R evaluates the error on over its own grid.
calib.vector.to.parameters <- function(x, params) {
  pars <- default.parameters()
  i <- 1
  for (param in params) {
    for (stratum in CALIB.STRATA) {
      pars[[param]][[stratum]] <- x[[i]]
      i <- i + 1
    }
  }
  return(pars)
}

# Constraints for the constrained scheme, in the c(x) <= 0 convention. They are
# evaluated on the raw calibration vector, so they read it in the layout
# described above: one block per parameter in the order of CALIB.PARAMS.FULL,
# each block holding one value per stratum.
# They encode what is known about the natural history independently of the
# targets: lesion onset and malignant progression do not decrease with age, and
# the fast pathway stays a minority of all lesions. Together they rule out
# the age-jagged and fast-pathway-dominated fits that match the targets equally well
# but are not clinically credible.
#
# `partial` says the vector only covers the strata a stepwise calibration has
# got through, which it takes one at a time: the same layout over the first
# strata, the rest not searched yet. What that decides is the monotonicity of
# every block over those strata; the fast-pathway share is a ratio of two whole
# blocks and waits for every stratum. The leading constant is satisfied by
# construction and gives the first step, which holds a single stratum and no
# pair to compare, a constraint to model.
calibration.constraints <- function(x, partial=FALSE) {
  x <- as.numeric(x)
  n.strata <- length(x) %/% length(CALIB.PARAMS.FULL)
  block <- function(param) {
    j <- match(param, CALIB.PARAMS.FULL)
    x[intersect(((j - 1) * n.strata + 1):(j * n.strata), seq_along(x))]
  }
  monotonic <- function(param) {
    values <- block(param)
    if (length(values) < 2) return(numeric(0))
    head(values, -1) - tail(values, -1)
  }
  slow <- block('p.lgl.onset')
  fast <- block('p.fpl.onset')
  fast.share <- if (length(slow) == length(CALIB.STRATA) && length(fast) == length(CALIB.STRATA))
    sum(fast) - 0.35 * sum(slow) else numeric(0)
  c(
    if (partial) -1,
    monotonic('p.lgl.onset'),
    monotonic('p.fpl.onset'),
    monotonic('p.hgl.progress'),
    fast.share
  )
}

# The same constraints in words: the app's agentic_constrained wrapper calls
# this with no arguments and appends the string to the agent's system prompt.
calibration.constraints.llm <- function() {
  n <- length(CALIB.STRATA)
  paste0(
    'The vector holds ', paste(CALIB.PARAMS.FULL, collapse=', '), ', in that ',
    'order, each as a block of ', n, ' values running over the age groups ',
    paste(CALIB.STRATA, collapse=', '), '.\n',
    'Every block is non-decreasing: no value is smaller than the one before it.\n',
    'Total p.fpl.onset is at most 35% of total p.lgl.onset.'
  )
}

# Relative squared error, averaged within each target series and summed over
# series. Series are on very different scales (incidence is per mil, the stage
# distribution is a proportion), so the error has to be scale free for all three
# to matter.
calibration.error <- function(pars, target) {
  calibration.strategy <- 'no_screening'
  result <- tryCatch({
    results <- run.simulation(calibration.strategy, pars)
    simulated <- results$outputs[[calibration.strategy]]

    error <- 0
    output <- list()
    for (series in names(target)) {
      target.values <- unlist(target[[series]])
      simulated.values <- simulated[[series]][names(target.values)]
      error <- error + mean(((simulated.values - target.values) / target.values)^2)
      output[[series]] <- simulated.values
    }
    list(error=error, output=output)
  }, error=function(e) {
    output <- lapply(target, function(series) {
      values <- rep(NA_real_, length(series))
      names(values) <- names(series)
      values
    })
    list(error=Inf, output=output)
  })
  return(result)
}

# Training set for the latent space methods: the initial guess scaled by a
# uniform random factor per parameter, kept inside [0, 1] as every calibrated
# parameter is an annual probability.
generate.training.dataset <- function(initial_guess, n, ...) {
  f.pars <- list(...)
  variation <- f.pars$variation
  if (is.null(variation)) variation <- 1

  n_params <- length(initial_guess)
  dataset <- matrix(NA_real_, nrow=n, ncol=n_params)
  for (i in 1:n) {
    factors <- runif(n_params, min=max(0, 1-variation), max=1+variation)
    dataset[i,] <- pmin(1, pmax(0, initial_guess * factors))
  }
  dataset <- dataset[sample(nrow(dataset)),]
  return(dataset)
}

# Same, but only feasible points: the constraints are what the embedding is
# meant to learn, so every row is made monotonic in age and is rejected unless
# the fast pathway stays a minority.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  f.pars <- list(...)
  variation <- f.pars$variation
  if (is.null(variation)) variation <- 1

  n_params <- length(CALIB.PARAMS.FULL)
  n_strata <- length(CALIB.STRATA)
  dataset <- matrix(NA_real_, nrow=n, ncol=n_params*n_strata)

  block <- function(j) ((j-1)*n_strata + 1):(j*n_strata)

  i <- 1
  attempts <- 0
  while (i <= n && attempts < 100*n) {
    attempts <- attempts + 1
    factors <- runif(n_params*n_strata, min=max(0, 1-variation), max=1+variation)
    candidate <- pmin(1, pmax(0, initial_guess * factors))
    # Sorting each parameter over the strata makes the row monotonic in age
    # without discarding it.
    for (j in seq_len(n_params)) candidate[block(j)] <- sort(candidate[block(j)])
    slow <- candidate[block(1)]
    fast <- candidate[block(2)]
    if (sum(fast) > 0.35 * sum(slow)) next
    dataset[i,] <- candidate
    i <- i + 1
  }
  if (i <= n) dataset <- dataset[seq_len(i-1),, drop=FALSE]

  dataset <- dataset[sample(nrow(dataset)),, drop=FALSE]
  return(dataset)
}

# Extra calibration plot: the two lesion onset curves the calibration is
# choosing between, which is where the ambiguity of this model lives.
plot.pathway.split <- function(data) {
  series <- function(param) {
    values <- data[startsWith(names(data), paste0(param, '.'))]
    if (length(values) == 0) return(NULL)
    data.frame(stratum=sub(paste0(param, '.'), '', names(values), fixed=TRUE),
               value=as.numeric(values),
               parameter=param)
  }
  plt.df <- do.call(rbind, Filter(Negate(is.null),
                                  list(series('p.lgl.onset'), series('p.fpl.onset'))))
  if (is.null(plt.df)) return(NULL)

  # The `text` aesthetic is what the calibration tab shows on hover: it renders
  # this plot with ggplotly(tooltip='text'), so a plot that does not map it gets
  # points with an empty hover label. ggplot2 warns about it as an unknown
  # aesthetic when the plot is drawn outside plotly, which is harmless.
  #
  # Significant digits rather than decimal places, one value at a time: fast-pathway
  # onset is an order of magnitude below low-grade lesion onset, so a fixed number of
  # decimals rounds it away and a shared format pads every other label out to the
  # decimals it needs.
  label <- function(x) format(signif(x, 4), scientific=FALSE)
  plt.df$text <- sprintf('%s\n%s: %s', plt.df$stratum, plt.df$parameter,
                         vapply(plt.df$value, label, character(1)))

  # The legend goes on top: ggplotly puts a bottom legend at a fixed fraction of
  # the plot height, which on a plot this short lands on the x axis title.
  plt <- ggplot(plt.df, aes(x=stratum, y=value, group=parameter, colour=parameter,
                            text=text)) +
    geom_point() +
    geom_line() +
    labs(x='Age group', y='Annual probability of lesion onset', colour='') +
    theme_minimal() +
    theme(legend.position='top')
  return(plt)
}

# Scheme descriptions. The app shows these as the label of the scheme, so they
# are short titles; the modelling detail behind each one is in the comments.

# Onset of the two lesion types that drive cancer here, one value per
# ten-year age group. Both count towards lesion prevalence and both end up as
# cancer, but a fast-pathway lesion gets there several times faster, so incidence
# constrains a weighted sum of the onsets and prevalence their unweighted sum:
# the targets leave a family of pathway mixes that fit almost equally well, and
# only the advanced-stage share tells them apart, weakly.
CALIB.DESCRIPTION.CORE <- 'Two-pathway lesion onset'

# Adds the annual probability that a high-grade lesion becomes cancer. Onset and
# progression compensate each other along a ridge in incidence while prevalence
# depends mostly on onset, and faster progression brings diagnoses forward, so
# raising it can improve one age group and worsen the next.
CALIB.DESCRIPTION.FULL <- 'Lesion onset and progression'

# Same parameters as the full scheme, restricted to fits that are non-decreasing
# in age and keep fast-pathway onset a minority (at most 35% of low-grade lesion
# onset), which excludes the age-jagged and fast-pathway-dominated fits.
CALIB.DESCRIPTION.CONSTRAINED <- 'Constrained onset and progression'

# ### TEST
#
# strategies <- sapply(get.strategies(), function(s) s$name)
# pars <- default.parameters()
# results <- run.simulation(strategies, pars)
# print(results$summary)
# print(calibration.error(pars, CALIB.TARGET)$error)
