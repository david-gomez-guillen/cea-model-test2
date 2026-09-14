# Cancer natural history model with two competing pathways: a slow pathway
# through two lesion grades and a fast pathway through a single lesion type.
#
# The model is deliberately harder to calibrate than a plain three-state Markov
# chain: the two pathways are only weakly distinguishable from the calibration
# targets, which makes the error surface multimodal. See overview.md.

MODEL.AGE.START <- 20
MODEL.AGE.END <- 89
MODEL.STRATA <- c('20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89')

# Health states. The slow and the fast pathway have their own preclinical
# (undiagnosed) cancer states because fast-pathway cancers progress faster and
# stay asymptomatic for longer, which is what makes the stage distribution at
# diagnosis informative about the pathway mix.
MODEL.STATES <- c(
  'normal',          #  1 no lesion
  'lgl',             #  2 low-grade lesion
  'hgl',             #  3 high-grade lesion
  'fpl',             #  4 fast-pathway lesion
  'pre.early.slow',  #  5 preclinical localized cancer, slow pathway
  'pre.late.slow',   #  6 preclinical advanced cancer, slow pathway
  'pre.early.fast',  #  7 preclinical localized cancer, fast pathway
  'pre.late.fast',   #  8 preclinical advanced cancer, fast pathway
  'clin.early',      #  9 diagnosed localized cancer, under treatment
  'clin.late',       # 10 diagnosed advanced cancer, under treatment
  'survivor',        # 11 post-treatment survivor
  'dead.cancer',     # 12 death from cancer
  'dead.other'       # 13 death from other causes
)

S.NORMAL <- 1; S.LGL <- 2; S.HGL <- 3; S.FPL <- 4
S.PE.SLOW <- 5; S.PL.SLOW <- 6; S.PE.FAST <- 7; S.PL.FAST <- 8
S.CE <- 9; S.CL <- 10; S.SURV <- 11; S.DCAN <- 12; S.DOTH <- 13

# States without a cancer diagnosis, i.e. the population at risk of a first
# cancer. Used as the denominator of incidence and lesion prevalence.
S.AT.RISK <- c(S.NORMAL, S.LGL, S.HGL, S.FPL, S.PE.SLOW, S.PL.SLOW, S.PE.FAST, S.PL.FAST)
S.LESION <- c(S.LGL, S.HGL, S.FPL)

# Artificial delay, in seconds, applied once per simulate() call. The model
# itself runs in a few milliseconds, which is unrealistically fast for anything
# that has to cope with a slow model (progress reporting, calibration budgets);
# this makes a run take a plausible amount of time. Set to 0 to disable.
MODEL.SIMULATION.DELAY <- 0

# Screening schedules, as the ages at which a round is offered.
SCREENING.SCHEDULES <- list(
  no_screening = list(modality = NULL, ages = integer(0)),
  test_biennial = list(modality = 'test', ages = seq(50, 74, by = 2)),
  procedure_10y = list(modality = 'procedure', ages = seq(50, 79, by = 10)),
  procedure_45 = list(modality = 'procedure', ages = seq(45, 79, by = 10))
)

# Base values of the three age-dependent (and calibrated) transition
# probabilities, one per stratum of MODEL.STRATA.
BASE.LGL.ONSET <- c(0.0050, 0.0100, 0.0180, 0.0270, 0.0320, 0.0340, 0.0370)
BASE.FPL.ONSET <- c(0.00018, 0.00035, 0.00056, 0.00077, 0.00088, 0.00095, 0.00102)
BASE.HGL.PROGRESS <- c(0.0053, 0.0064, 0.0094, 0.0139, 0.0188, 0.0225, 0.0251)

# Descriptors of every model parameter, in the shape the shiny app expects from
# get.parameters(). Kept here so that model.R is the single source of truth for
# the base case and can be run standalone.
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

stratum.of.age <- function(age) {
  idx <- min(length(MODEL.STRATA), max(1, (age - MODEL.AGE.START) %/% 10 + 1))
  MODEL.STRATA[idx]
}

# Parameters reach the model either as a single number or, when they are
# stratified (declared as such by the model, split by the user in the Parameters
# tab, or produced by a calibration), as one value per age stratum. This
# resolves both shapes to the value in force at a given age.
par.at.age <- function(value, age) {
  if (!is.list(value) && length(value) == 1) return(as.numeric(value))
  stratum <- stratum.of.age(age)
  nms <- names(value)
  if (!is.null(nms) && stratum %in% nms) return(as.numeric(value[[stratum]]))
  idx <- min(length(value), max(1, (age - MODEL.AGE.START) %/% 10 + 1))
  as.numeric(value[[idx]])
}

# Every parameter of the model is allowed to arrive stratified, because the user
# can split any of them by stratum in the Parameters tab. Resolving the whole
# list at the start of each cycle keeps the rest of the model working on plain
# numbers.
parameters.at.age <- function(pars, age) {
  lapply(pars, par.at.age, age = age)
}

# Annual probability of dying from something other than the modelled cancer,
# Gompertz in age.
other.cause.mortality <- function(age, base, rate) {
  min(1, base * exp(rate * (age - MODEL.AGE.START)))
}

# Build the one-cycle natural history transition matrix at a given age.
# Off-diagonal probabilities are clamped to [0, 1] and rows are rescaled when
# they would exceed 1, so that extreme parameter sets proposed by an optimizer
# still yield a valid Markov chain instead of an error.
transition.matrix <- function(age, p) {
  n <- length(MODEL.STATES)
  tp <- matrix(0, nrow = n, ncol = n)

  mo <- other.cause.mortality(age, p$mortality.other.base, p$mortality.other.rate)

  tp[S.NORMAL, S.LGL] <- p$p.lgl.onset
  tp[S.NORMAL, S.FPL] <- p$p.fpl.onset

  tp[S.LGL, S.HGL] <- p$p.lgl.progress
  tp[S.LGL, S.NORMAL] <- p$p.lgl.regress

  tp[S.HGL, S.PE.SLOW] <- p$p.hgl.progress
  tp[S.HGL, S.LGL] <- p$p.hgl.regress

  tp[S.FPL, S.PE.FAST] <- p$p.fpl.progress
  tp[S.FPL, S.NORMAL] <- p$p.fpl.regress

  tp[S.PE.SLOW, S.CE] <- p$p.symptomatic.early.slow
  tp[S.PE.SLOW, S.PL.SLOW] <- p$p.stage.progress.slow
  tp[S.PL.SLOW, S.CL] <- p$p.symptomatic.late.slow

  tp[S.PE.FAST, S.CE] <- p$p.symptomatic.early.fast
  tp[S.PE.FAST, S.PL.FAST] <- p$p.stage.progress.fast
  tp[S.PL.FAST, S.CL] <- p$p.symptomatic.late.fast

  tp[S.CE, S.SURV] <- p$p.cure.early
  tp[S.CE, S.DCAN] <- p$p.cancer.death.early
  tp[S.CL, S.SURV] <- p$p.cure.late
  tp[S.CL, S.DCAN] <- p$p.cancer.death.late
  tp[S.SURV, S.DCAN] <- p$p.survivor.death

  tp <- pmax(tp, 0)

  # Everyone alive is exposed to the same competing risk of other-cause death.
  alive <- setdiff(seq_len(n), c(S.DCAN, S.DOTH))
  tp[alive, S.DOTH] <- mo

  off <- rowSums(tp)
  excess <- off > 1
  if (any(excess)) tp[excess, ] <- tp[excess, ] / off[excess]

  diag(tp) <- 1 - rowSums(tp)
  tp[S.DCAN, S.DCAN] <- 1
  tp[S.DOTH, S.DOTH] <- 1

  return(tp)
}

# Per-state probability that a screening round detects the lesion or cancer a
# person is carrying, and the number of procedures the round generates.
# A test round sends positives to the procedure, so detection is the product of
# the two sensitivities; a procedure round only needs its own sensitivity. Late
# preclinical cancers are somewhat easier to detect than localized ones, and
# fast-pathway lesions are harder to detect, which is why they are missed more
# often.
screening.detection <- function(modality, p) {
  detect <- rep(0, length(MODEL.STATES))
  proc <- rep(0, length(MODEL.STATES))

  sens.proc <- rep(0, length(MODEL.STATES))
  sens.proc[S.LGL] <- p$sens.procedure.lgl
  sens.proc[S.HGL] <- p$sens.procedure.hgl
  sens.proc[S.FPL] <- p$sens.procedure.fpl
  sens.proc[c(S.PE.SLOW, S.PE.FAST)] <- p$sens.procedure.cancer
  sens.proc[c(S.PL.SLOW, S.PL.FAST)] <- min(1, p$sens.procedure.cancer * 1.05)
  # Fast-pathway cancers are harder to detect, so the procedure misses them more
  # often than slow-pathway ones.
  sens.proc[c(S.PE.FAST, S.PL.FAST)] <- sens.proc[c(S.PE.FAST, S.PL.FAST)] * p$rr.detection.fast

  if (identical(modality, 'test')) {
    sens.test <- rep(0, length(MODEL.STATES))
    sens.test[S.LGL] <- p$sens.test.lgl
    sens.test[S.HGL] <- p$sens.test.hgl
    sens.test[S.FPL] <- p$sens.test.fpl
    sens.test[c(S.PE.SLOW, S.PE.FAST)] <- p$sens.test.cancer
    sens.test[c(S.PL.SLOW, S.PL.FAST)] <- min(1, p$sens.test.cancer * 1.15)
    sens.test[c(S.PE.FAST, S.PL.FAST)] <- sens.test[c(S.PE.FAST, S.PL.FAST)] * p$rr.detection.fast

    positive <- sens.test
    # People without a lesion test positive at the complement of specificity.
    positive[S.NORMAL] <- 1 - p$spec.test
    positive[S.SURV] <- 1 - p$spec.test

    detect <- p$adherence.test * positive * sens.proc
    proc <- p$adherence.test * positive
  } else if (identical(modality, 'procedure')) {
    detect <- p$adherence.procedure * sens.proc
    proc <- rep(p$adherence.procedure, length(MODEL.STATES))
    proc[c(S.CE, S.CL, S.DCAN, S.DOTH)] <- 0
  }

  list(detect = pmin(1, detect), procedures = pmin(1, proc))
}

# Apply one screening round to the cohort. Detected lesions are removed (back to
# normal) and detected preclinical cancers are diagnosed at their current stage.
apply.screening <- function(cohort, modality, p) {
  sc <- screening.detection(modality, p)
  detect <- sc$detect

  removed <- cohort[S.LESION] * detect[S.LESION]
  found.early <- cohort[S.PE.SLOW] * detect[S.PE.SLOW] + cohort[S.PE.FAST] * detect[S.PE.FAST]
  found.late <- cohort[S.PL.SLOW] * detect[S.PL.SLOW] + cohort[S.PL.FAST] * detect[S.PL.FAST]

  cohort[S.LESION] <- cohort[S.LESION] - removed
  cohort[S.NORMAL] <- cohort[S.NORMAL] + sum(removed)
  cohort[S.PE.SLOW] <- cohort[S.PE.SLOW] * (1 - detect[S.PE.SLOW])
  cohort[S.PE.FAST] <- cohort[S.PE.FAST] * (1 - detect[S.PE.FAST])
  cohort[S.PL.SLOW] <- cohort[S.PL.SLOW] * (1 - detect[S.PL.SLOW])
  cohort[S.PL.FAST] <- cohort[S.PL.FAST] * (1 - detect[S.PL.FAST])
  cohort[S.CE] <- cohort[S.CE] + found.early
  cohort[S.CL] <- cohort[S.CL] + found.late

  n.procedures <- sum(cohort * sc$procedures)
  n.removals <- sum(removed)
  # Only people who are alive and not currently being treated for a cancer are
  # invited to the round.
  n.tests <- if (identical(modality, 'test')) {
    p$adherence.test * sum(cohort[c(S.AT.RISK, S.SURV)])
  } else 0

  list(cohort = cohort,
       n.tests = n.tests,
       n.procedures = n.procedures,
       n.removals = n.removals,
       found.early = found.early,
       found.late = found.late)
}

simulate <- function(strategies, pars, delay = MODEL.SIMULATION.DELAY) {
  if (is.numeric(delay) && length(delay) == 1 && !is.na(delay) && delay > 0) {
    Sys.sleep(delay)
  }

  ages <- seq(MODEL.AGE.START, MODEL.AGE.END)
  n.states <- length(MODEL.STATES)

  # Resolved once so that the discount curve stays well defined even if the user
  # splits the rate by stratum; everything else is resolved cycle by cycle.
  discount <- par.at.age(pars$discount, MODEL.AGE.START)

  summary.df <- data.frame()
  outputs <- list()
  cohort.states <- list()

  for (strategy in strategies) {
    schedule <- SCREENING.SCHEDULES[[strategy]]
    if (is.null(schedule)) stop(sprintf("Unknown strategy '%s'", strategy))

    cohort <- rep(0, n.states)
    cohort[S.NORMAL] <- 1
    trace <- matrix(NA_real_, nrow = length(ages), ncol = n.states,
                    dimnames = list(as.character(ages), MODEL.STATES))

    costs <- numeric(length(ages))
    qalys <- numeric(length(ages))
    at.risk <- numeric(length(ages))
    lesion.prevalence <- numeric(length(ages))
    new.early <- numeric(length(ages))
    new.late <- numeric(length(ages))

    for (i in seq_along(ages)) {
      age <- ages[i]
      p <- parameters.at.age(pars, age)
      cycle.cost <- 0
      cycle.disutility <- 0
      screen.early <- 0
      screen.late <- 0

      if (age %in% schedule$ages) {
        sc <- apply.screening(cohort, schedule$modality, p)
        cohort <- sc$cohort
        cycle.cost <- cycle.cost +
          sc$n.tests * p$cost.test +
          sc$n.procedures * p$cost.procedure +
          sc$n.removals * p$cost.removal +
          sc$n.procedures * p$p.procedure.complication * p$cost.complication
        cycle.disutility <- cycle.disutility +
          sc$n.procedures * p$p.procedure.complication * p$disutility.complication
        screen.early <- sc$found.early
        screen.late <- sc$found.late
      }

      trace[i, ] <- cohort

      state.costs <- rep(0, n.states)
      state.costs[S.CE] <- p$cost.treatment.early
      state.costs[S.CL] <- p$cost.treatment.late
      state.costs[S.SURV] <- p$cost.followup

      state.utilities <- rep(0, n.states)
      state.utilities[c(S.NORMAL, S.LGL, S.HGL, S.FPL,
                        S.PE.SLOW, S.PL.SLOW, S.PE.FAST, S.PL.FAST)] <- 1
      state.utilities[S.CE] <- p$utility.cancer.early
      state.utilities[S.CL] <- p$utility.cancer.late
      state.utilities[S.SURV] <- p$utility.survivor

      df <- 1 / (1 + discount)^(age - MODEL.AGE.START)
      costs[i] <- (sum(state.costs * cohort) + cycle.cost) * df
      qalys[i] <- (sum(state.utilities * cohort) - cycle.disutility) * df

      tp <- transition.matrix(age, p)

      # Registry incidence counts every cancer diagnosed during the cycle,
      # whether it surfaced through symptoms or was found by screening.
      new.early[i] <- cohort[S.PE.SLOW] * tp[S.PE.SLOW, S.CE] +
                      cohort[S.PE.FAST] * tp[S.PE.FAST, S.CE] + screen.early
      new.late[i] <- cohort[S.PL.SLOW] * tp[S.PL.SLOW, S.CL] +
                     cohort[S.PL.FAST] * tp[S.PL.FAST, S.CL] + screen.late
      at.risk[i] <- sum(cohort[S.AT.RISK])
      lesion.prevalence[i] <- sum(cohort[S.LESION])

      cohort <- as.numeric(cohort %*% tp)
    }

    strata.of.age <- sapply(ages, stratum.of.age)
    incidence <- tapply(seq_along(ages), factor(strata.of.age, levels = MODEL.STRATA),
                        function(idx) {
                          denom <- at.risk[idx]
                          mean(ifelse(denom > 1e-12, (new.early[idx] + new.late[idx]) / denom, 0))
                        })
    prevalence <- tapply(seq_along(ages), factor(strata.of.age, levels = MODEL.STRATA),
                         function(idx) {
                           denom <- at.risk[idx]
                           mean(ifelse(denom > 1e-12, lesion.prevalence[idx] / denom, 0))
                         })
    late.share <- tapply(seq_along(ages), factor(strata.of.age, levels = MODEL.STRATA),
                         function(idx) {
                           total <- sum(new.early[idx] + new.late[idx])
                           if (total > 1e-12) sum(new.late[idx]) / total else 0
                         })

    outputs[[strategy]] <- list(
      `Cancer incidence` = setNames(as.numeric(incidence), MODEL.STRATA),
      `Lesion prevalence` = setNames(as.numeric(prevalence), MODEL.STRATA),
      `Late-stage share` = setNames(as.numeric(late.share), MODEL.STRATA)
    )
    cohort.states[[strategy]] <- trace

    summary.df <- rbind(summary.df, data.frame(
      strategy = strategy,
      C = sum(costs),
      E = sum(qalys)
    ))
  }

  return(list(
    summary = summary.df,
    outputs = outputs,
    incidence = lapply(outputs, function(o) o$`Cancer incidence`),
    cohort.info = cohort.states
  ))
}
