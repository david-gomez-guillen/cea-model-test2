# Colorectal cancer natural history model with two competing carcinogenesis
# pathways (conventional adenoma-carcinoma sequence and serrated pathway).
#
# The model is deliberately harder to calibrate than a plain three-state Markov
# chain: the two pathways are only weakly distinguishable from the calibration
# targets, which makes the error surface multimodal. See overview.md.

MODEL.AGE.START <- 20
MODEL.AGE.END <- 89
MODEL.STRATA <- c('20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89')

# Health states. The conventional and the serrated route have their own
# preclinical (undiagnosed) cancer states because serrated cancers progress
# faster and stay asymptomatic for longer, which is what makes the stage
# distribution at diagnosis informative about the pathway mix.
MODEL.STATES <- c(
  'normal',        #  1 no lesion
  'lra',           #  2 low-risk (non-advanced) adenoma
  'aa',            #  3 advanced adenoma
  'ssl',           #  4 sessile serrated lesion
  'pre.early.c',   #  5 preclinical localized cancer, conventional pathway
  'pre.late.c',    #  6 preclinical regional/distant cancer, conventional pathway
  'pre.early.s',   #  7 preclinical localized cancer, serrated pathway
  'pre.late.s',    #  8 preclinical regional/distant cancer, serrated pathway
  'clin.early',    #  9 diagnosed localized cancer, under treatment
  'clin.late',     # 10 diagnosed regional/distant cancer, under treatment
  'survivor',      # 11 post-treatment survivor
  'dead.crc',      # 12 death from colorectal cancer
  'dead.other'     # 13 death from other causes
)

S.NORMAL <- 1; S.LRA <- 2; S.AA <- 3; S.SSL <- 4
S.PEC <- 5; S.PLC <- 6; S.PES <- 7; S.PLS <- 8
S.CE <- 9; S.CL <- 10; S.SURV <- 11; S.DCRC <- 12; S.DOTH <- 13

# States without a colorectal cancer diagnosis, i.e. the population at risk of a
# first cancer. Used as the denominator of incidence and lesion prevalence.
S.AT.RISK <- c(S.NORMAL, S.LRA, S.AA, S.SSL, S.PEC, S.PLC, S.PES, S.PLS)
S.LESION <- c(S.LRA, S.AA, S.SSL)

# Artificial delay, in seconds, applied once per simulate() call. The model
# itself runs in a few milliseconds, which is unrealistically fast for anything
# that has to cope with a slow model (progress reporting, calibration budgets);
# this makes a run take a plausible amount of time. Set to 0 to disable.
MODEL.SIMULATION.DELAY <- 0

# Screening schedules, as the ages at which a round is offered.
SCREENING.SCHEDULES <- list(
  no_screening = list(test = NULL, ages = integer(0)),
  fit_biennial = list(test = 'fit', ages = seq(50, 74, by = 2)),
  colonoscopy_10y = list(test = 'colonoscopy', ages = seq(50, 79, by = 10)),
  colonoscopy_45 = list(test = 'colonoscopy', ages = seq(45, 79, by = 10))
)

# Base values of the three age-dependent (and calibrated) transition
# probabilities, one per stratum of MODEL.STRATA.
BASE.ADENOMA.ONSET <- c(0.0050, 0.0100, 0.0180, 0.0270, 0.0320, 0.0340, 0.0370)
BASE.SSL.ONSET <- c(0.00018, 0.00035, 0.00056, 0.00077, 0.00088, 0.00095, 0.00102)
BASE.AA.PROGRESS <- c(0.0053, 0.0064, 0.0094, 0.0139, 0.0188, 0.0225, 0.0251)

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
  stratified.parameter('p.adenoma.onset',
                       'Annual probability of developing a low-risk adenoma',
                       BASE.ADENOMA.ONSET, 'Natural history (conventional)'),
  stratified.parameter('p.ssl.onset',
                       'Annual probability of developing a sessile serrated lesion',
                       BASE.SSL.ONSET, 'Natural history (serrated)'),
  stratified.parameter('p.aa.progress',
                       'Annual probability that an advanced adenoma becomes cancer',
                       BASE.AA.PROGRESS, 'Natural history (conventional)'),
  list(
    list(name = 'p.lra.progress',
         display.name = 'Annual probability that a low-risk adenoma becomes advanced',
         base.value = 0.015, class = 'Natural history (conventional)'),
    list(name = 'p.lra.regress',
         display.name = 'Annual probability that a low-risk adenoma regresses',
         base.value = 0.015, class = 'Natural history (conventional)'),
    list(name = 'p.aa.regress',
         display.name = 'Annual probability that an advanced adenoma regresses',
         base.value = 0.005, class = 'Natural history (conventional)'),
    list(name = 'p.ssl.progress',
         display.name = 'Annual probability that a serrated lesion becomes cancer',
         base.value = 0.040, class = 'Natural history (serrated)'),
    list(name = 'p.ssl.regress',
         display.name = 'Annual probability that a serrated lesion regresses',
         base.value = 0.020, class = 'Natural history (serrated)'),
    list(name = 'p.stage.progress.c',
         display.name = 'Annual probability of localized to advanced stage, conventional cancer',
         base.value = 0.300, class = 'Cancer progression'),
    list(name = 'p.stage.progress.s',
         display.name = 'Annual probability of localized to advanced stage, serrated cancer',
         base.value = 0.450, class = 'Cancer progression'),
    list(name = 'p.symptomatic.early.c',
         display.name = 'Annual probability of symptomatic diagnosis, localized conventional cancer',
         base.value = 0.270, class = 'Cancer progression'),
    list(name = 'p.symptomatic.late.c',
         display.name = 'Annual probability of symptomatic diagnosis, advanced conventional cancer',
         base.value = 0.550, class = 'Cancer progression'),
    list(name = 'p.symptomatic.early.s',
         display.name = 'Annual probability of symptomatic diagnosis, localized serrated cancer',
         base.value = 0.120, class = 'Cancer progression'),
    list(name = 'p.symptomatic.late.s',
         display.name = 'Annual probability of symptomatic diagnosis, advanced serrated cancer',
         base.value = 0.450, class = 'Cancer progression'),
    list(name = 'p.cure.early',
         display.name = 'Annual probability of cure after a localized cancer diagnosis',
         base.value = 0.550, class = 'Cancer survival'),
    list(name = 'p.crc.death.early',
         display.name = 'Annual probability of cancer death, localized cancer',
         base.value = 0.030, class = 'Cancer survival'),
    list(name = 'p.cure.late',
         display.name = 'Annual probability of cure after an advanced cancer diagnosis',
         base.value = 0.120, class = 'Cancer survival'),
    list(name = 'p.crc.death.late',
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
    list(name = 'adherence.fit',
         display.name = 'Proportion attending an invitation to FIT',
         base.value = 0.650, class = 'Screening (FIT)'),
    list(name = 'sens.fit.lra',
         display.name = 'FIT sensitivity for low-risk adenoma',
         base.value = 0.040, class = 'Screening (FIT)'),
    list(name = 'sens.fit.aa',
         display.name = 'FIT sensitivity for advanced adenoma',
         base.value = 0.240, class = 'Screening (FIT)'),
    list(name = 'sens.fit.ssl',
         display.name = 'FIT sensitivity for sessile serrated lesion',
         base.value = 0.080, class = 'Screening (FIT)'),
    list(name = 'sens.fit.cancer',
         display.name = 'FIT sensitivity for localized cancer',
         base.value = 0.550, class = 'Screening (FIT)'),
    list(name = 'spec.fit',
         display.name = 'FIT specificity',
         base.value = 0.950, class = 'Screening (FIT)'),
    list(name = 'adherence.colonoscopy',
         display.name = 'Proportion attending an invitation to colonoscopy',
         base.value = 0.550, class = 'Screening (colonoscopy)'),
    list(name = 'sens.colonoscopy.lra',
         display.name = 'Colonoscopy sensitivity for low-risk adenoma',
         base.value = 0.750, class = 'Screening (colonoscopy)'),
    list(name = 'sens.colonoscopy.aa',
         display.name = 'Colonoscopy sensitivity for advanced adenoma',
         base.value = 0.920, class = 'Screening (colonoscopy)'),
    list(name = 'sens.colonoscopy.ssl',
         display.name = 'Colonoscopy sensitivity for sessile serrated lesion',
         base.value = 0.550, class = 'Screening (colonoscopy)'),
    list(name = 'sens.colonoscopy.cancer',
         display.name = 'Colonoscopy sensitivity for localized cancer',
         base.value = 0.950, class = 'Screening (colonoscopy)'),
    list(name = 'rr.detection.serrated',
         display.name = 'Relative detection of serrated versus conventional cancer',
         base.value = 0.850, class = 'Screening (colonoscopy)'),
    list(name = 'p.colonoscopy.complication',
         display.name = 'Probability of a serious complication per colonoscopy',
         base.value = 0.002, class = 'Screening (colonoscopy)'),
    list(name = 'cost.fit',
         display.name = 'Cost per FIT performed',
         base.value = 30, class = 'Costs'),
    list(name = 'cost.colonoscopy',
         display.name = 'Cost per colonoscopy performed',
         base.value = 800, class = 'Costs'),
    list(name = 'cost.polypectomy',
         display.name = 'Additional cost per lesion removed',
         base.value = 250, class = 'Costs'),
    list(name = 'cost.complication',
         display.name = 'Cost per serious colonoscopy complication',
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
         display.name = 'Disutility of a serious colonoscopy complication',
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

# Annual probability of dying from something other than colorectal cancer,
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

  tp[S.NORMAL, S.LRA] <- p$p.adenoma.onset
  tp[S.NORMAL, S.SSL] <- p$p.ssl.onset

  tp[S.LRA, S.AA] <- p$p.lra.progress
  tp[S.LRA, S.NORMAL] <- p$p.lra.regress

  tp[S.AA, S.PEC] <- p$p.aa.progress
  tp[S.AA, S.LRA] <- p$p.aa.regress

  tp[S.SSL, S.PES] <- p$p.ssl.progress
  tp[S.SSL, S.NORMAL] <- p$p.ssl.regress

  tp[S.PEC, S.CE] <- p$p.symptomatic.early.c
  tp[S.PEC, S.PLC] <- p$p.stage.progress.c
  tp[S.PLC, S.CL] <- p$p.symptomatic.late.c

  tp[S.PES, S.CE] <- p$p.symptomatic.early.s
  tp[S.PES, S.PLS] <- p$p.stage.progress.s
  tp[S.PLS, S.CL] <- p$p.symptomatic.late.s

  tp[S.CE, S.SURV] <- p$p.cure.early
  tp[S.CE, S.DCRC] <- p$p.crc.death.early
  tp[S.CL, S.SURV] <- p$p.cure.late
  tp[S.CL, S.DCRC] <- p$p.crc.death.late
  tp[S.SURV, S.DCRC] <- p$p.survivor.death

  tp <- pmax(tp, 0)

  # Everyone alive is exposed to the same competing risk of other-cause death.
  alive <- setdiff(seq_len(n), c(S.DCRC, S.DOTH))
  tp[alive, S.DOTH] <- mo

  off <- rowSums(tp)
  excess <- off > 1
  if (any(excess)) tp[excess, ] <- tp[excess, ] / off[excess]

  diag(tp) <- 1 - rowSums(tp)
  tp[S.DCRC, S.DCRC] <- 1
  tp[S.DOTH, S.DOTH] <- 1

  return(tp)
}

# Per-state probability that a screening round detects the lesion or cancer a
# person is carrying, and the number of colonoscopies the round generates.
# A FIT round sends positives to colonoscopy, so detection is the product of the
# two sensitivities; a colonoscopy round only needs its own sensitivity. Late
# preclinical cancers are somewhat easier to detect than localized ones, and
# serrated lesions are flat, which is why they are missed more often.
screening.detection <- function(test, p) {
  detect <- rep(0, length(MODEL.STATES))
  colo <- rep(0, length(MODEL.STATES))

  sens.colo <- rep(0, length(MODEL.STATES))
  sens.colo[S.LRA] <- p$sens.colonoscopy.lra
  sens.colo[S.AA] <- p$sens.colonoscopy.aa
  sens.colo[S.SSL] <- p$sens.colonoscopy.ssl
  sens.colo[c(S.PEC, S.PES)] <- p$sens.colonoscopy.cancer
  sens.colo[c(S.PLC, S.PLS)] <- min(1, p$sens.colonoscopy.cancer * 1.05)
  # Serrated cancers are right-sided and flat-based, so a colonoscopy misses
  # them more often than conventional ones.
  sens.colo[c(S.PES, S.PLS)] <- sens.colo[c(S.PES, S.PLS)] * p$rr.detection.serrated

  if (identical(test, 'fit')) {
    sens.fit <- rep(0, length(MODEL.STATES))
    sens.fit[S.LRA] <- p$sens.fit.lra
    sens.fit[S.AA] <- p$sens.fit.aa
    sens.fit[S.SSL] <- p$sens.fit.ssl
    sens.fit[c(S.PEC, S.PES)] <- p$sens.fit.cancer
    sens.fit[c(S.PLC, S.PLS)] <- min(1, p$sens.fit.cancer * 1.15)
    sens.fit[c(S.PES, S.PLS)] <- sens.fit[c(S.PES, S.PLS)] * p$rr.detection.serrated

    positive <- sens.fit
    # People without a lesion test positive at the complement of specificity.
    positive[S.NORMAL] <- 1 - p$spec.fit
    positive[S.SURV] <- 1 - p$spec.fit

    detect <- p$adherence.fit * positive * sens.colo
    colo <- p$adherence.fit * positive
  } else if (identical(test, 'colonoscopy')) {
    detect <- p$adherence.colonoscopy * sens.colo
    colo <- rep(p$adherence.colonoscopy, length(MODEL.STATES))
    colo[c(S.CE, S.CL, S.DCRC, S.DOTH)] <- 0
  }

  list(detect = pmin(1, detect), colonoscopies = pmin(1, colo))
}

# Apply one screening round to the cohort. Detected lesions are removed
# (polypectomy, back to normal) and detected preclinical cancers are diagnosed
# at their current stage.
apply.screening <- function(cohort, test, p) {
  sc <- screening.detection(test, p)
  detect <- sc$detect

  removed <- cohort[S.LESION] * detect[S.LESION]
  found.early <- cohort[S.PEC] * detect[S.PEC] + cohort[S.PES] * detect[S.PES]
  found.late <- cohort[S.PLC] * detect[S.PLC] + cohort[S.PLS] * detect[S.PLS]

  cohort[S.LESION] <- cohort[S.LESION] - removed
  cohort[S.NORMAL] <- cohort[S.NORMAL] + sum(removed)
  cohort[S.PEC] <- cohort[S.PEC] * (1 - detect[S.PEC])
  cohort[S.PES] <- cohort[S.PES] * (1 - detect[S.PES])
  cohort[S.PLC] <- cohort[S.PLC] * (1 - detect[S.PLC])
  cohort[S.PLS] <- cohort[S.PLS] * (1 - detect[S.PLS])
  cohort[S.CE] <- cohort[S.CE] + found.early
  cohort[S.CL] <- cohort[S.CL] + found.late

  n.colonoscopies <- sum(cohort * sc$colonoscopies)
  n.polypectomies <- sum(removed)
  # Only people who are alive and not currently being treated for a cancer are
  # invited to the round.
  n.tests <- if (identical(test, 'fit')) {
    p$adherence.fit * sum(cohort[c(S.AT.RISK, S.SURV)])
  } else 0

  list(cohort = cohort,
       n.tests = n.tests,
       n.colonoscopies = n.colonoscopies,
       n.polypectomies = n.polypectomies,
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
        sc <- apply.screening(cohort, schedule$test, p)
        cohort <- sc$cohort
        cycle.cost <- cycle.cost +
          sc$n.tests * p$cost.fit +
          sc$n.colonoscopies * p$cost.colonoscopy +
          sc$n.polypectomies * p$cost.polypectomy +
          sc$n.colonoscopies * p$p.colonoscopy.complication * p$cost.complication
        cycle.disutility <- cycle.disutility +
          sc$n.colonoscopies * p$p.colonoscopy.complication * p$disutility.complication
        screen.early <- sc$found.early
        screen.late <- sc$found.late
      }

      trace[i, ] <- cohort

      state.costs <- rep(0, n.states)
      state.costs[S.CE] <- p$cost.treatment.early
      state.costs[S.CL] <- p$cost.treatment.late
      state.costs[S.SURV] <- p$cost.followup

      state.utilities <- rep(0, n.states)
      state.utilities[c(S.NORMAL, S.LRA, S.AA, S.SSL,
                        S.PEC, S.PLC, S.PES, S.PLS)] <- 1
      state.utilities[S.CE] <- p$utility.cancer.early
      state.utilities[S.CL] <- p$utility.cancer.late
      state.utilities[S.SURV] <- p$utility.survivor

      df <- 1 / (1 + discount)^(age - MODEL.AGE.START)
      costs[i] <- (sum(state.costs * cohort) + cycle.cost) * df
      qalys[i] <- (sum(state.utilities * cohort) - cycle.disutility) * df

      tp <- transition.matrix(age, p)

      # Registry incidence counts every cancer diagnosed during the cycle,
      # whether it surfaced through symptoms or was found by screening.
      new.early[i] <- cohort[S.PEC] * tp[S.PEC, S.CE] +
                      cohort[S.PES] * tp[S.PES, S.CE] + screen.early
      new.late[i] <- cohort[S.PLC] * tp[S.PLC, S.CL] +
                     cohort[S.PLS] * tp[S.PLS, S.CL] + screen.late
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
      `CRC incidence` = setNames(as.numeric(incidence), MODEL.STRATA),
      `Adenoma prevalence` = setNames(as.numeric(prevalence), MODEL.STRATA),
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
    incidence = lapply(outputs, function(o) o$`CRC incidence`),
    cohort.info = cohort.states
  ))
}
