# Colorectal cancer, two-pathway model

A cost-effectiveness model of colorectal cancer screening built as a harder
calibration test bed than `cea-model-test`. Cancer arises through **two parallel
routes**, which is what the clinical literature describes and what makes the
model difficult to calibrate: the targets constrain the *sum* of what the two
routes produce much better than they constrain the *split* between them, so the
error surface has several distinct local minima instead of one basin.

## Model structure

A Markov cohort model with annual cycles, followed from age 20 to age 89 (70
cycles), with no new entries.

**Conventional adenoma-carcinoma sequence** (slow, most cancers):

    Normal -> low-risk adenoma -> advanced adenoma -> preclinical cancer (localized -> advanced)

**Serrated pathway** (fast, a minority of cancers):

    Normal -> sessile serrated lesion -> preclinical cancer (localized -> advanced)

Low-risk adenomas and serrated lesions regress at a low annual rate; advanced
adenomas regress back to low-risk ones. Each route has its own pair of
preclinical (undiagnosed) cancer states, because serrated cancers progress to an
advanced stage faster and stay asymptomatic longer than conventional ones. This
is what makes the stage distribution at diagnosis carry information about the
pathway mix.

Preclinical cancers surface either through symptoms or through screening, and
become diagnosed cancer at the stage they had reached. Diagnosed patients are
treated and either recover (with a small annual risk of late recurrence and
death) or die of cancer. Everyone alive is exposed to a Gompertz other-cause
mortality that competes with all of the above.

The 13 states are `normal`, `lra`, `aa`, `ssl`, `pre.early.c`, `pre.late.c`,
`pre.early.s`, `pre.late.s`, `clin.early`, `clin.late`, `survivor`, `dead.crc`
and `dead.other`.

Within a cycle, screening is applied first, then costs and utilities are
accrued, then the natural history transition. Rows of the transition matrix that
extreme parameter sets would push above 1 are rescaled, so an optimizer probing
the corners of the search box gets a valid Markov chain rather than an error.

## Strategies

| Strategy | Description |
|---|---|
| `no_screening` | Natural history only. This is the strategy the calibration is run on. |
| `fit_biennial` | Faecal immunochemical test every two years from 50 to 74; positives go on to colonoscopy. |
| `colonoscopy_10y` | Colonoscopy every ten years from 50 to 79. |
| `colonoscopy_45` | The same, starting at 45. |

Detected lesions are removed and the person returns to `normal`; detected
preclinical cancers are diagnosed at their current stage. Serrated lesions are
flat, so both tests miss them far more often than conventional adenomas, and
serrated cancers are detected somewhat less often than conventional ones
(`rr.detection.serrated`). Colonoscopies carry a risk of a serious complication,
with its own cost and disutility. Surveillance after a positive finding is not
modelled.

## Parameters

Age-dependent parameters take one value per stratum. Everything else is a single
value. The `class` field only groups the parameters in the interface.

| Parameter | Base value | Class |
|---|---|---|
| `p.adenoma.onset` | 0.0050 - 0.0370 by age | Natural history (conventional) |
| `p.lra.progress` | 0.015 | Natural history (conventional) |
| `p.lra.regress` | 0.015 | Natural history (conventional) |
| `p.aa.progress` | 0.0053 - 0.0251 by age | Natural history (conventional) |
| `p.aa.regress` | 0.005 | Natural history (conventional) |
| `p.ssl.onset` | 0.00018 - 0.00102 by age | Natural history (serrated) |
| `p.ssl.progress` | 0.040 | Natural history (serrated) |
| `p.ssl.regress` | 0.020 | Natural history (serrated) |
| `p.stage.progress.c`, `p.stage.progress.s` | 0.300, 0.450 | Cancer progression |
| `p.symptomatic.early.c`, `p.symptomatic.late.c` | 0.270, 0.550 | Cancer progression |
| `p.symptomatic.early.s`, `p.symptomatic.late.s` | 0.120, 0.450 | Cancer progression |
| `p.cure.early`, `p.crc.death.early` | 0.550, 0.030 | Cancer survival |
| `p.cure.late`, `p.crc.death.late` | 0.120, 0.250 | Cancer survival |
| `p.survivor.death` | 0.010 | Cancer survival |
| `mortality.other.base`, `mortality.other.rate` | 0.0005, 0.085 | General |
| `adherence.fit`, `spec.fit` | 0.650, 0.950 | Screening (FIT) |
| `sens.fit.lra`, `sens.fit.aa`, `sens.fit.ssl`, `sens.fit.cancer` | 0.040, 0.240, 0.080, 0.550 | Screening (FIT) |
| `adherence.colonoscopy` | 0.550 | Screening (colonoscopy) |
| `sens.colonoscopy.lra`, `.aa`, `.ssl`, `.cancer` | 0.750, 0.920, 0.550, 0.950 | Screening (colonoscopy) |
| `rr.detection.serrated` | 0.850 | Screening (colonoscopy) |
| `p.colonoscopy.complication` | 0.002 | Screening (colonoscopy) |
| `cost.fit`, `cost.colonoscopy`, `cost.polypectomy`, `cost.complication` | 30, 800, 250, 3000 | Costs |
| `cost.treatment.early`, `cost.treatment.late`, `cost.followup` | 18000, 40000, 1200 | Costs |
| `utility.cancer.early`, `utility.cancer.late`, `utility.survivor` | 0.720, 0.500, 0.920 | Utilities |
| `disutility.complication` | 0.020 | Utilities |
| `discount` | 0.030 | General |

Sensitivity for advanced preclinical cancer is derived from the localized one
(times 1.15 for FIT, 1.05 for colonoscopy, capped at 1) rather than being a
parameter of its own.

## Strata

Ten-year age groups from `20-29` to `80-89`.

## Outputs

- `summary`: one row per strategy with the total discounted cost (`C`) and the
  total discounted quality-adjusted life years (`E`) of the cohort.
- `outputs`: per strategy, the three series the calibration uses, each by
  stratum:
  - `CRC incidence`, diagnosed cancers per person-year among those without a
    previous cancer diagnosis (screen-detected cancers included);
  - `Adenoma prevalence`, the proportion of the same population carrying any
    lesion, of either type, as a colonoscopy study would find it;
  - `Late-stage share`, the proportion of the cancers diagnosed in the stratum
    that were regional or distant.
- `incidence`: the incidence series alone, per strategy.
- `cohort.info`: the full cohort trace per strategy, one row per age.

## Calibration

All three schemes calibrate on `no_screening` against the same three target
series, which cover the six age groups from `30-39` to `80-89`. The `20-29`
group is burn-in: the cohort starts with no lesions, so its incidence is
essentially zero whatever the parameters. It is left out of the *targets* rather
than out of the scheme's strata, because the app builds the initial guess from
the base values over every stratum, so a scheme covering only some of them would
be handed a longer vector than it consumes.

| Scheme | Calibrated parameters | Dimension |
|---|---|---|
| `natural_history` | `p.adenoma.onset`, `p.ssl.onset`, `p.aa.progress` | 21 |
| `pathway_mix` | `p.adenoma.onset`, `p.ssl.onset` | 14 |
| `constrained` | the same three, with constraints | 21 |

The calibration vector is laid out parameter by parameter and, within a
parameter, stratum by stratum: the order the app builds its initial guess in
(`unlist(calib.pars[scheme$parameters])`), reports the calibrated values in and
shows them in the Parameters tab. Note that `calib.vector.to.pars()` in
`R/shiny_calibration.R` reads it the other way round, stratum by stratum. The
two orders coincide for a scheme that calibrates a single parameter, which is
what the app has been used with so far, but not for the schemes here: without
the loops in that function swapped, the values land on the wrong parameters and
the search box on each of them comes from the wrong base value.

The error is the mean squared *relative* deviation within each series, summed
over the three series. The series are on very different scales (incidence is per
mil, the other two are proportions), so a scale-free error is what makes all
three count.

The `constrained` scheme adds `vector_to_pars` and `constraints` for the
constrained BO wrappers. Feasible parameter sets are non-decreasing in age for
all three parameters and keep total serrated onset below 35% of total adenoma
onset; its training set generator only produces feasible rows.

The targets were generated by running the model on a known parameter set and
rounding to three significant digits, so the global optimum is worth about 1e-5
and an optimizer that stalls above that has been trapped rather than run out of
model. The reference solution is in `tests/reference_solution.R`, which lives
under `tests/` so that the agentic calibrator, which reads the model source to
build its prompt and skips that directory, does not get handed the answer.

### Why this is hard

- **The two routes are nearly interchangeable.** Both lesion types count towards
  the observed adenoma prevalence and both end up as cancer, so incidence
  constrains a weighted sum of the two onsets and prevalence constrains their
  unweighted sum. Only the share of cancers diagnosed at an advanced stage
  distinguishes them, and it moves by a few points across the whole plausible
  range of mixes. The result is a family of very different natural histories
  that fit almost equally well, with distinct local minima rather than one
  basin.
- **Onset and progression compensate.** Incidence through the conventional route
  depends on the product of `p.adenoma.onset` and `p.aa.progress`, while
  prevalence depends mostly on the first. That is a long curved valley in which
  a simplex crawls.
- **Responses are not monotone.** Serrated lesions drain the pool of people who
  would otherwise develop adenomas, and faster progression empties the pool of
  advanced adenomas and brings diagnoses forward. Raising a parameter can
  therefore improve one age group and worsen the next, which turns per-stratum
  residuals that would be convex into residuals with two roots.
- **Age groups are coupled.** Lesions take years to become cancer, so a
  parameter in one age group moves the observed incidence in the following ones;
  the problem is far from the near-separable one-parameter-per-stratum structure
  of `cea-model-test`.

### Measured behaviour

`tests/calibration_landscape.R` measures this: it runs Nelder-Mead from the base
case and from random points in the search box, reports where each run ends and
how much of the lesion onset it assigns to the serrated route, and compares the
result with a random search given the same number of model evaluations. On the
`natural_history` scheme, with a scope of 100% and 16 starting points, Nelder-
Mead is given about 5000 evaluations per run and ends between 0.0011 and 0.035,
a spread of 31x across starting points and never closer than three orders of
magnitude to the 1.8e-6 of the reference solution. The solutions assign between
2.8% and 6.6% of lesion onset to the serrated route, against 5.3% at the
reference solution: they are different natural histories, not the same one found
with different precision.

Run through the app's own algorithms on the same scheme, with an identity
mapping, 10 initial points and 100 iterations (three seeds):

| Algorithm | Best error found |
|---|---|
| `nm` | 0.0168 - 0.0179 |
| `bo_regular` | 0.0166 - 0.0184 |
| `bo_turbo` | 0.0014 - 0.0043 |
| `bo_baxus` | 0.051 - 0.069 |

At that budget TuRBO is an order of magnitude better than Nelder-Mead, while
plain EI over 21 dimensions has not yet paid off its exploration and BAxUS,
which assumes a low effective dimension, does worse than doing nothing.
