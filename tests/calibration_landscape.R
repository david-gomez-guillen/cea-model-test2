# Probe of the calibration landscape.
#
# Runs Nelder-Mead from the base case and from random starting points inside the
# box the app searches (the initial guess scaled by 1 +/- the calibration scope)
# and compares the results with a random search given the same number of model
# evaluations. Prints, for every solution found, the share of lesion onset that
# the solution assigns to the serrated pathway, which is what separates the
# basins of this model.
#
# Run from the repository root:
#   Rscript tests/calibration_landscape.R [scheme] [n_starts] [scope]
#
# e.g. Rscript tests/calibration_landscape.R natural_history 24 1

source('shiny_interface.R')
source('tests/reference_solution.R')

args <- commandArgs(trailingOnly=TRUE)
scheme.name <- if (length(args) >= 1) args[1] else 'natural_history'
n.starts <- if (length(args) >= 2) as.integer(args[2]) else 16
scope <- if (length(args) >= 3) as.numeric(args[3]) else 1

set.seed(20250813)

scheme <- get.calibration.schemes()[[scheme.name]]
params <- scheme$parameters
initial.guess <- scheme$initial_guess
lo <- pmax(0, initial.guess * (1 - scope))
hi <- pmin(1, initial.guess * (1 + scope))

n.evaluations <- 0
objective <- function(x) {
  n.evaluations <<- n.evaluations + 1
  x <- pmin(hi, pmax(lo, x))
  scheme$error_function(calib.vector.to.parameters(x, params), scheme$target)$error
}

# Share of total lesion onset assigned to the serrated pathway, which is where
# the ambiguity of this model lives.
serrated.share <- function(x) {
  block <- function(param) {
    j <- which(params == param)
    ((j-1)*length(CALIB.STRATA) + 1):(j*length(CALIB.STRATA))
  }
  adenoma <- sum(x[block('p.adenoma.onset')])
  ssl <- sum(x[block('p.ssl.onset')])
  ssl / (ssl + adenoma)
}

cat(sprintf('scheme: %s (%d parameters, scope %.0f%%)\n', scheme.name, length(initial.guess), 100*scope))
cat(sprintf('error at the reference solution: %.3g\n',
            objective(reference.vector(params))))
cat(sprintf('error at the base case (the initial guess): %.4g\n', objective(initial.guess)))
cat(sprintf('serrated onset share: reference %.3f, base case %.3f\n\n',
            serrated.share(reference.vector(params)), serrated.share(initial.guess)))

cat(sprintf('Nelder-Mead from %d starting points (the first is the base case)\n', n.starts))
solutions <- list()
for (i in seq_len(n.starts)) {
  x0 <- if (i == 1) initial.guess else runif(length(initial.guess), lo, hi)
  before <- n.evaluations
  fit <- optim(x0, objective, method='Nelder-Mead',
               control=list(maxit=5000, reltol=1e-12))
  solutions[[i]] <- list(error=fit$value,
                         x=pmin(hi, pmax(lo, fit$par)),
                         evaluations=n.evaluations - before)
}

errors <- sapply(solutions, function(s) s$error)
shares <- sapply(solutions, function(s) serrated.share(s$x))
budget <- round(mean(sapply(solutions, function(s) s$evaluations)))

report <- data.frame(start=ifelse(seq_along(errors) == 1, 'base case', 'random'),
                     error=signif(errors, 4),
                     serrated.share=round(shares, 3))
print(report[order(report$error),], row.names=FALSE)

cat(sprintf('\nbest %.4g / median %.4g / worst %.4g, %d evaluations per run on average\n',
            min(errors), median(errors), max(errors), budget))
cat(sprintf('spread between the best and the median local optimum: %.1fx\n',
            median(errors)/min(errors)))

random.errors <- replicate(budget, objective(runif(length(initial.guess), lo, hi)))
cat(sprintf('best of a random search with the same budget: %.4g\n', min(random.errors)))
cat(sprintf('total model evaluations: %d\n', n.evaluations))
