################################################################################
# Project: Suicide rates by sex and district in Greenland, 1983-2023
# Author: Anonymized for peer review
# Last updated: 30.10.26
################################################################################

# Setup
rm(list = ls())
library(brms)

setwd("C:/Users/beadl/Documents/Github/glsts")
set.seed(82)

# Reading StatBank files, where "-" means zero
rd <- function(f, nm) {
  d <- read.csv(f, fileEncoding = "latin1", na.strings = "-", col.names = nm)
  d[is.na(d)] <- 0
  d
}

# Loading raw data sets
d_old <- rd("BEXLD1_20261022-190132.csv",   # Suicide, 1983-2010
            c("year", "district", "cause", "sex", "age", "n"))
d_new <- rd("SUXLDA1A_20261022-190231.csv", # Suicide, 2002-2023
            c("age", "cause", "district", "sex", "year", "n"))   
d_pop <- rd("BEXSTC_20261022-190009.csv",   # Population by district, sex, age
            c("age", "district", "sex", "year", "pop"))           

# District labels (Ivittuut excluded)
districts <- c("Nanortalik", "Qaqortoq", "Narsaq", "Paamiut", "Nuuk",
               "Maniitsoq", "Sisimiut", "Kangaatsiaq", "Aasiaat",
               "Qasigiannguit", "Ilulissat", "Qeqertarsuaq", "Uummannaq",
               "Upernavik", "Qaanaaq", "Tasiilaq", "Ittoqqortoormiit")

# Districts in East Greenland
east <- c("Tasiilaq", "Ittoqqortoormiit")

# Defining the age breaks and labels
age_breaks <- c(-1, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, Inf)
age_labels <- c("0-14", "15-19", "20-24", "25-29", "30-34", "35-39",
                "40-44", "45-49", "50-54", "55-59", "60+")

# Posterior median and 90% credible interval
qs <- function(x) setNames(round(quantile(x, c(0.5, 0.05, 0.95)), 2), 
                           c("median", "lo90", "hi90"))

# Old table uses former names for Qaanaaq and Ittoqqortoormiit
d_old$district <- sub("Avanersuaq", "Qaanaaq", 
                      sub("Illoqqortoormiut", "Ittoqqortoormiit", 
                          d_old$district))
# New table codes sex as Man/Woman; recoded to Men/Women as in the old table
d_new$sex <- ifelse(d_new$sex == "Man", "Men", "Women")

# Stacking both tables with a source label
d <- rbind(cbind(d_old[c("year", "district", "sex", "age", "n")], src = "old"),
           cbind(d_new[c("year", "district", "sex", "age", "n")], src = "new"))
d <- subset(d, district %in% districts)
d$agegrp <- cut(d$age, age_breaks, labels = age_labels)

# Overlap 2002-2010: totals by sex and % of district-years that agree exactly
ov <- subset(d, year %in% 2002:2010)
print(xtabs(n ~ sex + src, ov))
dyx <- xtabs(n ~ year + district + src, ov)
print(round(100 * mean(dyx[, , "new"] == dyx[, , "old"]), 1))

# Old table to 2001, new table from 2002
d <- subset(d, (src == "old" & year <= 2001) | (src == "new" & year >= 2002))

# Mid-year population by district, sex and age group
pop <- subset(d_pop, district %in% districts & sex != "Total")
pop$agegrp <- cut(pop$age, age_breaks, labels = age_labels)
pop <- aggregate(pop ~ year + district + sex + agegrp, pop, sum)
# Matching each cell to the same cell one year later
key <- function(yr) paste(pop$district, pop$sex, pop$agegrp, yr)
pop$py <- (pop$pop + pop$pop[match(key(pop$year + 1), key(pop$year))]) / 2

# Deaths and population by district, sex, age group and year
panel <- merge(subset(pop, year %in% 1983:2023),
               aggregate(n ~ year + district + sex + agegrp, d, sum), 
               all.x = TRUE)
panel$n[is.na(panel$n)] <- 0
panel$east <- panel$district %in% east

# National rates by sex and age, 1983-2023, and expected deaths
std <- aggregate(cbind(n, py) ~ sex + agegrp, panel, sum)
std$rate <- std$n / std$py
print(transform(std, py = round(py), rate = round(1e5 * rate, 1)))
panel$E <- panel$py * std$rate[match(paste(panel$sex, panel$agegrp), 
                                     paste(std$sex, std$agegrp))]
crude <- tapply(std$n, std$sex, sum) / tapply(std$py, std$sex, sum)

# East Greenland ratio under direct and indirect standardization
std$w <- std$py / ave(std$py, std$sex, FUN = sum)
panel$period <- cut(panel$year, c(1982, 1996, 2010, 2023),
                    labels = c("1983-1996", "1997-2010", "2011-2023"))
a <- aggregate(cbind(n, py, E) ~ sex + period + east + agegrp, panel, sum)
a$dsr <- a$n / a$py * std$w[match(paste(a$sex, a$agegrp), 
                                  paste(std$sex, std$agegrp))]
a <- aggregate(cbind(dsr, n, E) ~ sex + period + east, a, sum)
a <- merge(subset(a, east), subset(a, !east), by = c("sex", "period"))
print(data.frame(a[c("sex", "period")], direct = round(a$dsr.x / a$dsr.y, 2),
                 indirect = round((a$n.x / a$E.x) / (a$n.y / a$E.y), 2)))

# Observed and expected suicides by district, sex and year
dy <- aggregate(cbind(O = n, E) ~ district + sex + year, panel, sum)
dy$t <- (dy$year - 2003) / 10
dy$district <- factor(dy$district, levels = districts)
print(tapply(dy$O, dy$sex, sum))

################################################################################
# Model formulas: linear and smooth district trends
################################################################################

# Defining models
f_lin <- bf(O ~ t + (1 + t | district) + offset(log(E)))
f_k5 <- bf(O ~ s(t, k = 5) + s(t, district, bs = "fs", k = 5) + offset(log(E)))
f_k8 <- bf(O ~ s(t, k = 8) + s(t, district, bs = "fs", k = 8) + offset(log(E)))

# Specifying priors
pri_lin <- c(prior(normal(0, 1), class = Intercept),
             prior(normal(0, 1), class = b),
             prior(normal(0, 0.5), class = sd),
             prior(lkj(2), class = cor))
pri_gam <- c(prior(normal(0, 1), class = Intercept),
             prior(normal(0, 1), class = b),
             prior(normal(0, 1), class = sds))
pri_nb <- c(pri_gam, prior(gamma(2, 0.1), class = shape)) 

# Fitting the models
fit <- function(f, s, family, prior) {
  brm(f, data = subset(dy, sex == s), family = family, prior = prior,
      chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 2026, init = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12),
      backend = "cmdstanr", refresh = 0)
}

fits <- list(women_lin = fit(f_lin, "Women", poisson(),     pri_lin),
             men_lin   = fit(f_lin, "Men",   poisson(),     pri_lin),
             women_k5  = fit(f_k5,  "Women", poisson(),     pri_gam),
             men_k5    = fit(f_k5,  "Men",   poisson(),     pri_gam),
             women_k8  = fit(f_k8,  "Women", poisson(),     pri_gam),
             men_k8    = fit(f_k8,  "Men",   poisson(),     pri_gam),
             women_nb  = fit(f_k5,  "Women", negbinomial(), pri_nb),
             men_nb    = fit(f_k5,  "Men",   negbinomial(), pri_nb))

################################################################################
# Results
################################################################################

# Convergence diagnostics
print(t(sapply(fits, function(m) c(
  max_rhat = round(max(rhat(m), na.rm = TRUE), 3),
  min_ess_ratio = round(min(neff_ratio(m), na.rm = TRUE), 2),
  divergent = sum(subset(nuts_params(m), Parameter == "divergent__")$Value)))))

# Fitted rates per 100,000 in 1983, 2003 and 2023
rates <- function(m, s) {
  nd <- expand.grid(district = factor(districts, levels = districts), 
                    year = c(1983, 2003, 2023), E = 1)
  nd$t <- (nd$year - 2003) / 10
  r <- posterior_epred(m, newdata = nd) * 1e5 * crude[[s]]
  colnames(r) <- nd$district
  lapply(split(seq_len(ncol(r)), nd$year), function(j) r[, j])
}

# Stacking weights and stacked posterior
R <- list(); R_models <- list(); I <- list()
for (s in c("Women", "Men")) {
  v <- paste0(tolower(s), c("_lin", "_k5", "_k8"))
  l <- lapply(fits[v], loo)
  w <- loo_model_weights(l, method = "stacking")
  print(loo_compare(l))
  print(round(w, 2))
  print(sapply(l, function(x) max(loo::pareto_k_values(x))))
  R_models[[s]] <- lapply(fits[v], rates, s = s)
  n <- ndraws(fits[[v[1]]])
  # Draws from each model in proportion to its weight, same draws for all years
  I[[s]] <- lapply(round(as.numeric(w) * n), function(nj) sample(n, nj))
  R[[s]] <- sapply(c("1983", "2003", "2023"), function(y)
    do.call(rbind, Map(function(r, i) r[[y]][i, , drop = FALSE], 
                       R_models[[s]], I[[s]])), simplify = FALSE)
}

# Summary measures, giving equal weight to each district
contrasts <- function(R) {
  ch <- log(R[["2023"]] / R[["1983"]])
  ie <- districts %in% east
  x <- list(change_typical = exp(rowMeans(ch)),
            change_typical_outside_east = exp(rowMeans(ch[, !ie])),
            east_level_2003 = exp(rowMeans(log(R[["2003"]][, ie])) -
                                    rowMeans(log(R[["2003"]][, !ie]))),
            east_change_ratio = exp(rowMeans(ch[, ie]) - rowMeans(ch[, !ie])))
  cbind(t(sapply(x, qs)), 
        pr_gt_1 = round(sapply(x, function(v) mean(v > 1)), 2))
}

for (s in c("Women", "Men")) print(contrasts(R[[s]]))

# District estimates (Table 1)
est <- do.call(rbind, lapply(c("Women", "Men"), function(s) {
  i <- dy$sex == s
  obs <- tapply(dy$O[i], dy$district[i], sum)
  ex <- tapply(dy$E[i], dy$district[i], sum)
  data.frame(sex = s, district = districts, obs = c(obs), exp = round(c(ex), 1),
             rate = round(c(1e5 * crude[[s]] * obs / ex)),
             r1983 = t(apply(R[[s]][["1983"]], 2, qs)),
             r2023 = t(apply(R[[s]][["2023"]], 2, qs)),
             change = t(apply(R[[s]][["2023"]] / R[[s]][["1983"]], 2, qs)),
             pr_up = round(colMeans(R[[s]][["2023"]] > R[[s]][["1983"]]), 2), 
             row.names = NULL)
}))

print(est[order(est$sex != "Women", -est$rate), ])

# National change, weighting districts by expected suicides
nat <- function(s, y) {
  i <- dy$sex == s & dy$year == y
  w <- tapply(dy$E[i], dy$district[i], sum)
  drop(R[[s]][[as.character(y)]] %*% w) / sum(w)
}

for (s in c("Women", "Men")) {
  ch <- nat(s, 2023) / nat(s, 1983)
  print(c(qs(ch), pr_up = round(mean(ch > 1), 2)))
}

# Probability that the female-to-male ratio rose
print(round(sapply(districts, function(d)
  mean(R$Women[["2023"]][, d] / R$Men[["2023"]][, d] >
         R$Women[["1983"]][, d] / R$Men[["1983"]][, d])), 2))

# Probability that Tasiilaq has the highest female rate in 2023
print(mean(districts[max.col(R$Women[["2023"]])] == "Tasiilaq"))

# Sensitivity to the choice of model
for (s in c("Women", "Men")) {
  print(lapply(R_models[[s]], contrasts))
  print(contrasts(rates(fits[[paste0(tolower(s), "_nb")]], s)))
}

# Posterior predictive check 
stat <- list(zero_share = function(y) mean(y == 0), max_count = max)
for (s in c("Women", "Men")) {
  v <- paste0(tolower(s), c("_lin", "_k5", "_k8"))
  Y <- do.call(rbind, Map(function(m, i) posterior_predict(m)[i, , drop = FALSE], 
                          fits[v], I[[s]]))
  o <- dy$O[dy$sex == s]
  print(cbind(observed = sapply(stat, function(f) round(f(o), 2)),
              t(sapply(stat, function(f) qs(apply(Y, 1, f))))))
}

################################################################################
# Session info
################################################################################

sessionInfo() 
cmdstanr::cmdstan_version()