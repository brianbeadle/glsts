# Sex-specific suicide trends in Greenland's small districts, 1983–2023

Replication material for the Short Communication "Sex-specific suicide trends in Greenland's small districts, 1983–2023" (*Mental Health & Prevention*). The tables below give the details that the article refers to as "see the online repository".

## Contents

- `analysis.R`: reproduces every number in the article and in the tables below.
- `BEXLD1_20261022-190132.csv`: suicide deaths by year, district, sex, and single year of age, 1983–2010 (StatBank Greenland, archived table BEXLD1, cause "E48 Suicide").
- `SUXLDA1A_20261022-190231.csv`: deaths with underlying cause "Suicide and attempted suicide" by district, sex, and single year of age, 2002–2023 (StatBank Greenland; displayed as table SUELDA1A, downloaded as SUXLDA1A). The category contains deaths only.
- `BEXSTC_20261022-190009.csv`: population on 1 January by district, sex, and single year of age, 1977–2026 (StatBank Greenland, table BEXSTC).

## Running the script

The script needs R with `brms` and the `cmdstanr` backend. Set the working directory in the script to the folder with the three data files. The script prints all results to the console; it saves no files. Fitting the eight models takes a few minutes on four cores.

Results were produced with R 4.5.1, brms 2.23.0, cmdstanr 0.8.0, CmdStan 2.38.0, loo 2.8.0, posterior 1.6.1, and mgcv 1.9-3 on Windows 11. The script ends by printing `sessionInfo()`.

## Data preparation

- BEXLD1 is used for 1983–2001 and SUXLDA1A for 2002–2023 (Table R1).
- Two district names in BEXLD1 were harmonized: Avanersuaq to Qaanaaq and Illoqqortoormiut to Ittoqqortoormiit.
- Missing values ("-") occur only in the count columns and denote zero.
- Ivittuut is excluded: it recorded no suicides, and its population fell from 233 in 1983 to zero by 2020.
- Three suicides without a district are excluded: two recorded outside Greenland (BEXLD1) and one in "other localities" (SUXLDA1A).
- All 17 districts are present in every year of all three tables. No district's population changes by more than 10.4% between consecutive years in 1983–2024, so there is no indication of boundary changes.
- The 17 districts contain at least 98.5% of Greenland's population in every year.
- Mid-year population is the average of the populations on 1 January of the year and the next.
- The analysis includes 1,972 suicides (1,508 men, 464 women) and 697 district-years for each sex.

## Model specification

All models are fitted separately for women and men, with the number of suicides in each district and year as a Poisson outcome and the log of expected suicides as an offset. Time is measured in decades from 2003.

| Model | Formula (brms) | Priors |
| --- | --- | --- |
| Linear | `O ~ t + (1 + t \| district) + offset(log(E))` | Intercept and slope normal(0, 1); district standard deviations half-normal(0, 0.5); correlation LKJ(2) |
| Smooth, K = 5 | `O ~ s(t, k = 5) + s(t, district, bs = "fs", k = 5) + offset(log(E))` | Intercept and linear part of the national smooth normal(0, 1); smoothing standard deviations half-normal(0, 1) |
| Smooth, K = 8 | as above with `k = 8` | as above |
| Negative binomial check | smooth model with K = 5 | as above; shape gamma(2, 0.1) |

Sampling: four chains of 4,000 iterations (2,000 warm-up), giving 8,000 draws; `adapt_delta = 0.99`, `max_treedepth = 12`, initial values zero, seed 2026.

The three Poisson models for each sex are combined by stacking on Pareto-smoothed leave-one-out predictive performance (`loo::loo_model_weights`, method "stacking"). The stacked posterior samples draws from each model in proportion to its weight. The same draws are used for 1983, 2003, and 2023 and for the posterior predictive check.

Following Richardson et al. (2004), a posterior probability of an increase above 0.8 or below 0.2 is treated as evidence of an increase or a decrease. Intervals are equal-tailed 90% credible intervals.

## Additional tables

### Table R1. Suicides in the overlap years 2002–2010, BEXLD1 and SUXLDA1A

| Sex | BEXLD1 | SUXLDA1A |
| --- | ---: | ---: |
| Women | 104 | 108 |
| Men | 310 | 314 |

85.6% of annual district counts (both sexes) agree exactly.

### Table R2. National suicide rates by sex and age, 17 districts, 1983–2023 (per 100,000)

| Age | Women, suicides | Women, rate | Men, suicides | Men, rate |
| --- | ---: | ---: | ---: | ---: |
| 0–14 | 16 | 5.9 | 23 | 8.1 |
| 15–19 | 106 | 127.0 | 261 | 299.4 |
| 20–24 | 87 | 98.4 | 366 | 391.3 |
| 25–29 | 65 | 70.9 | 236 | 234.1 |
| 30–34 | 45 | 50.5 | 146 | 143.0 |
| 35–39 | 30 | 37.4 | 107 | 108.8 |
| 40–44 | 35 | 47.9 | 109 | 116.1 |
| 45–49 | 21 | 31.2 | 82 | 94.1 |
| 50–54 | 21 | 34.1 | 69 | 86.9 |
| 55–59 | 17 | 33.0 | 45 | 67.9 |
| 60+ | 21 | 19.6 | 64 | 55.8 |

All-age rates: women 43.6, men 125.0 per 100,000. These are the rates used to compute expected suicides.

### Table R3. Ratio of East Greenland to other districts under direct and indirect age standardization

| Sex | Period | Direct | Indirect |
| --- | --- | ---: | ---: |
| Women | 1983–1996 | 2.25 | 2.44 |
| Women | 1997–2010 | 3.95 | 3.92 |
| Women | 2011–2023 | 5.13 | 5.02 |
| Men | 1983–1996 | 2.26 | 2.47 |
| Men | 1997–2010 | 3.21 | 3.25 |
| Men | 2011–2023 | 3.35 | 3.20 |

Direct standardization uses the population of the 17 districts over 1983–2023 as the standard. The periods serve only for this check.

### Table R4. MCMC diagnostics

| Model | Max R-hat | Min ESS ratio | Divergent transitions |
| --- | ---: | ---: | ---: |
| Women, linear | 1.002 | 0.24 | 0 |
| Men, linear | 1.003 | 0.18 | 0 |
| Women, smooth (K = 5) | 1.003 | 0.22 | 0 |
| Men, smooth (K = 5) | 1.002 | 0.26 | 0 |
| Women, smooth (K = 8) | 1.003 | 0.17 | 0 |
| Men, smooth (K = 8) | 1.002 | 0.28 | 0 |
| Women, negative binomial | 1.002 | 0.26 | 0 |
| Men, negative binomial | 1.004 | 0.29 | 0 |

### Table R5. Leave-one-out comparison and stacking weights

| Sex | Model | elpd difference | SE | Stacking weight | Max Pareto k |
| --- | --- | ---: | ---: | ---: | ---: |
| Women | Linear | 0.0 | 0.0 | 0.85 | 0.36 |
| Women | Smooth (K = 8) | -2.5 | 2.6 | 0.04 | 0.40 |
| Women | Smooth (K = 5) | -3.9 | 3.3 | 0.11 | 0.42 |
| Men | Smooth (K = 8) | 0.0 | 0.0 | 0.58 | 0.55 |
| Men | Linear | -1.2 | 3.6 | 0.42 | 0.35 |
| Men | Smooth (K = 5) | -2.2 | 1.9 | 0.00 | 0.57 |

The elpd difference is the difference in expected log predictive density from the best model for each sex. All Pareto k values are below 0.7.

### Table R6. Summary measures, stacked models

| Sex | Measure | Median | 90% interval | Pr(> 1) |
| --- | --- | ---: | ---: | ---: |
| Women | Change, typical district | 1.42 | 1.07–1.87 | 0.98 |
| Women | Change, typical district outside East Greenland | 1.33 | 0.98–1.78 | 0.94 |
| Women | East Greenland level, 2003 | 2.53 | 1.80–3.43 | 1.00 |
| Women | East Greenland change ratio | 1.66 | 0.96–3.69 | 0.93 |
| Women | National change | 1.28 | 0.99–1.66 | 0.94 |
| Men | Change, typical district | 0.97 | 0.80–1.20 | 0.40 |
| Men | Change, typical district outside East Greenland | 0.96 | 0.79–1.21 | 0.38 |
| Men | East Greenland level, 2003 | 2.48 | 2.06–3.01 | 1.00 |
| Men | East Greenland change ratio | 1.05 | 0.66–1.67 | 0.57 |
| Men | National change | 0.78 | 0.66–0.96 | 0.02 |

Changes are ratios of the 2023 to the 1983 rate. The typical district is the geometric mean over districts, with each district weighted equally. The East Greenland level is the ratio of geometric mean rates in East Greenland and the other districts in 2003, and the change ratio is the corresponding ratio of geometric mean changes. The national change weights districts by their expected suicides.

### Table R7. Suicides, rates, and change by district

| Sex | District | Suicides | Expected | Observed rate | Fitted 1983 | Fitted 2023 | Change | Pr(rise) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Women | Tasiilaq | 83 | 25.3 | 143 | 82 (46–120) | 195 (134–272) | 2.38 (1.32–4.85) | 0.99 |
| Women | Ittoqqortoormiit | 11 | 4.0 | 118 | 57 (27–103) | 119 (57–250) | 2.02 (0.99–6.00) | 0.95 |
| Women | Paamiut | 26 | 16.0 | 71 | 55 (37–87) | 72 (35–118) | 1.32 (0.53–2.39) | 0.74 |
| Women | Qaanaaq | 10 | 6.7 | 65 | 41 (19–69) | 78 (41–152) | 1.83 (0.96–5.29) | 0.94 |
| Women | Upernavik | 31 | 22.2 | 61 | 45 (24–68) | 67 (34–102) | 1.48 (0.77–2.84) | 0.86 |
| Women | Nanortalik | 20 | 18.3 | 48 | 38 (23–58) | 58 (33–97) | 1.52 (0.78–3.18) | 0.87 |
| Women | Narsaq | 17 | 15.6 | 48 | 37 (22–57) | 56 (30–96) | 1.51 (0.75–3.15) | 0.86 |
| Women | Qaqortoq | 31 | 29.7 | 46 | 36 (22–52) | 55 (35–85) | 1.54 (0.86–3.02) | 0.90 |
| Women | Ilulissat | 39 | 39.9 | 43 | 33 (19–47) | 53 (36–78) | 1.61 (0.96–3.29) | 0.94 |
| Women | Aasiaat | 26 | 30.0 | 38 | 35 (23–57) | 42 (25–67) | 1.20 (0.57–2.13) | 0.69 |
| Women | Uummannaq | 16 | 20.1 | 35 | 31 (18–49) | 41 (22–67) | 1.30 (0.63–2.55) | 0.75 |
| Women | Maniitsoq | 22 | 28.7 | 33 | 30 (18–45) | 42 (24–68) | 1.38 (0.73–2.85) | 0.81 |
| Women | Qasigiannguit | 9 | 12.0 | 33 | 34 (20–60) | 37 (17–69) | 1.11 (0.40–2.20) | 0.59 |
| Women | Sisimiut | 32 | 48.8 | 29 | 26 (15–38) | 35 (22–52) | 1.35 (0.74–2.70) | 0.80 |
| Women | Nuuk | 80 | 127.3 | 27 | 32 (22–49) | 25 (17–36) | 0.80 (0.41–1.42) | 0.27 |
| Women | Qeqertarsuaq | 5 | 8.2 | 26 | 31 (17–56) | 35 (15–69) | 1.17 (0.43–2.43) | 0.63 |
| Women | Kangaatsiaq | 6 | 11.1 | 24 | 28 (15–48) | 33 (16–63) | 1.22 (0.51–2.72) | 0.68 |
| Men | Ittoqqortoormiit | 37 | 12.6 | 367 | 278 (160–438) | 272 (141–468) | 0.98 (0.47–1.99) | 0.48 |
| Men | Tasiilaq | 193 | 74.7 | 323 | 257 (147–359) | 270 (164–377) | 1.06 (0.69–1.67) | 0.58 |
| Men | Qaanaaq | 34 | 21.3 | 200 | 161 (92–252) | 186 (109–299) | 1.15 (0.61–2.27) | 0.64 |
| Men | Qasigiannguit | 48 | 38.4 | 156 | 141 (89–203) | 152 (94–235) | 1.08 (0.60–1.97) | 0.58 |
| Men | Paamiut | 66 | 53.8 | 153 | 147 (101–203) | 150 (97–227) | 1.02 (0.60–1.80) | 0.52 |
| Men | Nanortalik | 73 | 60.0 | 152 | 135 (84–188) | 141 (86–209) | 1.05 (0.61–1.85) | 0.56 |
| Men | Narsaq | 60 | 50.6 | 148 | 148 (94–208) | 125 (77–188) | 0.85 (0.48–1.51) | 0.31 |
| Men | Uummannaq | 80 | 68.4 | 146 | 128 (86–178) | 161 (111–235) | 1.26 (0.76–2.23) | 0.76 |
| Men | Ilulissat | 149 | 128.3 | 145 | 202 (153–263) | 103 (72–155) | 0.51 (0.32–0.85) | 0.02 |
| Men | Upernavik | 86 | 74.4 | 144 | 140 (92–197) | 137 (95–194) | 0.98 (0.59–1.69) | 0.48 |
| Men | Qeqertarsuaq | 28 | 26.3 | 133 | 135 (85–209) | 140 (84–239) | 1.04 (0.56–2.04) | 0.54 |
| Men | Aasiaat | 92 | 94.1 | 122 | 123 (88–177) | 142 (99–221) | 1.16 (0.72–1.95) | 0.69 |
| Men | Kangaatsiaq | 33 | 35.4 | 117 | 98 (55–154) | 131 (81–211) | 1.34 (0.70–2.86) | 0.76 |
| Men | Qaqortoq | 82 | 93.2 | 110 | 123 (87–173) | 104 (72–159) | 0.85 (0.50–1.46) | 0.30 |
| Men | Maniitsoq | 86 | 99.8 | 108 | 113 (80–155) | 114 (78–180) | 1.02 (0.61–1.78) | 0.52 |
| Men | Nuuk | 257 | 408.6 | 79 | 117 (91–147) | 52 (39–67) | 0.44 (0.30–0.66) | <0.01 |
| Men | Sisimiut | 104 | 168.0 | 77 | 72 (48–100) | 87 (63–124) | 1.22 (0.73–2.13) | 0.73 |

Suicides and expected suicides are totals for 1983–2023. The observed rate is the indirectly age-standardized rate per 100,000 for 1983–2023 (observed divided by expected, times the national all-age rate). Fitted rates are posterior medians per 100,000 from the stacked models, with 90% credible intervals. Change is the ratio of the 2023 to the 1983 fitted rate. Pr(rise) is the posterior probability that the rate rose; to three decimals, Sisimiut (women) is 0.803. Districts are ordered by the observed rate within each sex.

### Table R8. Posterior probability that the female-to-male ratio of rates rose, 1983–2023

| District | Pr(rise) | District | Pr(rise) |
| --- | ---: | --- | ---: |
| Ilulissat | 0.99 | Maniitsoq | 0.73 |
| Tasiilaq | 0.96 | Paamiut | 0.67 |
| Qaqortoq | 0.90 | Sisimiut | 0.59 |
| Nuuk | 0.89 | Qeqertarsuaq | 0.56 |
| Ittoqqortoormiit | 0.88 | Aasiaat | 0.52 |
| Narsaq | 0.86 | Uummannaq | 0.52 |
| Upernavik | 0.80 | Qasigiannguit | 0.51 |
| Qaanaaq | 0.80 | Kangaatsiaq | 0.44 |
| Nanortalik | 0.77 | | |

To three decimals, Qaanaaq is 0.804 and Upernavik 0.796. Computed by combining draws from the stacked posteriors for women and men, which share no data or parameters. The posterior probability that Tasiilaq had the highest female rate in 2023 is 0.85.

### Table R9. Summary measures under each candidate model

| Sex | Model | Change, typical | Change, outside East | East level, 2003 | East change ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| Women | Linear | 1.42 (1.08–1.86) [0.98] | 1.34 (1.00–1.78) [0.95] | 2.60 (1.95–3.48) | 1.61 (0.98–3.20) [0.94] |
| Women | Smooth (K = 5) | 1.41 (1.03–1.93) [0.96] | 1.25 (0.90–1.75) [0.87] | 2.04 (1.51–2.73) | 2.67 (1.12–6.34) [0.97] |
| Women | Smooth (K = 8) | 1.42 (1.04–1.90) [0.97] | 1.38 (1.00–1.86) [0.95] | 2.10 (1.53–2.87) | 1.22 (0.69–2.53) [0.69] |
| Women | Negative binomial | 1.40 (1.02–1.93) [0.96] | 1.24 (0.89–1.75) [0.86] | 1.98 (1.49–2.65) | 2.70 (1.10–6.60) [0.96] |
| Men | Linear | 0.93 (0.79–1.09) [0.22] | 0.90 (0.76–1.08) [0.17] | 2.44 (2.07–2.86) | 1.21 (0.80–1.82) [0.78] |
| Men | Smooth (K = 5) | 1.05 (0.87–1.30) [0.67] | 1.04 (0.85–1.28) [0.62] | 2.42 (1.99–2.99) | 1.12 (0.63–2.06) [0.63] |
| Men | Smooth (K = 8) | 1.01 (0.82–1.24) [0.52] | 1.01 (0.83–1.25) [0.53] | 2.51 (2.05–3.08) | 0.95 (0.62–1.47) [0.42] |
| Men | Negative binomial | 1.04 (0.85–1.27) [0.61] | 1.02 (0.83–1.26) [0.56] | 2.32 (1.90–2.88) | 1.13 (0.62–2.10) [0.63] |

Posterior medians with 90% credible intervals and, in brackets, the posterior probability that the ratio exceeds one. The East Greenland level exceeds one with probability 1.00 under every model. The negative binomial model is the smooth model with K = 5.

### Table R10. Posterior predictive check, stacked models

| Sex | Statistic | Observed | Replicated median | 90% interval |
| --- | --- | ---: | ---: | ---: |
| Women | Share of district-years with no suicide | 0.59 | 0.58 | 0.54–0.61 |
| Women | Largest annual district count | 8 | 6 | 5–8 |
| Men | Share of district-years with no suicide | 0.23 | 0.22 | 0.19–0.25 |
| Men | Largest annual district count | 18 | 13 | 11–17 |

The largest male count, 18 suicides in Nuuk in 1984, lies above the replicated interval.

## Reference

Richardson, S., Thomson, A., Best, N., and Elliott, P. (2004). Interpreting posterior relative risk estimates in disease-mapping studies. *Environmental Health Perspectives* 112(9): 1016–1025. doi:10.1289/ehp.6740.
