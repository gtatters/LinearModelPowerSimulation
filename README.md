# Sampling Variability and Statistical Power

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

Every time you run an experiment you get a slightly different answer, even if
the true effect stays the same. This app makes that variability visible by
running 100 simulated experiments at once and displaying how the estimated
effect size and p-values vary across them.

You can switch between a **categorical predictor** (two-group comparison,
                                                    like control vs. treatment) and a **continuous predictor** (simple linear
                                                                                                                regression). Adjust the true difference, the standard deviation, and the
sample size, then press **Run Simulation** to see:
  
  - The population distributions and one random sample
- The distribution of effect-size estimates across 100 trials
- The distribution of p-values across 100 trials, with the chosen α level marked

## How to use

1. Choose a predictor type (categorical or continuous).
2. Set the true difference, standard deviation, and sample size.
3. Choose a significance level (α).
4. Press **Run Simulation** — try the same settings multiple times to see
how much the results fluctuate.
5. Increase sample size and observe the p-value distribution shifting toward zero.

## Learning goals

- Understand that a single experiment is just one draw from a distribution
of possible results
- See how statistical power depends on sample size, effect size, and variability
- Understand why small studies are unreliable even when the effect is real

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).