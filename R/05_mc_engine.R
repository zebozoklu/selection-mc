# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario table, number of reps R, estimator functions
# Outputs:
#   - raw draws (estimates per rep) OR
#   - aggregated summaries per scenario/estimator
# Sanity checks: reproducible with seed; progress logging; saves intermediate results
