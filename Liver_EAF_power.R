# --------------------------------------------------------------------#

# Title: Early Allograft Failure in Liver Transplants
# Developer: Fernando Palluzzi, PhD
# PI: AW Avolio, MD
# Centre: Agostino Gemelli
# Since: 04/08/2021

# --------------------------------------------------------------------#




## Priors

# Incidence of EAF at 90 days: 11.1% (derivation set)         [1]

# Incidence of EAF stratified by donor:
# - Liver transplantation from standard deceased donors:       7%
# - Liver transplantation from high-risk deceased donors:  10-12%
# - Liver transplantation from living donors:                3-4%

# Direct comparison EASE vs. L-GrAFT10:
# -      AUC(EASE) = 0.87 (0.83, 0.91)                        [2,3]
# - AUC(L-GrAFT10) = 0.72 (0.65, 0.78)                        [2,3]

# L-GrAFT7 C stat. = 0.78 (0.75, 0.82)                          [4]

# [1] Agopian VG et al. 2017; JAMA Surgery.
# [2] Avolio AW et al. 2021; Journal of Hepatology (Letter to the Editor).
# [3] Avolio AW et al. 2020; JAMA Surgery.
# [4] Agopian VG et al. 2021; Journal of Hepatology.

# R

library(pwr)
library(pROC)

randSwitch <- function(x, p) {
	y <- vector()
	for (j in 1:length(x)) {
		if (runif(1) >= p) {
			y <- c(y, ifelse(x[j] == 1, 0, 1))
		} else {
			y <- c(y, x[j])
		}
	}
	return(y)
}



#-----------------------------------------------------#
#  Power analysis A: EASE score algorithm validation  #
#-----------------------------------------------------#


# Current target n: 4000-5000 subjects (overall)
#                   1600-2000 subjects (prospective)
#                   2400-3000 subjects (retrospective)


### Simulation A-1. EASE vs. L-GrAFT10.

# - Incidence of EAF (p0) = 11.1%   [1]
# - L-GrAFT10 score AUC = 0.72      [2]
# - EASE score AUC = 0.87           [2]

# Goal.
# To maintain (or improve) EASE score performance (AUC) at day 7 
# of follow-up, rather than day 10, considering L-GrAFT10 as a baseline.

set.seed(123)

n0 <- 5000                         # Overall sample size (= highest target sample size)
p0 <- 0.111                        # EAF incidence (= 11.1%, from [1])

y <- rbinom(n0, 1, p0)             # EAF occurrence is sampled from a Binomial 
                                   # distribution with p0 = 0.111 and n0 = 5000.
                                   
# Score performance is simulated using p = P(x = k | y = k), for a given AUC.
# In other words, 1 - p is the prediction error (= (FP + FN)/n), such that 
# the estimated AUC for L-GrAFT10 and EASE scores are equal to the observed 
# AUC value.

# Generating ROC objects for L-GrAFT10 and EASE scores

x0 <- randSwitch(y, p = 0.712)     # AUC.L-GrAFT10 = 0.72 (baseline) [2]
roc0 <- roc(y, x0)

x1 <- randSwitch(y, p = 0.872)     # AUC.EASE = 0.87 [2]
roc1 <- roc(y, x1)

# Sample size calculation (based on deLong test)

power.roc.test(roc0, roc1, sig.level = 0.05, power = 0.8)
#         ncases = 36.33309
#      ncontrols = 310.3567
#           auc1 = 0.7232552
#           auc2 = 0.8713717
#      sig.level = 0.05
#          power = 0.8
#    alternative = two.sided
#              n = 36.33309 + 310.3567 = 346.6898   (~ 350 subjects)


### Simulation A-2. EASE vs. L-GrAFT7.

# - Incidence of EAF (p0) = 11.1%   [1]
# - L-GrAFT7 score AUC = 0.78       [1,3]
# - EASE score AUC = 0.87           [2]

# Goal.
# To maintain (or improve) EASE score performance (AUC) at day 7 
# of follow-up, rather than day 10, considering L-GrAFT7 as a baseline.
# This study takes into consideration the performances of the baseline 
# score to have a reference AUC at day 7.

set.seed(123)

n0 <- 5000                         # Overall target sample size (n0)
p0 <- 0.111                        # EAF incidence (= 11.1%, from [1])

y <- rbinom(n0, 1, p0)             # EAF occurrence is sampled from a Binomial 
                                   # distribution with p0 = 0.111 and n0 = 5000.

# Generating ROC objects for L-GrAFT7 and EASE scores

x0 <- randSwitch(y, p = 0.765)     # AUC0.L-GrAFT7 = 0.78 (baseline) [2]
roc0 <- roc(y, x0)

x1 <- randSwitch(y, p = 0.872)     # AUC.EASE = 0.87 [2]
roc1 <- roc(y, x1)

# Sample size calculation (based on deLong test)

power.roc.test(roc0, roc1, sig.level = 0.05, power = 0.8)
#         ncases = 82.52132
#      ncontrols = 704.8958
#           auc1 = 0.7817162
#           auc2 = 0.8713717
#      sig.level = 0.05
#          power = 0.8
#    alternative = two.sided
#              n = 82.52132 + 704.8958 = 787.4171   (~ 790 subjects)



#---------------------------------------------------------#
#  Power analysis B: Explore different incidences of EAF  #
#---------------------------------------------------------#

# Current target sizes (n):
#    800 subjects (living donors; 3% EAF incidence)
#   1000 subjects (standard deceased donors; 7% EAF incidence)
#    200 subjects (high-risk deceased donors; 10% EAF incidence)

# Goal
# To provide an estimation of the minimum sample size to achieve baseline 
# AUC, considering a different EAF incidence for each stratum.

# Stratum 1 (living donors), 3% EAF
power.roc.test(auc = 0.72, sig.level = 0.05, power = 0.8, kappa = (1 - 0.03)/0.03)
#         ncases = 13.26893
#      ncontrols = 429.0286
#            auc = 0.72
#      sig.level = 0.05
#          power = 0.8
#             n1 = 13.26893 + 429.0286 =  442.2975    (~ 442 subjects)

# Stratum 2 (standard deceased donors), 7% EAF
power.roc.test(auc = 0.72, sig.level = 0.05, power = 0.8, kappa = (1 - 0.07)/0.07)
#         ncases = 13.79953
#      ncontrols = 183.3367
#            auc = 0.72
#      sig.level = 0.05
#          power = 0.8
#             n2 = 13.79953 + 183.3367 = 197.1362    (~ 200 subjects)

# Stratum 3 (high-risk deceased donors), 10% EAF
power.roc.test(auc = 0.72, sig.level = 0.05, power = 0.8, kappa = (1 - 0.10)/0.10)
#         ncases = 14.22835
#      ncontrols = 128.0552
#            auc = 0.72
#      sig.level = 0.05
#          power = 0.8
#             n3 = 14.22835 + 128.0552 = 142.2836    (~ 142 subjects)



#----------------------------------------------------------#
#  Power analysis C: Class 5 EASE score cutoff validation  #
#----------------------------------------------------------#

# 3.36% of cases belong to class 5
p5 <- 0.0336

# Considering a minimum required AUC of 0.70 to be achieved for every class
# (including class 5):

power.roc.test(auc = 0.70, sig.level = 0.05, power = 0.8, kappa = (1 - p5)/p5)
#         ncases = 16.14696
#      ncontrols = 464.4174
#            auc = 0.70
#      sig.level = 0.05
#          power = 0.8
#              n = 16.14696 + 464.4174 = 480.5644   (> 480 subjects)


