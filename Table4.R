library(dplyr)
library(ggplot2)
library(latex2exp)
source("QAPpro.R")


# ---- load data ----

load("results/DepressionData/DP.rdata")


# ---- synthetic data sanity check ----

set.seed(2024)
n = 2.5e2
rho = 0.5

R1 = rnorm(n)
R2 = rho * R1 + sqrt(1 - rho^2) * rnorm(n)
S = R1 * rnorm(n)

E  = outer(1:n, 1:n, Vectorize(function(x, y) {(S[x]  + S[y]) /sqrt(2)}))
B1 = outer(1:n, 1:n, Vectorize(function(x, y) {(R1[x] + R1[y])/sqrt(2)}))
B2 = outer(1:n, 1:n, Vectorize(function(x, y) {(R2[x] + R2[y])/sqrt(2)}))

A = B2 + E

diag(A)  = 0
diag(B1) = 0
diag(B2) = 0

data_list = list(A = A, B1 = B1, B2 = B2)
res = dyadicLM("A ~ B1 + B2", data_list)
coefs   = res$coefs
var_mat = res$var_mat

coefs
var_mat

QAPpro(form = "A ~ B1 + B2",
       data_list = data_list,
       mode = "permute-e",
       num_perms = 100,
       var_scheme = "CR0",
       plot.flag = F,
       nameX.target = c("B1")
       )


# ---- variable key ----
# 1  dummy: at least one female
# 2  dummy: both female
# 3  age mean (centered)
# 4  age similarity (random, anonymised)
# 5  one student organisation (random, anonymised)
# 6  same student status
# 7  being friends
# 8  depression mean        -> depression-isolation hypothesis
# 9  depression similarity  -> depression-homophily hypothesis
# 10 depression mean * depression similarity
# 11 depression mean * being friends -> depression-friendship hypothesis


# ---- main analysis: permute-X (network 1) ----

data_list = c(dvs[1], ivs[[1]][2:12])
names(data_list) = c("Y", paste0("X", 1:11))

QAP_res = QAPpro(form = "Y ~ X1 + X2 + X3 + X4 + X7 + X8 + X9 + X10 + X11",
       data_list = data_list,
       mode = "permute-X",
       num_perms = 2000,
       var_scheme = "CR0",
       plot.flag = F,
       nameX.target = c("X8", "X9", "X10", "X11"),
       user.seed = 2024
       )

QAP_res$report

# Save the numerical results.  The first two rows reproduce Table 4; the
# remaining rows are the interaction terms discussed in the case study.
table4_results = cbind(
  term = c(
    "Depression mean",
    "Depression similarity",
    "Depression mean x depression similarity",
    "Depression mean x friendship"
  ),
  QAP_res$report
)
write.csv(table4_results,
          "results/DepressionData/Table4_results.csv",
          row.names = FALSE)


# ---- Depression level randomisation plot ----

MC = 2000
position_x = seq(-4, 4, length.out = MC)
line_normal = data.frame(x = position_x, density = dnorm(position_x))

record_level = data.frame(stat_S = QAP_res$record_perm_S[, 1])

hist_plot1 = record_level %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.3) +
  geom_line(data = line_normal, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}W_{Mean}^\pi$ )")) +
  theme_classic(base_size = 20)
hist_plot1

ggsave("results/DepressionData/DPLevel.png", hist_plot1, device = "png", width = 6, height = 6)


# ---- Depression similarity randomisation plot ----

record_sim = data.frame(stat_S = QAP_res$record_perm_S[, 2])

hist_plot2 = record_sim %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.30) +
  geom_line(data = line_normal, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}W_{Similarity}^\pi$ )")) +
  theme_classic(base_size = 20)
hist_plot2

ggsave("results/DepressionData/DPSimilarity.png", hist_plot2, device = "png", width = 6, height = 6)


# ---- permute-e check (network 2) ----

data_list = c(dvs[2], ivs[[2]][2:12])
names(data_list) = c("Y", paste0("X", 1:11))

QAPpro(form = "Y ~ X1 + X2 + X3 + X4 + X7 + X8 + X9 + X10 + X11",
       data_list = data_list,
       mode = "permute-e",
       num_perms = 1000,
       var_scheme = "CR0",
       plot.flag = F,
       nameX.target = c("X8"),
       user.seed = 2024
       )


# ---- model 2: permute-Y (network 1) ----

data_list = c(dvs[1], ivs[[1]][2:12])
names(data_list) = c("Y", paste0("X", 1:11))

QAP_res = QAPpro(form = "Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9",
       data_list = data_list,
       mode = "permute-Y",
       num_perms = 2000,
       var_scheme = "CR0",
       plot.flag = F,
       nameX.target = paste0("X", 8:9),
       user.seed = 2024
       )

QAP_res$report

record_m2 = data.frame(stat_S = QAP_res$record_perm_S[, 1])

record_m2 %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.25) +
  geom_line(data = line_normal, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\vartheta}$ )")) +
  theme_classic(base_size = 20)
