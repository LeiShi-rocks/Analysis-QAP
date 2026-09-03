library(dplyr)
library(ggplot2)
library(latex2exp)
source("QAPpro.R")


# ---- population setup ----

set.seed(2024)
n = 2.5e2

R = rnorm(n)
S = R * rnorm(n)

A = outer(1:n, 1:n, Vectorize(function(x, y) {(R[x] + R[y])/sqrt(2)}))
B = outer(1:n, 1:n, Vectorize(function(x, y) {(S[x] + S[y])/sqrt(2)}))

diag(A) = 0
diag(B) = 0


fit.lm = lm(c(B) ~ c(A))
cluster0 = factor(rep(1:n, each = n))
res = vcovCR(
  fit.lm,
  cluster = cluster0,
  type = "CR0"
)


# ---- Experiments for sampling distribution ----

set.seed(2024)
MC = 2e3
n = 2.5e2
rho = 0.5
record_stat = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

message("Sampling distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){

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

  record_stat$stat_NS[iter] = sqrt(n) * coefs[2]
  record_stat$stat_S[iter]  = coefs[2]/sqrt(var_mat[2,2])
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_stat, file = "results/Fig2/record_stat_MRQAP.rds")

# ---- sampling distribution plots ----

position_x    = seq(-4, 4, length.out = MC)
normal_density = dnorm(position_x)
line_normal    = data.frame(x = position_x, density = normal_density)

sampling_NS = record_stat %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.25) +
  geom_line(data = line_normal, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\vartheta}$ )")) +
  theme_classic(base_size = 20)

sampling_S = record_stat %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.25) +
  geom_line(data = line_normal, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $W$ )")) +
  theme_classic(base_size = 20)

ggsave("results/Fig2/sampling_NS_MAQAP.png", sampling_NS, device = "png", width = 6, height = 6)
ggsave("results/Fig2/sampling_S_MAQAP.png",  sampling_S,  device = "png", width = 6, height = 6)


# ---- Experiments for permutation distribution ----

set.seed(2024)
n = 2.5e2

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

MC = 2e3
record_perm = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

message("Permutation distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){
  permInd = sample(1:n)
  data_list = list(A = A, B1 = B1[permInd, permInd], B2 = B2)
  res = dyadicLM("A ~ B1 + B2", data_list)
  coefs   = res$coefs
  var_mat = res$var_mat

  record_perm$stat_NS[iter] = sqrt(n) * coefs[2]
  record_perm$stat_S[iter]  = coefs[2]/sqrt(var_mat[2,2])
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_perm, file = "results/Fig2/record_perm_MRQAP.rds")

# ---- permutation distribution plots ----

sd_super_NS = sqrt(var(record_stat$stat_NS))
position_x  = seq(-4, 4, length.out = MC)
line_NS     = data.frame(x = position_x, density = dnorm(position_x, sd = sd_super_NS))

position_w  = seq(0, 15, length.out = MC)
line_W      = data.frame(x = position_w, density = dchisq(position_w, df = 1))

perm_non_student = record_perm %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.25) +
  geom_line(data = line_NS, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\vartheta}^\pi$ )")) +
  theme_classic(base_size = 20)

perm_student = record_perm %>% ggplot(aes(x = stat_S^2)) +
  geom_histogram(aes(y = after_stat(density)), col = "white", fill = "grey60", binwidth = 0.15) +
  geom_line(data = line_W, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $W^\pi$ )")) +
  scale_x_continuous(limits = c(0, 7.5), breaks = seq(0, 7.5, by = 2.5)) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25)) +
  theme_classic(base_size = 20)

perm_non_student
perm_student

ggsave("results/Fig2/perm_NS_MAQAP.png", perm_non_student, device = "png", width = 6, height = 6)
ggsave("results/Fig2/perm_S_MAQAP.png",  perm_student,     device = "png", width = 6, height = 6)
