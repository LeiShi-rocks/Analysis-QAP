library(dplyr)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(mvtnorm)
library(latex2exp)



# ---- Preliminary ----

super.pop = function(A, B) {
  n = nrow(A)
  
  # Analysis based on U-statistics
  abar = 1/(n*(n-1)) * sum(A)
  bbar = 1/(n*(n-1)) * sum(B)
  chat = 1/(n*(n-1)) * sum((A -abar) * (B - bbar))
  vAhat = 1/(n*(n-1)) * sum((A - abar)^2)
  vBhat = 1/(n*(n-1)) * sum((B - bbar)^2)
  rhohat = chat/(sqrt(vAhat * vBhat))
  
  # Variance estimation
  Acenter = A  - abar
  Bcenter = B  - bbar
  AB = Acenter * Bcenter
  diag(AB) = 0
  vPhat = 4 * (n-1)/((n-2)*(n-4)) * sum((rowSums(AB, dim = 1L) / (n-1))^2) / (vAhat * vBhat)
  
  # report results
  list(
    abar = abar,
    bbar = bbar,
    chat = chat,
    vAhat = vAhat,
    vBhat = vBhat,
    rhohat = rhohat,
    vPhat = vPhat
  )
}




# ---- Walsh averages kernel ----

# ---- Conservative permutation test ----

set.seed(2024)
MC = 2e3
n = 5e2

record_super = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

message("Conservative sampling distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){
  
  U = runif(n, min = 0, max = 2*pi)
  R = sqrt(2) * sin(U)
  S = sqrt(2) * cos(U)
  
  A = outer(1:n, 1:n, Vectorize(function(x,y) {(R[x] + R[y])/sqrt(2)}))
  B = outer(1:n, 1:n, Vectorize(function(x,y) {(S[x] + S[y])/sqrt(2)}))
  diag(A) = 0
  diag(B) = 0
  
  res = super.pop(A, B)
  
  record_super$stat_NS[iter] = sqrt(n) * res$rhohat
  record_super$stat_S[iter] = sqrt(n) * res$rhohat/sqrt(res$vPhat)
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_super, file = "results/Fig1/record_super_conservative.rds")

# visualization
record_super = record_super[,1:2]
position_x = rep(seq(-4, 4, length.out = MC), 6)
normal_density = dnorm(position_x)
record_super = cbind(record_super, position_x, normal_density)

super_non_student = record_super %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(aes(x = position_x, y = normal_density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}$ )")) + 
  theme_classic(base_size = 20)

super_student = record_super %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(aes(x = position_x, y = normal_density), lty = 2) + 
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}/\hat{v}_P$ )")) + 
  theme_classic(base_size = 20)

super_non_student
super_student

ggsave("results/Fig1/super_non_student.png", super_non_student, device = "png")
ggsave("results/Fig1/super_student.png", super_student, device = "png")


set.seed(2024)

U = runif(n, min = 0, max = 2*pi)
R = sqrt(2) * sin(U)
S = sqrt(2) * cos(U)

A = outer(1:n, 1:n, Vectorize(function(x,y) {(R[x] + R[y])/sqrt(2)}))
B = outer(1:n, 1:n, Vectorize(function(x,y) {(S[x] + S[y])/sqrt(2)}))
diag(A) = 0
diag(B) = 0

MC = 2e3

record_perm = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

# DIPS permutation
message("Conservative permutation distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){
  permInd = sample(1:n)
  res = super.pop(A[permInd, permInd], B)
  record_perm$stat_NS[iter] = sqrt(n) * res$rhohat
  record_perm$stat_S[iter] = sqrt(n) * res$rhohat/sqrt(res$vPhat)
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_perm, file = "results/Fig1/record_perm_conservative.rds")

# visualization
record_perm = record_perm[,1:2]
sd_super_NS = sqrt(var(record_super$stat_NS))
position_x = seq(-4, 4, length.out = MC)
line_NS = data.frame(x = position_x, density = dnorm(position_x, sd = sd_super_NS))
line_S  = data.frame(x = position_x, density = dnorm(position_x))

perm_non_student = record_perm %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(data = line_NS, aes(x = x, y = density), lty = 2) + 
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}^\pi$ )")) + 
  theme_classic(base_size = 20)

perm_student = record_perm %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(data = line_S, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}^\pi/\hat{v}_P^\pi$ )")) + 
  theme_classic(base_size = 20)

perm_non_student
perm_student

ggsave("results/Fig1/perm_non_student.png", perm_non_student, device = "png")
ggsave("results/Fig1/perm_student.png", perm_student, device = "png")


# ---- Anti-conservative permutation test ----

set.seed(2024)
k = 19
# integrand <- function(x) {x^(2*k)}
integrand <- function(x) {sinh(3*x)^2}
res_sq <- integrate(integrand, lower = -2*pi, upper = 2*pi)

MC = 2e3
n = 5e4

U = runif(n, min = -2*pi, max = 2*pi)
# R = U^k / sqrt(res$value / (4*pi))
R = sinh(3*U) / sqrt(res_sq$value / (4*pi))
S = sqrt(2) * cos(U)

mean(R^2*S^2)

MC = 2e3
n = 5e2

integrand <- function(x) {sinh(3*x)^2}
res_sq <- integrate(integrand, lower = -2*pi, upper = 2*pi)

record_super = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

message("Anti-conservative sampling distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){
  
  U = runif(n, min = -2*pi, max = 2*pi)
  R = sinh(3*U) / sqrt(res_sq$value / (4*pi))
  S = sqrt(2) * cos(U)
  
  A = outer(1:n, 1:n, Vectorize(function(x,y) {(R[x] + R[y])/sqrt(2)}))
  B = outer(1:n, 1:n, Vectorize(function(x,y) {(S[x] + S[y])/sqrt(2)}))
  diag(A) = 0
  diag(B) = 0
  
  res = super.pop(A, B)
  
  record_super$stat_NS[iter] = sqrt(n) * res$rhohat
  record_super$stat_S[iter] = sqrt(n) * res$rhohat/sqrt(res$vPhat)
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_super, file = "results/Fig1/record_super_anti.rds")

# visualization
record_super = record_super[,1:2]
position_x = rep(seq(-4, 4, length.out = MC), 6)
normal_density = dnorm(position_x)
record_super = cbind(record_super, position_x, normal_density)

super_non_student = record_super %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(aes(x = position_x, y = normal_density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}$ )")) + 
  theme_classic(base_size = 20)

super_student = record_super %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(aes(x = position_x, y = normal_density), lty = 2) + 
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}/\hat{v}_P$ )")) + 
  theme_classic(base_size = 20)

super_non_student
super_student

ggsave("results/Fig1/super_non_student_anti.png", super_non_student, device = "png")
ggsave("results/Fig1/super_student_anti.png", super_student, device = "png")


set.seed(2024)

integrand <- function(x) {sinh(3*x)^2}
res_sq <- integrate(integrand, lower = -2*pi, upper = 2*pi)

U = runif(n, min = -2*pi, max = 2*pi)
R = sinh(3*U) / sqrt(res_sq$value / (4*pi))
S = sqrt(2) * cos(U)

A = outer(1:n, 1:n, Vectorize(function(x,y) {(R[x] + R[y])/sqrt(2)}))
B = outer(1:n, 1:n, Vectorize(function(x,y) {(S[x] + S[y])/sqrt(2)}))
diag(A) = 0
diag(B) = 0

MC = 2e3

record_perm = data.frame(
  stat_NS = rep(0, MC),
  stat_S = rep(0, MC)
)

# DIPS permutation
message("Anti-conservative permutation distribution")
pb = txtProgressBar(min = 0, max = MC, style = 3)
for (iter in 1:MC){
  permInd = sample(1:n)
  res = super.pop(A[permInd, permInd], B)
  record_perm$stat_NS[iter] = sqrt(n) * res$rhohat
  record_perm$stat_S[iter] = sqrt(n) * res$rhohat/sqrt(res$vPhat)
  setTxtProgressBar(pb, iter)
}
close(pb)

saveRDS(record_perm, file = "results/Fig1/record_perm_anti.rds")

# visualization
record_perm = record_perm[,1:2]
sd_super_NS = sqrt(var(record_super$stat_NS))
position_x = seq(-4, 4, length.out = MC)
line_NS = data.frame(x = position_x, density = dnorm(position_x, sd = sd_super_NS))
line_S  = data.frame(x = position_x, density = dnorm(position_x))

perm_non_student = record_perm %>% ggplot(aes(x = stat_NS)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(data = line_NS, aes(x = x, y = density), lty = 2) + 
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}^\pi$ )")) + 
  theme_classic(base_size = 20)

perm_student = record_perm %>% ggplot(aes(x = stat_S)) +
  geom_histogram(aes(y = ..density..), col = "white", fill = "grey60", binwidth = 0.20) +
  geom_line(data = line_S, aes(x = x, y = density), lty = 2) +
  xlab(TeX(r"( $\sqrt{n}\hat{\rho}^\pi/\hat{v}_P^\pi$ )")) + 
  theme_classic(base_size = 20)

perm_non_student
perm_student

ggsave("results/Fig1/perm_non_student_anti.png", perm_non_student, device = "png")
ggsave("results/Fig1/perm_student_anti.png", perm_student, device = "png")
