library(rstan)
library(tidyverse)
library(driver)
library(grid)
library(tidybayes)
library(latex2exp)

setwd("~/Research/zero_types/results/2018-04-11_simulation/")

set.seed(4)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Datasets ----------------------------------------------------------------

dat <- list()

# Universal priors
lambda_prior <- c(5, 3)
theta_prior <- c(.5, .5)
gamma_prior <- c(1, 1)

# Type I Zeroes
lambda1 <- .5
type1 <- c(0, 1, 0, 2, 0, 0)
dat$type1 <- within(list(), {
  y <- type1
  z <- rep(1, length(y))
  x <- rep(1, length(y))
  N <- length(y)
  N_person <- 1
  N_batch <- 1
  
  lambda_prior <- lambda_prior
  theta_prior <- theta_prior
  gamma_prior <- gamma_prior
})

# Type IIa Zeroes
lambdas2a <- c(1.4, 0.6, 3.2)
batch <- rep(1:3, each=5)
type2a <- c(0,1,2,3,1,0,1,0,2,0,8,1,2,0,5)
dat$type2a <- within(list(), {
  y <- type2a
  x <- batch
  z <- rep(1, length(y))
  N <- length(y)
  N_person <- length(unique(z))
  N_batch <- length(unique(x))
  lambda_prior <- lambda_prior
  theta_prior <- theta_prior
  gamma_prior <- gamma_prior
})

# Type IIb Zeroes
lambda2b <- 1
type2b <- c(2, 0, 0, 0, 0, 1, 0, 3, 2, 3, 1, 1, 1, 2, 0, 2, 2, 1, 2, 1)
# Sim from bernouli with p=0.3
quench <- c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0)
type2b[quench==1] <- 0
dat$type2b <- within(list(), {
  y <- type2b
  x <- rep(1, length(y))
  z <- rep(1, length(y))
  N <- length(y)
  N_person <- 1
  N_batch <- 1
  
  lambda_prior <- lambda_prior
  theta_prior <- theta_prior
  gamma_prior <- gamma_prior
})

# # Type IIb Zeroes (Batch Simulation)
lambda2bBATCH <- 1
batch <- rep(1:3, each=5)
type2bBATCH <- c(0,1,2,1,1,0,0,0,0,0,0,1,3,1,0)
type2bBATCH[batch==2] <- 0
dat$type2bBATCH <- within(list(), {
  y <- type2bBATCH
  x <- batch
  z <- rep(1, length(y))
  N <- length(y)
  N_person <- 1
  N_batch <- length(unique(x))
  
  lambda_prior <- lambda_prior
  theta_prior <- theta_prior
  gamma_prior <- gamma_prior
})

# Type III Zeroes
lambdas3 <- c(1.4, 0, 3.2)
person <- rep(1:3, each=5)
type3 <- c(0,1,2,3,1,0,0,0,0,0,8,1,2,0,5)
dat$type3 <- within(list(), {
  y <- type3
  x <- person
  z <- person
  N <- length(y)
  N_person <- length(unique(z))
  N_batch <- N_person
  lambda_prior <- lambda_prior
  theta_prior <- theta_prior
  gamma_prior <- gamma_prior
})


# Models ------------------------------------------------------------------

fit.t1m1 <- stan("m1.stan", data=dat$type1)
fit.t2bm1 <- stan("m1.stan", data=dat$type2b)
fit.t2bBATCHm1 <- stan("m1.stan", data=dat$type2bBATCH)
fit.t3m1 <- stan("m1.stan", data=dat$type3)

fit.t1m2b <- stan("m2b.stan", data=dat$type1)
fit.t2bm2b <- stan("m2b.stan", data=dat$type2b)
fit.t2bBATCHm2b <- stan("m2b.stan", data=dat$type2bBATCH)
fit.t3m2b <- stan("m2b.stan", data=dat$type3)

fit.t1m3 <- stan("m3.stan", data=dat$type1, control=list(adapt_delta=0.9))
fit.t2bm3 <- stan("m3.stan", data=dat$type2b, control=list(adapt_delta=0.9))
fit.t2bBATCHm3 <- stan("m3.stan", data=dat$type2bBATCH, control=list(adapt_delta=0.9))
fit.t3m3 <- stan("m3.stan", data=dat$type3, control=list(adapt_delta=0.9))

## Special Models ##
# Type 0 
fit.t1m0a <- stan("m0.stan", data=c(dat$type1, k=.05))
fit.t1m0b <- stan("m0.stan", data=c(dat$type1, k=.5))
fit.t1m0c <- stan("m0.stan", data=c(dat$type1, k=1))
fit.t3m0 <- stan("m0.stan", data=c(dat$type3, k=.5))

# Type IIa
fit.t2am1 <- stan("m1.stan", data=dat$type2a)
fit.t2am2a <- stan("m2a.stan", data=dat$type2a)
#fit.t2am2b <- stan("m2b_for_t2a.stan", data=dat$type2a) # here theta_1 = 0 to mimic m2a (not used in manuscript)
fit.t2am2b <- stan("m2b.stan", data=dat$type2a)
fit.t2am3 <- stan("m3.stan", data=dat$type2a)

# Other application of Type IIa (for Type 2b)
fit.t2bm2a <- stan("m2a.stan", data=dat$type2b)
fit.t2bBATCHm2a <- stan("m2a.stan", data=dat$type2bBATCH)



fits <- list()
fits$t1 <- list(`Model 0\n(k=0.05)` = fit.t1m0a, 
                `Model 0\n(k=0.5)` = fit.t1m0b,
                `Model 0\n(k=1.0)` = fit.t1m0c, 
                `Model I` = fit.t1m1, 
                `Model IIb` = fit.t1m2b, 
                `Model III` = fit.t1m3)
fits$t2a <- list(`Model I` = fit.t2am1, 
                 `Model IIa` = fit.t2am2a,
                 `Model IIb` = fit.t2am2b, 
                 `Model III` = fit.t2am3)
fits$t2b <- list(`Model I` = fit.t2bm1,
                 `Model IIa` = fit.t2bm2a,
                 `Model IIb` = fit.t2bm2b, 
                 `Model III` = fit.t2bm3)
fits$t2bBATCH <- list(`Model I` = fit.t2bBATCHm1,
                      `Model IIa` = fit.t2bBATCHm2a,
                      `Model IIb` = fit.t2bBATCHm2b, 
                      `Model III` = fit.t2bBATCHm3)
fits$t3 <- list(`Model 0\n(k=1.0)` = fit.t3m0,
                `Model I` = fit.t3m1, 
                `Model IIb` = fit.t3m2b, 
                `Model III` = fit.t3m3)

# Plots for Type I --------------------------------------------------------

llambda <- lambda1
pt1 <- fits$t1 %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  #map(rstan::extract, pars=c("lambda", "lambda0")) %>% 
  map(data.frame) %>% 
  map(rownames_to_column, "iter") %>% 
  map(gather, param, val, -iter) %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>%
  mutate(model = factor(model, levels=c("Model III", "Model IIb", "Model I", 
                                        "Model 0\n(k=0.05)", "Model 0\n(k=0.5)", 
                                        "Model 0\n(k=1.0)"))) %>% 
  select(-param) %>% 
  ggplot(aes(x=val, y=model)) +
  geom_halfeyeh(.prob=c(.95, .50)) + 
  geom_segment(x = llambda, xend = llambda, y = .7, yend=1.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 1.7, yend=2.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 2.7, yend=3.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 3.7, yend=4.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 4.7, yend=5.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 5.7, yend=6.3, color="darkred") +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12))+
  coord_cartesian(xlim=c(0,3))
ggsave("type1.pdf",plot=pt1, height=2.8, width=4, units="in")



# Plots for Type IIa ------------------------------------------------------

llambda <- lambdas2a[1]
pt2a <- fits$t2a %>% 
  #map(rstan::extract, pars=c("lambda", "lambda0")) %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  map(data.frame) %>% 
  map(rownames_to_column, "iter") %>% 
  map(gather, param, val, -iter) %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>% 
  select(-param) %>% 
  mutate(model = factor(model, levels=c("Model III", "Model IIb",  "Model IIa", "Model I"))) %>% 
  ggplot(aes(x=val, y=model)) +
  geom_halfeyeh(.prob=c(.95, .50)) + 
  geom_segment(x = llambda, xend = llambda, y = .7, yend=1.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 1.7, yend=2.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 2.7, yend=3.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 3.7, yend=4.3, color="darkred") +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12)) +
  coord_cartesian(xlim=c(0,4))
ggsave("type2a.pdf", plot=pt2a, height=1.8, width=4, units="in")

# Plots for Type IIb -------------------------------------------------------

llambda <- lambda2b
pt2b <- fits$t2b %>% 
  #map(rstan::extract, pars=c("lambda", "lambda0")) %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  map(data.frame) %>% 
  map(rownames_to_column, "iter") %>% 
  map(gather, param, val, -iter) %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>% 
  select(-param) %>% 
  mutate(model = factor(model, levels=c("Model III", "Model IIb", "Model IIa", "Model I"))) %>% 
  ggplot(aes(x=val, y=model)) +
  geom_halfeyeh(.prob=c(.95, .50)) + 
  geom_segment(x = llambda, xend = llambda, y = .7, yend=1.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 1.7, yend=2.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 2.7, yend=3.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 3.7, yend=4.3, color="darkred") +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12))
ggsave("type2b.pdf", plot=pt2b, height=1.8, width=4, units="in")


# Plots for Type IIbBATCH -------------------------------------------------------

llambda <- lambda2bBATCH
pt2bBATCH <- fits$t2bBATCH %>% 
  #map(rstan::extract, pars=c("lambda", "lambda0")) %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  map(data.frame) %>% 
  map(rownames_to_column, "iter") %>% 
  map(gather, param, val, -iter) %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>% 
  select(-param) %>% 
  mutate(model = factor(model, levels=c("Model III", "Model IIb", "Model IIa", "Model I"))) %>% 
  ggplot(aes(x=val, y=model)) +
  geom_halfeyeh(.prob=c(.95, .50)) + 
  geom_segment(x = llambda, xend = llambda, y = .7, yend=1.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 1.7, yend=2.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 2.7, yend=3.3, color="darkred") +
  geom_segment(x = llambda, xend = llambda, y = 3.7, yend=4.3, color="darkred") +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12))
ggsave("type2bBATCH.pdf", plot=pt2bBATCH, height=2.5, width=4, units="in")


# Plots for Type III ------------------------------------------------------

# Annotations
llambdas <- data.frame(lambda =lambdas3, person=factor(1:length(lambdas3))) %>% 
  bind_rows(., ., ., .id="model") %>% 
  mutate(model = as.double(model)) %>% 
  mutate(person = factor(paste("Person",person))) %>% 
  rename(val = lambda)
llambdas.log <- data.frame(lambda =lambdas3, person=factor(1:length(lambdas3))) %>% 
  bind_rows(., ., ., .id="model") %>% 
  mutate(model = as.double(model)) %>% 
  filter(person!=2) %>% 
  mutate(person = factor(paste("Person",person))) %>% 
  rename(val = lambda)
arrows <- data.frame(val = 0, person="Person 2", model=c(1,2,3,4))

# Main Plot (Not Log Scale)
pt3 <- fits$t3 %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  modify_depth(2, gather_array,val, iter, person) %>% 
  map(bind_rows, .id="param") %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>% 
  select(-param) %>% 
  mutate(model = factor(model, levels=c("Model III", "Model IIb", "Model I", "Model 0\n(k=1.0)"))) %>% 
  mutate(person = factor(paste("Person",person))) %>% 
  ggplot(aes(x=val, y = model)) +
  geom_halfeyeh() +
  facet_grid(person~.) + 
  geom_segment(data=llambdas, aes(xend=val,y=model-.3, yend=model+.3), 
               color="darkred")+
  coord_cartesian(xlim=c(0, 7)) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        strip.text.y = element_text(size=12, angle = 0))
ggsave("type3.pdf", plot=pt3, height=5, width=5, units="in")


# Focused Log Scale Plot
pt3log <- fits$t3 %>% 
  map(rstan::extract, pars=c("lambda")) %>% 
  modify_depth(2, gather_array,val, iter, person) %>% 
  map(bind_rows, .id="param") %>% 
  bind_rows(.id="model") %>% 
  filter(param != "lambda0") %>% 
  select(-param) %>% 
  mutate(model = factor(model, levels=c("Model III", "Model IIb", "Model I", "Model 0\n(k=1.0)"))) %>% 
  filter(person==2) %>% 
  mutate(person = factor(paste("Person",person))) %>% 
  ggplot(aes(x=val, y = model)) +
  geom_halfeyeh() +
  facet_grid(person~.) + 
  #geom_segment(data=llambdas.log, aes(xend=val,y=model-.3, yend=model+.3), 
  #             color="darkred")+
  geom_segment(data=arrows, aes(x=1e-6, xend=1e-7, y=model, yend=model), 
               arrow=arrow(length=unit(0.1, "inches")), color="darkred") +
  theme_minimal() +
  scale_x_log10(breaks=10^c(seq(-7,4, by=2)), limits=c(1e-7, 1e4)) +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        strip.text.y = element_text(size=12, angle = 0))
ggsave("type3.log.pdf", plot=pt3log, height=2, width=5, units="in")



# How many Samples to Estimate Dataset 1 with Model IIb? ------------------


# Type I Zeroes

make_dat1 <- function(n){
  lambda1 <- .5
  type1 <- rpois(n, lambda1)
  # return dataframe list
  type1 <- within(list(), {
    y <- type1
    z <- rep(1, length(y))
    x <- rep(1, length(y))
    N <- length(y)
    N_person <- 1
    N_batch <- 1
    
    lambda_prior <- lambda_prior
    theta_prior <- theta_prior
    gamma_prior <- gamma_prior
  })
  return(type1)
}

posteriors <- list()
# n1 <- 100
n2 <- 30
# ns <-  seq(5, n1, by=5)
ns <- c(5, 10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320, 640, 960, 1280)
k <- 1
for (i in ns){
  for (j in 1:n2){
    print(i)
    fit <- stan("m2b.stan", data=make_dat1(i), iter = 4000, chains=1)
    l <- rstan::extract(fit, pars="lambda")$lambda[,1]
    p <- c(i,quantile(l, c(0.025, 0.975)), mean(l))
    names(p) <- c("n","p2.5", "p97.5", "mean")
    posteriors[[k]] <- p
    k <- k+1
  }
}    
posteriors_tidy <- posteriors %>% 
  do.call(rbind,.) %>% 
  as.tibble() %>% 
  mutate(model="m2b")

posteriors <- list()
k <- 1
for (i in ns){
  for (j in 1:n2){
    print(i)
    fit <- stan("m1.stan", data=make_dat1(i), iter = 4000, chains=1)
    l <- rstan::extract(fit, pars="lambda")$lambda[,1]
    p <- c(i,quantile(l, c(0.025, 0.975)), mean(l))
    names(p) <- c("n","p2.5", "p97.5", "mean")
    posteriors[[k]] <- p  
    k <- k+1
  }
}    
posteriors_tidy <- posteriors %>% 
  do.call(rbind,.) %>% 
  as.tibble() %>% 
  mutate(model="m1") %>% 
  bind_rows(posteriors_tidy)

save(posteriors_tidy, file="posteriors_tidy.RData")
# load("posteriors_tidy.RData")

posteriors_tidy %>% 
  filter(mean < 1000) %>% 
  filter(n %in% ns) %>% 
  mutate(model = ifelse(model=="m1", "Model I", "Model IIb")) %>% 
  rename(Model = model) %>% 
  ggplot(aes(x = factor(n), y = mean, fill= Model)) +
  geom_hline(yintercept=0.5, color="darkred") + 
  geom_boxplot() + 
  theme_minimal()+
  xlab("Number of Samples") +
  theme(axis.text = element_text(size=12), 
        axis.title.y = element_blank()) +
  scale_fill_brewer(palette = "Set1")
ggsave("posterior_consistency.pdf", height=4, width=7)


# Posterior correlation of theta and lambda in model IIb ------------------

logit <- function(x){log(x/(1-x))}

# for t1 simulation
tmp <- rstan::extract(fit.t1m2b, pars=c("lambda", "theta", "lp__"))
tmp <- cbind(tmp$lambda, tmp$theta,tmp$lp__)
colnames(tmp) <- c("Lambda", "Theta", "LogPosterior")
cor.test(tmp[,1], tmp[,2], method="spearman")

find_hull <- function(df) df[chull(df[,1], df[,2]), ]
hull90 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .90)) %>% 
  find_hull()

hull80 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .80)) %>% 
  find_hull()

hull95 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .95)) %>% 
  find_hull()

hullt1 <- bind_rows(hull80, hull90, hull95, .id="hull")

tmpt1 <- as.data.frame(tmp)

# For t3 simulation
tmp <- rstan::extract(fit.t3m2b, pars=c("lambda", "theta", "lp__"))
tmp <- cbind(tmp$lambda[,2], tmp$theta[,2], tmp$lp__)
colnames(tmp) <- c("Lambda", "Theta", "LogPosterior")
cor.test(tmp[,1], tmp[,2], method="spearman")

find_hull <- function(df) df[chull(df[,1], df[,2]), ]
hull90 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .90)) %>% 
  find_hull()

hull80 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .80)) %>% 
  find_hull()

hull95 <- as.data.frame(tmp) %>% 
  arrange(Lambda, Theta) %>% 
  filter(LogPosterior > quantile(LogPosterior, .95)) %>% 
  find_hull()

hullt3 <- bind_rows(hull80, hull90, hull95, .id="hull")

tmpt3 <- as.data.frame(tmp)

hull <- bind_rows(hullt1, hullt3, .id="simulation")
tmp <- bind_rows(tmpt1, tmpt3, .id="simulation")

hull <- mutate(hull, simulation=ifelse(simulation==1, "Simulation 1", "Simulation 5"))
tmp <- mutate(tmp, simulation=ifelse(simulation==1, "Simulation 1", "Simulation 5"))

ggplot(tmp, aes(x = Lambda, y = Theta)) +
  geom_point(alpha=0.7) +
  geom_polygon(data = hull, aes(group=hull), fill=NA, color="red", size=1,  alpha=0.5) +
  scale_color_gradient2(low = "black", mid="darkred", high="red", midpoint=-11) +
  facet_grid(~simulation, scales="free") +
  theme_minimal() +
  scale_x_log10() +
  xlab(TeX("$\\lambda$")) +
  ylab(TeX("$\\theta$")) +
  theme(axis.title = element_text(size=13), 
        axis.text = element_text(size=10), 
        strip.text = element_text(size=13))
ggsave("lambda_theta_posterior.pdf", height=4, width=7, units="in")
