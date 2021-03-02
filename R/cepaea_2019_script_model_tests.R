library(cmdstanr) ## Stan backend (another)
library(brms) ## the interface we are using
library(tidyverse)
library(bayesplot)
library(tidybayes)
library(matrixStats)
library(here)

## some useful default settings
options(mc.cores = 4) ## reduce/increase depending on cores available
N_chains <- 4
N_warmup <- 2000 
N_iter <- N_warmup + 2000 


raw_behaviour <- read_csv(here("data","cepaea_2019_behaviour.csv"))

data <- raw_behaviour %>% 
  filter(dead==0) %>% 
  #remove the dead individual from the dataset
  #we then standardize covariates
  mutate(is.shaded = as.numeric(population=="shaded")-mean(as.numeric(population=="is.shaded")),
         scale_size = scale(shell_diameter)[,1])

# Testing model specs

## the movement model

#we then define our movement variable
data_test <- data %>% 
  mutate(has.moved = FPT1 < 1200) %>% 
  mutate(mvt = FPT3-FPT1)

#let's try to fit a simple model with this movement variable. 
# movement data is time, so a lognormal model should be right-ish, no?:

mod_mvt_test <- brm(
  bf(mvt|subset(has.moved==TRUE) + ## we exclude non-moving snails 
         cens(FPT3>=1200) ~ is.shaded * (band_number + fusion) + 
    scale_size + (1|series)+(1|fullID)), family=lognormal,
    data=data_test, chains=4,iter=2000,warmup=1000,
    prior=c(set_prior("normal(0,1)",class="Intercept"), # priors are probably bad but this is just for testing
                set_prior("normal(0,1)",class="b"),
                set_prior("normal(0,1)", class="sd"),
                set_prior("normal(0,1)",class="sigma")),
        backend="cmdstanr",seed=42,control=list(adapt_delta=0.99,max_treedepth=15))

newdata <- subset(data_test, has.moved ==TRUE)
ppc_ribbon(
  yrep = (predict(mod_mvt_test, newdata=newdata, summary = FALSE)),
  x = rank(predict(mod_mvt_test, newdata=newdata)[, 1]),
  y = newdata$mvt,
  prob = 0.5, prob_outer = 0.95
)

# fairly easy to see that this model is rubbish 
# and doesn't even come close to capturing the very slow individuals, 
# even with individual-level random effects and fixed covariates. 
# Changing distributions for e.g. gamma or others (not shown) does not solve the problem. 
# Let's try something else, by inverting the movement variable, so that it becomes something closer to a speed:


data_test$mvt2 <- 1/data_test$mvt
data_test$mvt3 <- (data_test$mvt2 - mean(subset(data_test$mvt2,data_test$has.moved==TRUE)))/
  sd(subset(data_test$mvt2,data_test$has.moved==TRUE))
## scaling is obviously weird because of the censoring,
## but let's do it anyway to stabilise model


mod_mvt_test2 <- brm(
  bf(mvt3|subset(has.moved==TRUE) + ## we exclude non-moving snails 
         cens(- 1*(FPT3>=1200)) ~ # because the mvt metric is inverted, the censure is too (right- to left-) 
       is.shaded * (band_number + fusion) + 
    scale_size + (1|series)+(1|fullID)),
    data=data_test, chains=4,iter=8000,warmup=4000,
    prior=c(set_prior("normal(0,1)",class="Intercept"),
                set_prior("normal(0,1)",class="b"),
                set_prior("normal(0,1)", class="sd"),
                set_prior("normal(0,1)",class="sigma")),
        backend="cmdstanr",seed=42,control=list(adapt_delta=0.999,max_treedepth=10))

newdata <- subset(data_test, has.moved ==TRUE)
ppc_ribbon(
  yrep = (predict(mod_mvt_test2, newdata=newdata, summary = FALSE)),
  x = rank(predict(mod_mvt_test2, newdata=newdata)[, 1]),
  y = newdata$mvt3,
  prob = 0.5, prob_outer = 0.95
)

# seems to behave much much better. We'll thus analyse movement as 1/time in the final model.

## the food model

## let's try to do something similar with the food data
## food data is bounded between 0 and 1.5 g, 
# so a Beta model can be appropriate,
# and to go around the beta restrictions on 0 and 100% values, 
# we can reasonably say any 0 is merely "below the detection limit"
# with a left-censoring indicator

# let's try that

data_test <- data_test %>% 
  mutate(food = replace(Food_intake, which(Food_intake ==0), 0.0001)) %>% # replace 0s by "at the detection limit" 
  mutate(food = food/1.5) # convert to proportion of total food eaten
  
mod_food_test <- brm(
    bf(food|cens(- 1*(food<=0.0001)) ~  
         is.shaded * (band_number + fusion) + 
         scale_size + (1|series)+(1|fullID)), family=Beta,
    data=data_test, chains=4,iter=2000,warmup=1000,
    prior=c(set_prior("normal(0,1)",class="Intercept"),
            set_prior("normal(0,1)",class="b"),
            set_prior("normal(0,1)", class="sd")),
    backend="cmdstanr",seed=42,control=list(adapt_delta=0.99,max_treedepth=15))

ppc_ribbon(
  yrep = (predict(mod_food_test, summary = FALSE)),
  x = rank(predict(mod_food_test)[, 1]),
  y = data_test$food,
  prob = 0.5, prob_outer = 0.95
)
## looks good?
pp_check(mod_food_test) # ah nope looks bad, predicts an excess of low high values

# so let's just try to scale the raw values and put that in a gaussian?
data_test <- data_test %>% 
  mutate(food2 = scale(Food_intake)[,1])

mod_food_test2 <- brm(
  bf(food2 ~  
       is.shaded * (band_number + fusion) + 
       scale_size + (1|series)+(1|fullID)),
  data=data_test, chains=4,iter=2000,warmup=1000,
  prior=c(set_prior("normal(0,1)",class="Intercept"),
          set_prior("normal(0,1)",class="b"),
          set_prior("normal(0,1)", class="sd")),
  backend="cmdstanr",seed=42,control=list(adapt_delta=0.99,max_treedepth=15))

ppc_ribbon(
  yrep = (predict(mod_food_test2, summary = FALSE)),
  x = rank(predict(mod_food_test2)[, 1]),
  y = data_test$food2,
  prob = 0.5, prob_outer = 0.95
)
## looks good?
pp_check(mod_food_test2) # ah nope looks bad, predicts an excess of low high values

## seems much much better, even with disregarding the lower boundary
