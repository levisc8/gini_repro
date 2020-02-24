# IBM Simulations. Each is pulled directyl from Data-driven Modelling of Structured
# Populations, Ellner, Childs & Rees 2016. For now, using the chapter 2 versions
# which are the simplest possible forms of the model.

# This chapter contains IBMs for two species - one monocarpic thistle and
# and one mammal. below are parameters that define each vital rate for the thistle.

# Carlina parameters -------------

library(tibble)
library(reldist) #For Gini coefficient estimation
library(tidyverse)

car_par_true <- c(

  ## survival
  surv_int  =  -0.65,
  surv_z    =   0.75,

  ## flowering
  flow_int  = -18.00,
  flow_z    =   6.9,

  ## growth
  grow_int  =   0.96,
  grow_z    =   0.59,
  grow_sd   =   0.67,

  ## recruit size
  rcsz_int  =   -.08,
  rcsz_sd   =   0.76,

  ## seed size
  seed_int  =   1.00,
  seed_z    =   2.20,

  ## recruitment probability
  p_r       =   0.007
)


##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time

G_z1z <- function(z1, z, m_par) {

  mu <- m.par["grow_int"] + m.par["grow_z"] * z           # mean size next year
  sig <- m.par["grow_sd"]                                 # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)            # pdf that you are size z1 given you were size z

  return(p_den_grow)
}

## Survival function, logistic regression

s_z <- function(z, m_par) {

  linear_p <- m_par["surv_int"] + m_par["surv_z"] * z  # linear predictor
  p        <- 1 / (1 + exp( -linear_p ))               # logistic transformation to probability
  return(p)

}

## Probability of flowering function, logistic regression

p_bz <- function(z, m_par) {

  linear_p <- m_par["flow_int"] + m_par["flow_z"] * z      # linear predictor
  p        <- 1 / (1 + exp( -linear_p ))                   # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(z, m_par) {

  N <- exp(m_par["seed_int"] + m_par["seed_z"] * z)    # seed production of a size z plant

  return(N)
}

## Recruit size pdf

c_0z1 <- function(z1, m_par) {

  mu       <- m_par["rcsz_int"]
  sig      <- m_par["rcsz_sd"]
  p_deRecr <- dnorm(z1, mean = mu, sd = sig)              # pdf of a size z1 recruit

  return(p_deRecr)
}

# Carlina IBM simulation -----------------

yr            <- 1
n_yrs         <- 200
init_pop_size <- 50

# initial population sizes and ages

z   <- rnorm(init_pop_size,
             mean = car_par_true["rcsz_int"],
             sd   = car_par_true["rcsz_sd"])

id  <- seq_len(length(z))

# calculate initial pop size and mean size

pop_size_t <- init_pop_size

# I add types to the NAs so it's easier to remember data types when I'm translating
# to another language (e.g. c++). They function exactly the same NAs.

sim_data <- tibble(
  id    = id,
  z     = z,
  repr  = rep(NA_integer_, length(z)),
  seeds = rep(NA_integer_, length(z)),
  surv  = rep(NA_integer_, length(z)),
  z1    = rep(NA_real_,    length(z)),
  age   = rep(0,           length(z)),
  alive = rep(TRUE,        length(z)),
  yr    = rep(1L,          length(z))
)

pop_size <- dim(sim_data)[1]
id <- sim_data$id
age <- sim_data$age
z <- sim_data$z

# iterate the model using the 'true' parameters and store data in a
# list

y <- 1

set.seed(121)

while(y < n_yrs && dim(sim_data)[1] < 150000) {

  # Temporary data list to store data from a particular year

  temp <- list(
    id    = id,
    z     = z,
    repr  = rep(NA_integer_, pop_size),
    seeds = rep(NA_integer_, pop_size),
    surv  = rep(NA_integer_, pop_size),
    z1    = rep(NA_real_, pop_size),
    age   = age,
    alive = rep(NA, pop_size),
    yr    = rep(y, pop_size)
  )

  # generate binomial random number for the probability of flowering,
  # where the probability of flowering depends on your size z, this is a
  # vector of 0's and 1's, you get a 1 if you flower

  temp$repr <- rbinom(n    = pop_size, #Number of outputs
                      prob = p_bz(temp$z, car_par_true), #Individual probability for each output
                      size = 1) #Number of trials to give a success (just 1)


  # number of plants that flowered

  repr_n <- sum(temp$repr)

  # we'll assume plant make a Poisson distributed number of seeds

  temp$seeds[temp$repr == 1] <- rpois(repr_n, #Number of entries to generate
                                      b_z(temp$z[temp$repr == 1],
                                          car_par_true)) #Believe this is the Poisson-mean(and var) to generate each entry

  # generate the number of recruits

  recr_n <- ifelse(repr_n == 0, #If none managed to reproduce
                   0, #Then equals 0
                   rbinom(n    = 1, #Else generate a single number
                          size = sum(temp$seeds, na.rm = TRUE), #That tells us the sum over all seeds
                          prob = car_par_true["p_r"])) #Of the probability of germinating (identical)
  #^So, note, number of recruits is in fact a function of number of seeds produced

  # generate new recruit sizes and get an index to start their IDs with.
  # Index gets used later

  rc_sz <- rnorm(recr_n,
                 mean = car_par_true["rcsz_int"], #Start their sizes(z1) at the intercept (state 0)
                 sd   = car_par_true["rcsz_sd"])

  rc_id_start <- max(sim_data$id) + 1


  # Vector of expected sizes z1. Gets used later in
  E_z1 <- car_par_true['grow_int'] + car_par_true['grow_z'] * temp$z

  for(i in seq_len(pop_size)) {

    # for the non-reproductive plants generate random number for survival

    temp$surv[i]       <- ifelse(
      as.logical(temp$repr[i]), #If it bred
      0, #Then oops, no survival
      rbinom(n    = 1, #Otherwise, give me one number
             size = 1, #That gives either a success or not (one trial)
             prob = s_z(temp$z[i], car_par_true)) #Based on a probability given by the individual's size
    )

    # let them grow (if they're still alive!)

    temp$z1[i] <- ifelse(
      as.logical(temp$surv[i]), #If it survived
      rnorm(n = 1, #Give me one number
            mean = E_z1[i], #With a mean based on the expected size for that size
            sd = car_par_true["grow_sd"]), #And some randomness (previously established)
      NA_real_ #Otherwise, it gets an NA
    )

    temp$alive[i] <- as.logical(temp$surv[i])

    temp$age[i]   <- temp$age[i] + 1

  }

  # Next, if we have any recruits, append them to temp.
  if(recr_n > 0) {

    rc_id <- seq(rc_id_start,
                 (rc_id_start + recr_n - 1),
                 by = 1)
    temp$id <- c(temp$id, rc_id)
    temp$z  <- c(temp$z, rep(NA_real_,
                             recr_n))
    temp$repr <- c(temp$repr, rep(NA_integer_,
                                  recr_n))
    temp$seeds <- c(temp$seeds, rep(NA_integer_,
                                    recr_n))
    temp$surv  <- c(temp$surv, rep(NA_integer_,
                                   recr_n))
    temp$z1    <- c(temp$z1, rc_sz)
    temp$age   <- c(temp$age,
                    rep(0, recr_n))
    temp$alive <- c(temp$alive, rep(TRUE,
                                    recr_n))
    temp$yr    <- c(temp$yr, rep(y, recr_n))


  }

  pop_size <- sum(temp$alive)

  if(y == 1) {

    sim_data$repr  <- temp$repr
    sim_data$seeds <- temp$seeds
    sim_data$surv  <- temp$surv
    sim_data$z1    <- temp$z1
    sim_data$age   <- temp$age
    sim_data$alive <- temp$alive
    sim_data$yr    <- y

  } else {
    sim_data   <- rbind(sim_data,
                        data.frame(temp))
  }

  id <- sim_data$id[sim_data$yr == y & !is.na(sim_data$z1)]
  z <- sim_data$z1[sim_data$yr == y & !is.na(sim_data$z1)]
  age <- sim_data$age[sim_data$yr == y & !is.na(sim_data$z1)]

  y <- y + 1
}

# quick sanity checks for combinations of variables that shouldn't be possible

stopifnot(! any(sim_data$repr & is.na(sim_data$z), na.rm = TRUE))

# Reproduction is fatal, cannot be repr and have size_next

stopifnot(! any(sim_data$repr & ! is.na(sim_data$z1), na.rm = TRUE))
stopifnot(! any(sim_data$repr & is.na(sim_data$seeds), na.rm = TRUE))

stopifnot(
  ! any(
    is.na(sim_data$z) &
      ! is.na(sim_data$z1) &
      sim_data$age > 0
  )
)

# store the output data. You can read it back in using the snippet commented out
# beneath the saveRDS() call

#saveRDS(sim_data, file = 'data/monocarp_ibm_data.rds')

#monocarp_ibm_data <- readRDS('monocarp_ibm_data.rds')

## Some simple checks ##

sim_data <- sim_data %>%
  group_by(yr) %>%
  mutate(pop_size = n())

# Soay IBM -----

# Soay sheep parameters -----------------

soa_par_true <- c(

  ## survival
  surv_int  = -9.65e+0,
  surv_z    =  3.77e+0,

  ## growth
  grow_int  =  1.41e+0,
  grow_z    =  5.57e-1,
  grow_sd   =  7.99e-2,

  ## reproduce or not
  repr_int  = -7.23e+0,
  repr_z    =  2.60e+0,

  ## recruit or not
  recr_int  =  1.93e+0,

  ## recruit size
  rcsz_int  =  3.62e-1,
  rcsz_z    =  7.09e-1,
  rcsz_sd   =  1.59e-1
)

# Growth function
G_z1z <- function(z1, z, m_par) {

  mu <- m.par["grow_int"] + m.par["grow_z"] * z           # mean size next year
  sig <- m.par["grow_sd"]                                 # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)            # pdf that you are size z1 given you were size z

  return(p_den_grow)
}

## Survival function, logistic regression
s_z <- function(z, m_par) {

  linear_p <- m_par["surv_int"] + m_par["surv_z"] * z  # linear predictor
  p        <- 1 / (1 + exp( -linear_p ))               # logistic transformation to probability
  return(p)

}

## Probability of reproducing, logistic regression
p_bz <- function(z, m_par) {

  linear_p <- m_par["repr_int"] + m_par["repr_z"] * z      # linear predictor
  p        <- 1 / (1 + exp( -linear_p ))                   # logistic transformation to probability
  return(p)
}

# Recruitment function (N_B - from birth in spring to first summer), logistic regression
pr_z <- function(m_par) {

  linear_p <- m_par["recr_int"]           # linear predictor
  p        <- 1 / (1 + exp( -linear_p ))  # logistic transformation to probability

  return(p)
}

# Recruit size function

c_z1z <- function(z1, z, m_par) {

  mean       <- m_par["rcsz_int"] + m_par["rcsz_z"] * z   # mean size of recuits next year
  sd         <- m_par["rcsz_sd"]                          # sd about mean
  p_den_rcsz <- dnorm(z1, mean = mean, sd = sd)           # pdf for offspring

  return(p_den_rcsz)
}

master_data <- data.frame(matrix(ncol = 16, nrow = 0))
names(master_data) <- c("Iteration", "id", "z", "repr", "surv", "z1", "age",
                        "sex", "alive", "yr", "off1", "off2", "gini1", "gini2",
                        "pop_size", 'total_offspring')

set.seed(121)

for (iterate in 1:50) {
  # This one differs slightly in that there is a maternal effect of size
  # on recruit size.

  yr            <- 1
  n_yrs         <- 200
  init_pop_size <- 50

  # initial population sizes and ages. The 3.2 is pulled from the script in
  # the IPM book. I am not entirely sure why they picked that size for the
  # parent of them all.

  z   <- rnorm(init_pop_size,
               mean = soa_par_true["rcsz_int"] + soa_par_true['rcsz_z'] * 3.2,
               sd   = soa_par_true["rcsz_sd"])

  id  <- seq_len(length(z))

  # calculate initial pop size and mean size

  pop_size_t <- init_pop_size

  sim_data <- tibble(
    id    = id,
    z     = z,
    repr  = rep(NA_integer_, length(z)),
    surv  = rep(NA_integer_, length(z)),
    z1    = rep(NA_real_,    length(z)),
    age   = rep(0,           length(z)),
    sex   = rbinom(length(z), prob = 1/2, size = 1),
    alive = rep(TRUE,        length(z)),
    yr    = rep(1L,          length(z)),
    off1  = rep(NA_real_, length(z)),
    off2  = rep(NA_real_, length(z)),
    gini1 = rep(NA_real_, length(z)),
    gini2 = rep(NA_real_, length(z))
  )

  pop_size <- dim(sim_data)[1]
  id       <- sim_data$id
  age      <- sim_data$age
  z        <- sim_data$z
  sex      <- sim_data$sex

  # A constant to define the minimum age an individual has to be to reproduce.
  # this is to keep the model from generating 0 year old individuals
  # who are reproductive (which I think is nonsense). Can/should be adjusted
  # to explore implications

  min_repr_age <- 1

  # iterate the model using the 'true' parameters and store data in a
  # list

  y <- 1

  while(y < n_yrs && dim(sim_data)[1] < 15000) {

    # Temporary data list to store data from a particular year

    temp <- list(
      id    = id,
      z     = z,
      repr  = rep(NA_integer_, pop_size),
      surv  = rep(NA_integer_, pop_size),
      z1    = rep(NA_real_, pop_size),
      age   = age,
      sex   = sex,
      alive = rep(NA, pop_size),
      yr    = rep(y, pop_size),
      off1  = rep(NA_real_, pop_size),
      off2  = rep(NA_real_, pop_size),
      gini1 = rep(NA_real_, pop_size),
      gini2 = rep(NA_real_, pop_size)
    )

    # The IBM here looks a little different because the census timing is post-reproductive.
    # thus, survival and growth happen first, and then we figure out who reproduces.
    # It is helpful that reproduction is not fatal in this example.

    # However, for the purposes of our simulation, I'm going to switch things up a little bit:
    # instead, we're going to calculate their growth, then *potential* reproductive output,
    # followed by a skewing of that reproductive output. Then we're going to kill 'em off,
    # and only allow the reproducers that survived to produce offspring. Make sense?

    temp$surv <- rbinom(n = pop_size,
                        prob = s_z(temp$z,
                                   soa_par_true),
                        size = 1)

    # let them grow (if they're still alive!)

    E_z1                    <- soa_par_true['grow_int'] + soa_par_true['grow_z'] * temp$z

    temp$z1                 <- rnorm(length(z), mean = E_z1, sd = soa_par_true['grow_sd'])

    temp$z1[temp$surv == 0] <- NA_real_

    # generate binomial random number for the probability of reproducing,
    # where the probability of reproducing depends on your size z and sex.

    repr_ind <- which(temp$surv == 1 &
                        temp$sex  == 1 &
                        temp$age  >= min_repr_age)

    temp$repr[repr_ind] <- rbinom(n    = length(repr_ind),
                                  prob = p_bz(temp$z[repr_ind],
                                              soa_par_true),
                                  size = 1)

    # now we add in some skewed reproductive output to fit a Chi-squared distribution. This lets us
    # skew reproduction towards fewer individuals

    if (sum(temp$repr, na.rm = TRUE) > 0) { #Conditioning it on there being any offspring at all
      original    <- temp$repr[which(temp$repr == 1)] #Taking reproductive individuals as producing 1 offspring
      sumOriginal <- sum(original)

      chisqDist <- rchisq(n = length(original), df = 15) #Should increase inequality a bit
      sumChisq  <- sum(chisqDist)
      final     <- chisqDist*sumOriginal/sumChisq
      #^This standardizes by the sum so that the total offspring are equal to the original (no growth difference)

      temp$off1[which(temp$repr == 1)]  <- final
      temp$gini1[which(temp$repr == 1)] <- gini(final)

      chisqDist <- rchisq(n = length(original), df = 5) #Should increase inequality by more
      sumChisq  <- sum(chisqDist)
      final     <- chisqDist*sumOriginal/sumChisq

      temp$off2[which(temp$repr == 1)] <- final
      temp$gini2[which(temp$repr == 1)] <- gini(final)

      temp$repr <- temp$off1
    }

    # Time for a stochastic mortality event. Roughly every ten years, we'll get a random ~10% mortality event

    diceRoll <- sample(seq(1, 10), size = 1)

    if (diceRoll == 10) {

      randomDead   <- rep(0, length(temp$surv))
      aliveThusFar <- which(temp$surv == 1)

      randomDead[aliveThusFar] <- rbinom(n = length(aliveThusFar),
                                         size = 1,
                                         prob = 0.9) #Kills off ~10% of the remaining alive

      differences <- which(temp$surv == 1 & randomDead == 0) #Ones that need correcting

      temp$surv[differences] <- 0
      temp$z1[differences] <- NA_real_
      temp$repr[differences] <- NA_integer_
    }

    # number of fertile adults

    sex_ind  <- which(temp$repr > 0)
    #^Changed this from "== 1" because of non-integer values in post-transform reproduction
    recr_sex <- rep(NA_integer_, pop_size)

    recr_sex[sex_ind] <- rbinom(n    = length(sex_ind),
                                size = 1,
                                prob = 0.5)

    # generate recruits, assign size. Another key difference
    # from this model vs the book model is that they discard males, while this
    # model retains them.

    real_recr_ind <- which(!is.na(recr_sex))
    realized_recr <- rep(NA_integer_, pop_size)

    realized_recr[real_recr_ind] <- temp$repr[real_recr_ind]*rbinom(n    = length(real_recr_ind),
                                                                    size = 1,
                                                                    prob = pr_z(soa_par_true))
    #^Here we now multiply the 0/1 of survival by the total offspring (for post-standardized cases)

    use_z_ind     <- which(temp$sex  == 1 & temp$surv     == 1 &
                             temp$repr > 0 & realized_recr > 0)

    total_realized_offspring <- round(sum(temp$repr[use_z_ind], na.rm = TRUE))

    # generate new recruit sizes and get an index to start their IDs with.
    # Index gets used later

    rc_sz <- c()

    if (length(use_z_ind) > 1 & total_realized_offspring > 0) {
      for (offspring in 1:total_realized_offspring) {
        sampledSize <- sample(temp$z[use_z_ind],
                              size = 1,
                              prob = temp$repr[use_z_ind])
        #^Sample a reproductive adult size with probability proportional to its allocation of reproduction
        rc_sz       <- c(rc_sz, rnorm(n = 1,
                                      mean = soa_par_true["rcsz_int"] +
                                        soa_par_true["rcsz_z"]*sampledSize,
                                      sd = soa_par_true["rcsz_sd"]))
        #^Generate offspring size based off that adult
      }
    } else if (length(use_z_ind) == 1 & total_realized_offspring > 0) {
      for (offspring in 1:total_realized_offspring) {
        sampledSize <- temp$z[use_z_ind]
        rc_sz       <- c(rc_sz, rnorm(n = 1,
                                      mean = soa_par_true["rcsz_int"] +
                                        soa_par_true["rcsz_z"]*sampledSize,
                                      sd = soa_par_true["rcsz_sd"]))
        #^This is specifically when there's only one parent; "sample" yells at you for providing one probability
      }

    }

    # They lived! Time to increment the age variable.
    temp$age[temp$surv == 1] <- temp$age[temp$surv == 1] + 1

    if(length(rc_sz) > 0) {

      rc_id_start <- max(sim_data$id) + 1
      recr_n      <- length(rc_sz)

      # Now, append all the info to temp
      rc_id <- seq(rc_id_start,
                   (rc_id_start + length(rc_sz) - 1),
                   by = 1)
      temp$id <- c(temp$id, rc_id)
      temp$z  <- c(temp$z, rep(NA_real_, recr_n))
      temp$repr <- c(temp$repr,
                     rep(NA_integer_,
                         recr_n))
      temp$surv  <- c(temp$surv, rep(NA_integer_,
                                     recr_n))
      temp$z1    <- c(temp$z1, rc_sz)
      temp$age   <- c(temp$age,
                      rep(0, recr_n))
      temp$alive <- c(temp$alive,
                      rep(TRUE,
                          recr_n))
      if (total_realized_offspring <= length(use_z_ind)) {
        temp$sex <- c(temp$sex,
                      recr_sex[sample(use_z_ind, size = total_realized_offspring)])
      } else {
        temp$sex <- c(temp$sex,
                      recr_sex[sample(use_z_ind, size = total_realized_offspring, replace = TRUE)])
      }
      temp$yr    <- c(temp$yr,
                      rep(y, recr_n))
      temp$off1  <- c(temp$off1,
                      rep(NA_real_, recr_n))
      temp$off2  <- c(temp$off2,
                      rep(NA_real_, recr_n))
      temp$gini1 <- c(temp$gini1,
                      rep(NA_real_, recr_n))
      temp$gini2 <- c(temp$gini2,
                      rep(NA_real_, recr_n))



    }


    temp$alive             <- as.logical(temp$surv)

    if(y == 1) {

      sim_data$repr  <- temp$repr
      sim_data$surv  <- temp$surv
      sim_data$z1    <- temp$z1
      sim_data$age   <- temp$age
      sim_data$alive <- temp$alive
      sim_data$sex   <- temp$sex
      sim_data$yr    <- y

    } else {
      sim_data   <- rbind(sim_data,
                          data.frame(temp))
    }

    # Set up next year's simulation

    id       <- sim_data$id[sim_data$yr == y & !is.na(sim_data$z1)]
    z        <- sim_data$z1[sim_data$yr == y & !is.na(sim_data$z1)]
    age      <- sim_data$age[sim_data$yr == y & !is.na(sim_data$z1)]
    sex      <- sim_data$sex[sim_data$yr == y & !is.na(sim_data$z1)]
    pop_size <- length(sim_data$z[sim_data$yr == y & sim_data$alive])

    y <- y + 1

  }

  # Series of sanity checks for the data. These are quick and dirty - you will
  # want to check them a bit more thoroughly I'd imagine.
  stopifnot(! any(sim_data$repr == 1 & sim_data$age < min_repr_age, na.rm = TRUE))
  stopifnot(! any(sim_data$repr == 1 & sim_data$sex == 0, na.rm = TRUE))
  stopifnot(! any(sim_data$repr == 1 & is.na(sim_data$z1), na.rm = TRUE))
  stopifnot(! any(sim_data$repr == 1 & sim_data$surv == 0, na.rm = TRUE))

  Iteration <- rep(iterate, nrow(sim_data))
  sim_data <- cbind(Iteration, sim_data)

  master_data <- rbind(master_data, sim_data)
}

master_data <- master_data %>%
  group_by(Iteration, yr) %>%
  mutate(pop_size = n())

ggplot(master_data, aes(x = yr, y = pop_size, group = Iteration)) + geom_line()

## Here is the original code, pre-iteration ##

# This one differs slightly in that there is a maternal effect of size
# on recruit size.

yr            <- 1
n_yrs         <- 200
init_pop_size <- 50

# initial population sizes and ages. The 3.2 is pulled from the script in
# the IPM book. I am not entirely sure why they picked that size for the
# parent of them all.

set.seed(121)

z   <- rnorm(init_pop_size,
             mean = soa_par_true["rcsz_int"] + soa_par_true['rcsz_z'] * 3.2,
             sd   = soa_par_true["rcsz_sd"])

id  <- seq_len(length(z))

# calculate initial pop size and mean size

pop_size_t <- init_pop_size

sim_data <- tibble(
  id    = id,
  z     = z,
  repr  = rep(NA_integer_, length(z)),
  surv  = rep(NA_integer_, length(z)),
  z1    = rep(NA_real_,    length(z)),
  age   = rep(0,           length(z)),
  sex   = rbinom(length(z), prob = 1/2, size = 1),
  alive = rep(TRUE,        length(z)),
  yr    = rep(1L,          length(z)),
  off1  = rep(NA_real_, length(z)),
  off2  = rep(NA_real_, length(z)),
  gini1 = rep(NA_real_, length(z)),
  gini2 = rep(NA_real_, length(z))
)

pop_size <- dim(sim_data)[1]
id       <- sim_data$id
age      <- sim_data$age
z        <- sim_data$z
sex      <- sim_data$sex

# A constant to define the minimum age an individual has to be to reproduce.
# this is to keep the model from generating 0 year old individuals
# who are reproductive (which I think is nonsense). Can/should be adjusted
# to explore implications

min_repr_age <- 1

# iterate the model using the 'true' parameters and store data in a
# list

y <- 1

set.seed(121)

while(y < n_yrs && dim(sim_data)[1] < 15000) {

  # Temporary data list to store data from a particular year

  temp <- list(
    id    = id,
    z     = z,
    repr  = rep(NA_integer_, pop_size),
    surv  = rep(NA_integer_, pop_size),
    z1    = rep(NA_real_, pop_size),
    age   = age,
    sex   = sex,
    alive = rep(NA, pop_size),
    yr    = rep(y, pop_size),
    off1  = rep(NA_real_, pop_size),
    off2  = rep(NA_real_, pop_size),
    gini1 = rep(NA_real_, pop_size),
    gini2 = rep(NA_real_, pop_size)
  )

  # The IBM here looks a little different because the census timing is post-reproductive.
  # thus, survival and growth happen first, and then we figure out who reproduces.
  # It is helpful that reproduction is not fatal in this example.

  # However, for the purposes of our simulation, I'm going to switch things up a little bit:
  # instead, we're going to calculate their growth, then *potential* reproductive output,
  # followed by a skewing of that reproductive output. Then we're going to kill 'em off,
  # and only allow the reproducers that survived to produce offspring. Make sense?

  temp$surv <- rbinom(n = pop_size,
                      prob = s_z(temp$z,
                                 soa_par_true),
                      size = 1)

  # let them grow (if they're still alive!)

  E_z1                    <- soa_par_true['grow_int'] + soa_par_true['grow_z'] * temp$z

  temp$z1                 <- rnorm(length(z), mean = E_z1, sd = soa_par_true['grow_sd'])

  temp$z1[temp$surv == 0] <- NA_real_

  # generate binomial random number for the probability of reproducing,
  # where the probability of reproducing depends on your size z and sex.

  repr_ind <- which(temp$surv == 1 &
                      temp$sex  == 1 &
                      temp$age  >= min_repr_age)

  temp$repr[repr_ind] <- rbinom(n    = length(repr_ind),
                                prob = p_bz(temp$z[repr_ind],
                                            soa_par_true),
                                size = 1)

  # now we add in some skewed reproductive output to fit a Chi-squared distribution. This lets us
  # skew reproduction towards fewer individuals

  if (sum(temp$repr, na.rm = TRUE) > 0) { #Conditioning it on there being any offspring at all
    original    <- temp$repr[which(temp$repr == 1)] #Taking reproductive individuals as producing 1 offspring
    sumOriginal <- sum(original)

    chisqDist <- rchisq(n = length(original), df = 15) #Should increase inequality a bit
    sumChisq  <- sum(chisqDist)
    final     <- chisqDist*sumOriginal/sumChisq
    #^This standardizes by the sum so that the total offspring are equal to the original (no growth difference)

    temp$off1[which(temp$repr == 1)]  <- final
    temp$gini1[which(temp$repr == 1)] <- gini(final)

    chisqDist <- rchisq(n = length(original), df = 5) #Should increase inequality by more
    sumChisq  <- sum(chisqDist)
    final     <- chisqDist*sumOriginal/sumChisq

    temp$off2[which(temp$repr == 1)] <- final
    temp$gini2[which(temp$repr == 1)] <- gini(final)
  }

  # Time for a stochastic mortality event. Roughly every ten years, we'll get a random ~10% mortality event

  diceRoll <- sample(seq(1, 10), size = 1)

  if (diceRoll == 10) {

    randomDead   <- rep(0, length(temp$surv))
    aliveThusFar <- which(temp$surv == 1)

    randomDead[aliveThusFar] <- rbinom(n = length(aliveThusFar),
                                       size = 1,
                                       prob = 0.9) #Kills off ~10% of the remaining alive

    differences <- which(temp$surv == 1 & randomDead == 0) #Ones that need correcting

    temp$surv[differences] <- 0
    temp$z1[differences] <- NA_real_
    temp$repr[differences] <- NA_integer_
  }

  # number of fertile adults

  sex_ind  <- which(temp$repr > 0)
  #^Changed this from "== 1" because of non-integer values in post-transform reproduction
  recr_sex <- rep(NA_integer_, pop_size)

  recr_sex[sex_ind] <- rbinom(n    = length(sex_ind),
                              size = 1,
                              prob = 0.5)

  # generate recruits, assign size. Another key difference
  # from this model vs the book model is that they discard males, while this
  # model retains them.

  real_recr_ind <- which(!is.na(recr_sex))
  realized_recr <- rep(NA_integer_, pop_size)

  realized_recr[real_recr_ind] <- temp$repr[real_recr_ind]*rbinom(n    = length(real_recr_ind),
                                                                  size = 1,
                                                                  prob = pr_z(soa_par_true))
  #^Here we now multiply the 0/1 of survival by the total offspring (for post-standardized cases)

  use_z_ind     <- which(temp$sex  == 1 & temp$surv     == 1 &
                           temp$repr > 0 & realized_recr > 0)

  total_realized_offspring <- round(sum(temp$repr[use_z_ind], na.rm = TRUE))

  # generate new recruit sizes and get an index to start their IDs with.
  # Index gets used later

  rc_sz <- c()

  if (total_realized_offspring > 1) {
    for (offspring in 1:total_realized_offspring) {
      sampledSize <- sample(temp$z[use_z_ind],
                            size = 1,
                            prob = temp$repr[use_z_ind])
      #^Sample a reproductive adult size with probability proportional to its allocation of reproduction
      rc_sz       <- c(rc_sz, rnorm(n = 1,
                                    mean = soa_par_true["rcsz_int"] +
                                      soa_par_true["rcsz_z"]*sampledSize,
                                    sd = soa_par_true["rcsz_sd"]))
      #^Generate offspring size based off that adult
    }
  } else if (total_realized_offspring == 1) {
    sampledSize <- temp$z[use_z_ind]
    rc_sz       <- c(rc_sz, rnorm(n = 1,
                                  mean = soa_par_true["rcsz_int"] +
                                    soa_par_true["rcsz_z"]*sampledSize,
                                  sd = soa_par_true["rcsz_sd"]))
    #^This is specifically when there's only one parent; "sample" yells at you for providing one probability
  }

  # They lived! Time to increment the age variable.
  temp$age[temp$surv == 1] <- temp$age[temp$surv == 1] + 1

  if(length(rc_sz) > 0) {

    rc_id_start <- max(sim_data$id) + 1
    recr_n      <- length(rc_sz)

    # Now, append all the info to temp
    rc_id <- seq(rc_id_start,
                 (rc_id_start + length(rc_sz) - 1),
                 by = 1)
    temp$id <- c(temp$id, rc_id)
    temp$z  <- c(temp$z, rep(NA_real_, recr_n))
    temp$repr <- c(temp$repr,
                   rep(NA_integer_,
                       recr_n))
    temp$surv  <- c(temp$surv, rep(NA_integer_,
                                   recr_n))
    temp$z1    <- c(temp$z1, rc_sz)
    temp$age   <- c(temp$age,
                    rep(0, recr_n))
    temp$alive <- c(temp$alive,
                    rep(TRUE,
                        recr_n))
    temp$sex   <- c(temp$sex,
                    recr_sex[sample(use_z_ind, size = total_realized_offspring)])
    temp$yr    <- c(temp$yr,
                    rep(y, recr_n))
    temp$off1  <- c(temp$off1,
                    rep(NA_real_, recr_n))
    temp$off2  <- c(temp$off2,
                    rep(NA_real_, recr_n))
    temp$gini1 <- c(temp$gini1,
                    rep(NA_real_, recr_n))
    temp$gini2 <- c(temp$gini2,
                    rep(NA_real_, recr_n))



  }


  temp$alive             <- as.logical(temp$surv)

  if(y == 1) {

    sim_data$repr  <- temp$repr
    sim_data$surv  <- temp$surv
    sim_data$z1    <- temp$z1
    sim_data$age   <- temp$age
    sim_data$alive <- temp$alive
    sim_data$sex   <- temp$sex
    sim_data$yr    <- y

  } else {
    sim_data   <- rbind(sim_data,
                        data.frame(temp))
  }

  # Set up next year's simulation

  id       <- sim_data$id[sim_data$yr == y & !is.na(sim_data$z1)]
  z        <- sim_data$z1[sim_data$yr == y & !is.na(sim_data$z1)]
  age      <- sim_data$age[sim_data$yr == y & !is.na(sim_data$z1)]
  sex      <- sim_data$sex[sim_data$yr == y & !is.na(sim_data$z1)]
  pop_size <- length(sim_data$z[sim_data$yr == y & sim_data$alive])

  y <- y + 1

}

# Series of sanity checks for the data. These are quick and dirty - you will
# want to check them a bit more thoroughly I'd imagine.
stopifnot(! any(sim_data$repr == 1 & sim_data$age < min_repr_age, na.rm = TRUE))
stopifnot(! any(sim_data$repr == 1 & sim_data$sex == 0, na.rm = TRUE))
stopifnot(! any(sim_data$repr == 1 & is.na(sim_data$z1), na.rm = TRUE))
stopifnot(! any(sim_data$repr == 1 & sim_data$surv == 0, na.rm = TRUE))

saveRDS(sim_data, file = 'data/soay_ibm_data.rds')

## Simple checks ##

sim_data <- sim_data %>%
  group_by(yr) %>%
  mutate(pop_size = n(),
         total_offspring = sum(off1, na.rm = T))

ggplot(sim_data, aes(x = yr, y = pop_size)) + geom_point()

## Converting a series of numbers into an equal-sum chisq-distribution ##

original <- rep(1, 100)
sumOriginal <- sum(original)

chisqDist <- rchisq(n = length(original), df = 5)
sumChisq <- sum(chisqDist)
final <- chisqDist*sumOriginal/sumChisq

gini(final)

# Iceplant IBM -----------
# Carpobrotus spp from Israeli drone data

set.seed(121)

ice_par_true <- c(
  grow_int   = 0.01566,
  grow_z     = 0.89849,
  grow_sd_z  = -0.1754,
  surv_int   = 2.7897,
  surv_z     = 0.7312,
  flow_int   = 1.021,
  flow_z     = 1.103,
  flow_n_int = 1.664,
  flow_n_z   = 0.753,
  est_p      = 0.0084,
  recr_mean  = -3.103,
  recr_sd    = 1.064
)

# some helper functions

s_z <- function(z, par) {

  lin_p <- par['surv_int'] + par['surv_z'] * z
  return(1/(1 + exp( - lin_p )))

}

flow_z <- function(z, par) {

  lin_p <- par['flow_int'] + par['flow_z'] * z
  return(1/(1 + exp( - lin_p )))

}

flow_n_z <- function(z, par) {

  exp(par['flow_n_int'] + par['flow_n_z'] * z)

  lin_p <- par['surv_int'] + par['surv_z'] * z
  return(1/(1 + exp( - lin_p )))

}

flow_z <- function(z, par) {

  lin_p <- par['flow_int'] + par['flow_z'] * z
  return(1/(1 + exp( - lin_p )))

}

flow_n_z <- function(z, par) {

  exp(par['flow_n_int'] + par['flow_n_z'] * z)

}

# Set up the simulation parameters

y             <- 1
n_yrs         <- 200
init_pop_size <- 50

# initial population sizes and ages. We assume they all come from seeds

z   <- rnorm(init_pop_size,
             mean = ice_par_true["recr_mean"],
             sd   = ice_par_true["recr_sd"])

id  <- seq_len(length(z))

# calculate initial pop size and mean size

pop_size_t <- init_pop_size

sim_data <- tibble(
  id     = id,
  z      = z,
  repr   = rep(NA_integer_, length(z)),
  surv   = rep(NA_integer_, length(z)),
  flow_n = rep(NA_integer_, length(z)),
  z1     = rep(NA_real_,    length(z)),
  age    = rep(0,           length(z)),
  alive  = rep(TRUE,        length(z)),
  yr     = rep(1L,          length(z))
)

pop_size <- dim(sim_data)[1]
id       <- sim_data$id
age      <- sim_data$age
z        <- sim_data$z


while(y < n_yrs & dim(sim_data)[1] < 100000) {

  # Hold the data somewhere temporarily
  temp <- list(
    id     = id,
    z      = z,
    repr   = rep(NA_integer_, length(z)),
    surv   = rep(NA_integer_, length(z)),
    flow_n = rep(NA_integer_, length(z)),
    z1     = rep(NA_integer_, length(z)),
    age    = age,
    alive  = rep(TRUE,        length(z)),
    yr     = rep(y,           length(z))
  )

  # Reproduction happens first in this model

  temp$repr <- rbinom(n    = pop_size,
                      prob = flow_z(z, ice_par_true),
                      size = 1)

  flow_ind  <- which(temp$repr == 1)
  repr_n    <- length(flow_ind)

  # Get the number of recruits. We'll assign them ages, etc in a minute.

  temp$flow_n[flow_ind] <- rpois(repr_n,
                                 flow_n_z(temp$z[flow_ind],
                                          ice_par_true))

  flow_tot              <- sum(temp$flow_n, na.rm = TRUE)

  n_recr                <- rbinom(n    = 1,
                                  size = flow_tot,
                                  prob = ice_par_true['est_p'])
  # The plants must now survive and grow.

  temp$surv             <- rbinom(n    = pop_size,
                                  prob = s_z(z, ice_par_true),
                                  size = 1)

  surv_ind              <- which(temp$surv == 1)

  temp$alive            <- as.logical(temp$surv)

  # Ice plant differs from the others in that growth becomes more determinstic
  # as the plants get bigger. The SD is a function of size in this model, and the
  # coefficient is negative

  E_z1                  <- ice_par_true['grow_int'] +
    ice_par_true['grow_z'] *
    temp$z[surv_ind]
  sd_z1                 <- sqrt(exp(2 * ice_par_true['grow_sd_z'] * temp$z[surv_ind]))

  temp$z1[surv_ind]     <- rnorm(n    = length(surv_ind),
                                 mean = E_z1,
                                 sd   = sd_z1)

  # Now, if there are any, give recruits a size and package them up
  # to fit with the surviving adults

  if(n_recr > 0) {

    recr_z      <- rnorm(n_recr,
                         mean = ice_par_true['recr_mean'],
                         sd   = ice_par_true['recr_sd'])

    rc_id_start <- max(sim_data$id)

    rc_id       <- seq(rc_id_start,
                       (rc_id_start + n_recr - 1),
                       by = 1)

    temp$id     <- c(temp$id, rc_id)
    temp$z      <- c(temp$z, rep(NA_real_,
                                 n_recr))
    temp$repr   <- c(temp$repr, rep(NA_integer_,
                                    n_recr))
    temp$flow_n <- c(temp$flow_n, rep(NA_integer_,
                                      n_recr))
    temp$surv   <- c(temp$surv, rep(NA_integer_,
                                    n_recr))
    temp$z1     <- c(temp$z1, recr_z)
    temp$age    <- c(temp$age,
                     rep(0, n_recr))
    temp$alive  <- c(temp$alive, rep(TRUE,
                                     n_recr))
    temp$yr     <- c(temp$yr, rep(y, n_recr))


  }


  pop_size <- sum(temp$alive)

  if(y == 1) {

    sim_data$repr   <- temp$repr
    sim_data$surv   <- temp$surv
    sim_data$flow_n <- temp$flow_n
    sim_data$z1     <- temp$z1
    sim_data$age    <- temp$age
    sim_data$alive  <- temp$alive
    sim_data$yr     <- y

  } else {
    sim_data   <- rbind(sim_data,
                        data.frame(temp))
  }

  id  <- sim_data$id[sim_data$yr == y & !is.na(sim_data$z1)]
  z   <- sim_data$z1[sim_data$yr == y & !is.na(sim_data$z1)]
  age <- sim_data$age[sim_data$yr == y & !is.na(sim_data$z1)] + 1

  y <- y + 1

}

# Plants can die after reproducing, so we don't have nearly as many things to
# check (I don't think, but may regret those words later...)
stopifnot(! any(sim_data$repr == 1 & is.na(sim_data$flow_n), na.rm = TRUE))
stopifnot(! any(sim_data$repr == 1 & is.na(sim_data$z), na.rm = TRUE))

saveRDS(sim_data, 'data/carpobrotus_ibm_data.rds')
