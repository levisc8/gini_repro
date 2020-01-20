# IBM Simulations. Each is pulled directyl from Data-driven Modelling of Structured
# Populations, Ellner, Childs & Rees 2016. For now, using the chapter 2 versions
# which are the simplest possible forms of the model.

# This chapter contains IBMs for two species - one monocarpic thistle and
# and one mammal. below are parameters that define each vital rate for the thistle.

# Carlina parameters -------------

library(tibble)

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

  temp$repr <- rbinom(n    = pop_size,
                      prob = p_bz(temp$z, car_par_true),
                      size = 1)


  # number of plants that flowered

  repr_n <- sum(temp$repr)

  # we'll assume plant make a Poisson distributed number of seeds

  temp$seeds[temp$repr == 1] <- rpois(repr_n,
                                      b_z(temp$z[temp$repr == 1],
                                          car_par_true))

  # generate the number of recruits

  recr_n <- ifelse(repr_n == 0,
                   0,
                   rbinom(n    = 1,
                          size = sum(temp$seeds, na.rm = TRUE),
                          prob = car_par_true["p_r"]))

  # generate new recruit sizes and get an index to start their IDs with.
  # Index gets used later

  rc_sz <- rnorm(recr_n,
                 mean = car_par_true["rcsz_int"],
                 sd   = car_par_true["rcsz_sd"])

  rc_id_start <- max(sim_data$id) + 1


  # Vector of expected sizes z1. Gets used later in
  E_z1 <- car_par_true['grow_int'] + car_par_true['grow_z'] * temp$z

  for(i in seq_len(pop_size)) {

    # for the non-reproductive plants generate random number for survival

    temp$surv[i]       <- ifelse(
      as.logical(temp$repr[i]),
      0,
      rbinom(n    = 1,
             size = 1,
             prob = s_z(temp$z[i], car_par_true))
    )

    # let them grow (if they're still alive!)

    temp$z1[i] <- ifelse(
      as.logical(temp$surv[i]),
      rnorm(n = 1,
            mean = E_z1[i],
            sd = car_par_true["grow_sd"]),
      NA_real_
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

saveRDS(sim_data, file = 'data/monocarp_ibm_data.rds')

# monocarp_ibm_data <- readRDS('data/monocarp_ibm_data.rds')

# soay IBM -----

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

# This one differs slightly in that there is a maternal affect of size
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
  yr    = rep(1L,          length(z))
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
    yr    = rep(y, pop_size)
  )

  # The IBM here looks a little different because the census timing is post-reproductive.
  # thus, survival and growth happen first, and then we figure out who reproduces.
  # It is helpful that reproduction is not fatal in this example.

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

  # number of fertile adults

  repr_n <- sum(temp$repr, na.rm = TRUE)

  sex_ind  <- which(temp$repr == 1)
  recr_sex <- rep(NA_integer_, pop_size)

  recr_sex[sex_ind] <- rbinom(n    = length(sex_ind),
                              size = 1,
                              prob = 0.5)

  # generate recruits, assign size. Another key difference
  # from this model vs the book model is that they discard males, while this
  # model retains them.

  real_recr_ind <- which(!is.na(recr_sex))
  realized_recr <- rep(NA_integer_, pop_size)

  realized_recr[real_recr_ind] <- rbinom(n    = length(real_recr_ind),
                                         size = 1,
                                         prob = pr_z(soa_par_true))

  use_z_ind     <- which(temp$sex  == 1 & temp$surv     == 1 &
                         temp$repr == 1 & realized_recr == 1)

  # generate new recruit sizes and get an index to start their IDs with.
  # Index gets used later

  rc_sz <- rnorm(length(use_z_ind),
                 mean = soa_par_true["rcsz_int"] +
                   soa_par_true['rcsz_z'] * temp$z[use_z_ind],
                 sd   = soa_par_true["rcsz_sd"])

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
                    recr_sex[use_z_ind])
    temp$yr    <- c(temp$yr,
                    rep(y, recr_n))


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
