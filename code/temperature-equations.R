
# convert degrees C to Kelvin for Arrhenius equation
C_to_K <- function(C) C + 273.15 

# Euler's number
euler <- exp(1) 

# Boltzmann constant, units = eV (electron volt)
k <- 8.6173303 * 10^(-5) 

# Arrhenius equation
arrhenius_eqn <- function(E_a,   # activation energy of metabolic processes, units = eV/K
                          T)     # temperature in Kelvin
{
  euler^(-E_a/(k*T))
}

# initial intrinsic growth rate, r (can think of as the intercept)
r0 <- function(r_base, # baseline r for a given temperature 
               E_B,    # activation energy for metabolism in units of eV 
               T){
  # r_base = r0 * arrhenius_eqn to get observed r value, so:
  r_base / arrhenius_eqn(E_a = E_B, T = T) # = r0
}

# temperature scaling of intrinsic growth rate
r_scaling <- function(r0, E_B, T){
  r0 * arrhenius_eqn(E_a = E_B, T = T)
}

# initial carrying capacity, K
K0 <- function(K_base, # baseline K for a given temperature
               E_B,
               E_S,    # activation energy for resource supply rate in units of eV
               T){
  K_base / euler^(-1*(E_S-E_B)/(k*T))
}

# temperature scaling of carrying capacity
K_scaling <- function(K0, E_B, E_S, T){
  K0 * euler^(-1*(E_S-E_B)/(k*T))
}

# initial mortality rate, m
m0 <- function(m_base, # baseline m for a given temperature
               E_m,    # activation energy for mortality in units of eV
               T){
  # m_base = m0 * arrenius_eqn to get observed m value, so:
  m_base / arrhenius_eqn(E_a = E_m, T = T)
}

# temperature scaling of mortality rate
m_scaling <- function(m0, E_m, T){
  m0 * arrhenius_eqn(E_a = E_m, T = T)
}

# conversion efficiency (e) is considered to be temperature dependent,
# which is attributed to data from Peters 1983

# initial attack rate, a
a0 <- function(a_base, # baseline a for a given temperature   
               v0_C,   # initial body velocity of consumer
               v0_R,   # initial body velocity of resource
               E_vC,   # activation energy of consumer body velocity
               E_vR,   # activation energy of resource body velocity
               T_C,    # consumer body temperature in Kelvin
               T_R){   # resource body temperature in Kelvin
  # a_base = a0 * sqrt(v0_C^2 * arrhenius_eqn for C + v0_R^2 * arrhenius_eqn for R)
  a_base / sqrt((v0_C^2)*arrhenius_eqn(E_a = 2*E_vC, T = T_C) + (v0_R^2)*arrhenius_eqn(E_a = 2*E_vR, T = T_R))
}

# temperature scaling of attack rate
# Osmond et al. 2017, Am Nat set v0 = 1 for both consumer and resource
# and E_v = 0.46 for both the consumer and resource
# this assumes they are both active foragers
a_scaling <- function(a0, v0_C, v0_R, E_vC, E_vR, T_C, T_R){
  a0 * sqrt((v0_C^2)*arrhenius_eqn(E_a = 2*E_vC, T = T_C) + (v0_R^2)*arrhenius_eqn(E_a = 2*E_vR, T = T_R))
}
