####################################################
#### Create Figure 3 from Falla & Cunniffe 2025 ####
####################################################

## This script contains all code needed to re-create Figure 3 from Falla & Cunniffe 2025,
## including saving it to a pdf

rm(list = ls())

# packages
library(deSolve)
library(ggplot2)
library(latex2exp)
library(dplyr)
library(gridExtra)

#### FUNCTIONS ####
companion_plant_ode <- function(times, y, parms) {
  
  ### ODE function describing the model state equations in Falla & Cunniffe 2025
  
  omega <- parms[["omega"]]
  eta <- parms[["eta"]]
  H <- parms[["H"]]
  R <- parms [["R"]]
  b <- parms[["b"]]
  a <- parms[["a"]]
  gamma <- parms[["gamma"]]
  rho <- parms[["rho"]]
  
  # aphid migration
  q <- parms[["q"]]
  pi <- parms[["pi"]]
  lambda <- parms[["lambda"]]
  
  # plant attractiveness
  v_i <- parms[["v_i"]]
  v_r <- parms[["v_r"]]
  v_b <- parms[["v_b"]]
  
  # plant acceptability
  e_i <- parms[["e_i"]]
  e_r <- parms[["e_r"]]
  e_b <- parms[["e_b"]]
  
  # states
  I <- y[["I"]]
  Z <- y[["Z"]]
  X <- y[["X"]]
  
  S <- H - I - R
  
  # R plant-induced susceptible plant attractiveness (v_rs) + acceptability (e_rs)
  v_rs <- 1 - v_b*(R/H)
  e_rs <- 1 - e_b*(R/H)
  
  # number of plants weighted by plant attractiveness
  H_hat <- v_rs*(S + v_i*I) + v_r*R
  
  # aphid dispersal rate
  phi <- H_hat / (eta*(q*H_hat + (1 - q)*omega*(e_rs*v_rs*(v_i*e_i*I + S) + v_r*e_r*R))) 
  
  # aphid infectivity loss rate
  inf_loss_s <- rho + (1 - rho)*e_rs*omega
  inf_loss_i <- rho*(1 - a*(1 - e_rs*e_i*omega)) + (1 - rho)*e_rs*e_i*omega
  inf_loss_r <- rho + (1 - rho)*e_r*omega
  tau <- (1 - q)*phi*(inf_loss_s*v_rs*S + inf_loss_i*v_i*v_rs*I + inf_loss_r*v_r*R) / H_hat
  
  # state equations
  dI <- (1 - q)*phi*Z*b*v_rs*S / H_hat - gamma*I
  dX <- (1 - pi)*lambda*(H_hat/H) + tau*Z - (1 - q)*phi*a*X*(1 - e_i*e_rs*omega)*v_i*v_rs*I / H_hat - q*phi*X
  dZ <- pi*lambda*(H_hat/H) + (1 - q)*phi*a*X*(1 - e_i*e_rs*omega)*v_i*v_rs*I / H_hat - tau*Z - q*phi*Z
  
  return(list(c(dI, dX, dZ)))
}

calculate_R0 <- function(parms) {
  
  ### Function that for a given named vector of parameters, returns the value R0
  
  omega <- parms[["omega"]]
  eta <- parms[["eta"]]
  H <- parms[["H"]]
  R <- parms [["R"]]
  b <- parms[["b"]]
  a <- parms[["a"]]
  gamma <- parms[["gamma"]]
  rho <- parms[["rho"]] # assume always lose infectivity after probing
  
  # aphid migration
  q <- parms[["q"]]
  pi <- parms[["pi"]] # not needed, must be 0 for R0 to be possible (otherwise will always be some infection)
  lambda <- parms[["lambda"]] 
  
  if (pi != 0) {
    warning("R0 won't be accurate unless pi=0")
  }
  
  # plant attractiveness
  v_i <- parms[["v_i"]]
  v_r <- parms[["v_r"]]
  v_b <- parms[["v_b"]]
  
  # plant acceptability
  e_i <- parms[["e_i"]]
  e_r <- parms[["e_r"]]
  e_b <- parms[["e_b"]]
  
  v_rs <- 1 - v_b*(R/H)
  e_rs <- 1 - e_b*(R/H)
  lose_inf_s <- rho + (1 - rho)*e_rs*omega
  lose_inf_r <- rho + (1 - rho)*e_r*omega
  
  numerator <- (1 - q)^2 * v_rs^2 * v_i * (H - R) * a * b * (1 - e_i*e_rs*omega) * lambda
  denominator <- gamma * q * H * ((1 - q)*(lose_inf_s*v_rs*(H - R) + lose_inf_r*v_r*R) + q*(v_rs*(H - R) + v_r*R))
  
  return(numerator / denominator)
}

calculate_eqX <- function(parms) {
  
  ### Function that for a given named vector of parameters, returns the 
  ### disease-free (I=0, Z=0, pi=0) equilibrium for the number of aphids (X)
  
  v_rs <- 1 - parms[["v_b"]]*parms[["R"]]/parms[["H"]]
  e_rs <- 1 - parms[["e_b"]]*parms[["R"]]/parms[["H"]]
  
  S <- parms[["H"]] - parms[["R"]] # no I
  H_hat <- v_rs*S + parms[["v_r"]]*parms[["R"]]
  
  return(parms[["lambda"]] * parms[["eta"]] * 
           (parms[["q"]]*H_hat + (1 - parms[["q"]])*parms[["omega"]]*
              (v_rs*e_rs*S + parms[["v_r"]]*parms[["e_r"]]*parms[["R"]])) /
           (parms[["q"]]*parms[["H"]])
  )
}

covary_parms_get_R0 <- function(parm1_vals, parm1_name, parm2_vals, parm2_name, parms) {
  
  ### Function that, for all combinations of given values of 2 model parameters (parm1 and parm2),
  ### calculates R0 and returns data frame with columns: parm1 values, parm2 values, R0. 
  ### "parm1_vals" and "parm2_vals" are vectors of values for the parameters with names of "parm1_name" 
  ### and "parm2_name" (given as strings)
  
  parms_df <- expand.grid(parm1_vals, parm2_vals) # dataframe with all combinations of parm1 and parm2 values
  names(parms_df) <- c(parm1_name, parm2_name)
  
  R0_vals <- apply(parms_df, 1, function(row) {
    parms[[parm1_name]] <- row[1]
    parms[[parm2_name]] <- row[2]
    calculate_R0(parms)
  })

  return(data.frame(parms_df, R0 = R0_vals))
  
}

#### PARAMETERS ####

parms_basic <- c(
  # APHID PROBING AND FEEDING
  omega = 0.2, # base probability of feeding on a plant after probing
  rho = 0.7, # probability of aphid losing infectivity from probing 
  
  # PLANT POPULATION SIZES
  R = 0.25*600, # number of virus-resistant companion plants
  H = 600, # total number of plants (S + I + R)
  
  # VIRUS ACQUISITION/INOCULATION PROBIBILITY
  b = 0.8, # probability of virus inoculation when aphid probes S plant
  a = 0.8, # probability of virus acquisition when aphid probes I plant
  
  # PLANT DEATH/REPLANTING RATE
  gamma = 0.02, # plant death rate (/day)
  
  # VECTOR PREFERENCE (ATTRACTIVENESS + ACCEPTABILITY)
  v_i = 1, # infected plant virus-induced attractiveness
  v_r = 1, # companion plant attractiveness
  v_b = 0, # degree that companion plants reduce susceptible (S,I) plant attractiveness (0-1)
  e_i = 1, # infected plant virus-induced acceptability
  e_r = 1, # companion plant acceptability
  e_b = 0, # degree that companion plants reduce susceptible (S,I) plant acceptability (0-1)
  
  # APHID FLIGHT
  eta = 0.8, # aphid feed length (days)
  
  # APHID MIGRATION
  q = 0.2, # probability aphid emigrates from field per flight
  pi = 0, # proportion of immigrating aphids that arrive infective (Z) - default 0 so R0 can be used
  lambda = 40 # number of immigrating aphids entering field (/day)
)

## REPELLENT COMPANION PARAMETERISATION
parms_repel <- parms_basic
parms_repel[["v_b"]] <- 0.8
parms_repel[["e_b"]] <- 0.8
parms_repel[["v_r"]] <- 0.8
parms_repel[["e_r"]] <- 0.4

## TRAP COMPANION PARAMETERISATION
parms_trap <- parms_basic
parms_trap[["v_b"]] <- 0
parms_trap[["e_b"]] <- 0
parms_trap[["v_r"]] <- 3 
parms_trap[["e_r"]] <- 3

#### INITIAL STATES ####

eqX_basic <- calculate_eqX(parms_basic) # X(0) is at disease-free equilibrium
init_states_basic <- c(I = 0, X = eqX_basic, Z = 1)

eqX_trap <- calculate_eqX(parms_trap) 
init_states_trap <- c(I = 0, X = eqX_trap, Z = 1)

eqX_repel <- calculate_eqX(parms_repel) 
init_states_repel <- c(I = 0, X = eqX_repel, Z = 1)


#### CREATE FIG 3a-c: Epidemic trajectories under different R plant parameterisations ####

tmax <- 300 # days

# create parameterisation for when R=0
parms_noR <- parms_basic
parms_noR[["R"]] <- 0
eqX_noR <- calculate_eqX(parms_noR)

# run trajectories
trajec_basic <- data.frame(deSolve::ode(y = init_states_basic,
                                       times = 0:tmax, 
                                       func = companion_plant_ode,
                                       parms = parms_basic))
trajec_trap <- data.frame(deSolve::ode(y = init_states_trap,
                                        times = 0:tmax, 
                                        func = companion_plant_ode,
                                        parms = parms_trap))
trajec_repel <- data.frame(deSolve::ode(y = init_states_repel,
                                       times = 0:tmax, 
                                       func = companion_plant_ode,
                                       parms = parms_repel))
trajec_no_R <- data.frame(deSolve::ode(y = c(I = 0, X = eqX_noR, Z = 1),
                                       times = 0:tmax, 
                                       func = companion_plant_ode,
                                       parms = parms_noR))

all_trajecs_df <- data.frame(rbind(trajec_no_R, trajec_basic, trajec_trap, trajec_repel),
                             parms = rep(c("No R plants", "Basic R plants", "Trap R plants", "Repellent R plants"),
                                         each = nrow(trajec_no_R)))
all_trajecs_df$propI <- all_trajecs_df$I / parms_basic[["H"]] # get proportion of infected plants (I/H)

## PLOTTING
general_format_a_to_c <- list( # list of formatting for plots a-c
  labs(x = "Time (days)", col = NULL),
  theme_bw(),
  theme(legend.position = "none",
        text = element_text(size = 15)),
  scale_color_manual(values = c("No R plants" = "#7CAE00",
                                "Basic R plants" = "#C77CFF",
                                "Trap R plants" = "#00BFC4",
                                "Repellent R plants" = "#F8766D"),
                     breaks = c("No R plants",
                                "Basic R plants",
                                "Trap R plants",
                                "Repellent R plants")))
# a) I trajectories
all_trajecs_Iplot_w_leg <- ggplot(data = all_trajecs_df, aes(x = time, y = propI, col = parms)) +
  geom_line() +
  general_format_a_to_c +
  labs(y = TeX("Proportion of infected plants $(I/H)$"),
       title = "(a)") +
  theme(legend.position = "bottom") # show legend for extraction

trajec_leg <- cowplot::get_legend(all_trajecs_Iplot_w_leg) # extract legend
all_trajecs_Iplot <- all_trajecs_Iplot_w_leg + theme(legend.position = "none") # remove legend

# b) X trajectories
all_trajecs_Xplot <- all_trajecs_Iplot + aes(y = X) +
  labs(y = TeX("Number of un-infective aphids $(X)$"),
       title = "(b)")

# b) Z trajectories
all_trajecs_Zplot <- all_trajecs_Iplot + aes(y = Z) +
  labs(y = TeX("Number of infective aphids $(Z)$"),
       title = "(c)")

#### CREATE FIG 3d-i: Proportion R plants versus R0 for varying R plant attractiveness (v_r), for TRAP vs REPELLENT ####

R_vals <- seq(0, 1, length.out = 50) * parms_basic[["H"]]
omega_vals <- c(0.2, 0.6, 0.9)
vr_vals_trap <- c(1, 1.2, 1.5, 2, 3)
vr_vals_repel <- c(0.01, 0.1, 0.5, 0.8)

all_vr_R0_df <- data.frame() # initialise empty data frame

for (omega in omega_vals) { # loop over omega values, calculating R0 for each v_r value (for trap and repel)
  parms_trap[["omega"]] <- omega
  R_vr_TRAP <- covary_parms_get_R0(parm1_vals = vr_vals_trap, 
                                                parm1_name = "v_r", 
                                                parm2_vals = R_vals, 
                                                parm2_name = "R", 
                                                parms = parms_trap)
  parms_repel[["omega"]] <- omega
  R_vr_REPEL <- covary_parms_get_R0(parm1_vals = vr_vals_repel, 
                                            parm1_name = "v_r", 
                                            parm2_vals = R_vals, 
                                            parm2_name = "R", 
                                            parms = parms_repel)
  
  # append to data frame
  all_vr_R0_df <- rbind(all_vr_R0_df,
                        data.frame(rbind(R_vr_TRAP, R_vr_REPEL),
                                   omega_val = paste("omega =", omega),
                                   R_type = rep(c("trap", "repel"), 
                                                times = c(nrow(R_vr_TRAP), nrow(R_vr_REPEL)))
                                   ))
  
}

## PLOTTING (plot individually for ease of combination with panels a-c)

# general format for plots d-f
general_format_d_to_f <- list( 
  ylim(0, max(all_vr_R0_df %>% filter(R_type == "repel") %>% pull(R0))),
    geom_hline(yintercept = 1, lty = 2, col = "grey46"),
    labs(x =  TeX(r"(Proportion of companion plants ($R/H$))"),
         y = expression(R[0]),
         col = TeX("R plant attractiveness $(v_r)$")),
    theme_bw(),
    scale_colour_viridis_d(),
    theme(legend.position = "none",
          text = element_text(size = 15)))

# d) Repellent plants, omega = 0.2 (low)
R0_low_om_REPEL_plot_w_leg <- ggplot(all_vr_R0_df %>% filter(R_type == "repel" & omega_val == "omega = 0.2"), 
                                 aes(x = R/parms_repel[["H"]], y = R0, col = as.factor(v_r))) +
  geom_line() +
  general_format_d_to_f +
  labs(title = TeX(r"((d) $\omega = 0.2$)")) +
  theme(legend.position = "bottom")

R0_repel_leg <- cowplot::get_legend(R0_low_om_REPEL_plot_w_leg)
R0_low_om_REPEL_plot <- R0_low_om_REPEL_plot_w_leg + theme(legend.position = "none")

# e) Repellent plants, omega = 0.6 (mid)
R0_mid_om_REPEL_plot <- R0_low_om_REPEL_plot %+% (all_vr_R0_df %>% filter(R_type == "repel" & omega_val == "omega = 0.6")) +
  labs(title = TeX(r"((e) $\omega = 0.6$)"))

# f) Repellent plants, omega = 0.9 (high)
R0_high_om_REPEL_plot <- R0_low_om_REPEL_plot %+% (all_vr_R0_df %>% filter(R_type == "repel" & omega_val == "omega = 0.9")) +
  labs(title = TeX(r"((f) $\omega = 0.9$)"))

# general format for plots g-i
general_format_g_to_i <- list( 
  ylim(0, max(all_vr_R0_df %>% filter(R_type == "trap") %>% pull(R0))),
  geom_hline(yintercept = 1, lty = 2, col = "grey46"),
  labs(x =  TeX(r"(Proportion of companion plants ($R/H$))"),
       y = expression(R[0]),
       col = TeX("R plant attractiveness $(v_r)$")),
  theme_bw(),
  scale_colour_brewer(palette = "YlOrRd"),
  theme(legend.position = "none",
        text = element_text(size = 15)))

# g) Trap plants, omega = 0.2 (low)
R0_low_om_TRAP_plot_w_leg <- ggplot(all_vr_R0_df %>% filter(R_type == "trap" & omega_val == "omega = 0.2"), 
                                     aes(x = R/parms_repel[["H"]], y = R0, col = as.factor(v_r))) +
  geom_line() +
  general_format_g_to_i +
  labs(title = TeX(r"((g) $\omega = 0.2$)")) +
  theme(legend.position = "bottom")

R0_trap_leg <- cowplot::get_legend(R0_low_om_TRAP_plot_w_leg) # extract legend
R0_low_om_TRAP_plot <- R0_low_om_TRAP_plot_w_leg + theme(legend.position = "none")

# h) Trap plants, omega = 0.6 (mid)
R0_mid_om_TRAP_plot <- R0_low_om_TRAP_plot %+% (all_vr_R0_df %>% filter(R_type == "trap" & omega_val == "omega = 0.6")) +
  labs(title = TeX(r"((h) $\omega = 0.6$)"))

# i) Trap plants, omega = 0.9 (high)
R0_high_om_TRAP_plot <- R0_low_om_TRAP_plot %+% (all_vr_R0_df %>% filter(R_type == "trap" & omega_val == "omega = 0.9")) +
  labs(title = TeX(r"((i) $\omega = 0.9$)"))

#### FORMAT ALL FIGURE PANELS TOGETHER ####

# layout matrix for grid.arrange
layout <- matrix(c(1,  2, 3,
                NA, 4, NA,
                5,  6, 7,
                NA, 8, NA,
                9, 10, 11,
                NA, 12, NA),
              ncol = 3, byrow = T)


pdf("Figure_3.pdf", height = 13, width = 12) # save final figure to pdf in current directory
gridExtra::grid.arrange(all_trajecs_Iplot, all_trajecs_Xplot, all_trajecs_Zplot,
                        trajec_leg,
                        R0_low_om_REPEL_plot, R0_mid_om_REPEL_plot, R0_high_om_REPEL_plot,
                        R0_repel_leg,
                        R0_low_om_TRAP_plot, R0_mid_om_TRAP_plot, R0_high_om_TRAP_plot,
                        R0_trap_leg,
                        heights = c(1, 0.2, 1, 0.2, 1, 0.2),
                        layout_matrix = layout)
dev.off()

