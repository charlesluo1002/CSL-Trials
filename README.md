# CSL Clinical Trials Design and Analysis

Clinical trial design and analysis procedures covering phase 2a, phase 2b and phase 3 trials of a kidney failure therapy. 

## Phase 2a

Simulation_of_Bronze_Study.R - 
From Posterior distributions for the parameters that are generated using MCMC given previous studies data, simulate the unconditional probability success of the phase 2a Bronze study.

Shiny - 
Interactive Cloud app that allows for calculation and visualization of assurance given different input values.

## Phase 2b

phase2b_3_joint_endpoints_power_calculation.R - 
Calculates table of upper prediction bounds for hazard ratio given 2 assumed surrogate endpoints. It then approximates the boundary and estimate the power.

Linear_and_quadratic_decision_boundary_plotting.R - 
Identifies the boundary, plot it together with the linear or quadratic fitted curve.

H0_H1_power_PoS_plots.R - 
Produces 2 plots. First one gives distribution under 2 hypotheses H0 and H1, with alpha and beta highlighted. Second one gives distribution for power calculation vs distribution for PoS calculation, critical value highlighted.

## Phase 3

Phase3_PoS_SGLT2i_included.R - 
Estimates the probability of success through simulation. Outputs tables of PoS with rows and columns as the targeted endpoints for different groups in phase 2a.

plot for epsilon elicitation.R - 
Plot of the distribution of the ratio of endpoints between 2 subgroups for different epsilons, extreme tails are highlighted for better expert elicitation.

Overall_power.R - 
Calculates phase 3 power where success requires both overall and SGLT2i subgroup outcomes to be significant.
