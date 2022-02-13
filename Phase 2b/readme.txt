R codes for phase 2b/3 power and decision boundary estimation using 2 surrogate endpoints. SAS codes for SEs of slope

1. phase2b_3_joint_endpoints_power_calculation.R
	Program that calculates table of upper prediction bounds for hazard ratio given 2 assumed surrogate endpoints. It then approximates the boundary and estimate the power.

2. Linear_and_quadratic_decision_boundary_plotting.R
	Identifies the boundary, plot it together with the linear or quadratic fitted curve.

3. H0_H1_power_PoS_plots.R
	Produces 2 plots. First one gives distribution under H0 and H1, with alpha and beta highlighted. Second one gives distribution for power calculation vs distribution for PoS calculation, critical value highlighted.

4. SAS code for finding SEs
	Contains 2 SAS code, the one called 'joint_endpoint_SE_charles.sas' is the important one that calculates SEs of eGFR slope for various settings in phase 2b/3
