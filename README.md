# sign_size
 Contains code to replicate the simulations in Section 5 and the empirical application in section 6. The only required change for the code to run is to change the "path" string at the beginning for whatever the correct path is in your machine.

In particular: 

/simulation: 
 - "simulation_goncalves_dgp.m" runs and saves the simulations
 - "simulation_goncalves_dgp_create_table.m" loads the simulation results and reports them in a Latex code that we then use for the paper.

/main_oil_irfs: 
-"main_oil_irfs.m" loads the data and runs the local projections presented in the paper. Several estimation options are commented in the code.

/_auxiliary_functions: contains auxiliary functions used by the main scripts.

/_data : contains all data used in the empirical application. In particular, it contains the data from the replication files in Kaenzig (2021,AER), as well as a different database not used in the paper.

