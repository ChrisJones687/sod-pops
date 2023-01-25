## Setup multirun
# remotes::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)
library(terra)
library(folderfun)
library(doParallel)
# library(plyr)

setff("In", "C:/Users/cmjone25/Desktop/SOD_OR/")


infected_file <- ffIn("Infections/inf_2015_eu1.tif")
host_file <- ffIn("Hosts/hosts_2016.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temp <- TRUE
temperature_coefficient_file <- ffIn("Weather/weather_coef_2016.tif")
means <- read.table("parameters/eu1_2020_posteriors.csv", header = F)
parameter_means <- t(means)
parameter_means <- parameter_means[1,]
parameter_means[1] <- parameter_means[1]
parameter_cov_matrix <- read.table("parameters/eu1_2020_cov_posteriors.csv", header = F)
precip <- FALSE
precipitation_coefficient_file <- ""
model_type <- "SI"
latency_period <- 0
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2020-01-01'
end_date <- '2020-12-31'
use_survival_rates <- FALSE
survival_rate_month <- 3
survival_rate_day <- 15
survival_rates_file <- ""
use_lethal_temperature <- FALSE
temperature_file <- ''
lethal_temperature <- -15
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
mortality_frequency = "year"
mortality_frequency_n = 1
management <- FALSE
treatment_dates <- c('2020-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"
number_of_iterations <- 200
number_of_cores <- 10
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
random_seed <- NULL
output_frequency <- "year"
output_frequency_n <- 1
movements_file <- ""
use_movements <- FALSE
start_exposed <- FALSE
generate_stochasticity <- TRUE
establishment_stochasticity <- TRUE
movement_stochasticity <- TRUE
dispersal_stochasticity <- TRUE
establishment_probability <- 0.5
dispersal_percentage <- 0.99
quarantine_areas_file <- ""
use_quarantine <- FALSE
use_spreadrates <- FALSE
use_overpopulation_movements <- FALSE
overpopulation_percentage <- 0.75
leaving_percentage <- .50
leaving_scale_coefficient <- 5
exposed_file <- ""
mask <- NULL
write_outputs <- "None"
output_folder_path <- ""
network_filename <- ""
network_movement <- "walk"
use_initial_condition_uncertainty <- FALSE
use_host_uncertainty <- FALSE

eu1_2016_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)

sum(values(eu1_2016_mr$simulation_mean))
sum(values(eu1_2016_mr$simulation_sd))

sim_inf_eu1_2016 <- c(eu1_2016_mr$simulation_mean, eu1_2016_mr$simulation_sd)
writeRaster(sim_inf_eu1_2016, ffIn("sim_infections/eu1_2016_end_nt.tif"), overwrite = TRUE)

infected_file <- ffIn("sim_infections/eu1_2016_end_nt.tif")
host_file <- ffIn("Hosts/hosts_2017.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temperature_coefficient_file <- ffIn("Weather/weather_coef_2017.tif")
use_initial_condition_uncertainty <- TRUE


eu1_2017_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)


sum(values(eu1_2017_mr$simulation_mean))
sum(values(eu1_2017_mr$simulation_sd))

sim_inf_eu1_2017 <- c(eu1_2017_mr$simulation_mean, eu1_2017_mr$simulation_sd)
writeRaster(sim_inf_eu1_2017, ffIn("sim_infections/eu1_2017_end_nt.tif"), overwrite = TRUE)

infected_file <- ffIn("sim_infections/eu1_2017_end_nt.tif")
host_file <- ffIn("Hosts/hosts_2018.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temperature_coefficient_file <- ffIn("Weather/weather_coef_2018.tif")
use_initial_condition_uncertainty <- TRUE
# 
# infected_file <- ffIn("End of Year Infections/end_inf_2017_eu1.tif")
# host_file <- ffIn("Hosts/hosts_2018.tif")
# total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
# temp <- TRUE
# temperature_coefficient_file <- ffIn("Weather/weather_coef_2018.tif")
# means <- read.table("parameters/eu1_2020_posteriors.csv", header = F)
# parameter_means <- t(means)
# parameter_means <- parameter_means[1,]
# parameter_means[1] <- parameter_means[1]
# parameter_cov_matrix <- read.table("parameters/eu1_2020_cov_posteriors.csv", header = F)
# precip <- FALSE
# precipitation_coefficient_file <- ""
# model_type <- "SI"
# latency_period <- 0
# time_step <- "week"
# season_month_start <- 1
# season_month_end <- 12
# start_date <- '2020-01-01'
# end_date <- '2020-12-31'
# use_survival_rates <- FALSE
# survival_rate_month <- 3
# survival_rate_day <- 15
# survival_rates_file <- ""
# use_lethal_temperature <- FALSE
# temperature_file <- ''
# lethal_temperature <- -15
# lethal_temperature_month <- 1
# mortality_on <- FALSE
# mortality_rate <- 0
# mortality_time_lag <- 0
# mortality_frequency = "year"
# mortality_frequency_n = 1
# management <- FALSE
# treatment_dates <- c('2020-12-24')
# treatments_file <- ""
# treatment_method <- "ratio"
# natural_kernel_type <- "cauchy"
# anthropogenic_kernel_type <- "cauchy"
# natural_dir <- "NONE"
# anthropogenic_dir <- "NONE"
# number_of_iterations <- 200
# number_of_cores <- 10
# pesticide_duration <- c(0)
# pesticide_efficacy <- 1.0
# random_seed <- NULL
# output_frequency <- "year"
# output_frequency_n <- 1
# movements_file <- ""
# use_movements <- FALSE
# start_exposed <- FALSE
# generate_stochasticity <- TRUE
# establishment_stochasticity <- TRUE
# movement_stochasticity <- TRUE
# dispersal_stochasticity <- TRUE
# establishment_probability <- 0.5
# dispersal_percentage <- 0.99
# quarantine_areas_file <- ""
# use_quarantine <- FALSE
# use_spreadrates <- FALSE
# use_overpopulation_movements <- FALSE
# overpopulation_percentage <- 0.75
# leaving_percentage <- .50
# leaving_scale_coefficient <- 5
# exposed_file <- ""
# mask <- NULL
# write_outputs <- "None"
# output_folder_path <- ""
# network_filename <- ""
# network_movement <- "walk"
# use_initial_condition_uncertainty <- FALSE
# use_host_uncertainty <- FALSE

eu1_2018_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)


sum(values(eu1_2018_mr$simulation_mean))

sim_inf_eu1_2018 <- c(eu1_2018_mr$simulation_mean, eu1_2018_mr$simulation_sd)
writeRaster(sim_inf_eu1_2018, ffIn("sim_infections/eu1_2018_end_nt.tif"), overwrite = TRUE)

infected_file <- ffIn("sim_infections/eu1_2018_end_nt.tif")
host_file <- ffIn("Hosts/hosts_2019.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temperature_coefficient_file <- ffIn("Weather/weather_coef_2019.tif")
use_initial_condition_uncertainty <- TRUE


eu1_2019_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)


eu1_sim_2019 <- eu1_2019_mr$simulation_mean
eu1_2019 <- rast(ffIn("End of Year Infections/end_inf_2019_eu1.tif"))
eu1_sim_2019[eu1_2019 > 0]
eu1_2018 <- rast(ffIn("End of Year Infections/end_inf_2018_eu1.tif"))
sum(values(eu1_2018))
sum(values(eu1_sim_2019))
sum(values(eu1_2019_mr$simulation_mean))
sum(values(eu1_2019))

sim_inf_eu1_2019 <- c(eu1_2019_mr$simulation_mean, eu1_2019_mr$simulation_sd * 2)
writeRaster(sim_inf_eu1_2019, ffIn("sim_infections/eu1_2019_end_nt.tif"), overwrite = TRUE)

infected_file <- ffIn("sim_infections/eu1_2019_end_nt.tif")
host_file <- ffIn("Hosts/hosts_2020.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temperature_coefficient_file <- ffIn("Weather/weather_coef_2020.tif")

eu1_2020_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)


eu1_sim_2020 <- eu1_2020_mr$simulation_mean
eu1_2020 <- rast(ffIn("End of Year Infections/end_inf_2020_eu1.tif"))
eu1_sim_2020[eu1_2020 > 0]
eu1_2019 <- rast(ffIn("End of Year Infections/end_inf_2019_eu1.tif"))
sum(values(eu1_2019))
sum(values(eu1_sim_2020))
sum(values(eu1_2020))

sim_inf_eu1_2020 <- c(eu1_2020_mr$simulation_mean, eu1_2020_mr$simulation_sd *2)
writeRaster(sim_inf_eu1_2020, ffIn("sim_infections/eu1_2020_end_nt.tif"))


infected_file <- ffIn("sim_infections/eu1_2020_end_nt.tif")
host_file <- ffIn("Hosts/hosts_2021.tif")
total_populations_file <- ffIn("Hosts/lemma_max100m.tif")
temperature_coefficient_file <- ffIn("Weather/weather_coef_2021.tif")

eu1_2021_mr <- pops_multirun(infected_file,
                             host_file,
                             total_populations_file,
                             parameter_means,
                             parameter_cov_matrix,
                             temp,
                             temperature_coefficient_file,
                             precip,
                             precipitation_coefficient_file,
                             model_type,
                             latency_period,
                             time_step,
                             season_month_start,
                             season_month_end,
                             start_date,
                             end_date,
                             use_survival_rates,
                             survival_rate_month,
                             survival_rate_day,
                             survival_rates_file,
                             use_lethal_temperature,
                             temperature_file,
                             lethal_temperature,
                             lethal_temperature_month,
                             mortality_on,
                             mortality_rate,
                             mortality_time_lag,
                             mortality_frequency,
                             mortality_frequency_n,
                             management,
                             treatment_dates,
                             treatments_file,
                             treatment_method,
                             natural_kernel_type,
                             anthropogenic_kernel_type,
                             natural_dir,
                             anthropogenic_dir,
                             number_of_iterations,
                             number_of_cores,
                             pesticide_duration,
                             pesticide_efficacy,
                             random_seed,
                             output_frequency,
                             output_frequency_n,
                             movements_file,
                             use_movements,
                             start_exposed,
                             generate_stochasticity,
                             establishment_stochasticity,
                             movement_stochasticity,
                             dispersal_stochasticity,
                             establishment_probability,
                             dispersal_percentage,
                             quarantine_areas_file,
                             use_quarantine,
                             use_spreadrates,
                             use_overpopulation_movements,
                             overpopulation_percentage,
                             leaving_percentage,
                             leaving_scale_coefficient,
                             exposed_file,
                             mask,
                             write_outputs,
                             output_folder_path,
                             network_filename,
                             network_movement,
                             use_initial_condition_uncertainty,
                             use_host_uncertainty)


eu1_sim_2021 <- eu1_2021_mr$simulation_mean
eu1_2021 <- rast(ffIn("End of Year Infections/end_inf_2021_eu1.tif"))
sum(values(eu1_2020))
sum(values(eu1_sim_2020))
sum(values(eu1_2021))

sim_inf_eu1_2021 <- c(eu1_2021_mr$simulation_mean, eu1_2021_mr$simulation_sd *2)
writeRaster(sim_inf_eu1_2021, ffIn("sim_infections/eu1_2021_end_nt.tif"))
