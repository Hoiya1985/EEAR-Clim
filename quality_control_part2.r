library(stringr)
library(glue)
library(magrittr)
library(dplyr)
library(data.table)
source("./Dataset_v2/Scripts/quality_control_functions.r")

work_dir <- "" #working directory
provider <- "" #provider's name
diag_head <- "" #logfile header

data_meta <- NA #metadata
data_origin <- "" #starting directory

diag_df <- data.frame(Code = character(), 
					  Date = as.Date(character()), 
					  T_clim = character(), 
					  Tmin_test = character(), 
					  Tmax_test = character(), 
					  Minmax_test = character(),
					  Hydro_test = character(), 
					  Dew_point_test = character(), 
					  Wind_clim = character(), 
					  Clear_sky_radiation = character(), 
					  Snow_clim = character(), 
					  Snowdrift_test = character(), 
					  Variable = character()
					 )
					 
test_names <- c("T Clim", "Tmin Clim", "Tmax Clim", "Minmax", "Dew Point", 
		        "Clear Sky Radiation", "Snow Clim", "Snowdrift")
stats <- data_meta$Code

for (i in stats) {

	print(glue("Quality Control Station {data_meta$Name[data_meta$Code == i]}"))
	
	file_name <- Sys.glob(glue("{data_origin}{i}_*.txt"))[1]
	
	data_stat <- fread(file_name)
	
	data_stat$T_QC <- as.numeric(data_stat$T)
	data_stat$TMIN_QC <- as.numeric(data_stat$TMIN)
	data_stat$TMAX_QC <- as.numeric(data_stat$TMAX)
	data_stat$P_QC <- as.numeric(data_stat$P)
	data_stat$SP_QC <- as.numeric(data_stat$SP)
	data_stat$RH_QC <- as.numeric(data_stat$RH)
	data_stat$WS_QC <- as.numeric(data_stat$WS)
	data_stat$WD_QC <- as.numeric(data_stat$WD)
	data_stat$G_QC <- as.numeric(data_stat$G)
	data_stat$HS_QC <- as.numeric(data_stat$HS)
	data_stat$HN_QC <- as.numeric(data_stat$HN)
	
	#Mean Temperature
	if (sum(is.na(data_stat$T_QC)) != nrow(data_stat)) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
	
		flag_mat[, 1] <- clim_temperature(data_stat$T_QC, data_stat$TMIN_QC, data_stat$TMAX_QC, test = "complete")
	
		data_stat$T_QC[flag_mat[, 1] == 1] <- NaN
		data_stat$TMIN_QC[flag_mat[, 1] == 1] <- NaN
		data_stat$TMAX_QC[flag_mat[, 1] == 1] <- NaN
	
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "T", i) %>% rbind(diag_df, .)
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMIN", i) %>% rbind(diag_df, .)
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMAX", i) %>% rbind(diag_df, .)
						  
		rm(flag_mat)
		
	}
	
	#Minimum and Maximum Temperature
	if ((sum(is.na(data_stat$TMIN_QC)) != nrow(data_stat)) | (sum(is.na(data_stat$TMIN_QC)) != nrow(data_stat))) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
	
		flag_mat[, 2] <- clim_temperature(data_min = data_stat$TMIN_QC, data_max = data_stat$TMAX_QC, test = "partial_min")
		flag_mat[, 3] <- clim_temperature(data_min = data_stat$TMIN_QC, data_max = data_stat$TMAX_QC, test = "partial_max")
		flag_mat[, 4] <- clim_temperature(data_min = data_stat$TMIN_QC, data_max = data_stat$TMAX_QC, test = "minmax")
	
		data_stat$TMIN_QC[(flag_mat[, 2] == 1) | (flag_mat[, 3] == 1) | (flag_mat[, 4] == 1)] <- NaN
		data_stat$TMAX_QC[(flag_mat[, 2] == 1) | (flag_mat[, 3] == 1) | (flag_mat[, 4] == 1)] <- NaN
	
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMIN", i) %>% rbind(diag_df, .)
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMAX", i) %>% rbind(diag_df, .)
						  
		rm(flag_mat)
		
	}
	
	#Relative Humidity
	if (sum(is.na(data_stat$RH_QC)) != nrow(data_stat)) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
	
		flag_mat[, 5] <- clim_rh(data_stat$RH_QC, data_stat$T_QC)
	
		data_stat$RH_QC[flag_mat[, 5] == 1] <- NaN
	
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "RH", i) %>% rbind(diag_df, .)
						  
		rm(flag_mat)
		
	}
	
	#Global Radiation
	if (sum(is.na(data_stat$G_QC)) != nrow(data_stat)) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
	
		flag_mat[, 6] <- clear_sky_rad(data_stat$G_QC, data_stat$SP_QC, data_stat$Date, lat = data_meta$Latitude[data_meta$Code == i])
	
		data_stat$G_QC[flag_mat[, 6] == 1] <- NaN
	
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "G", i) %>% rbind(diag_df, .)
						  
		rm(flag_mat)
		
	}
	
	#Snow Depth and Snowfall
	if ((sum(is.na(data_stat$HS_QC)) != nrow(data_stat)) | (sum(is.na(data_stat$HN_QC)) != nrow(data_stat))) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
	
		flag_mat[, 7] <- clim_snow(data_stat$HS_QC, data_stat$HN_QC)
		flag_mat[, 8] <- snowdrift(data_stat$HS_QC, data_stat$T_QC, data_stat$P_QC, data_stat$WS_QC, data_stat$HN_QC)
	
		data_stat$HS_QC[(flag_mat[, 7] == 1) | (flag_mat[, 8] == 1)] <- NaN
		data_stat$HN_QC[flag_mat[, 7] == 1] <- NaN
		
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "HS", i) %>% rbind(diag_df, .)
		diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "HN", i) %>% rbind(diag_df, .)
	
		rm(flag_mat)
		
	}
	
	fwrite(data_stat, file_name)
	
}

vars <- c("T", "TMIN", "TMAX", "P", "SP", "RH", "WS", "WD", "G", "HS", "HN")
for (i in vars) {

	diag_out <- diag_df[diag_df$Variable == i, ] %>% dplyr::select(., !all_of("Variable"))
	diag_file <- glue("{work_dir}{provider}_{i}_log_part2.txt")
	fwrite(diag_out, diag_file)
	
}
