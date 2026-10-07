library(stringr)
library(glue)
library(magrittr)
library(dplyr)
library(data.table)
library(geosphere)
library(heatwaveR)
library(tidyverse)
library(zoo)
source("./spatial_qc_functions.r")
source("./ref_search.r")
source("./Dataset_v2/Scripts/quality_control_functions.r")

work_dir <- "./Dataset_v2/"

provider <- "MeteoSwiss"
diag_head <- glue("{provider} Quality Spatial Control Log")

data_meta <- fread(glue("{work_dir}Metadata/{provider}_meta.txt"))
meta_all <- data.frame(Code = character(), 
					   Name = character(), 
					   Longitude = double(), 
					   Latitude = double(), 
					   Elevation = integer(), 
					   T = character(),
					   TMIN = character(), 
					   TMAX = character(), 
					   P = character(),
					   SP = character(),
					   RH = character(), 
					   WS = character(), 
					   WD = character(),
					   G = character(), 
					   HS = character(), 
					   HN = character(),
					   Date_Start = integer(), 
					   Date_End = integer(),
					   Provider = character()
					  )
meta_in <- Sys.glob(glue("{work_dir}Metadata/*_meta.txt"))

for (file in meta_in) {

	meta_tmp <- fread(file)
	meta_tmp$Provider <- str_split(file, glue("{work_dir}Metadata/"))[[1]][2] %>% str_split(., "_meta.txt") %>% .[[1]] %>% .[1]
	meta_all <- rbind(meta_all, meta_tmp)
	
}

data_origin <- glue("{work_dir}Quality_Control/{provider}/")
origin_search <- glue("{work_dir}Quality_Control/")

diag_df <- data.frame(Code = character(),
					  Date = as.Date(character()),
					  Corroboration_test = character(),
					  Anomaly_Re_test = character(),
					  Wet_test = character(),
					  Dry_test = character(),
					  Wind_test = character(), 
					  Rad_test = character(),
					  Snow_test = character(),
					  Sp_test = character(),
					  Variable = character()
					 )
					 
test_names <- c("Corroboration Test", "Anomaly Re Test", "Wet Test", "Dry Test",
				"Wind Test", "Rad Test", "Snow Test", "Sp Test")
stats <- data_meta$Code

for (i in stats[1420:length(stats)]) {

	print(glue("Quality Control Station {data_meta$Name[data_meta$Code == i]}"))
	
	file_name <- Sys.glob(glue("{data_origin}{i}_*.txt"))[1]
	
	data_stat <- fread(file_name)
	
	data_stat$T_H <- as.numeric(data_stat$T_QC)
	data_stat$TMIN_H <- as.numeric(data_stat$TMIN_QC)
	data_stat$TMAX_H <- as.numeric(data_stat$TMAX_QC)
	data_stat$P_H <- as.numeric(data_stat$P_QC)
	data_stat$SP_H <- as.numeric(data_stat$SP_QC)
	data_stat$RH_H <- as.numeric(data_stat$RH_QC)
	data_stat$WS_H <- as.numeric(data_stat$WS_QC)
	data_stat$WD_H <- as.numeric(data_stat$WD_QC)
	data_stat$G_H <- as.numeric(data_stat$G_QC)
	data_stat$HS_H <- as.numeric(data_stat$HS_QC)
	data_stat$HN_H <- as.numeric(data_stat$HN_QC)
	
	if (nrow(data_stat) < (366*3)) {
	
		fwrite(data_stat, file_name)
		next
	}
	
	#Mean Temperature
	
	if (!is.na(data_meta$T[data_meta$Code == i]) & (data_meta$T[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$T == "x", ]
	
		data_ref <- ref_search(data_stat$T_QC, data_stat$Date, "T_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "")
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 1] <- corroboration_test(data_stat$T_QC, data_ref, data_stat$Date, th = 10)
			flag_mat[, 2] <- anomaly_re_test(data_stat$T_QC, data_ref, data_stat$Date, th = 10, var_type = "")
	
			data_stat$T_H[(flag_mat[, 1] == 1) & (flag_mat[, 2] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "T", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Minimum Temperature
	
	if (!is.na(data_meta$TMIN[data_meta$Code == i]) & (data_meta$TMIN[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$TMIN == "x", ]
	
		data_ref <- ref_search(data_stat$TMIN_QC, data_stat$Date, "TMIN_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "")
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 1] <- corroboration_test(data_stat$TMIN_QC, data_ref, data_stat$Date, th = 10)
			flag_mat[, 2] <- anomaly_re_test(data_stat$TMIN_QC, data_ref, data_stat$Date, th = 10, var_type = "")
	
			data_stat$TMIN_H[(flag_mat[, 1] == 1) & (flag_mat[, 2] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMIN", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Maximum Temperature
	
	if (!is.na(data_meta$TMAX[data_meta$Code == i]) & (data_meta$TMAX[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$TMAX == "x", ]
	
		data_ref <- ref_search(data_stat$TMAX_QC, data_stat$Date, "TMAX_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "")
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 1] <- corroboration_test(data_stat$TMAX_QC, data_ref, data_stat$Date, th = 10)
			flag_mat[, 2] <- anomaly_re_test(data_stat$TMAX_QC, data_ref, data_stat$Date, th = 10, var_type = "")
	
			data_stat$TMAX_H[(flag_mat[, 1] == 1) & (flag_mat[, 2] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "TMAX", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Precipitation
	
	if (!is.na(data_meta$P[data_meta$Code == i]) & (data_meta$P[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$P == "x", ]
	
		data_ref <- ref_search(data_stat$P_QC, data_stat$Date, "P_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "", h = NA, r = 0.7, wd = TRUE)
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
		
		if (ncol(as.matrix(data_ref)) >= 3) {
		
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
		
			ref_wd <- fread("wetdry_stations.txt") %>% as.matrix(.)
			dmin <- min(ref_wd, na.rm = TRUE)
	
			flag_mat[, 3] <- wet_dry_test(data_stat$P_QC, data_ref, data_stat$Date, "wet", 0.3, dmin, ref_wd)
			flag_mat[, 4] <- wet_dry_test(data_stat$P_QC, data_ref, data_stat$Date, "dry", 0.8, dmin, ref_wd)
	
			data_stat$P_H[(flag_mat[, 3] == 1) | (flag_mat[, 4] == 1)] <- NaN
	
			flag_mat[, 1] <- corroboration_test(data_stat$P_H, data_ref, data_stat$Date, win_clim = 14, th = 50, zeros = TRUE)
			flag_mat[, 2] <- anomaly_re_test(data_stat$P_H, data_ref, data_stat$Date, win_clim = 14, th = 50, var_type = "ratio", zeros = TRUE)
	
			data_stat$P_H[(flag_mat[, 1] == 1) & (flag_mat[, 2] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "P", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Pressure
	
	if (!is.na(data_meta$SP[data_meta$Code == i]) & (data_meta$SP[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$SP == "x",]
	
		data_ref <- ref_search(data_stat$SP_QC, data_stat$Date, "SP_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "", p = TRUE)
		
		if (is.null(data_ref)) {data_ref <- array(NA, dim = c(1, 1, 1))}
	
		if (ncol(as.matrix(data_ref[, , 1])) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			
			flag_mat[, 8] <- sp_spatial_test(cbind(data_stat$SP_QC, data_stat$T_QC, data_stat$RH_QC), data_ref[, , 1:3], data_meta$Elevation[data_meta$Code == i], data_ref[1, , 4])
	
			data_stat$SP_H[(flag_mat[, 8] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "SP", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Relative Humidity
	
	if (!is.na(data_meta$RH[data_meta$Code == i]) & (data_meta$RH[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$RH == "x",]
	
		data_ref <- ref_search(data_stat$RH_QC, data_stat$Date, "RH_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "")
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 1] <- corroboration_test(data_stat$RH_QC, data_ref, data_stat$Date, th = 20)
			flag_mat[, 2] <- anomaly_re_test(data_stat$RH_QC, data_ref, data_stat$Date, th = 20, var_type = "")
	
			data_stat$RH_H[(flag_mat[, 1] == 1) & (flag_mat[, 2] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "RH", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
		
	}
	
	#Wind Speed
	
	if (!is.na(data_meta$WS[data_meta$Code == i]) & (data_meta$WS[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$WS == "x",]
	
		data_ref <- ref_search(data_stat$WS_QC, data_stat$Date, "WS_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "", d = 30000, h = 250, r = 0.5)
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 5] <- wind_spatial_test(data_stat$WS_QC, data_stat$Date, data_ref)
			data_stat$WS_H[flag_mat[, 5] == 1] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "WS", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
	}
	
	#Global Radiation
	
	if (!is.na(data_meta$G[data_meta$Code == i]) & (data_meta$G[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$G == "x",]
	
		data_ref <- ref_search(data_stat$G_QC, data_stat$Date, "G_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "")
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 6] <- rad_spatial_test(data_stat$G_QC, data_stat$Date, data_ref)
			data_stat$G_H[flag_mat[, 6] == 1] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "G", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
	}
	
	#Snow Depth and Snowfall
	
	if (!is.na(data_meta$HS[data_meta$Code == i]) & (data_meta$HS[data_meta$Code == i] == "x")) {
	
		meta_test <- meta_all[meta_all$HS == "x",]
	
		data_ref <- ref_search(data_stat$HS_QC, data_stat$Date, "HS_QC", meta_test, data_meta[data_meta$Code == i], origin_search, "", d = 30000, h = 250, r = 0.7)
		
		if (is.null(data_ref)) {data_ref <- matrix(NA, nrow = 1, ncol = 1)}
	
		if (ncol(as.matrix(data_ref)) >= 3) {
	
			flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names))
			flag_mat[, 7] <- snow_test(data_stat$HS_QC, data_stat$T_H, data_stat$Date, data_ref, win = 14)
			data_stat$HS_H[(flag_mat[, 7] == 1) | (flag_mat[, 8] == 1)] <- NaN
			data_stat$HN_H[(flag_mat[, 7] == 1) | (flag_mat[, 8] == 1)] <- NaN
	
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "HS", code = i) %>% rbind(diag_df, .)
			diag_df <- diag_fill(flag_mat, data_stat$Date, test_names, "HN", code = i) %>% rbind(diag_df, .)
	
			rm(flag_mat)
		
		}
	}
	
	fwrite(data_stat, file_name)
	
}

vars <- c("T", "TMIN", "TMAX", "P", "SP", "RH", "WS", "WD", "G", "HS", "HN")
for (i in vars) {

	diag_out <- diag_df[diag_df$Variable == i, ] %>% dplyr::select(., !all_of("Variable"))
	if (i %in% c("T", "TMIN", "TMAX", "P", "RH")) {
	
		fg <- which(((as.data.frame(diag_out)[3] == "no") & (as.data.frame(diag_out)[4] == "yes")) | ((as.data.frame(diag_out)[3] == "yes") & (as.data.frame(diag_out)[4] == "no")))
		diag_out <- diag_out[-fg, ]
	
	}
	diag_file <- glue("{work_dir}Logs/{provider}/{provider}_{i}_log_space.txt")
	fwrite(diag_out, diag_file)
	
}