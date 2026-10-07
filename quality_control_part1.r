library(stringr)
library(glue)
library(magrittr)
library(dplyr)
library(data.table)
source("./Dataset_v2/Scripts/quality_control_functions.r")

work_dir <- "" #working directory
provider <- "" #provider's name
diag_head <- "" #Log-file header

data_meta <- NA #read metadata
data_origin <- "" #starting path
origin_out <- "" #out path
if (!dir.exists(origin_out)) {system(glue("mkdir {origin_out}"))}

diag_df <- data.frame(Code = character(), 
					  Date = as.Date(character()), 
					  Phy_test = character(), 
					  Rep_test = character(), 
					  Time_step = character(), 
					  Outliers = character(), 
					  Variable = character(),
					  Value = double()
					 )
test_names <- c("Physical Test", "Repeated Values", "Time Step", "Outliers")
stats <- data_meta$Code

for (i in stats) {

	print(glue("Quality Control Station {data_meta$Name[data_meta$Code == i]}"))
	
	file_name <- Sys.glob(glue("{data_origin}{i}_*.txt"))[1]
	
	data_stat <- fread(file_name)
	
	#Mean Temperature
	data_stat$T <- as.numeric(data_stat$T)
	if (sum(data_stat$T[!is.na(data_stat$T)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$T, "T")
		flag_mat[, 2] <- repeat_test(data_stat$T, "T")
	
		data_stat$T[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$T, "T")
		flag_mat[, 4] <- outlier(data_stat$T, "T")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "T", i) 
		diag_tmp$Value <- data_stat$T[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Minimum Temperature
	data_stat$TMIN <- as.numeric(data_stat$TMIN)
	if (sum(data_stat$TMIN[!is.na(data_stat$TMIN)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$TMIN, "TMIN")
		flag_mat[, 2] <- repeat_test(data_stat$TMIN, "TMIN")
	
		data_stat$TMIN[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$TMIN, "TMIN")
		flag_mat[, 4] <- outlier(data_stat$TMIN, "TMIN")
		
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "TMIN", i) 
		diag_tmp$Value <- data_stat$TMIN[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Maximum Temperature
	data_stat$TMAX <- as.numeric(data_stat$TMAX)
	if (sum(data_stat$TMAX[!is.na(data_stat$TMAX)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$TMAX, "TMAX")
		flag_mat[, 2] <- repeat_test(data_stat$TMAX, "TMAX")
	
		data_stat$TMAX[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$TMAX, "TMAX")
		flag_mat[, 4] <- outlier(data_stat$TMAX, "TMAX")
		
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "TMAX", i) 
		diag_tmp$Value <- data_stat$TMAX[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
	
		rm(flag_mat)
		
	}
	
	#Precipitation
	data_stat$P <- as.numeric(data_stat$P)
	if (sum(data_stat$P[!is.na(data_stat$P)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$P, "P")
		flag_mat[, 2] <- repeat_test(data_stat$P, "P")
	
		data_stat$P[flag_mat[, 2] == 1] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$P, "P")
		flag_mat[, 4] <- clim_precipitation(data_stat$P, data_stat$T)
		
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "P", i) 
		diag_tmp$Value <- data_stat$P[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
	
		rm(flag_mat)
		
	}
	
	#Pressure
	data_stat$SP <- as.numeric(data_stat$SP)
	if (sum(data_stat$SP[!is.na(data_stat$SP)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$SP, "SP")
		flag_mat[, 2] <- repeat_test(data_stat$SP, "SP")
	
		data_stat$SP[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 4] <- outlier(data_stat$SP, "SP")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "SP", i) 
		diag_tmp$Value <- data_stat$SP[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Relative Humidity
	data_stat$RH <- as.numeric(data_stat$RH)
	if (sum(data_stat$RH[!is.na(data_stat$RH)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$RH, "RH")
		flag_mat[, 2] <- repeat_test(data_stat$RH, "RH")
	
		data_stat$RH[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 4] <- outlier(data_stat$RH, "RH")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "RH", i) 
		diag_tmp$Value <- data_stat$RH[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Wind Speed
	data_stat$WS <- as.numeric(data_stat$WS)
	if (sum(data_stat$WS[!is.na(data_stat$WS)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$WS, "WS")
		flag_mat[, 2] <- repeat_test(data_stat$WS, "WS")
	
		data_stat$WS[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$WS, "WS")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "WS", i) 
		diag_tmp$Value <- data_stat$WS[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Wind Direction
	data_stat$WD <- as.numeric(data_stat$WD)
	if (sum(data_stat$WD[!is.na(data_stat$WD)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$WD, "WD")
		flag_mat[, 2] <- repeat_test(data_stat$WD, "WD")
	
		data_stat$WD[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$WD, "WD")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "WD", i) 
		diag_tmp$Value <- data_stat$WD[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Global Radiation
	data_stat$G <- as.numeric(data_stat$G)
	if (sum(data_stat$G[!is.na(data_stat$G)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$G, "G")
		flag_mat[, 2] <- repeat_test(data_stat$G, "G")
	
		data_stat$G[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 4] <- outlier(data_stat$G, "G")
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "G", i) 
		diag_tmp$Value <- data_stat$G[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	#Snow Depth
	data_stat$HS <- as.numeric(data_stat$HS)
	if (sum(data_stat$HS[!is.na(data_stat$HS)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$HS, "HS")
	
		data_stat$HS[flag_mat[, 1] == 1] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$HS, "HS")
		flag_mat[, 4] <- outlier(data_stat$HS, "HS", out_th = 9, out99 = 5)
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "HS", i) 
		diag_tmp$Value <- data_stat$HS[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
	
		rm(flag_mat)
		
	}
	
	#Snowfall
	data_stat$HN <- as.numeric(data_stat$HN)
	if (sum(data_stat$HN[!is.na(data_stat$HN)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$HN, "HN")
		flag_mat[, 2] <- repeat_test(data_stat$HN, "HN")
	
		data_stat$HN[(flag_mat[, 1] == 1) | (flag_mat[, 2] == 1)] <- NaN
	
		flag_mat[, 3] <- time_step(data_stat$HN, "HN")
		flag_mat[, 4] <- outlier(data_stat$HN, "HN", out_th = 9, out99 = 5)
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "HN", i) 
		diag_tmp$Value <- data_stat$HN[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
	
		rm(flag_mat)
		
	}
	
	data_stat$CC <- as.numeric(data_stat$CC)
	if (sum(data_stat$CC[!is.na(data_stat$CC)]) != 0) {
	
		flag_mat <- matrix(0, nrow = nrow(data_stat), ncol = length(test_names)) #quality control matrix
	
		flag_mat[, 1] <- phi_test(data_stat$CC, "CC")
	
		data_stat$CC[flag_mat[, 1] == 1] <- NaN
	
		diag_tmp <- diag_fill(flag_mat, data_stat$Date, test_names, "CC", i) 
		diag_tmp$Value <- data_stat$CC[data_stat$Date %in% diag_tmp$Date]
		diag_df <- rbind(diag_df, diag_tmp)
						  
		rm(flag_mat)
		
	}
	
	file_out <- str_replace_all(file_name, data_origin, origin_out)
	fwrite(data_stat, file_out)
	
}

vars <- c("T", "TMIN", "TMAX", "P", "SP", "RH", "WS", "WD", "G", "HS", "HN", "CC")
work_out <- glue("{work_dir}Logs/{provider}/")
if (!dir.exists(work_out)) {system(glue("mkdir {work_out}"))}

for (i in vars) {

	diag_out <- diag_df[diag_df$Variable == i, ] %>% dplyr::select(., !all_of("Variable"))
	if (nrow(diag_out) != 0) {diag_out$Error_N <- sprintf("%04d",1:nrow(diag_out))}
	diag_file <- glue("{work_out}/{provider}_{i}_log.txt")
	fwrite(diag_out, diag_file)
	
	if (nrow(diag_out) == 0) {next}
	if (i == "P") {
		diag_out <- diag_out %>% dplyr::select(., c("Code", "Date", "Physical Test", "Time Step", "Outliers", "Value", "Error_N"))
		diag_out <- diag_out[(as.data.frame(diag_out)[,3] == "yes") | (as.data.frame(diag_out)[,4] == "yes") | (as.data.frame(diag_out)[,5] == "yes"),]
	} else {
		diag_out <- diag_out %>% dplyr::select(., c("Code", "Date", "Time Step", "Outliers", "Value", "Error_N"))
		diag_out <- diag_out[(as.data.frame(diag_out)[,3] == "yes") | (as.data.frame(diag_out)[,4] == "yes"),]
	}
	if (nrow(diag_out) != 0) {
		diag_file <- glue("{work_dir}Logs/{provider}_{i}_log.txt")
		fwrite(diag_out, diag_file)
	}
	
}

for (i in vars) {

	if (!file.exists(glue("{work_dir}Logs/{provider}_{i}_log.txt"))) {next}
	diag_test <- fread(glue("{work_dir}Logs/{provider}_{i}_log.txt"))
	diag_in <- fread(glue("{work_out}/{provider}_{i}_log.txt"))
	
	if (i == "P") {
	
		diag_out <- diag_in[(as.data.frame(diag_in)[, 4] == "yes") | (diag_in$Error_N %in% diag_test$Error_N) ,]
		
	} else {
	
		diag_out <- diag_in[(as.data.frame(diag_in)[, 3] == "yes") | (as.data.frame(diag_in)[, 4] == "yes") | (diag_in$Error_N %in% diag_test$Error_N) ,]
		
	}
	
	fwrite(diag_out, glue("{work_out}/{provider}_{i}_log.txt"))
	
}