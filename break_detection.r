library(glue)
library(magrittr)
library(lubridate)
library(data.table)
library(dplyr)
library(fields)
library(stringr)
library(geosphere)
library(tidyverse)
source("./Dataset_v2/Scripts/ref_search.r")
source("./Dataset_v2/Scripts/homogenization_functions.r")

options(dplyr.summarise.inform = FALSE)

work_dir <- "" #Working directory
main_dir <- "" #Starting directory
out_dir <- "" #Output directory
ac_dir <- "" #ACMANT Program directory

variable <- ""
var_ac <- "" #TT,RR" #variable name for Acmant, FF wind, RS radiation, RR snow depth-snowfall, HH humidity, PP pressure
name_ac <- "" #TM, TN, TX, PP, SP, RH, WS, WD, GG, HS, HN
distance <- 100000
hgt <- 300
r_cor <- NA #r=0.9 (tas,tmax,tmin),0.8 (rad,pres,rh), 0.7 (prec,hs,hn), 0.5 wind
clim_type <- NA #1 tas,tmax,tmin,rh,pres 2 prec,wind,hs,hn,rad 
snht_th <- NA #snht=35 hs,hn; 25 otherwise	
ratio <- FALSE

region <- ""
regs <- "" #Providers

meta_all <- NA #All metadata

data_meta <- NA #Read metadata of selected provider
stats <- data_meta$Code
break_df <- data.frame(Code = character(), Variable = character(), Homo = integer())
ac_batch <- c()

for (i in stats) {
	
	data_yn <- as.data.frame(data_meta[data_meta$Code == i, ])[colnames(data_meta) == variable]
	if (is.na(data_yn) | (data_yn != "x")) {next}
	data_in <- NA #Read data
	data_in <- data.frame(Date = data_in$Date, Value = as.data.frame(data_in)[which(colnames(data_in) == glue("{variable}_QC"))])
	colnames(data_in) <- c("Date", "Value")
	nas <- which(!is.na(data_in$Value))
	if (length(nas) == 0) {next}
	n <- interval(as.Date(glue("{data_in$Date[nas[1]]}")), as.Date(glue("{data_in$Date[nas[length(nas)]]}"))) %>% time_length(., "years")
		
	if (n < 30) {
		
		break_df <- data.frame(Code = i, Variable = variable, Homo = "N < 30") %>% rbind(break_df, .)
		next
			
	}
	
	print(glue("Breakpoint detection station {data_meta$Name[data_meta$Code == i]}"))
	
	meta_test <- meta_all[which(as.data.frame(meta_all)[which(colnames(meta_all) == variable)] == "x"), ]
	#Search for references
	data_ref <- ref_search(data_in$Value, data_in$Date, glue("{variable}_QC"), meta_all, data_meta[data_meta$Code == i], 
						   work_dir, out_dir, d = distance, h = hgt, r = r_cor, na_th = 70, sel_mode = "homo", n_min = 4, n_max = 25) 
	#No references case					   	   
	if (is.null(data_ref)) {
		
		break_df <- data.frame(Code = i, Variable = variable, Homo = "Insufficient number of reference stations") %>% 
						rbind(break_df, .)
		next
			
	} 
	
	ref_stat <- fread(glue("{out_dir}reference_stations.txt"))
	
	if (ratio) {
	
		data_mon <- data_in %>% 
				mutate(Month = month(Date), Year = year(Date)) %>% 
				group_by(Month, Year) %>% 
				summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * sum(Value, na.rm = TRUE)) %>% 
				arrange(., Year)
				
	} else {
	
		data_mon <- data_in %>% 
					mutate(Month = month(Date), Year = year(Date)) %>% 
					group_by(Month, Year) %>% 
					summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * mean(Value, na.rm = TRUE)) %>% 
					arrange(., Year)
					
	}
					
	y_start <- data_mon$Year[1]
	y_end <- data_mon$Year[nrow(data_mon)]
		
	#drop years with no data in every station
	ref_na <- apply(as.matrix(data_ref), 1, function(x) {sum(!is.na(x))})
	if (length(ref_na[ref_na == 0]) != 0) {
	
		nas <- which(ref_na == 0)
		y_start <- data_mon$Year[nas[length(nas)]] + 1
			
		data_stat <- data_mon$Mon_mean[(data_mon$Year >= y_start) & (data_mon$Year <= y_end)]
		data_ref <- data_ref[(data_mon$Year >= y_start) & (data_mon$Year <= y_end), ]
			
	} else {
		
		data_stat <- data_mon$Mon_mean
			
	}
		
	if (y_start > y_end) {
		
		break_df <- data.frame(Code = i, Variable = variable, Homo = "N < 30") %>% rbind(break_df, .)
		next
			
	}
		
	if ((y_end - y_start + 1) < 30) {
		
		break_df <- data.frame(Code = i, Variable = variable, Homo = "N < 30") %>% rbind(break_df, .)
		next
			
	}
			
	na_all <- apply(data_ref, 2, function(x) {sum(!is.na(x))}) %>% length(.)
	
	if (na_all == 0) {
			
		break_df <- data.frame(Code = i, Variable = variable, Homo = "No valid data in reference stations") %>% 
							rbind(break_df, .)
		next
			
	}
		
	print(glue("References n.: {nrow(ref_stat)}"))
	
	break_df <- data.frame(Code = i, Variable = variable, Homo = glue("Yes (N. Ref: {nrow(ref_stat)})")) %>% 
							rbind(break_df, .)
	
	data_mon <- data_mon[data_mon$Year %in% c(y_start:y_end), ]
	mat_ref <- rbind(t(as.matrix(data_stat)), t(data_ref))
	
	if (data_mon$Month[1] != 1) {
	
		y_start <- y_start + 1
		n_start <- which(data_mon$Month == 1)
		data_mon <- data_mon[n_start[1]:nrow(data_mon), ]
		mat_ref <- mat_ref[, n_start[1]:ncol(mat_ref)]
		
	}
	
	if (data_mon$Month[nrow(data_mon)] != 12) {
	
		y_end <- y_end - 1
		n_end <- which(data_mon$Month == 12)
		data_mon <- data_mon[1:n_end[length(n_end)], ]
		mat_ref <- mat_ref[, 1:n_end[length(n_end)]]
		
	}
	
	#----------------------------------------------------------------------------------------------
	#ACMANT
	
	if (!dir.exists(glue("{out_dir}Acmant/{variable}"))) {system(glue("mkdir {out_dir}Acmant/{variable}"))}
	ac_vec <- c()
	if (which(regs == region) <= length(letters)) {
		
		name_start <- letters[regs == region]
		
	} else {
		
		name_start <- which(regs == region) - 26
		
	}
	
	name_end <- paste(str_split(i, "")[[1]][4:7], collapse= "" )
	dataset_code <- glue("{str_to_title(name_start)}{name_end}")
	ac_opt <- acmant_options(variable, dataset_code)
	
	meta_sel <- rbind(meta_all[meta_all$Code == i], ref_stat)
			
	if (length(y_start:y_end) > 200) {
			
		y_start_ac <- (y_end + 1) - 200
		m_start <- ((12 * length(y_start_ac:y_end)) - (200 * 12)) + 1
		mat_ref_ac <- mat_ref[, m_start:ncol(mat_ref)]
			
	} else {
		
		y_start_ac <- y_start 
		mat_ref_ac <- mat_ref
		
	}
			
	acmant_prepare(mat_ref_ac, meta_all, meta_sel, y_start_ac, y_end, glue("{out_dir}Acmant/{variable}/Input/"), dataset_code, -999.9)

	ac_vec <- c(ac_vec, glue("cd {ac_dir}\""))

	ac_vec <- c(ac_vec, glue("copy \"C:\\Users\\Giulio Bongiovanni\\OneDrive\\Desktop\\Dataset_v2\\Homogenization\\Acmant\\{variable}\\Input\\{ac_opt$Name}*.txt\" {ac_dir}\\Input\\\""))
			
	if (variable %in% c("P", "HS", "HN")) {
				
		ac_vec <- c(ac_vec, glue("(echo {ac_opt$Code} && echo {ac_opt$Res} && echo {dataset_code} && echo {nrow(mat_ref_ac)} && echo {y_start_ac} && echo {length(y_start_ac:y_end)} && echo {ac_opt$Snow_Start} && echo {ac_opt$Snow_End} && echo {ac_opt$Outlier} && echo {ac_opt$Output}) | ACMANT4run.BAT"))
			
	} else {
			
		if (variable %in% c("WS", "WD")) {
			
			ac_vec <- c(ac_vec, glue("(echo {ac_opt$Code} && echo {ac_opt$Res} && echo {dataset_code} && echo {nrow(mat_ref_ac)} && echo {y_start_ac} && echo {length(y_start_ac:y_end)} && echo {ac_opt$Outlier} && echo {ac_opt$Output}) | ACMANT4run.BAT"))
			
		} else {
				
			ac_vec <- c(ac_vec, glue("(echo {ac_opt$Code} && echo {ac_opt$Res} && echo {dataset_code} && echo {nrow(mat_ref_ac)} && echo {y_start_ac} && echo {length(y_start_ac:y_end)} && echo {ac_opt$Season} && echo {ac_opt$Outlier} && echo {ac_opt$Output}) | ACMANT4run.BAT"))
			
		}
			
	}
			
	ac_vec <- c(ac_vec, glue("del {ac_dir}\\Input\\*.txt\""))
	ac_vec <- c(ac_vec, glue("rename {ac_dir}\\Output\\{dataset_code}_breaks.txt\" \"{i}_breaks.txt\""))
	ac_vec <- c(ac_vec, glue("copy {ac_dir}\\Output\\{i}_breaks.txt\" \"C:\\Users\\Giulio Bongiovanni\\OneDrive\\Desktop\\Dataset_v2\\Homogenization\\Acmant\\{variable}\\"))
	ac_vec <- c(ac_vec, glue("del {ac_dir}\\Output\\{i}_breaks.txt\""))
	if (length(ac_vec) != 0) {ac_batch <- c(ac_batch, ac_vec)}
	
	#----------------------------------------------------------------------------------------------
	#RHTest
	
	if (!dir.exists(glue("{out_dir}RH_Test/{variable}"))) {system(glue("mkdir {out_dir}RH_Test/{variable}"))}
	source("./Dataset_v2/Scripts/RHtestsV4.r")
	base_out <- glue("{out_dir}RH_Test/{variable}/{i}_base.dat")
	ref_out <- glue("{out_dir}RH_Test/{variable}/{i}_ref.dat")
	rh_prepare(mat_ref[1, ], t(mat_ref[2:nrow(mat_ref), ]), base_out, ref_out, -999, y_start, y_end)
	outfile <- glue("{out_dir}RH_Test/{variable}/{i}_out")
	rh_run(base_out, ref_out, outfile, -999)
			
	#----------------------------------------------------------------------------------------------
	#Climatol
	
	if (!dir.exists(glue("{out_dir}Climatol/{variable}"))) {system(glue("mkdir {out_dir}Climatol/{variable}"))}
	library(climatol)
	clim_dir <- glue("{out_dir}Climatol/{variable}/")
	climatol_out <- climatol_prepare(mat_ref, ref_stat, data_meta[data_meta$Code == i], clim_dir, variable, y_start, y_end, i, -9999)
	warn = getOption("warn")
	options(warn = -1)
	climatol_run(meta_stat = i, var_name = variable, y_start = y_start, y_end = y_end, missing_flag = -9999, var_type = clim_type, snht = snht_th, 
				 out = 35, clim_dir = clim_dir, main_dir = main_dir) 
	options(warn = warn)
			
	#----------------------------------------------------------------------------------------------
		
	if (exists("mat_ref")) {rm(mat_ref)}
	if (file.exists(glue("{out_dir}reference_stations.txt"))) {system(glue("rm {out_dir}reference_stations.txt"))}
		
}

if (length(ac_batch) != 0) {writeLines(ac_batch, glue("{out_dir}Acmant/{variable}/Acmant_{region}.bat"))}
fwrite(break_df, glue("{out_dir}breakpoint_summary_{region}_{variable}.txt"))