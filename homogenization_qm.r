library(glue)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyverse)
library(fields)
library(stringr)
library(geosphere)
source("./Dataset_v2/Scripts/ref_search.r")
source("./Dataset_v2/Scripts/homogenization_functions.r")
options(dplyr.summarise.inform = FALSE)

work_dir <- "" #working directory
variable <- "" #variable to homogenize
out_dir <- "" #output directory
method <- "" #add for T,TMIN,TMAX,SP,RH ratio for P,WS,HS,HN,G

zeros <- FALSE #TRUE for precipitation, wind speed, snowfall, snow depth, global radiation, relative humidity, cloud cover
pcent <- FALSE #TRUE for relative humidity, cloud cover
y_overlap <- 3 #minimum number of overlapping years
y_homo <- 5 #minimum length (in years) of segments to be homogenized
distance <- 100000
hgt <- 300
r_cor <- NA #r=0.9 (tas,tmax,tmin),0.8 (rad,pres,rh), 0.7 (prec,hs,hn), 0.6 (cc), 0.5 wind
r_over <- NA #correlation threshold between segments:0.7 temperature, 0.6 precipitation-rh-hs-rad-pres-hn, 0.5 cc, 0.4 wind
qp <- seq(0.05, 0.95, by = 0.05) #quantile probabilities
seg_min <- 3 #minimum number of references segments
th_low <- NA #minimum accepted value #-50 temperature,0 precipitation-rad-rh-hs-wind-cc, 500 pres
th_high <- NA #maximum accepted value #50 temperature, 75 wind, 100 rh-cc, 1100 pres, 300 hs, 250 hn, 1000 prec, 138.24 rad
#p_low <- 10
#p_high <- 90
dec <- 1

breakpoints <- NA #read breakpoint statistics
meta <- NA àread breakpoints metadata

meta_all <- NA #all metadata

for (i in meta$Code) {

	out_file <- Sys.glob(glue("{out_dir}{meta$Provider[meta$Code == i]}/{i}_*.txt"))[1]
	data_out <- fread(out_file) %>% as.data.frame(.)
	move_in <- which(colnames(data_out) == glue("{variable}_QC"))
	move_out <- which(colnames(data_out) == glue("{variable}_H"))

	if (meta$Homo[meta$Code == i] != "NH") {
	
		data_out[, move_out] <- data_out[, move_in]
		fwrite(data_out, out_file)
		next
		
	}
	
	print(glue("Homogenization station {meta$Name[meta$Code == i]}"))
	
	data_in <- data.frame(Date = data_out$Date, Value = data_out[, move_in])
	
	meta_ref <- meta_all[which(as.data.frame(meta_all)[which(colnames(meta_all) == variable)] == "x"), ]
	data_ref <- ref_search(data_in$Value, data_in$Date, glue("{variable}_QC"), meta_ref, meta[meta$Code == i, ], 
						   out_dir, work_dir, d = distance, h = hgt, r = r_cor, na_th = 70, sel_mode = "homo", n_min = 4, n_max = 100) 
	ref_stat <- fread(glue("{work_dir}reference_stations.txt"))
	rm(data_ref)
	
	print(glue("References n.: {nrow(ref_stat)}"))
	
	data_break <- breakpoints[breakpoints$Code == i, ]
	data_break <- data_break[order(Date), ] #order breakpoints by date
	
	if (zeros) {data_in$Value[data_in$Value < 0.1] <- NA}
	
	init <- 1 #starting index in time series
	for (j in 1:(nrow(data_break) - 1)) {
	
		init <- c(init, which(data_in$Date == data_break$Date[j]))
			
	}
	
	for (j in nrow(data_break):1) {
			
		break_start <- which(data_in$Date == data_break$Date[j]) - 1 #index of breakpoints in time series
		if (length(break_start) == 0) {next}
		#candidate station splitting
		#segment before break
		data_sub <- data_in[init[j]:break_start, ] #sub-segment of data for candidate station
		sub_length <- interval(data_sub$Date[nrow(data_sub)], data_sub$Date[1]) %>% time_length(., "years") %>% abs(.)
		
		if (sub_length < y_homo) {
		
			print(glue("Period {j} too short: no correction applied"))
			next
		
		}
		
		#segment after the break
		#in case of only one break, the time series is splitted in only two segments
		if (j == nrow(data_break)) {
					
			data_next <- data_in[(break_start + 1):nrow(data_in), ]
				
		} else {
					
			break_end <- which(data_in$Date == data_break$Date[j + 1]) - 1
			data_next <- data_in[(break_start + 1):break_end, ]
				
		}
		
		q_seg <- qm_quantile(data_sub)
		q_base <- qm_quantile(data_next)
		vc <- c()	 #corrections matrix
		
		#Loop in references
		for (k in ref_stat$Code) {
					
			#read reference data
			ref_file <- Sys.glob(glue("{out_dir}{ref_stat$Provider[ref_stat$Code == k]}/{k}_*.txt"))[1]
			data_ref <- fread(ref_file, select = c("Date", glue("{variable}_QC")), col.names = c("Date", "Val")) %>% 
						merge(data.frame(Date = data_in$Date, Value = data_in$Value), . ,all.x = TRUE) %>% 
						select(., !Value)
		
			if (zeros) {data_ref$Val[data_ref$Val < 0.1] <- NA}
			#breakpoints file for reference
			break_ref <- breakpoints[breakpoints$Code == k, ]
			
			if (nrow(break_ref) != 0) {
						
				break_ref <- data_break[order(Date), ] #order breakpoints by date #order dates on breakpoint file
				#loop on breakpoints
						
				for (k0 in 1:nrow(break_ref)) {
					
					#starting date of segments
					if (k0 == 1) {
						
						date_ref <- as.Date(glue("{ref_stat$Date_Start[ref_stat$Code == k]}"), "%Y%m%d")
							
					} else {
								
						date_ref <- as.Date(break_ref$Date[k0 - 1])
							
					}
							
					if ((break_ref$Date[k0] - 1) < date_ref) {next} #no valid starting date
					ref_seq <- seq(date_ref, break_ref$Date[k0] - 1, by = "1 day") #date sequence of reference segment
					sub_mat <- qm_overlap(data_sub, data_next, data_ref, data_break[j, ], break_ref[k0, ], ref_seq, y_overlap, "before", r_over)
					next_mat <- qm_overlap(data_sub, data_next, data_ref, data_break[j, ], break_ref[k0, ], ref_seq, y_overlap, "after", r_over)
					if (!is.null(sub_mat) & !is.null(next_mat)) {break}
						
				}
					
			} else { #Case with no breaks
						
				ref_seq <- as.Date(data_ref$Date)
				break_no <- data.frame(Code = k, Date = NA, Value = NA, Variable = variable)
				sub_mat <- qm_overlap(data_sub, data_next, data_ref, data_break[j, ], break_no, ref_seq, y_overlap, "before", r_over)
				next_mat <- qm_overlap(data_sub, data_next, data_ref, data_break[j, ], break_no, ref_seq, y_overlap, "after", r_over)
					
			}
			
			if (!exists("sub_mat") | !exists("next_mat")) {next}
			
			#Compute quantiles
			if (!is.null(sub_mat) & !is.null(next_mat)) {
			
				colnames(sub_mat) <- c("Date", "Value")
				colnames(next_mat) <- c("Date", "Value")

				if (sum(is.na(sub_mat$Value)) != nrow(sub_mat)) {
						
					rq_seg <- qm_quantile(sub_mat)
					rq_base <- qm_quantile(next_mat)

				} else {
				
					next#no valid correlations
					
				} 
					
			} else {
			
				next #no overlapping period data
			
			} 
			
			#First correction computation
			
			if (method == "ratio") {
			
				adj <- (q_base / rq_base) / (q_seg / rq_seg)
					
			} else {
						
				adj <- (q_base - rq_base) - (q_seg - rq_seg)
					
			}
			
			adj[which(adj == Inf, arr.ind = TRUE)] <- NaN
			
			if (length(adj[!is.na(adj)]) >= length(qp)) {
						
				#Smoothing
				adj_sm <- qm_smooth(adj) #smoothed corrections
				adj_def <- matrix(NA, nrow = 12, ncol = length(qp)) #matrix of final corrections
				#loop on months to compute final corrections
				for (m in 1:12) {
					
					adj_mon <- adj_sm[m, ] #corrections of selected month
					
					if (length(adj_mon[!is.na(adj_mon)]) > 0) {
						
						adj_def[m, ] <- qm_slope(adj_mon, q_base[m, ], method) #check negative slopes
						
					} else {
								
						adj_def[m, ] <- rep(NA, length(qp))
							
					}
						
				}
					
			} else {
				
				adj_def <- matrix(NA, nrow = 12, ncol = length(qp))
				
			}
			
			#data adjusting
			vc <- cbind(vc, qm_correct(adj_def, data_sub, q_base, method))
			rm(sub_mat, next_mat)
			
		}
		
		#final homogeneous series
		if (is.null(vc)) {
					
			print("No correlated segments")
			next
				
		} else {
					
			if (ncol(vc) >= seg_min) {
						
				hom <- apply(vc, 1, median, na.rm = TRUE)
					
			} else {
						
				print("Correlated segments below threshold")
				next
					
			}
				
		}
		
		data_out[, move_out] <- data_out[, move_in]
		data_tmp <- data_out[, move_out][init[j]:break_start]
		
		if (sum(is.na(hom)) == length(hom)) {
					
			print(glue("No correction break {j}: no valid adjustment"))
				
		} else {
					
			print(glue("Break {j} Corrected"))
			if (zeros) {data_tmp[data_tmp < 0.1] <- NA}
			if (pcent) {data_tmp[(data_tmp > p_low) & (hom < 0) & !is.na(data_tmp)] <- NA}
			if (pcent) {data_tmp[(data_tmp < p_high) & (hom > 100) & !is.na(data_tmp)] <- NA}
			data_tmp[!is.na(data_tmp)] <- hom[!is.na(data_tmp)]
			#if (!is.na(th_low)) {data_tmp[data_tmp < th_low] <- th_low}
			#if (!is.na(th_high)) {data_tmp[data_tmp > th_high] <- th_high}
			data_out[, move_out][init[j]:break_start][!is.na(data_tmp)] <- data_tmp[!is.na(data_tmp)]
				
		}
		
	}
			
	data_out[, move_out] <- round(as.numeric(data_out[, move_out]), dec)
	fwrite(data_out, out_file)

}