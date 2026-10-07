ref_search <- function(data, datetime, variable, meta, meta_stat, data_dir, homo_dir = "", d = 50000, h = 100, r = 0.8, na_th = 80, wd = FALSE, sel_mode = "qc", n_min = 3, n_max = 10, p = FALSE) {
	
	stat_geo <- distGeo(cbind(meta_stat$Longitude, meta_stat$Latitude), cbind(meta$Longitude[meta$Code != meta_stat$Code], meta$Latitude[meta$Code != meta_stat$Code]))
	ref_stat <- which(stat_geo <= d)
	if (!is.na(h)) {
	
		stat_hgt <- abs(meta_stat$Elevation - meta$Elevation[meta$Code != meta_stat$Code])
		ref_stat <- ref_stat[stat_hgt[ref_stat] <= h]
		
	}
	
	if (length(ref_stat) > (n_max * 2)) {
		
		ref_tmp <- ref_stat[stat_geo[ref_stat] <= (0.5 * d)]
		
		if (length(ref_tmp) > n_max) {
		
			ref_stat <- ref_tmp
			
		} else {
		
			ref_tmp <- ref_stat[stat_geo[ref_stat] <= (0.75 * d)]
			
			if (length(ref_tmp) > n_max) {
		
				ref_stat <- ref_tmp
				
			} else {
			
				ref_stat <- ref_stat[stat_geo[ref_stat] <= (0.5 * d)]
			
			}
			
		}

	}
	
	meta_ref <- meta[meta$Code != meta_stat$Code]
	if (sel_mode == "homo") {
		
		data_ref <- data.frame(Date = as.Date(character()), Value = double(), Ref = character())
		
	} else {
	
		data_ref <- c()
		
	}
	for (i in ref_stat) {
	
		ref_file <- glue("{data_dir}{meta_ref$Provider[i]}/{meta_ref$Code[i]}_{meta_ref$Elevation[i]}_{meta_ref$Longitude[i]}_{meta_ref$Latitude[i]}_{meta_ref$Date_Start[i]}-{meta_ref$Date_End[i]}.txt")
		data_tmp <- fread(ref_file, select = c("Date", variable), col.names = c("Date", "Value")) %>% merge(data.frame(Date = datetime, Val = data), . ,all.x = TRUE) %>% select(., !Val)
		if (sel_mode == "homo") {
		
			data_ref <- data.frame(Date = data_tmp$Date, Value = data_tmp$Value, Ref = i) %>% rbind(data_ref, .)
			
		} else {
		
			data_ref <- cbind(data_ref, data_tmp$Value)
			
		}
		
	}
	
	if (length(ref_stat) >= n_min) {
	
		if (sel_mode == "homo") {
		
			if (variable %in% c("P_QC","G_QC","HN_QC")) {
			
				data_mon <- data.frame(Date = datetime, Value = data) %>% 
							mutate(Month = month(Date), Year = year(Date)) %>% 
							group_by(Month, Year) %>% 
							summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * sum(Value, na.rm = TRUE)) %>% 
							arrange(., Year)
							
				data_ref_mon <- data_ref %>% 
								mutate(Month = month(Date), Year = year(Date)) %>% 
								group_by(Month, Year, Ref) %>% 
								summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * sum(Value, na.rm = TRUE)) %>% 
								arrange(., Year, Ref)
			
			} else {
		
				data_mon <- data.frame(Date = datetime, Value = data) %>% 
							mutate(Month = month(Date), Year = year(Date)) %>% 
							group_by(Month, Year) %>% 
							summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * mean(Value, na.rm = TRUE)) %>% 
							arrange(., Year)
							
				data_ref_mon <- data_ref %>% 
								mutate(Month = month(Date), Year = year(Date)) %>% 
								group_by(Month, Year, Ref) %>% 
								summarise(Mon_mean = NA^(sum(!is.na(Value))/length(Value) < 0.8) * mean(Value, na.rm = TRUE)) %>% 
								arrange(., Year, Ref)
								
			}
								
			if (variable == "HS_QC") {
				
				data_diff <- shift(data_mon$Mon_mean[data_mon$Month %in% c(1:5, 11:12)], type = "lag")
				data_ref <- spread(data_ref_mon[data_ref_mon$Month %in% c(1:5, 11:12), ], Ref, Mon_mean) %>% arrange(., Year)
				data_ref <- as.data.frame(data_ref)[3:ncol(data_ref)] %>% as.matrix(.)
				ref_diff <- apply(data_ref, 2, function(x) {shift(x, type = "lag")})
				data_mon <- data_mon$Mon_mean[data_mon$Month %in% c(1:5, 11:12)]
		
			} else {
				
				data_ref <- spread(data_ref_mon, Ref, Mon_mean) %>% arrange(., Year)
				data_ref <- as.data.frame(data_ref)[3:ncol(data_ref)] %>% as.matrix(.)
				data_diff <- shift(data_mon$Mon_mean, type = "lag")
				ref_diff <- apply(data_ref, 2, function(x) {shift(x, type = "lag")})
				data_mon <- data_mon$Mon_mean
				
			}
		
				
		} else {
	
			if (variable == "HS_QC") {
		
				data_diff <- shift(data[month(datetime) %in% c(1:5, 11:12)], type = "lag")
				ref_diff <- apply(data_ref[month(datetime) %in% c(1:5, 11:12), ], 2, function(x) {shift(x, type = "lag")})
			
			} else {
	
				data_diff <- shift(data, type = "lag")
				ref_diff <- apply(data_ref, 2, function(x) {shift(x, type = "lag")})
			
			}
			
		}
		
		data_cor <- apply(ref_diff, 2, function(x) {cor(x, data_diff, method = "spearman", use = "pairwise.complete.obs")})
		if (sum((data_cor >= r) & !is.na(data_cor)) < n_min) {
		
			ref_stat <- ref_stat[(data_cor >= (r - 0.1)) & !is.na(data_cor)]
			data_ref <- data_ref[, (data_cor >= (r - 0.1)) & !is.na(data_cor)]
			
		} else {
			
			ref_stat <- ref_stat[(data_cor >= r) & !is.na(data_cor)]
			data_ref <- data_ref[, (data_cor >= r) & !is.na(data_cor)]
			
		}
		
	}
	
	if (length(ref_stat) >= n_min) {
	
		if (sel_mode == "homo") {
		
			na_data <- which(!is.na(data_mon))
			
		} else {
	
			na_data <- which(!is.na(data))
			
		}
		
		data_ref_tmp <- data_ref[na_data[1]:na_data[length(na_data)], ]
	
		nas <- apply(data_ref_tmp, 2, function(x) {100 * sum(!is.na(x)) / length(x)})
		ref_stat <- ref_stat[nas >= na_th]
		data_ref <- data_ref[, nas >= na_th]
	
	}
	
	if (length(ref_stat) < n_min) {
	
		print(glue("Warning: Station {meta_stat$Name} has not enough reference stations"))
		data_ref <- c()
		return(data_ref)
		
	} else {
	
		if (length(ref_stat) > n_max) {
		
			st_df <- data.frame(Ref = ref_stat, Dis = stat_geo[ref_stat], Id = 1:length(ref_stat))
			st_df <- st_df[order(st_df$Dis), ]
			data_ref <- data_ref[, st_df$Id[1:n_max]]
			ref_stat <- ref_stat[st_df$Id[1:n_max]]
			
		}
	
		if (wd) {
	
			fwrite(as.data.frame(stat_geo[ref_stat]), "./wetdry_stations.txt")
		
		}
	
		if (p) {
	
			data_tmp <- data_ref
			data_ref <- array(NA, dim = c(nrow(data_tmp), ncol(data_tmp), 4))
			data_ref[, , 1] <- data_tmp
		
			for (i in ref_stat) {
			
				ref_file <- glue("{data_dir}{meta_ref$Provider[i]}/{meta_ref$Code[i]}_{meta_ref$Elevation[i]}_{meta_ref$Longitude[i]}_{meta_ref$Latitude[i]}_{meta_ref$Date_Start[i]}-{meta_ref$Date_End[i]}.txt")
				data_t <- fread(ref_file, select = c("Date", "T_QC"), col.names = c("Date", "Value")) %>% merge(data.frame(Date = datetime, Val = data), . ,all.x = TRUE) %>% select(., !Val)
				data_ref[, which(ref_stat == i), 2] <- data_t$Value
			
				ref_file <- glue("{data_dir}{meta_ref$Provider[i]}/{meta_ref$Code[i]}_{meta_ref$Elevation[i]}_{meta_ref$Longitude[i]}_{meta_ref$Latitude[i]}_{meta_ref$Date_Start[i]}-{meta_ref$Date_End[i]}.txt")
				data_t <- fread(ref_file, select = c("Date", "RH_QC"), col.names = c("Date", "Value")) %>% merge(data.frame(Date = datetime, Val = data), . ,all.x = TRUE) %>% select(., !Val)
				data_ref[, which(ref_stat == i), 3] <- data_t$Value
			
				data_ref[, which(ref_stat == i), 4] <- meta_ref$Elevation[i]
			
			}
		
		}
		
		if (sel_mode != "homo") {
	
			print(glue("{variable} Reference stations selected: {length(ref_stat)}"))
			
		}
		
		if (sel_mode == "homo") {
		
			meta_df <- meta[meta$Code != meta_stat$Code, ][ref_stat, ]
			fwrite(meta_df, glue("{homo_dir}reference_stations.txt"))
			
		}
	
		return(data_ref)
		
	}
	
}