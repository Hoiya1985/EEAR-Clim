create_monthly <- function(work_dir, variable, date_start, date_end, homo_dir) {

	regions <- Sys.glob(glue("{work_dir}*"))
	date_long <- seq(as.Date(glue("{date_start}")), as.Date(glue("{date_end}")), by = "1 day")
	options(dplyr.summarise.inform = FALSE)
	meta_stat <- data.frame(Code = character(), Region = character(), Years = double())
	mat <- c()
	for (i in regions) {
	
		file_in <- Sys.glob(glue("{work_dir}{i}/*.txt"))
	
		for (file in file_in) {
		
			code <- str_split(file, glue("{work_dir}{i}/"))[[1]][2] %>% str_split(., "_") %>% .[[1]] %>% .[1]
			
			data_in <- fread(file)
			data_in <- data.frame(Date = data_in$Date, Value = as.data.frame(data_in)[which(colnames(data_in) == variable)])
			
			nas <- which(!is.na(data_in$Value))
			yr <- interval(as.Date(glue("{data_in$Date[nas[1]]}"), "%Y%m%d"), as.Date(glue("{data_in$Date[nas[length(nas)]]}"), "%Y%m%d")) %>% 
				  time_length(., "years")
			meta_stat <- data.frame(Code = code, Region = i, Years = yr) %>% rbind(meta_stat, .)
			
			day_date <- data.frame(Date = date_long)
			day_df <- merge(day_date, data_in, by.y = "Date", all.x = TRUE)
		
			if (variable %in% c("P", "G", "HN")) {
			
				day_df <- day_df %>% 
						  mutate(Year = year(Date), Month = month(Date)) %>% 
						  group_by(Month, Year) %>% 
						  summarise(Mon_mean = NA^(sum(!is.na(as.numeric(Value))) / length(Value) < 0.8) * sum(as.numeric(Value), na.rm = TRUE))
		
			} else {
				day_df <- day_df %>% 
						  mutate(Year = year(Date), Month = month(Date)) %>% 
						  group_by(Month, Year) %>% 
						  summarise(Mon_mean = NA^(sum(!is.na(as.numeric(Value))) / length(Value) < 0.8) * mean(as.numeric(Value), na.rm = TRUE))
			}
		
			day_df <- day_df[order(day_df$Year), ]
			day_df$Mon_mean <- round(as.numeric(day_df$Mon_mean), 1)
			day_df$Mon_mean[is.na(day_df$Mon_mean)] <- NaN
			mat <- rbind(mat, day_df$Mon_mean)
		}
	}

	out_dat <- glue("{homo_dir}{variable}_{date_start}-{date_end}.dat")
	fwrite(data.table(mat), out_dat, sep = " ", col.names = FALSE)
	fwrite(meta_stat, glue("{homo_dir}{variable}_lengths.txt"))

}

climatol_prepare <- function(mat_ref, meta, meta_stat, out_dir, variable, y_start, y_end, suffix, missing_flag) {
	
	in_df <- as.data.frame(cbind(meta_stat$Longitude, meta_stat$Latitude, meta_stat$Elevation, glue("{meta_stat$Code}"), glue("{meta_stat$Name}_{meta_stat$Code}")))
	in_df <- as.data.frame(cbind(meta$Longitude, meta$Latitude, meta$Elevation, glue("{meta$Code}"), glue("{meta$Name}_{meta$Code}"))) %>% rbind(in_df, .)
	out_df <- in_df[!duplicated(in_df), ]
	meta_out <- glue("{out_dir}{variable}-{suffix}_{y_start}-{y_end}.est")
	write.table(out_df, meta_out, row.names = FALSE, col.names = FALSE)
	fread(meta_out, sep = " ") %>% write.table(., meta_out, row.names = FALSE, col.names = FALSE)
	
	mat_ref[is.na(mat_ref)] <- missing_flag
	mat_ref <- mat_ref[!duplicated(in_df), ]
	data_out <- glue("{out_dir}{variable}-{suffix}_{y_start}-{y_end}.dat")
	fwrite(data.table(mat_ref), data_out, sep = " ", col.names = FALSE)
	
	file_out <- glue("{out_dir}{variable}-{suffix}_{y_start}-{y_end}")
	return(file_out)
	
}

rh_prepare <- function(data_stat, data_ref, base_out, ref_out, missing_flag, y_start, y_end, tau_cor = 0.3) {

	#computation of reference series
	data_cor <- apply(data_ref, 2, function(x) {cor(as.numeric(x), as.numeric(data_stat), method = "pearson", use = "pairwise.complete.obs")})
	weight <- exp(-(1 - data_cor^2) / (tau_cor^2/log(2)))
	ref_series <- apply(data_ref, 1, function(x) weighted.mean(as.numeric(x), weight, na.rm = TRUE))
	yy <- matrix(rep(seq(y_start, y_end, by = 1), 12), nrow = length(y_start:y_end), ncol = 12) %>% t(.) %>% as.vector(.)
	data_stat[is.na(data_stat)] <- missing_flag
	ref_series[is.na(ref_series)] <- missing_flag
	#write input files
	base_stat <- data.frame(Year = yy, Month = sprintf("%2d", rep(1:12, length(y_start:y_end))), Day = sprintf("%02d", rep(0, length(yy))), 
							Data = sprintf("%6.1f", round(as.numeric(data_stat), 1)))
	fwrite(base_stat, base_out, sep = "\t", col.names = FALSE)
	refs_stat <- data.frame(Year = yy, Month = sprintf("%2d", rep(1:12, length(y_start:y_end))), Day = sprintf("%02d", rep(0, length(yy))),
						    Data = sprintf("%6.1f", round(ref_series, 1)))
	fwrite(refs_stat, ref_out, sep = "\t", col.names = FALSE)
	
}

acmant_options <- function(variable, reg, tr = "m", of = "N", op = "Y", seas = "S", m1 = 11, m2 = 5) {

	vars <- c("T", "TMIN", "TMAX", "P", "SP", "RH", "WS", "WD", "G", "HS", "HN")
	v_code <- c("TT", "TT", "TT", "RR", "PP", "HH", "FF", "FF", "SS", "RR", "RR")
	
	if (variable %in% vars) {
	
		var_code <- v_code[vars == variable]
		
	} else {
	
		var_code <- "TT"
		
	}
	
	options_df <- data.frame(Name = reg, Code = var_code, Res = tr, Outlier = of, Output = op, Season = seas, Snow_Start = m1, Snow_End = m2)
	
	return(options_df)
	
}
	
acmant_prepare <- function(mat_ref, meta, meta_sel, y_start, y_end, out_dir, stat_code, missing_flag) {

	for (i in 1:nrow(mat_ref)) {
	
		num_code <- sprintf("%04.f", i)
		ac_txt <- glue("{out_dir}{stat_code}{num_code}.txt")
		mat_ac <- y_start:y_end
		mat_ac <-  matrix(mat_ref[i, ], nrow = 12) %>% t(.) %>% cbind(mat_ac, .)
		mat_ac[is.na(mat_ac)] <- missing_flag
		cat(meta_sel$Name[i], file = ac_txt, sep = "\n")
		write.table(mat_ac, ac_txt, append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
	}
	
}

climatol_run <- function(meta_stat, var_name, y_start, y_end, missing_flag, var_type, snht, out, clim_dir, main_dir) {

	setwd(clim_dir)
	h_saf <- safely(homogen)
	p <- h_saf(varcli = glue("{var_name}-{meta_stat}"), anyi = y_start, anyf = y_end, test = "snht", std = var_type, swa = 60, dz.max = out, inht = snht, 
			   na.strings = missing_flag, sufbrk = "", raway = 1, verb = FALSE)
	if (!is.null(p$error)) {
	
		closeAllConnections()
		print(p$error)

	}
	setwd(main_dir)

}

rh_run <- function(base_out, ref_out, outfile, missing_flag) {

	FindU.wRef(Bseries = base_out, Rseries = ref_out, output = outfile, MissingValueCode = glue("{missing_flag}"), p.lev = 0.95, Iadj = 10000, 
			   Mq = 10, GUI = F, Ny4a = 0)
	StepSize.wRef(Bseries = base_out, Rseries = ref_out, InCs = glue("{outfile}_mCs.txt"), output = outfile, 
				  MissingValueCode = glue("{missing_flag}"), p.lev = 0.95, Iadj = 10000, Mq = 10, GUI = F, Ny4a = 0)
	out_tmp <- readLines(glue("{outfile}_fCs.txt"), -1)
	out_tmp <- str_replace_all(out_tmp, "0.950","0.950 ")
	writeLines(out_tmp, glue("{outfile}_fCs.txt"))
	if (length(out_tmp) != 1) {
		out_sign <- fread(glue("{outfile}_fCs.txt"), skip = 1)
		v_test <- min(as.numeric(out_sign$V8), na.rm = TRUE)
		th_tmp <- str_split(out_sign$V10, "-")
		th_min <- th_tmp[[which(out_sign$V8 == v_test)]][1] %>% as.numeric(.)
		if (v_test < th_min) {
			out_sign <- out_sign[-which(out_sign$V8 == v_test), ]
			out_new <- readLines(glue("{outfile}_fCs.txt"), 1)
			for (i in 1:nrow(out_sign)) {
				out_new <- c(out_new, paste(out_sign[i, ], collapse = " "))
			}
			writeLines(out_new, glue("{outfile}_fCs.txt"))
		}
	}
}

y_border <- function(break_df, metay_start, metay_end, yth_low = 2, yth_high = 2) {

	y_diff_start <- interval(as.Date(break_df$Date), as.Date(glue("{metay_start}"), "%Y%m%d")) %>% time_length(., "years") %>% abs(.)
	y_diff_end <- interval(as.Date(glue("{metay_end}"), "%Y%m%d"), as.Date(break_df$Date)) %>% time_length(., "years") %>% abs(.)
	break_df <- break_df[(y_diff_start > yth_low) & (y_diff_end > yth_high), ]
	
	return(break_df)

}

y_between <- function(break_df, yth = 2) {

	brk_shift <- shift(break_df$Date, type = "lag")
	y_diff <- interval(as.Date(break_df$Date), as.Date(brk_shift)) %>% time_length(., "years") %>% abs(.)
	off_shift <- which(y_diff <= yth)
		
	if (length(off_shift) != 0) {
		
		off_del <- c()
		val_shift <- shift(break_df$Value, type = "lag")
		
		for (i in off_shift) {
			
			if (break_df$Value[i] <= val_shift[i]) {
				
				off_del <- c(off_del, i)
					
			} else {
				
				off_del <- c(off_del, i - 1)
					
			}
				
		}
			
		break_df <- break_df[-off_del, ]
			
	}
	
	return(break_df)

}

y_comp <- function(homo_tmp, homo_df, code, methods = c("Climatol", "RH_Test", "Acmant"), th = 2) {
	
	m_vec <- c()
	for (i in seq_along(methods)) {m_vec <- rbind(m_vec, c(methods[i], methods[-i]))}
	
	for (i in methods) {
	
		test_df <- homo_tmp[homo_tmp$Method == i, ]
		
		if (nrow(test_df) != 0) {
		
			c_remove <- c()
			
			im <- m_vec[methods == i, -1]
		
			for (j in 1:nrow(test_df)) {
			
				for (k in im) {
		
					ref_df <- homo_tmp[homo_tmp$Method == k, ]
					n_skip <- which(homo_tmp$Method == k)[1] - 1
			
					if (nrow(ref_df) != 0) {
			
						yc <- interval(test_df$Date[j], ref_df$Date) %>% time_length(., "years") %>% abs(.)
						n_pass <- which(yc <= th)
						pass_min <- which((yc <= th) & (yc == min(yc))) + n_skip
			
						if (length(n_pass) != 0) {
			
							homo_df <- data.frame(Code = code, 
											      Date = test_df$Date[j], 
											      Value = test_df$Value[j], 
											      Method = i
											     ) %>% rbind(homo_df, .)					 
							homo_tmp <- homo_tmp[-pass_min, ]
							c_remove <- c(c_remove, j)
							
						}
				
					}
					
				}
			
			}
		
			c_remove <- unique(c_remove)
			if (!is.null(c_remove)) {homo_tmp <- homo_tmp[-c_remove, ]}
			
		}
		
	}
	
	return(homo_df)

}

qm_overlap <- function(data_sub, data_next, data_ref, data_brk, ref_brk, dr, yrs, seg, r_th) {

	#count days of overlapping around breakpoints
	int_start <- length(intersect(data_sub$Date, dr))
	int_end <- length(intersect(data_next$Date, dr))
							
	if ((int_start >= (365 * yrs)) & (int_end >= (365 * yrs))) {
	
		if ((as.Date(data_brk$Date) - dr[1]) > 7300) { #maximum of overlapping 20 years considered
		
			cut_start <- which(data_ref$Date == (data_brk$Date - 7300))
			
		} else {
		
			cut_start <- which(data_ref$Date == dr[1])
			
		}
		
		if (!is.na(ref_brk$Date)) {
			
			ref <- ref_brk$Date - 1
			
		} else {
		
			ref <- data_ref$Date[nrow(data_ref)]
			
		}
		
		if ((ref - data_brk$Date) > 7300) {
		
			cut_end <- which(data_ref$Date == (data_brk$Date + 7300))
		
		} else {
			
			if (!is.na(ref_brk$Date)) {
				
				cut_end <- which(data_ref$Date == (ref_brk$Date - 1))
			
			} else {
				
				cut_end <- nrow(data_ref)
			
			}
			
		}
		
		if (length(cut_start) == 0) {cut_start <- 1}
		if (length(cut_end) == 0) {cut_end <- nrow(data_ref)}
		ref_sub <- data_ref[cut_start:cut_end, ] #cutting
		
		#computation of correlation
		ref_cor <- merge(data_next, ref_sub, all.x = TRUE)
		rc <- cor(as.numeric(ref_cor$Value), as.numeric(ref_cor$Val), method = "pearson", use = "pairwise.complete.obs")
		
		#if reference is enough correlated create overlaped matrices
		if ((rc >= r_th) & !is.na(rc)) {
		
			sub_mat <- merge(data_sub, ref_sub, all.x = TRUE) %>% select(., !Value) #segment reference matrix
			next_mat <- merge(data_next, ref_sub, all.x = TRUE) %>% select(., !Value) #basis reference matrix
			
		} else {
		
			sub_mat <- NULL
			next_mat <- NULL
			
		}
		
	} else {
	
		sub_mat <- NULL
		next_mat <- NULL
		
	}

	if (seg == "before") {
	
		return(sub_mat)
		
	}
	if (seg == "after") {
	
		return(next_mat)
		
	}
}

qm_quantile <- function(data_sub, qp = seq(0.05, 0.95, by = 0.05)) {

	#Month feature filling
	data_sub$Mon <- month(as.Date(data_sub$Date))
	#Quantiles matrix
	q_mat <- matrix(0, nrow = 12, ncol = length(qp))
	#every monthly quantile also consider data of adiacent months
	#data_q: data for quantiles
	for (i in 1:12) {
		#January
		if (i == 1) { #January
		
			data_q <- data_sub$Value[(data_sub$Mon > 11) | (data_sub$Mon < 3)]
			
		} else {
		
			if (i == 12) { #December
			
				data_q <- data_sub$Value[(data_sub$Mon > 10) | (data_sub$Mon < 2)]
				
			} else { #Other months
			
				data_q <- data_sub$Value[(data_sub$Mon >= (i - 1)) & (data_sub$Mon <= (i + 1))]
				
			}
			
		}
		
		q_mat[i,] <- quantile(x = as.numeric(data_q), probs = qp,na.rm = TRUE)
		#q_mat[i,][(q_mat[i,] > 0) & (q_mat[i,] < 0.1)] <- NaN
	}
	
	return(q_mat)
	
}

qm_smooth <- function(adj_mat, qp = seq(0.05, 0.95, by = 0.05)) {

	qpl <- length(seq(0.05, 0.95, by = 0.05))
	smooth_out <- matrix(NA, nrow = 12, ncol = qpl)
	
	for (i in 1:12) {
	
		if (i == 1) { #January
		
			for (j in 1:qpl) {
			
				if (j == 1) { #first percentile
				
					smooth_out[i, j] <- (adj_mat[12, j] + 2 * adj_mat[i, j] + adj_mat[i+1, j] + adj_mat[i, j+1]) / 5
					
				} else { #last percentile
				
					if (j == qpl) {
					
						smooth_out[i, j] <- (adj_mat[12, j] + 2 * adj_mat[i, j] + adj_mat[i+1, j] + adj_mat[i, j-1]) / 5
						
					} else { #other percentiles
					
						smooth_out[i, j] <- (adj_mat[12, j] + adj_mat[i, j] + adj_mat[i+1, j] + adj_mat[i, j-1] + adj_mat[i, j+1]) / 5
						
					}
					
				}
				
			}
			
		} else {
		
			if (i == 12) {#December
			
				for (j in 1:qpl) {
				
					if (j == 1) {#first percentile
					
						smooth_out[i, j] <- (adj_mat[i-1, j] + 2 * adj_mat[i, j] + adj_mat[1, j] + adj_mat[i, j+1]) / 5
						
					} else {
					
						if (j == qpl) { #last percentile
						
							smooth_out[i, j] <- (adj_mat[i-1, j] + 2 * adj_mat[i, j] + adj_mat[1, j] + adj_mat[i, j-1]) / 5
							
						} else {
						
							smooth_out[i, j] <- (adj_mat[i-1, j] + adj_mat[i, j] + adj_mat[1, j] + adj_mat[i, j-1] + adj_mat[i, j+1]) / 5
							
						}
						
					}
					
				}
				
			} else {
			
				for (j in 1:qpl) {
				
					if (j == 1) {
					
						smooth_out[i, j] <- (adj_mat[i-1, j] + 2 * adj_mat[i, j] + adj_mat[i+1, j] + adj_mat[i, j+1]) / 5
						
					} else {
					
						if (j == qpl) {
						
							smooth_out[i, j] <- (adj_mat[i-1, j] + 2 * adj_mat[i, j]+adj_mat[i+1, j] + adj_mat[i, j-1]) / 5
							
						} else {
						
							smooth_out[i, j] <- (adj_mat[i-1, j] + adj_mat[i, j] + adj_mat[i+1, j]+ adj_mat[i, j-1] + adj_mat[i, j+1]) / 5
							
						}
						
					}
					
				}
				
			}
			
		}
		
	}
	
	return(smooth_out)
	
}

qm_slope <- function(adj_vec, q_vec, method, qp = seq(0.05, 0.95, by = 0.05)) {

	qpl <- length(qp)

	if ((qpl %% 2) != 0) {
  
		q_inf <- (qpl - 1) / 2
		q_sup <- ((qpl + 1) / 2) + 1
	
		if (method == "ratio") {
	
			for (i in q_sup:qpl) {
			
				nas <- c(adj_vec[i], adj_vec[i-1], q_vec[i], q_vec[i-1])
			
				if (sum(is.na(nas)) == 0) {
				
					if ((adj_vec[i] / adj_vec[i-1]) < (q_vec[i-1] / q_vec[i])) {
				
						adj_vec[i] <- adj_vec[i-1] * q_vec[i-1] / q_vec[i]
					
					}
				
				}
			
			}
		
			for (i in q_inf:1) {
		
				nas <- c(adj_vec[i], adj_vec[i+1], q_vec[i], q_vec[i+1])
			
				if (sum(is.na(nas)) == 0) {
			
					if ((adj_vec[i] / adj_vec[i+1]) > (q_vec[i+1] / q_vec[i])) {
				
						adj_vec[i] <- adj_vec[i+1] * q_vec[i+1] / q_vec[i]
					
					}
				
				}
			
			}
		
		} else {
	
			for (i in q_sup:qpl) {
		
				nas <- c(adj_vec[i], adj_vec[i-1], q_vec[i], q_vec[i-1])
			
				if (sum(is.na(nas)) == 0) {
			
					if ((adj_vec[i] - adj_vec[i-1]) < (q_vec[i-1] - q_vec[i])) {
				
						adj_vec[i] <- adj_vec[i-1] + q_vec[i-1] - q_vec[i]
					
					}
				
				}
			
			}
		
			for (i in q_inf:1) {
		
				nas <- c(adj_vec[i], adj_vec[i+1], q_vec[i], q_vec[i+1])
			
				if (sum(is.na(nas)) == 0) {
			
					if ((adj_vec[i] - adj_vec[i+1]) > (q_vec[i+1] - q_vec[i])) {
				
						adj_vec[i] <- adj_vec[i+1] + q_vec[i+1] - q_vec[i]
					
					}
				
				}
			
			}
		
		}
	
	}
  
	return(adj_vec)
  
}

qm_correct <- function(adj_mat, data_seg, q_seg, method, qp = seq(0.05, 0.95, by = 0.05)) {

  seg_out <- rep(NA, nrow = data_seg)
  
	for (i in 1:nrow(data_seg)) {
	
		m <- as.numeric(month(as.Date(data_seg$Date[i])))
		
		if (!is.na(data_seg$Value[i]) & !is.na(q_seg[m, 1]) & (data_seg$Value[i] < q_seg[m, 1])) {
			
			adj_out <- adj_mat[m, 1]
		
		} else {
			
			if (!is.na(data_seg$Value[i]) & !is.na(q_seg[m, length(qp)]) & (data_seg$Value[i] > q_seg[m, length(qp)])) {
				
				adj_out <- adj_mat[m, length(qp)]
			
			} else {
				
				low <- qp[which(q_seg[m, ] <= data_seg$Value[i])]
				high <- qp[which(q_seg[m, ] >= data_seg$Value[i])]
				
				if (!length(low) | !length(high)) {
				
					seg_out[i] <- data_seg$Value[i]
					next
					
				}
				
				i_qp <- which(qp == (min(high) + max(low)) / 2)
				
				if ((max(low) < min(high)) & ((min(high) - max(low)) == 0.05) & (length(i_qp) != 0)) {
					
					adj_out <- adj_mat[m, i_qp]
				
				} else {
					
					if (max(low) == min(high)) {
						
						thv <- (abs(max(low) - 0.5) * sign(max(low) - 0.5)) + 0.5
						adj_out <- adj_mat[m, which.min(abs(qp - thv))]
					
					} else {
						
						if ((max(low) < 0.5) & (min(high) < 0.5)) {
							
							adj_out <- adj_mat[m, which(qp == max(low))]
						
						} else {
							
							if ((max(low) > 0.5) & (min(high) > 0.5)) {
								
								adj_out <- adj_mat[m, which(qp == min(high))]
							
							} else {
								
								adj_out <- adj_mat[m, which(qp == 0.5)]
							
							}
						
						}
					
					}
				
				}
			
			}
		
		}
		
		if (!length(adj_out)) {
		
			seg_out[i] <- data_seg$Value[i]
			
		} else {
		
			if (method == "ratio") {
		
				seg_out[i] <- data_seg$Value[i] * adj_out
		
			} else {
			
				seg_out[i] <- data_seg$Value[i] + adj_out
		
			}
			
		}
		#seg_out[i] <- adj_out
	}
	
	return(seg_out)
	
}