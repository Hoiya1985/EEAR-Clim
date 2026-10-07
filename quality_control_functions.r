phi_test <- function(data, variable, th_usr = c(NaN, NaN)) {

	flag_vec <- integer(length(data))
	th_def <- data.frame(Var = c("T", "TMIN", "TMAX", "P", "SP", "RH", "WS", "WD", "G", "HS", "HN", "CC"),
						 Lower = c(-50, -50, -50, 0, 500, 0, 0, 0, 0, 0, 0, 0),
						 Upper = c(50, 50, 50, 900, 1100, 100, 75, 360, 138.24, 800, 250, 100))
	
	if (is.na(th_usr[1]) & is.na(th_usr[2])) {
	
		th_min <- th_def$Lower[th_def$Var == variable]
		th_max <- th_def$Upper[th_def$Var == variable]
		
	} else {
	
		th_min <- th_usr[1]
		th_max <- th_usr[2]
		
	}
	
	flag_vec[((data > th_max) | (data < th_min)) & !is.na(data)] <- 1
	
	return(flag_vec)

}

repeat_test <- function(data, variable, rep_th = 5, rep0 = 180, wind_na = 90) {
	
	flag_vec <- integer(length(data))
	
	re_df <- with(rle(data), data.frame(number = values, length = lengths, start = cumsum(lengths) - lengths + 1, end = cumsum(lengths)))
	
	if ((variable != "P") & (variable != "HN")) {
		
		re <- re_df[(re_df$length >= rep_th) & !is.na(re_df$number), ]
		if (nrow(re) > 0) {
			flags <- mapply(x = re$start, y = re$end, function(x, y) {seq(x, y, by = 1)}) %>% unlist(.) %>% as.vector(.)
			
			if (variable %in% c("WS", "WD")) {
			
				for (i in flags) {
				
					if ((i > 3) & (i < (length(data) - 3))) {
					
						fg_na <- 100 * sum(is.na(data[(i - 3):(i + 3)])) / length(data[(i - 3):(i + 3)])
						if (fg_na >= wind_na) {
							flag_vec[i] <- 1
						}
						
					} else {
					
						flag_vec[i] <- 1
						
					}
				}
		
				
			} else {
			
				flag_vec[flags] <- 1
				
			}
			
		}

	} else {
	
		re <- re_df[(re_df$length >= rep_th) & (re_df$number != 0) & !is.na(re_df$number), ]
		if (nrow(re) > 0) {
		
			flags <- mapply(x = re$start, y = re$end, function(x, y) {seq(x, y, by = 1)}) %>% unlist(.) %>% as.vector(.)
			flag_vec[flags] <- 1
			
		}
		
		if (variable == "P") {
	
			re_df0 <- with(rle(data), data.frame(number = values, length = lengths, start = cumsum(lengths) - lengths + 1,end = cumsum(lengths)))
			re0 <- re_df0[(re_df0$length >= rep0) & (re_df0$number == 0) & !is.na(re_df0$number), ]
			if (nrow(re0) > 0) {
		
				zeros <- mapply(x = re0$start, y = re0$end, function(x, y) {seq(x, y, by = 1)}) %>% unlist(.) %>% as.vector(.)
				flag_vec[zeros] <- 1
			
			}
			
		}
	}
	
	return(flag_vec)					

}

time_step <- function(data, variable, day_th = NA) {

	day_df <- data.frame(Var = c("T", "TMIN", "TMAX"),
						 Value = c(20, 20, 20))
	hs_df <- c(150, -50, 50)
						 
	flag_vec <- integer(length(data))
	
	if (is.na(day_th)) {
		day_th <- day_df$Value[day_df$Var == variable]
	}
	
	if (variable %in% c("T", "TMIN", "TMAX")) {
	
		t_bef <- abs(data - shift(data, type="lag"))
		t_aft <- abs(data - shift(data, type="lead"))
		flag_vec[((t_bef > day_th) & (t_aft > day_th))] <- 1
		
	}
	
	if (variable == "WS") {
	
		dv <- abs(data - shift(data, type="lag"))
		tlim <- quantile(data, 0.95, na.rm = TRUE) + (8.9 * IQR(data, na.rm = TRUE))
		flag_vec[dv > tlim] <- 1
	}
	
	if (variable == "HS") {
	
		dv_bef <- data - shift(data, type="lag")
		dv_aft <- data - shift(data, type="lead")
		flag_vec[(dv_bef > hs_df[1]) | (dv_bef < hs_df[2]) | ((abs(dv_bef) >= hs_df[3]) & (abs(dv_aft) >= hs_df[3]))] <- 1
		
	}
	
	return(flag_vec)
	
}
		

outlier <- function(data, variable, out_th = 3, out99 = NA) {
	
	flag_vec <- integer(length(data))
	
	if (variable %in% c("T", "TMIN", "TMAX", "SST")) {
	
		out_mad <- abs(data - median(data, na.rm = TRUE)) / mad(data, na.rm = TRUE)
		flag_vec[(out_mad > out_th) & !is.na(out_mad)] <- 1
		
	} else {
	
		if (variable %in% c("HS", "HN")) {
		
			flag_vec[data >= (out_th * quantile(data[data != 0], 0.95, na.rm = TRUE))] <- 1
			if (!is.na(out99)) {
				flag_vec[data >= (out99 * quantile(data[data != 0], 0.99, na.rm = TRUE))] <- 1
			}
		} else {
	
			med <- median(data, na.rm = TRUE)
			w1 <- max(c(min(data, na.rm = TRUE),quantile(data, 0.25, na.rm = TRUE) -(1.5 * IQR(data, na.rm = TRUE))), na.rm = TRUE)
			w2 <- min(c(max(data, na.rm = TRUE),quantile(data, 0.75, na.rm = TRUE) +(1.5 * IQR(data, na.rm = TRUE))), na.rm = TRUE)
		
			mad1 <- (data[data > med] - med) / (w2 - med)
			mad2 <- (med - data[data < med]) / (med - w1)
		
			flag_vec[!is.na(mad1) & (mad1 > out_th)] <- 1
			flag_vec[!is.na(mad2) & (mad2 > out_th)] <- 1
			
		}
		
	}
	
	return(flag_vec)

}

clim_temperature <- function(data_t = NULL, data_min = NULL, data_max = NULL, test = "", th_min = 0.01, th_max = 30) {
	
	flag_vec <- integer(length(data))
	
	if (test == "complete") {
	
		flag_vec[((data_min > data_t) & (data_max < data_t)) & !is.na(data_t) & !is.na(data_min) & !is.na(data_max)] <- 1
		
	}
	
	if (test == "partial_min") {
	
		flag_vec[(data_min > data_max) & !is.na(data_max) & !is.na(data_min)] <- 1
		
	}
	
	if (test == "partial_max") {
	
		flag_vec[(data_max < data_min) & !is.na(data_min) & !is.na(data_max)] <- 1
		
	}
	
	if (test == "minmax") {
	
		flag_vec[((data_max-data_min) < th_min) | ((data_max-data_min) > th_max)] <- 1
		
	}
	
	return(flag_vec)
	
}

clim_precipitation <- function(data, data_t, th1 = 9, th2 = 5) {
	
	flag_vec <- integer(length(data))
	
	if (sum(is.na(data_t)) != 0) {
	
		q1 <- which((data >= (th1 * quantile(data[data != 0], 0.95, na.rm = TRUE))) & ((data_t > 0) | is.na(data_t)))
		q2 <- which((data >= (th2 * quantile(data[data != 0], 0.95, na.rm = TRUE))) & (data_t <= 0))
		
		q_out <- c(q1, q2)
		q_out <- q_out[which(data_t[q_out] < 0)]
		
	} else {
	
		q_out <- which(data >= (th1 * quantile(data[data != 0], 0.95, na.rm = TRUE)))
		
	}
	
	flag_vec[q_out] <- 1
	
	return(flag_vec)	

}

#lawrence 2005; not applied
clim_rh <- function(data, data_t, a = 17.625, b = 243.04, th_td = 0.2, th_min = -7, th_max = 12) {

	flag_vec <- integer(length(data))
	
	rat <- (a * data_t) / (b + data_t)
	td <- (b * (log(0.01 * data) + rat)) / (a - log(0.01 * data) - rat)
		
	flag_vec[((td - data_t) > th_td) & !is.na(td - data_t)] <- 1
		
	rh <- 100 * exp(((a * td) / (b + td)) - rat)
		
	flag_vec[(((rh - data) <= th_min) | ((rh - data) >= th_max)) & !is.na(rh - data)] <- 1
	
	return(flag_vec)

}

clim_wind <- function(data_speed, data_dir) {

	flag_vec <- integer(length(data_speed))

	wind_var <- which((data_speed > 0.1) & !is.na(data_speed))
	wind_const <- which((data_speed <= 0.1) & !is.na(data_speed))
	
	flag_vec[(data_speed > 0.1) & (data_dir == 0) & !is.na(data_speed) & !is.na(data_dir)] <- 1
	flag_vec[(data_speed <= 0.1) & (data_dir != 0) & !is.na(data_speed) & !is.na(data_dir)] <- 1
	
	return(flag_vec)

}

hypso_test <- function(data_x, data_y, z_x, z_y, e0 = 6.11, mr0 = 0.622, a = 17.625, b = 243.04, rd = 287.05, g = 9.81) {

	flag_vec <- integer(length(data))
	
	es <- e0 * exp((a * data_x[, 2]) / (data_x[, 2] + b))
	mr <- (0.01 * data_x[, 3] * 1000 * mr0 * es) / (data_x[, 1] - es)
	tv_x <- (273.15 + (data_x[, 2]) * (1 + (mr / mr0)) / (1 + mr))
	
	es <- e0 * exp((a * data_y[, 2]) / (data_y[, 2] + b))
	mr <- (0.01 * data_y[, 3] *1000 * mr0 * es) / (data_y[, 1] - es)
	tv_y <- (273.15 + (data_y[, 2]) * (1 + (mr / mr0)) / (1 + mr))
	
	tv <- apply(cbind(tv_x, tv_y), 1, mean, na.rm = TRUE)
	
	p_h <- data_y[, 1] * exp(- (g * (z_x - z_y) / (rd * tv)))

	return(p_h)
	
}

#allen et al 1996; not applied
clear_sky_rad <- function(data, data_p, dates, lat) {

	flag_vec <- integer(length(data))
	
	p <- 0.1 * data_p
	s <- 1367 * 0.0864
	dd <- lubridate::day(dates)
	mm <- month(dates)
	yy <- year(dates)
	
	doy <- as.integer((275 * mm / 9) - 30 + dd) - 2 #Day of the year
	doy[mm < 3] <- doy[mm < 3] + 2
	doy[(mm > 2) & lubridate::leap_year(doy)] <- doy[(mm > 2) & lubridate::leap_year(doy)] + 1
	delta <- 0.409 * sin((2 * pi * doy / 365) - 1.39)
	dr <- 1 + (0.033 * cos(2 * pi * doy / 365))
	phi <- lat * pi /180
	sun <- acos(-(tan(phi) * tan(delta)))
	rad_et <- s * dr * ((sun * sin(phi) * sin(delta)) + (cos(phi) * cos(delta) * sin(sun))) #extraterrestrial radiation
	ci <- data / rad_et
	ktb <- 0.5 + ((ci - 0.7) * 5) #clarity coefficient
	kt <- exp(-((0.0018 * p) / (ktb * sin(45)))) #clearness index
	rad_cs <- kt * rad_et #clear-sky radiation
	cs <- apply(cbind(data, rad_cs), 1, sum)
	if (!is.na(sum(cs, na.rm = TRUE)) & (sum(!is.na(cs)) > 3)) {
		
		r <- cor(as.numeric(data), as.numeric(rad_cs), use = "pairwise.complete.obs")
		if (!is.na(r)) {
			data_cor <- cor.test(as.numeric(data), as.numeric(rad_cs), method = "pearson", alternative = "two.sided")$p.value
			if (data_cor > 0.05) {
			
				flag_vec <- rep(1,length(data))
			
			} else {
		
				flag_vec <- integer(length(data))
			
			}
		} else {
	
			flag_vec <- integer(length(data))
		
		}
		
	}
	
	return(flag_vec)

}

#Bertolotto et al. 2018
clim_snow <- function(data, data2) {

	flag_vec <- integer(length(data))
	
	flag_vec[((data < 0) & (data2 >= 0) & !is.na(data) & !is.na(data2)) | ((data >= 0) & (data2 < 0) & !is.na(data) & !is.na(data2))] <- 1
	
	return(flag_vec)
	
}

snowdrift <- function(data, data_t, data_p, data_wind, data_hn, th_t = 0, th_p = 0, th_wind = 7, th_snow = 5) {

	flag_vec <- integer(length(data))
	delta_t <- shift(data_t, type = "lag")
	delta_snow <- data - shift(data, type = "lag")	
		
	flag_vec[((delta_t <= th_t) & !is.na(delta_t)) & 
			 (((data_p == th_p) & !is.na(data_p)) | ((data_p != th_p) & (data_hn == 0) & !is.na(data_p) & (!is.na(data_hn)))) & 
			 ((data_wind > th_wind) & !is.na(data_wind)) & !is.na(data) & (abs(delta_snow) > th_snow)] <- 1
			 
	return(flag_vec)
	
}

diag_fill <- function(flag_mat, datetime, test_vec, variable, code, flag = 1) {

	flag_sum <- apply(flag_mat, 1, sum, na.rm = TRUE)
	
	if (length(flag_sum[flag_sum != 0]) == 0) {
	
		return(NULL)
		
	} else {
	
		test_mat <- matrix("", nrow = length(which(flag_sum != 0)), ncol = length(test_vec))
		for (i in 1:length(test_vec)) {
	
			test_mat[, i] <- ifelse(flag_mat[flag_sum != 0, i] == flag, "yes", "no")
		
		}
	
		diag_df <- cbind(data.frame(Code = code, Date = datetime[flag_sum != 0]), 
						test_mat, 
						data.frame(Variable = variable)
						) %>% as.data.frame(.)
	
		colnames(diag_df) <- c("Code", "Date", test_vec, "Variable")
	
		return(diag_df)
		
	}
	
}