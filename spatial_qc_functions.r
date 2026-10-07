anomaly <- function(data, datetime, win_clim = 15, zeros = FALSE, wind = FALSE) {

	data_df <- data.frame(t = datetime, temp = data)
	if (zeros) {
	
		data_df$temp[data_df$temp == 0] <- NaN
		
	}
	
	clim <- ts2clm(data_df, climatologyPeriod = c(data_df$t[1], data_df$t[nrow(data_df)]), windowHalfWidth = win_clim)
	clim <- clim[!duplicated(clim$t), ]
	if (wind) {
		
		anom <- data / clim$seas[!is.na(clim$t)]
	
	} else {
		
		anom <- data - clim$seas[!is.na(clim$t)]
	
	}
	
	return(anom)
	
}

corroboration_test <- function(data, data_ref, datetime, win_clim = 7, time_win = 3, th, th_rel = 50, th_na = 3, zeros = FALSE) {

	flag_vec <- integer(length(data))
	anom <- anomaly(data, datetime, win_clim, zeros)
	anom_ref <- apply(data_ref, 2, function(x) {anomaly(x, datetime, win_clim, zeros)})
	
	if (zeros) {
	
		data_ma <- embed(rbind(matrix(NA, nrow = 1, ncol = ncol(data_ref)), anom_ref, matrix(NA, nrow = 1, ncol = ncol(data_ref))), time_win)
		data_out <- apply(data_ma, 2, function(x) {abs(anom - x)})
		data_rel <- apply(data_ma, 2, function(x) {abs(100 * (anom - x) / anom)})
		na_out <- apply(data_out, 1, function(x) sum(!is.na(x)))
		na_rel <- apply(data_rel, 1, function(x) sum(!is.na(x)))
		th_out <- apply(data_out, 1, function(x) sum(x < th, na.rm = TRUE))
		th_rel <- apply(data_rel, 1, function(x) sum(x < th_rel, na.rm = TRUE))
		
		flag_vec[!is.na(data) & (th_out == 0) & (na_out >= th_na) & (th_rel == 0) & (na_rel >= th_na)] <- 1
	
	} else {
	
		data_out <- embed(rbind(matrix(NA, nrow = 1, ncol = ncol(data_ref)), anom_ref, matrix(NA, nrow = 1, ncol = ncol(data_ref))), time_win) %>% 
						  apply(., 2, function(x) {abs(anom - x)})
						  
		na_out <- apply(data_out, 1, function(x) sum(!is.na(x)))
		th_out <- apply(data_out, 1, function(x) sum(x < th, na.rm = TRUE))
		
		flag_vec[!is.na(data) & (th_out == 0) & (na_out >= th_na)] <- 1
		
	}
	
	return(flag_vec)
	
}

corr_fact <- function(data, data_ref, cf) {

	if (cf == "ratio") {
	
		cf_out <- mean(data, na.rm = TRUE) / mean(data_ref[!is.na(data)], na.rm = TRUE)
		
	} else {
	
		cf_out <- mean(data, na.rm = TRUE) - mean(data_ref[!is.na(data)], na.rm = TRUE)
		
	}
	
	return(cf_out)
	
}

gap_fill <- function(data, data_ref, nas, tau_cor, var_type) {


	cor_data <- apply(data_ref, 2, function(x) {cor(x, data, method = "spearman", use = "pairwise.complete.obs")})
	weight <- exp(-(1 - cor_data^2) / (tau_cor^2/log(2)))
	cf <- apply(data_ref, 2, function(x) {corr_fact(data, x, var_type)})
	
	if (length(nas) > 1) {
	
		if (var_type == "ratio") {
		
			fill_value <- apply(data_ref[nas, ], 1, function(x) weighted.mean(x * cf, weight, na.rm = TRUE))
			
		} else {
		
			fill_value <- apply(data_ref[nas, ], 1, function(x) weighted.mean(x + cf, weight, na.rm = TRUE))
			
		}
		
	} else {
	
		if (var_type == "ratio") {
		
			fill_value <- weighted.mean(data_ref[nas, ] * cf, weight, na.rm = TRUE)
			
		} else {
		
			fill_value <- weighted.mean(data_ref[nas, ] + cf, weight, na.rm = TRUE)
			
		}
		
	}
	
	data[nas] <- round(fill_value, 1)
	
	return(data)

}

anomaly_re_test <- function(data, data_ref, datetime, th, th_rel = 50, th_nan = 3, win_clim = 15, tau_cor = 0.3, var_type, zeros = FALSE) {

	flag_vec <- integer(length(data))
	nas <- which(is.na(data))
	
	if (length(nas) != 0) {
	
		data_test <- gap_fill(data, data_ref, nas, tau_cor, var_type)
		
	} else {
	
		data_test <- data
		
	}
	
	anom <- anomaly(data, datetime, win_clim, zeros)
	anom_ref <- apply(data_ref, 2, function(x) {anomaly(x, datetime, win_clim, zeros)})
	data_nan <- apply(anom_ref, 1, function(x) sum(!is.na(x), na.rm = TRUE))
	data_new <- apply(data_test + anom_ref - anom, 1, mean, na.rm = TRUE)
	
	if (zeros) {
	
		flag_vec[(abs(100 * (data - data_new) / data) >= th_rel) & ((data - data_new) >= th) & (data_nan >= th_nan) & !is.na(data)] <- 1
		
	} else {
	
		flag_vec[(abs(data - data_new) >= th) & (data_nan >= th_nan)] <- 1
		
	}
	
	return(flag_vec)
	
}

wet_dry_test <- function(data, data_ref, datetime, test, fwd, dmin, ref_wd, cov_min = 5, cov_max = 9, fmin = 2.7, fmin_cov = 3.2, wd_th = 20000, dth = 15000, wet_th = 0.3, dry_th = 0.3, th_test = 0, th_nan = 3) {

	flag_vec <- integer(length(data))
	cvs <- which((month(datetime) >= cov_min) & (month(datetime) <= cov_max))
	
	if (test == "wet") {
	
		pth <- fwd + (fmin * dmin / 150000)
		data_tmp <-  embed(rbind(matrix(NA, nrow = 1, ncol = ncol(data_ref)), data_ref, matrix(NA, nrow = 1, ncol = ncol(data_ref))) ,3)
		wet_data <- apply(data_tmp, 1, function(x) sum(x >= wet_th, na.rm = TRUE))
		data_valid <- apply(data_ref, 1, function(x) sum(!is.na(x), na.rm = TRUE))
		flag_wet <- which((wet_data == th_test) & (data_valid >= th_nan))
		flag_vec[setdiff(flag_wet[(data[flag_wet] > pth) & !is.na(data[flag_wet])], cvs)] <- 1
		
		if (length(ref_wd[ref_wd <= 20000]) > 1) {
		
			pth <- fwd + (fmin_cov * dmin / 150000)
			data_tmp <-  embed(rbind(matrix(NA, nrow = 1, ncol = length(ref_wd[ref_wd <= wd_th])), data_ref[, ref_wd <= wd_th], matrix(NA, nrow = 1, ncol = length(ref_wd[ref_wd <= wd_th]))) ,3)
			wet_data <- apply(data_tmp, 1, function(x) sum(x >= wet_th, na.rm = TRUE))
			data_valid <- apply(data_ref[, ref_wd <= wd_th], 1, function(x) sum(!is.na(x), na.rm = TRUE))
			flag_wet <- which((wet_data == th_test) & (data_valid >= th_nan))
			flag_vec[setdiff(flag_wet[(data[flag_wet] > pth) & !is.na(data[flag_wet])], cvs)] <- 1
			
		}
		
	}
	
	if (test == "dry") {
	
		pth <- fwd + (fmin * dmin / 150000)
		st_test <- which((data < dry_th) & !is.na(data))
		if (length(st_test) != 0) {
		
			data_tmp <-  embed(rbind(matrix(NA, nrow = 1, ncol = ncol(data_ref)), data_ref, matrix(NA, nrow = 1, ncol = ncol(data_ref))) ,3)
			ref_test <- apply(data_tmp[st_test, ], 2, function(x) x - pth)
			data_test <- apply(data_ref[st_test, ], 2, function(x) x - pth)
			dry_data <- apply(ref_test, 1, function(x) sum(x <= th_test, na.rm = TRUE))
			data_valid <- apply(data_test, 1, function(x) sum(!is.na(x), na.rm = TRUE))
			flag_dry <- which((dry_data == th_test) & (data_valid >= th_nan))
			flag_vec[setdiff(st_test[flag_dry][!is.na(data[st_test[flag_dry]])], cvs)] <- 1
		}
		
		if (length(ref_wd[ref_wd <= 20000]) > 1) {
		
			pth <- fwd + (fmin_cov * dmin / 150000)
			data_tmp <-  embed(rbind(matrix(NA, nrow = 1, ncol = length(ref_wd[ref_wd <= wd_th])), data_ref[, ref_wd <= wd_th], matrix(NA, nrow = 1, ncol = length(ref_wd[ref_wd <= wd_th]))) ,3)
			ref_test <- apply(data_tmp[st_test, ], 2, function(x) x - pth)
			data_test <- apply(data_ref[st_test, ref_wd <= wd_th], 2, function(x) x - pth)
			dry_data <- apply(ref_test, 1, function(x) sum(x <= th_test, na.rm = TRUE))
			data_valid <- apply(data_test, 1, function(x) sum(!is.na(x), na.rm = TRUE))
			flag_dry <- which((dry_data == th_test) & (data_valid >= th_nan))
			flag_vec[setdiff(st_test[flag_dry][!is.na(data[st_test[flag_dry]])], cvs)] <- 1
		}
		
	}
	
	return(flag_vec)
	
}

snow_test <- function(data, data_t, datetime, data_ref, win = 7, hs_step = 50, rel_th = 60, snow_th = 30, th = 0, th_na = 3, t_th = 2.4) {

	flag_vec <- integer(length(data))
	
	flag_tmp <- corroboration_test(data, data_ref, datetime, th = snow_th, th_rel = rel_th, zeros = TRUE, win_clim = win)
	
	snow_bef <- data - shift(data, type = "lag")
	snow_aft <- shift(data, type = "lead") - data
	
	data_rel <- apply(abs(100 * (data_ref - data) / data), 1, function(x) sum(x < rel_th, na.rm = TRUE))
	data_valid <- apply(data_ref, 1, function(x) sum(!is.na(x), na.rm = TRUE))
	flag_vec[(flag_tmp == 1) & (flag_vec == 0) & (data_rel == th) & (data_valid == th_na)] <- 1
	flag_vec[(flag_vec == 1) & (snow_bef <= 0) & (snow_aft <= 0) & (data_t > t_th) & !is.na(data_t)] <- 0
	flag_vec[(flag_vec == 1) & (sign(snow_bef) == sign(snow_aft))] <- 0
	flag_vec[(flag_vec == 1) & (sign(snow_bef) != sign(snow_aft)) & (data_t < t_th) & !is.na(data_t)] <- 0
	
	flag_vec[(flag_tmp == 1) & ((snow_bef > hs_step) | (snow_aft > hs_step)) & (data_t > t_th) & !is.na(data_t)] <- 1
	flag_vec[(flag_tmp == 1) & (abs(snow_bef) > hs_step) & (abs(snow_aft) > hs_step) & (sign(snow_bef) != sign(snow_aft))] <- 1
	flag_vec[(flag_tmp == 1) & ((snow_bef < -hs_step) | (snow_aft < -hs_step)) & (data_t <= t_th) & !is.na(data_t)] <- 1
	
	return(flag_vec)
	
}

wind_spatial_test <- function(data, datetime, data_ref, win_clim = 7, th_ma = 3.7, th_ds = 7.5, th_cv = 7.5, ma_low = 0.5, ds_low = 0.044, cv_low = 0.22, th = 0, calm_th = 1, zeros = FALSE, wind = TRUE) {

	flag_vec <- integer(length(data))
	anom <- anomaly(data, datetime, win_clim, zeros, wind)
	anom_ref <- apply(data_ref, 2, function(x) {anomaly(x, datetime, win_clim, zeros, wind)})
	
	ma <- rollapply(data = anom, width = 1+(win_clim*2), FUN = mean, na.rm = TRUE)
	ds <- rollapply(data = anom, width = 1+(win_clim*2), FUN = sd, na.rm = TRUE)
	cv <- ds / ma
	
	ref_ts <- apply(anom_ref, 1, mean, na.rm = TRUE)
	p95 <- quantile(ref_ts, 0.95, na.rm = TRUE)
	iqr <- IQR(ref_ts, na.rm = TRUE)
	
	flag_test <- matrix(0, nrow = length(data), ncol = 3)
	
	ma_up <- p95 + (th_ma * iqr)
	ds_up <- p95 + (th_ds * iqr)
	cv_up <- p95 + (th_cv * iqr)
	
	flag_test[(ma <= ma_up) & (ma >= ma_low), 1] <- 1
	flag_test[(ds <= ds_up) & (ds >= ds_low), 1] <- 1
	flag_test[(cv <= cv_up) & (cv >= cv_low), 1] <- 1
	
	flags <- apply(flag_test, 1, sum)
	
	flag_vec[!is.na(data) & (data > calm_th) & (flags == th)] <- 1
	
	return(flag_vec)

}

rad_ci <- function(data, data_ref, datetime, n = 3, err_th = 1.4826) {

	data_ref_tbl <- as.data.frame(data_ref - data)
	colnames(data_ref_tbl) <- glue("St{1:ncol(data_ref)}")
	data_ref_df <- data_ref_tbl %>% gather(., "Station", "Val")
	delta <- data.frame(Date = rep(datetime, ncol(data_ref)), Value = data_ref_df$Val, Station = data_ref_df$Station)

	bias_mon <- delta %>% mutate(Month = month(Date), Year = year(Date)) %>% 
						 group_by(Month, Year, Station) %>% 
						 summarise(Bias_m = median(Value, na.rm = TRUE))
						 
	bias <- bias_mon %>% group_by(Year) %>% summarise(Bias = median(Bias_m, na.rm = TRUE))
	bias <- bias[order(bias$Year), ]
	err_mad <- bias_mon %>% group_by(Year) %>% summarise(Mad = err_th * median(abs(Bias_m), na.rm = TRUE))
	err_mad <- err_mad[order(err_mad$Year), ]
	
	conf_int <- data.frame(Time = c(bias$Year), Low = bias$Bias - (n * err_mad$Mad), High = bias$Bias + (n * err_mad$Mad))
	
	return(conf_int)
	
}

rad_win <- function(delta, ci_df, datetime, index, w, d_th, valid_th, max_out) {

	flag <- integer(length(delta))
	d_start <- index
	d_end <- index + w - 1
	ci <- ci_df[ci_df$Year == year(datetime[index]), ]
	data_test <- delta[d_start:d_end]
	data_valid <- 100 * sum(data_test, na.rm = TRUE) / w
	data_over <- 100 * sum(data_test[(data_test  > ci$High) & (abs(data_test) > d_th)])
	data_under <- 100 * sum(data_test[(data_test < ci$Low) & (abs(data_test) > d_th)])
			
	if ((data_valid >= valid_th) & (data_over >= max_out) | (data_under >= max_out)) {
				
		flag[d_start:d_end] <- 1
					
	}
	
	return(flag)
	
}

rad_spatial_test <- function(data, datetime, data_ref, k = 5, w = 15, d_th = 5, valid_th = 80, max_out = 80) {

	conf_int <- rad_ci(data, data_ref, datetime)
	
	flag_vec <- integer(length(data))
	flag_mat <- matrix(0, nrow = length(data), ncol = ncol(data_ref))
	
	for (i in 1:ncol(data_ref)) {
	
		delta <- data_ref[, i] - data

		for (j in seq(1, length(data), by = 5)) {
	
			flag_mat[, i] <- rad_win(delta, conf_int, datetime, j, w, d_th, valid_th, max_out)
			
		}
		
	}
	
	flag <- apply(flag_mat, 1, function (x) sum(x == 0, na.rm = TRUE))
	flag_vec[(flag == 0) & !is.na(data)] <- 1
	
	return(flag_vec)
	
}

sp_spatial_test <- function(data, data_ref, datetime, z1, z2, th_p = 2, th = 0, th_na = 3) {

	flag_vec <- integer(nrow(data))
	ref_p <- matrix(NA, nrow = nrow(data), ncol = ncol(as.matrix(data_ref[, , 1])))
	for (i in 1:ncol(ref_p)) {
	
		ref_p[, i] <- hypso_test(data, data_ref[, i, ], z1, z2[i])
		
	}
	
	anom <- anomaly(data[, 1], datetime, win_clim = 15, zeros = FALSE)
	anom_ref <- apply(ref_p, 2, function(x) {anomaly(x, datetime, win_clim = 15, zeros = FALSE)})
	
	data_rel <- abs(anom - anom_ref)
	data_test <- apply(data_rel, 1, function(x) sum(x <= th_p, na.rm = TRUE))
	data_valid <- apply(data_rel, 1, function(x) sum(!is.na(x), na.rm = TRUE))
	
	flag_vec[(data_test == th) & (data_valid >= th_na) & !is.na(data[, 1])] <- 1
	
	return(flag_vec)
	
}