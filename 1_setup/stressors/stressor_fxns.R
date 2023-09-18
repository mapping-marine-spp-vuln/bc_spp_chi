### Functions for dealing with stressor rasters


log_plus <- function(x, add_frac = 100) {
  ### Log(x+1) alternative that adds a small portion to the values, so zeros
  ### will be converted to a negative value.  After transforming, the log of  
  ### the added value is subtracted from all values, resulting in values of 
  ### zero returning to zero, and resulting in consistent scale-independent
  ### results (as long as add_frac remains the same).
  ### add_frac is the fraction of the max value to be added. Adding 1/100th
  ### of the max value makes a reasonable curve; 1/1000th makes a very steep
  ### for low values then quickly levels off.
  if(str_detect(tolower(class(x)), 'raster')) {
    
    x_rng <- minmax(x)
    if(x_rng[1] < 0) stop('Raster to be transformed cannot have negative values!')
    add <- x_rng[2] / add_frac
    y <- log(x + add)
    y <- y - log(add)
  } else {
    if(min(x, na.rm = TRUE) < 0) stop('Vector to be transformed cannot have negative values!')
    add <- max(x, na.rm = TRUE) / add_frac
    # if(add > 1) add <- 1
    y <- log(x + add)
    y <- y - log(add)
  }
  return(y)
}
### For rescaling a stressor to range from 0 to 1, based on min and max, with
### option for using a particular quantile as ref point (default 1.0, i.e., max)
rescale_stressor <- function(x, qtile = 1.0) {
  
  if(str_detect(tolower(class(x)), 'raster')) {
    y <- x ### temp raster to modify
    x_rng <- minmax(x)
    if(is.null(qtile)) {
      ref = x_rng[2]
    } else {
      ref = quantile(values(x), qtile, na.rm = TRUE)
    }
    values(y) <- (values(x) - x_rng[1]) / (ref - x_rng[1])
    ### in case of quantile, clip results to 1
    values(y)[values(y) > 1] <- 1
    return(y)    
  } else {
    stop('Stressor to be rescaled must be a raster or spatraster layer!')
  }
}

check_missing_cells <- function(x, warn = TRUE) {
  ### check raster x against the ocean area raster to make sure all cells have
  ### a non-NA value
  ocean_r <- rast(here('_spatial/ocean_bc_1km.tif'))
  
  tmp <- mask(ocean_r, x, inverse = TRUE)
  
  if(sum(!is.na(values(tmp))) > 0 & warn) warning('Raster contains NAs for some ocean cells!')
  return(tmp)
}

focal_gapfill <- function(x, iters = 1, method = 'mean') {
  ### use focal() to gapfill NAs using the mean of nearby cells
  w  <- matrix(1, 3, 3) #c(0,1,0,1,0,1,0,1,0), nrow=3)
  w2 <- matrix(1, 5, 5)
  w3 <- matrix(1, 7, 7)
  y <- x ### create copy of x
  for(i in 1:iters) {
    ### iterate over y, focal each time
    message('focal iteration ', i)
    if(i > 10) w <- w2 ### increase window for high iterations
    if(i > 25) w <- w3
    y <- focal(y, w, mean, na.rm = TRUE, NAonly = TRUE, pad = TRUE)
  }
  ### mask the result and return
  ocean_r <- rast(here('_spatial/ocean_bc_1km.tif'))
  y_masked <- mask(y, ocean_r)
  return(y_masked)
}
