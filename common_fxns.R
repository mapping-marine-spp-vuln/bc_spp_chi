### set options
options(readr.show_types = FALSE)

#####################################.
#### Helper functions in general ####
#####################################.
here_anx <- function(f = '', ...) { 
  ### create file path to git-annex dir for project
  f <- paste(f, ..., sep = '/')
  f <- stringr::str_replace_all(f, '\\/+', '/')
  f_anx <- sprintf('/home/shares/ohi/spp_vuln/bc_spp_chi/%s', f)
  return(f_anx)
}

smash2or <- function(x, sep = ' ') {
  ### collapse to OR string
  paste(x, sep = sep, collapse = '|')
}

### Helper functions for gathering species rangemaps
check_tryerror <- function(l) {
  x <- sapply(l, class) %>% 
    unlist() %>% as.vector()
  return(any(stringr::str_detect(tolower(x), 'error')))
}

get_one_map <- function(f) {
  if(file.exists(f)) {
    df <- data.table::fread(f)
    
    if('presence' %in% names(df)) {   ### IUCN spp map with presence col
      df <- df %>% filter(presence != 5) %>% select(-presence)
    }
    if('prob' %in% names(df)) {
      df <- df %>% filter(prob >= 0.5) %>% select(-prob)
    }
    return(df)
  } else {
    warning('No map found for ', basename(f))
    return(NULL)
  }
}

get_spp_info <- function() {
  data.table::fread(here('_data/bc_spp_vuln_info.csv'))
}

collect_spp_rangemaps <- function(spp_vec, file_vec, idcol = 'species', parallel = TRUE) {
  ### give a vector of species names (or IDs) and filenames; 
  ### default: read in using parallel::mclapply
  message('Collecting ', n_distinct(file_vec), ' maps for ', 
          n_distinct(spp_vec), ' species...')
  if(parallel == TRUE) {
    out_maps_list <- parallel::mclapply(file_vec, mc.cores = 40, FUN = get_one_map)
  } else {
    out_maps_list <- lapply(file_vec, FUN = get_one_map)
  }
  if(check_tryerror(out_maps_list)) {
    stop('Try-error found when assembling species rangemaps!')
  }
  message('... Binding maps...')
  out_maps_df <- out_maps_list %>%
    setNames(spp_vec) %>%
    purrr::compact() %>%
    data.table::rbindlist(idcol = idcol) %>%
    distinct()
  
  return(out_maps_df)
}


################################################.
####   Helper functions for spatial stuff   ####
################################################.

### Simple Features functions
create_bc_sf <- function() {
  bc_bbox <- st_bbox(c(xmin = -140, xmax = -114, ymax = 60, ymin = 46), crs = st_crs(4326))
  bc_bbox_sf <- st_as_sfc(bc_bbox)
  return(bc_bbox_sf)
}

clip_to_bc <- function(x) {
  ### for SF features, transform to wgs84, clip to +-180 and +-90
  id <- unique(x$iucn_sid)
  epsg <- st_crs(x)$epsg
  if(epsg != 4326 | is.na(epsg)) {
    message('  Original EPSG = ', epsg, '; Proj4 = ', st_crs(x)$proj4string,
            '\n...converting to EPSG:4326 WGS84 for clipping ', id)
    x <- st_transform(x, 4326)
  }
  
  ### check that bounding box is within BC bbox
  bc_bbox <- st_bbox(c(xmin = -140, xmax = -114, ymax = 60, ymin = 46), crs = st_crs(4326))
  bc_bbox_sf <- create_bc_sf()
  x_bbox <- st_bbox(x)
  if(x_bbox$xmin < bc_bbox$xmin | 
     x_bbox$xmax > bc_bbox$xmax |
     x_bbox$ymin < bc_bbox$ymin | 
     x_bbox$ymax > bc_bbox$ymax) {
    
    ### turn off spherical geometry temporarily!
    check_use_s2 <- sf::sf_use_s2()  ### store current value
    message('  Some bounds outside bc bounding box - clipping ', id)
    suppressMessages(sf::sf_use_s2(FALSE)) ### turn it off
    
    suppressMessages(suppressWarnings(z <- st_intersection(x, bc_bbox_sf)))
    
    suppressMessages(sf::sf_use_s2(check_use_s2)) ### turn spherical geom to orig value
    
    message('  Smoothing to half degree max')
    z <- smoothr::densify(z, max_distance = 0.5)
  } else {
    message('  All bounds OK, no clipping necessary for ', id)
    z <- x
  }
  return(z)
}

get_bc_rast <- function() {
  rast_base <- terra::rast(here('_spatial/ocean_bc_1km.tif')) %>%
    terra::setValues(1:terra::ncell(.))
  return(rast_base)
}

map_to_bc <- function(df, by = 'cell_id', which, xfm = NULL, ocean_mask = TRUE) {
  if(!by %in% names(df)) stop('Dataframe needs a valid column for "by" (e.g., cell_id)!')
  if(!which %in% names(df)) stop('Dataframe needs a valid column for "which" (e.g., n_spp)!')
  if(any(is.na(df[[by]]))) stop('Dataframe contains NA values for "by"!')
  if(length(df[[by]]) != length(unique(df[[by]]))) stop('Dataframe contains duplicate observations of ', by, '!')
  
  ### Instead of raster::subs (which is pretty slow), just make a 
  ### vector of all the new values by joining the dataframe to a
  ### new dataframe with complete cell_id column, then replace all 
  ### raster values at once, very fast!
  out_rast <- get_bc_rast()
  df1 <- data.frame(cell_id = 1:ncell(out_rast)) %>%
    dt_join(df, by = 'cell_id', type = 'left')
  values(out_rast) <- df1[[which]]
  
  ### this keeps the layer name from the base rast... swap with "which"
  names(out_rast) <- which

  if(!is.null(xfm)) {
    if(class(xfm) != 'function') stop('xfm argument must be a function!')
    out_rast <- xfm(out_rast)
  }
  
  if(ocean_mask) {
    out_rast <- out_rast %>%
      terra::mask(terra::rast(here('_spatial/ocean_bc_1km.tif')))
  }
  return(out_rast)
}


####################################################.
####         Pooled variance functions          ####
####################################################.

pooled_var <- function(x_bar, y_bar, s_x, s_y, n_x, n_y) {
  ### convert std dev to var
  var_x <- ifelse(is.na(s_x), 0, s_x^2)
  var_y <- ifelse(is.na(s_y), 0, s_y^2)
  
  var_xy_clean <- ((n_x - 1)*var_x + (n_y - 1)*var_y) / (n_x + n_y - 1)
  var_xy_error <- (n_x * n_y) * (x_bar - y_bar)^2 / ((n_x + n_y)*(n_x + n_y - 1))
  
  return(var_xy_clean + var_xy_error)
}

iterated_pooled_var <- function(mean_vec, sdev_vec, n_vec, flag = FALSE) {
  if(!all.equal(length(mean_vec), length(sdev_vec), length(n_vec))) {
    stop('Mean, std dev, and n vectors must all be equal length!')
  }
  if(length(mean_vec) == 1) {
    warning('Only one element - no need for pooled variance!')
    return(sdev_vec[1]^2)
  }
  ### initialize values for first in list
  mean_x <- mean_vec[1]; s_x <- sdev_vec[1]; n_x <- n_vec[1]
  for(i in 2:length(mean_vec)) { ## i <- 2
    if(flag) message('  ... processing iteration ', i - 1, '...')
    
    mean_y <- mean_vec[i]
    s_y    <- sdev_vec[i]
    n_y    <- n_vec[i]
    var_out <- pooled_var(x_bar = mean_x, y_bar = mean_y, 
                          n_x = n_x, n_y = n_y, 
                          s_x = s_x, s_y = s_y)
    
    ### set up values for next iteration
    mean_x <- (mean_x * n_x + mean_y * n_y) / (n_x + n_y)
    s_x <- sqrt(var_out)
    n_x <- n_x + n_y
  }
  return(var_out)
}

### test that the function returns the correct variance value
# set.seed(42)
# n_vec <- sample(1:20, size = 123, replace = TRUE)
# x_list <- list()
# for(x_i in seq_along(n_vec)) { ### x_i <- 1
#   x_list[[x_i]] <- rnorm(mean = sample(5 * c(1:10), size = 1),
#                        sd   = sample(.5 * c(1:10), size = 1),
#                        n = n_vec[x_i])
# }
# 
# ### initialize values for first term
# mean_vec <- sapply(x_list, mean)
# sdev_vec <- sapply(x_list, sd)
# n_vec    <- sapply(x_list, length)
# 
# var_out <- iterated_pooled_var(mean_vec, sdev_vec, n_vec)
# 
# var_check <- var(unlist(x_list))
# 
# abs(var_check - var_out) < 1e-10

####################################################.
####     Functional vulnerability functions     ####
####################################################.

calc_fe <- function(n_spp) {
  k <- n_spp - 1
  fv <- 0.5^k
} 

calc_spp_cell_fe <- function(spp_cells, spp_fe) {
  ### parallelize this across smaller chunks to keep group_by from crashing 
  ### everything - but not for every cell individually!  
  ### Set up 1000 different cell groups across the 100k(ish) cells in the chunk
  ### to divide work between dplyr and parallel...
  cell_id_df <- data.frame(cell_id = spp_cells$cell_id %>% unique()) %>%
    mutate(cell_gp = rep(1:250, length.out = n()))
  cell_gps <- cell_id_df$cell_gp %>% unique()
  
  fe_df <- spp_cells %>%
    oharac::dt_join(spp_fe %>% select(species, fe_id), 
                    by = 'species', type = 'left')
  
  ### use the number of observations to limit the number of cores... more
  ### observations, fewer cores!
  n_cores <- ceiling(25 / ceiling((nrow(fe_df) / 8e7)))
  
  fv_list <- parallel::mclapply(cell_gps, mc.cores = n_cores,
                                FUN = function(gp) { 
                                  ### gp <- 5
                                  cell_ids <- cell_id_df %>% filter(cell_gp == gp) %>% .$cell_id
                                  
                                  x <- fe_df %>%
                                    filter(cell_id %in% cell_ids) %>%
                                    data.table() %>% 
                                    ### data.table syntax for: group_by/mutate (see below)
                                    .[, ':='(n_spp = length(unique(species)),
                                             n_fe  = length(unique(fe_id))),
                                      by = .(cell_id)] %>%
                                    .[, ':='(n_spp_fe = length(unique(species))),
                                      by = .(cell_id, fe_id)] %>%
                                    .[, ':='(fv = calc_fe(n_spp_fe)),
                                      by = .(cell_id, fe_id)]
                                  
                                  return(x)
                                })
  if(check_tryerror(fv_list)) {
    stop('Try-error found when calculating species/cell functional vulnerability!')
  }
  fv_df <- data.table::rbindlist(fv_list)
  return(fv_df)
}