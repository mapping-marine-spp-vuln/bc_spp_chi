### Support functions for this project

valid_check <- function(spp_shp) {
  valid <- st_is_valid(spp_shp)
  ### can return a vector if multiple polygons with same ID
  
  if(all(valid)) {
    ### all good!
    return(spp_shp)
  } 
  
  ### at least one invalid geometry to deal with

  id <- unique(spp_shp$iucn_sid)
  message('  Found invalid geometries in ', id)

  ### first try st_make_valid?
  message('  trying to fix ', id, ' using st_make_valid()')
  spp_shp1 <- st_make_valid(spp_shp)
  if(all(st_is_valid(spp_shp1))) {
    message('  ...fixed invalid geoms in ', id, ' using st_make_valid()')
    return(spp_shp1)
  }
  
  ### idea from:
  ### https://stackoverflow.com/questions/68478179/how-to-resolve-spherical-geometry-failures-when-joining-spatial-data
  check_use_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE)) ### turn off spherical geometry
  
  ### try zero-distance buffer method
  message('  trying to fix ', id, ' using zero-distance buffer')
  suppressMessages(suppressWarnings({
    spp_shp2   <- st_buffer(spp_shp1, dist = 0) %>%
      st_cast('MULTIPOLYGON')
  }))
  
  if(all(st_is_valid(spp_shp2))) {
    ### check area before and after buffering to make sure no loss
    area_pre  <- st_area(spp_shp) %>% as.numeric() / 1e6
    area_post <- st_area(spp_shp2) %>% as.numeric() / 1e6
    
    area_check <- all.equal(area_pre, area_post)
    area_ratio <- max(sum(area_pre) / sum(area_post),
                      sum(area_post) / sum(area_pre))
    ### set sf_use_s2 back to prior value
    suppressMessages(sf::sf_use_s2(check_use_s2))
    
    
    ### error check to make sure the buffer didn't lose polygons
    if(area_check == FALSE | 
       (class(area_check) != 'logical' & area_ratio > 1.001)) {
      ### use all.equal() for near equality, and for comparing all 
      ### elements in case of a vector.  If a difference, choose an arbitrary
      ### threshold for "close enough".
      message('Error: ', id, ': area_pre = ', round(sum(area_pre), 3), 
              '; area_post = ', round(sum(area_post), 3), 
              '; area_ratio = ', round(area_ratio, 5), '; not equal!')
      stop('Area_pre and area_post not equal for ', id, '!')
    } else {
      message('  ...fixed invalid geoms in ', id, ' using zero-distance buffer!')
      return(spp_shp2)
    }
  } else {
    stop('  ...invalid geometry in ', id, ' not fixed with zero-distance buffer!')
  }
  ### set sf_use_s2 back to prior value
  suppressMessages(sf::sf_use_s2(check_use_s2))
  return(NULL)
}

fix_fieldnames <- function(poly_sf) {
  names(poly_sf)[names(poly_sf) %in% c('binomial', 'binomil')] <- 'sciname'
  names(poly_sf)[names(poly_sf) %in% c('id_no', 'sisid')] <- 'iucn_sid'
  names(poly_sf)[names(poly_sf) %in% c('presenc')] <- 'presence'
  
  if(!'subpop' %in% names(poly_sf)) {
    poly_sf$subpop <- NA_character_
    ### if shape doesn't have subpop column, add it as NA
  }

  if(!'presence' %in% names(poly_sf)) {
    poly_sf <- poly_sf %>%
      mutate(presence = 1)
  }
  return(poly_sf)
}

fix_turtle_polys <- function(shp) {
  message('Replacing TURTLES with buffered turtle polygons!!!')
  ### the turtle polys have issues with the western boundary - not quite at
  ### +180; some have issues on the east as well. Buffer problem spp by
  ### 0.15 degrees before clipping.  Near the equator this adds an error of
  ### ~ 15 km!  But clearly the shitty polygons are creating an error as well.
  ### Identify the problem ones:
  ### Caretta caretta 3897;          Dermochelys coriacea 6494; 
  ### Eretmochelys imbricata 8005;   Chelonia mydas 4615
  ### I'll just buffer all subpops, for ease.  Land gets masked out later.
  
  turtles_buffered_file <- here_anx('iucn_spp/turtles_shp_buffered.shp')
  
  if(!file.exists(turtles_buffered_file)) {
    message('...Creating a temporary buffered turtles shapefile...')
    polys_rept <- read_sf(shp, type = 6) %>%
      janitor::clean_names() %>%
      fix_fieldnames()
    
    polys_rept_buff <- polys_rept %>%
      filter(iucn_sid %in% c(3897, 6494, 8005, 4615)) # %>%
      # st_buffer(dist = 0.25) 
    polys_rept_non_buff <- polys_rept %>%
      filter(!iucn_sid %in% c(3897, 6494, 8005, 4615))
    polys_rept_fixed <- rbind(polys_rept_non_buff, polys_rept_buff)
    st_write(polys_rept_fixed, turtles_buffered_file, delete_layer = TRUE)
  }
  
  ### now read in as polys_all the fixed buffered file
  polys_all <- st_read(turtles_buffered_file)
  return(polys_all)
}

match_to_map <- function(poly_sf, maps_df) {
  ### Fix the shapefile IUCN id vs. subpop IUCN id.  The inner_join
  ### keeps only the polygon features of species still to be
  ### rasterized (from the id_fix dataframe).
  id_fix <- maps_df %>%
    filter(shp_file == shp) %>%
    select(shp_iucn_sid, iucn_sid, subpop, depth_zone, sciname) %>%
    distinct()
  
  polys_match <- poly_sf %>%
    select(shp_iucn_sid = iucn_sid, shp_sciname = sciname, subpop, presence, geometry) %>%
    mutate(presence = ifelse(presence == 0, 1, presence),
           subpop   = as.character(subpop)) %>%
    inner_join(id_fix, by = c('shp_iucn_sid', 'subpop')) 
  
  return(polys_match)
}

buffer_tiny_polys <- function(spp_shp) {
  ### Check that CRS is in meters.
  # poly_crs <- crs(spp_shp) %>%
  #   as.character()
  # if(!str_detect(poly_crs, '\\+units=m')) {
  #   stop('For buffer_tiny_polys(), expecting units in meters')
  # }
  id <- unique(spp_shp$iucn_sid)
  
  ### Identify small polygons by their total area.  With 1 km cells, look
  ### for polygons a little larger than 1e6 m^2
  thresh_m2 <- units::set_units(1.30 * 1e6, m^2)
  ### Cell resolution is 1000 meters.  A buffer of ~half
  ### that would nearly guarantee all polygons to result in at least
  ### one raster cell, unless the polygon is near a corner.
  poly_area_m2 <- st_area(spp_shp$geometry)
  buffer_dist <- 500 ### in meters
  
  if(any(poly_area_m2 < thresh_m2)) { ### first just check across all polygons
    message('  ...found tiny polygons for ', id)
    spp_shp_buffered <- spp_shp
    for(z in 1:nrow(spp_shp)) {
      ### z <- 1
      poly_area_z <- poly_area_m2[z]
      if(poly_area_z < thresh_m2) { ### then check poly by poly}
        msg_stem <- '  Tiny polygon in %s (%s km^2 < %s km2 threshold); buffering by %s km...'
        message(sprintf(msg_stem, id, round(poly_area_z/1e6), thresh_m2/1e6, buffer_dist/1e3))
        spp_shp_buffered[z, ] <- spp_shp[z, ] %>%
          st_buffer(dist = buffer_dist)
      }
    }
    return(spp_shp_buffered %>% st_cast('MULTIPOLYGON'))
  } else {
    # message('  ...all good - no tiny polygons for ', id)
    return(spp_shp)
  }
}

clip_to_depth <- function(spp_rast, depth_zone) {
  ### depth clip if necessary; otherwise clip to bathy raster (which previously
  ### was clipped to area raster - so cells with any marine area will be kept,
  ### and non-marine cells will be dropped).
  
  if(length(depth_zone) != 1) stop('Non-unique depth_zone field!')
  ### this shouldn't happen - each spp should have only one depth
  
  if(depth_zone == '< 20 m') {
    ### intertidal, very shallow spp
    spp_rast <- mask(spp_rast, rast_shallow)
  } else if(depth_zone == '< 200 m') {
    ### spp on the continental shelf
    spp_rast <- mask(spp_rast, rast_neritic)
  } else {
    ### rast_bathy covers the entire ocean - effectively masks out land
    spp_rast <- mask(spp_rast, rast_bathy)
  }
  
  return(spp_rast)
}

