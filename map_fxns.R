### load map layers

rgn_sf <- read_sf(here('_spatial/bc_eez/bc_eez_incl_inland.gpkg'))

land_sf <- read_sf(here('_spatial/bc_eez/bc_continent.gpkg')) %>%
  st_crop(rgn_sf)

meow_sf <- read_sf(here('_spatial/meow_rgns/meow_rgns.shp')) %>%
  janitor::clean_names() %>%
  st_transform(st_crs(rgn_sf)) %>%
  st_crop(rgn_sf) %>%
  mutate(ecoregion = str_replace(ecoregion, 'Pacific F', '\nPacific F'),
         ecoregion = str_replace(ecoregion, '/', '/\n'),
         ecoregion = str_replace(ecoregion, 'Vancouver Coast', 'Vancouver\nCoast'))

format_map <- function(p) {
  p <- p + geom_sf(data = meow_sf,  aes(geometry = geometry),  
                     color = 'grey80', fill = NA, size = .1) +
    geom_sf(data = rgn_sf, aes(geometry = geom), 
            color = 'yellow', fill = NA, alpha = .2, size = .5) +
    geom_sf(data = rgn_sf, aes(geometry = geom), 
            color = 'red', fill = NA, size = .2, linetype = 'dashed') +
    geom_sf(data = land_sf, aes(geometry = geom), 
            color = 'grey60', fill = 'grey90', size = .1) +
    geom_sf_text(data = meow_sf, aes(geometry = geometry, label = ecoregion),
                 size = 2.2, fontface = 'bold') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(panel.background = element_rect(color = NA, fill = 'grey70'))
  return(p)
}

str_clrs_df <- tribble(
  ~stressor,                             ~str_lbl,  ~str_fill, ~str_clr,
  'sst_rise',                          'SST rise',  '#4B0055', 'red',
  'sst_extremes',                  'SST extremes',  '#91689A', 'red',
  'ocean_acidification',    'Ocean acidification',  '#009B95', 'red',
  'sea_level_rise',              'Sea level rise',  '#94F8F1', 'red',
  'uv_radiation',                  'UV radiation',  '#5DC8C2', 'red',
  'shipping_large',                    'Shipping',  '#FDE333', NA,
  # 'benth_str',               'Benthic structures',  '#FFF9E3', NA,
  'biomass_removal',           'Targeted fishing',  '#00588B', NA,
  'bycatch',                    'Fishing bycatch',  '#C4DCFF', NA,
  'fishing_benthic_dest', 'Benthic dest. fishing',  '#6A98CA', NA,
  'direct_human',          'Direct human impacts',  '#53CC67', NA,
  'light',                      'Light pollution',  '#BFFFC6', NA,
  'nutrient',                'Nutrient pollution',  '#76EA87', NA
) %>%
  mutate(str_lbl = fct_inorder(str_lbl))

