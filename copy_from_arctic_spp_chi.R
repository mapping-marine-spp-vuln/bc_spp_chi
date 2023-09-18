library(tidyverse)

### copy directory structure
w <- list.dirs(here::here('../arctic_spp_chi'), recursive = TRUE, full.names = TRUE)
w <- str_remove(w, 'bc_spp_chi/../')
w <- w[!str_detect(w, '\\.git|.Rproj')]

z <- w %>% str_replace('arctic_spp_chi', 'bc_spp_chi')
lapply(z, dir.create) ### dir.create takes a single path

### copy necessary files
x <- list.files(here::here('../arctic_spp_chi'), recursive = TRUE, full.names = TRUE)
x <- str_remove(x, 'bc_spp_chi/../')
x <- x[!str_detect(x, '\\.git|.Rproj|tmp/|int/|_output/')]

y <- x %>% str_replace('arctic_spp_chi', 'bc_spp_chi')

zxcv <- file.copy(x, y, recursive = TRUE)

### Now for the git-annex, just copy the directory structure
ww <- list.dirs('/home/shares/ohi/spp_vuln/arctic_spp_chi', recursive = TRUE, full.names = TRUE)
ww <- ww[!str_detect(ww, '\\.git|.Rproj')]

zz <- ww %>% str_replace('arctic_spp_chi', 'bc_spp_chi')
lapply(zz, dir.create) ### dir.create takes a single path
