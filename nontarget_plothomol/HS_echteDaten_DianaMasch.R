#echte Daten von Michael
library(readxl)
#read excel sheet with corresponding compounds
forhomol_Michael <- read_excel("D:/HomologeReihen/homologereihen/HS_search_simulateLCMS_DianaMasch/Input_R/Kopie von forhomol.xlsx",
                           )
forhomol <- forhomol_Michael %>%
  select(
    `Mass`,
    `m/z`,
    `RT`,
    `Adduct`,
    `Row Names`,`TandemLCpart`,
  ) %>%
  filter(TandemLCpart == "RP"
    )
  )

forhomol_HS <- forhomol[, c(2,3)]

# add a dummy coloumn with intensities, needed for nontarget::homol.search() 
forhomol_HS$intensity <- '999'
forhomol_HS <- forhomol_HS %>% rename("m/z" = "mz")

# cmpds_Dino %>% 
#   filter(between(mz, 50, 1500)) %>% 
#   sample_frac(size = 0.1) %>% 
#   ggplot(aes(x = logp, y = mz)) +
#   geom_point()

# organizing
forhomol_HS$intensity <- as.numeric(forhomol_HS$intensity)
forhomol_HS <- forhomol_HS %>% 
  na.exclude()
head(forhomol_HS)
forhomol_HS$mz <- as.numeric(forhomol_HS$mz)
forhomol_HS$RT <- as.numeric(forhomol_HS$RT)
forhomol_HS$intensity <- as.numeric(forhomol_HS$intensity)
forhomol_HS <- as.data.frame(forhomol_HS)
col_order <- c("mz", "intensity","RT" )
forhomol_HS<- forhomol_HS[, col_order]

# forhomol_HS <- forhomol_HS[1:500,]

homol_echt <- homol.search(forhomol_HS,
                      isotopes,	elements = c("C", "H", "O"), use_C = T, # 
                      minmz=10, 	maxmz=35,
                      minrt=0.1,  maxrt=1, #0
                      ppm=T,
                      mztol=1,  rttol=0.1, # 0.5
                      minlength=3, #5
                      mzfilter=F,
                      spar=.45, 	R2=.98,
                      plotit=T, mat_size=3) # do the thing...


#plothomol(homol_echt)
plothomol_Diana(homol_echt)
# #randomized m/z
# forhomol_HS$mz <- sample(1200, size = nrow(forhomol_HS), replace = TRUE)
# #permuate
# permute()
# head(forhomol_HS)
# forhomol_HS_permute <- gtools::permute(forhomol_HS$RT)
# head(forhomol_HS_permute)
# forhomol_HS$RT <- forhomol_HS_permute
# homo_filter <- homol_echt[[1]] %>%
#   filter(homol_echt[[1]]$"HS IDs" != "0")
# 
# merged_HS <- merge(homo_filter, forhomol_Michael, by.x = c('mz'), by.y= "m/z")
# 
# write.csv(merged_HS, file ="homol_echt.csv")
