# load required libraries ------------------------------------------------------
library(tidyverse)
library(rcdk)
library(MetaboCoreUtils)
library(devtools)
library(nontarget)
data(adducts)
data(isotopes)

# read compounds and select only required columns ------------------------------
cmpds <- read_tsv("D:/homologereihen/homologereihen/Compound-SBtab.tsv", skip = 1) %>% 
  select(`!ID`,
         `!Notes:Name_neutral`,
         `!Notes:Formula_Neutral`,
         `!Notes:InChI_neutral`,
         `!Notes:InChIKey_neutral`,
         `!Notes:SMILES_neutral`) %>% 
  filter(!is.na(`!Notes:InChIKey_neutral`),
         !is.na(`!Notes:SMILES_neutral`))

# calculate the exact mass and logP value --------------------------------------
cmpds <- cmpds %>% 
  mutate(exactmass = unlist(lapply(.$`!Notes:Formula_Neutral`, function(x) {
    get.formula(x)@mass})),
    logp = unlist(lapply(.$`!Notes:SMILES_neutral`, function(x) {
      mol <- parse.smiles(x)[[1]]
      get.xlogp(mol)}))
  )

# calculate m/z values from exact mass, RT from RTI ladder ---------------------
# read rti ladder
rti <- read_tsv("D:/homologereihen/homologereihen/HS_search_simulateLCMS_DianaMasch/Input_R/rti.txt")

# calculate theoretical adduct masses
cmpds_adducts <- bind_cols(cmpds,
                           as_tibble(mass2mz(cmpds$exactmass, c("[M+H]+", "[M+Na]+"))))


cmpds_adducts_long <- pivot_longer(cmpds_adducts, cols = -colnames(cmpds)) %>% 
  dplyr::rename(adduct = name,
         mz = value) %>% 
  mutate(RT = approx(rti$logP, rti$RT, logp)$y)

cmpds_adducts_short <- cmpds_adducts_long[,c(7,8,11)]


cmpds_adducts_long_ritbind <-  merge(rti,cmpds_adducts_short, all = T)
cmpds_Dino <- cmpds_adducts_long[, c(11,10)]

# add a dummy coloumn with intensities, needed for nontarget::homol.search() 
cmpds_adducts_long_ritbind$intensity <- '999'

# cmpds_Dino %>% 
#   filter(between(mz, 50, 1500)) %>% 
#   sample_frac(size = 0.1) %>% 
#   ggplot(aes(x = logp, y = mz)) +
#   geom_point()

# organizing
cmpds_adducts_long_ritbind$intensity <- as.numeric(cmpds_adducts_long_ritbind$intensity)
cmpds_adducts <- cmpds_adducts %>% 
  na.exclude()
head(cmpds_adducts)

cmpds_adducts_long_ritbind$mz <- as.numeric(cmpds_adducts_long_ritbind$mz)
cmpds_adducts_long_ritbind$RT <- as.numeric(cmpds_adducts_long_ritbind$RT)
cmpds_adducts_long_ritbind$intensity <- as.numeric(cmpds_adducts_long_ritbind$intensity)
cmpds_adducts <- as.data.frame(cmpds_adducts)
col_order <- c("mz", "intensity","RT" )
cmpds_adducts<- cmpds_adducts[, col_order]

cmpds_Dino <- cmpds_Dino[1:500,]
 
 
# searching for homologue series in the 3-coloumn-df
homol <- homol.search(cmpds_Dino,
                      isotopes,	elements = c("C", "H", "O"), use_C = T, # 
                      minmz=10, 	maxmz=35,
                      minrt=0.1,  maxrt=1, #0
                      ppm=T,
                      mztol=1,  rttol=0.1, # 0.5
                      minlength=3, #5
                      mzfilter=F,
                      spar=.45, 	R2=.98,
                      plotit=T, mat_size=3) # do the thing...


plothomol(homol)
# 
###############################################################################

library(plotly)
library(tibble)
library(ggplot2)

###############################################################################
# plothomol_Diana <- function(homol) {
  homo_tibble <- tibble::as_tibble(homol[[1]]) 

homo_tibble <- homo_tibble %>%
  mutate(hsid=homo_tibble$`HS IDs`) %>%
  mutate(mzsplit=homo_tibble$'m/z increment') %>%
  tidyr::separate_rows(c('hsid', 'mzsplit'), sep="/") %>%
  mutate(hsid=as.factor(hsid)) 

summary(homo_tibble)
homo_filter <- homo_tibble %>%
  filter(homo_tibble$"HS IDs" != "0")
# write.csv(homo_filter, file = "simulate_lcms_homolseries.txt",sep = "\t", row.names = T, col.names = T, dec = ",")

# grouping mz increments for plot legend/interactive plotly 
homo_filter$mzsplit <- as.numeric(homo_filter$mzsplit)
homo_filter$delta_mz = NA
homo_filter$delta_mz[homo_filter$mzsplit < 14] = "x <14"
homo_filter$delta_mz[homo_filter$mzsplit < 15 & homo_filter$mzsplit > 14] = ">14 x <15"
homo_filter$delta_mz[homo_filter$mzsplit < 18 & homo_filter$mzsplit > 15] = ">15 x <18"
homo_filter$delta_mz[homo_filter$mzsplit < 21 & homo_filter$mzsplit > 18] = ">18 x <21"
homo_filter$delta_mz[homo_filter$mzsplit < 30 & homo_filter$mzsplit > 21] = ">21 x <30"
homo_filter$delta_mz[homo_filter$mzsplit < 32 & homo_filter$mzsplit > 30] = ">30 x <32"
homo_filter$delta_mz[homo_filter$mzsplit < 60 & homo_filter$mzsplit > 32] = ">32 x <60"
homo_filter$delta_mz[homo_filter$mzsplit < 90 & homo_filter$mzsplit > 60] = ">60 x <90"
homo_filter$delta_mz[homo_filter$mzsplit < 120 & homo_filter$mzsplit > 90] = ">90 x <120"
homo_filter$delta_mz[homo_filter$mzsplit < 150 & homo_filter$mzsplit > 120] = ">120 x <150"
homo_filter$delta_mz[homo_filter$mzsplit < 200 & homo_filter$mzsplit > 150] = ">150 x <200"

neworderLegend = c("x <14", ">14 x <15", ">15 x <18", ">18 x <21", ">21 x <30", ">30 x <32", ">32 x <60")
homo_filter$delta_mz <- factor(homo_filter$delta_mz, levels=neworderLegend)

homo_plot <- ggplot(data=homo_filter, aes(x = mz, y = RT, text = paste(
  "m/z increment ", homo_filter$`m/z increment`, "\n",
  "HS ID ", homo_filter$`HS IDs`, "\n",
  sep = ""))) +
  labs( title = "homologue series; subset" )+
  geom_line(data=homo_filter, aes(x=RT, y=mz,group = hsid, color=delta_mz, alpha = 1), lwd=1)+ #x=mz, y=rt
  geom_point(alpha = 1, color = "darkgrey", size = 0.7)+
  #theme(legend.title = element_text())+
  theme_bw()

ggplotly_homo_plot <- ggplotly(homo_plot) %>% plotly::add_annotations( text="m/z increment", xref="paper", yref="paper",
                                                                       x=1.02, xanchor="left",
                                                                       y=0.8, yanchor="bottom",    # Same y as legend below
                                                                       legendtitle=TRUE, showarrow=FALSE ) %>%
  plotly::layout( legend=list(y=0.8, yanchor="top" ) )

ggplotly_homo_plot 
plothomol_Diana(homol)
###############################################################################
#combining the ouput of homol.search with the annotations from the start
cmpds_reunited <- merge(homo_filter, cmpds_adducts_long, by.x = c('mz'), by.y= "mz")

cmpds_annotated <- cmpds_reunited %>%
  select(hsid, mzsplit, everything())

cmpds_annotated$delta_mz <- as.character(cmpds_annotated$delta_mz)

write.csv(cmpds_annotated, file ="cmpds_annotated_subsetworm.csv")
