# load required libraries ------------------------------------------------------
library(tidyverse)
library(rcdk)
library(MetaboCoreUtils)
library(devtools)
library(nontarget)
data(adducts)
data(isotopes)
library(circlize)
library(devtools)
library(gdata)
library(jsonlite)
library(metaMS)
library(plotrix)
library(qdapRegex)
library(rJava)
library(rcdk)
library(rinchi)
library(rromeo)
library(readxl)
library(CAMERA)
library(RMassBank)
library(stringr)
library(stringi)
library(squash)
library(tools)
library(webchem)
library(zeallot)
library(BBmisc)
library(purrr)
library(schoolmath)
library(classyfireR)
library(sen2r)
library(sjmisc)
library(OrgMassSpecR)
library(ChemmineOB)
library(rJava)

#read excel sheet with corresponding compounds
Wormandtomato <- read_excel("D:/HomologeReihen/homologereihen/HS_search_simulateLCMS_DianaMasch/Input_R/Wormandtomato.xlsx",
                            sheet = "slycopersium", skip = 1)
# View(Wormandtomato)
# changing from MG to Mg, and FE, SE, MO manually in Excel sheet


# filtering NAs and only relevant coloumns

tomato <- Wormandtomato %>%
  select(
    `!ID`,
    `!Notes:Name_neutral`,
    `!Notes:Formula_Neutral`,
    `!Notes:InChI_neutral`,
    `!Notes:InChIKey_neutral`,
    `!Notes:ChEBI_neutral`
  ) %>%
  filter(
    !is.na(`!Notes:InChIKey_neutral`),!is.na(`!Notes:ChEBI_neutral`),
    !is.na('!Notes:Formula_Neutral')
  )

# database request for InChIkeys to SMILES -from anusha
gsmiles = c()

for (i in 1:length(tomato$'!Notes:InChIKey_neutral')) {
  tes <-
    tryCatch({
      webchem::get_cid(stringr::str_trim(as.character(tomato$`!Notes:InChIKey_neutral`[i])), from = "inchikey")
    }, error = function(cond) {
      message("webchem not able to get cid from Inchikey")
    })
  tes1 <-
    tryCatch({
      tes$cid
    }, error = function(cond) {
      message("Inchikey to CID did not convert")
    })
  tes2 <-
    tryCatch({
      webchem::pc_prop(
        as.numeric(tes1[1]),
        properties = c(
          "MolecularFormula",
          "MolecularWeight",
          "CanonicalSMILES",
          "InChI",
          "InChIKey"
        )
      )
    }, error = function(cond) {
      message("Inchikey to CID did not convert so did not get properties")
    })
  tes2IN <-
    tryCatch({
      tes2$InChI
    }, error = function(cond) {
      message("Inchi to Inchikey failed because of CID not converting-1")
    })
  In <- as.character(tomato$`!Notes:InChI_neutral` [i])
  SM <-
    tryCatch({
      tes2$CanonicalSMILES
    }, error = function(cond) {
      message("Inchi to Inchikey failed because of CID not converting-2")
    })
  if (!sjmisc::is_empty(SM)) {
    gsmiles = c(gsmiles, SM)
  } else if (!sjmisc::is_empty(IN)) {
    mol <-
      tryCatch({
        rinchi::parse.inchi(IN)
      }, error = function(cond) {
        message("name is empty")
      })
    SM <-
      tryCatch({
        rcdk::get.smiles(mol[[1]])
      }, error = function(cond) {
        message("name is empty")
      })
    gsmiles = c(gsmiles, SM)
  }
  else{
    gsmiles = c(gsmiles, "NA")
  }
}

inchi_tomato <- tomato %>%
  mutate (SMILES = gsmiles) %>%
  filter(.$SMILES != "NA")

# calculate mass and logp
tomato <- inchi_tomato %>%
  mutate(exactmass = unlist(lapply(.$`!Notes:Formula_Neutral`, function(x) {
    get.formula(x)@mass
  })),
  logp = unlist(lapply(.$SMILES, function(x) {
    mol <- parse.smiles(x)[[1]]
    get.xlogp(mol)
  })))


# calculate m/z values from exact mass, RT from RTI ladder ---------------------
# read rti ladder
rti <- read_tsv("rti.txt")

# calculate theoretical adduct masses
tomato_adducts <- bind_cols(tomato,
                            as_tibble(mass2mz(tomato$exactmass, c("[M+H]+", "[M+Na]+"))))



tomato_adducts_long <-
  pivot_longer(tomato_adducts, cols = -colnames(tomato)) %>%
  rename(adduct = name,
         mz = value) %>%
  mutate(RT = approx(rti$logP, rti$RT, logp)$y)

tomato_Dino <- tomato_adducts_long[, c(11, 12)]

# add a dummy coloumn with intensities, needed for nontarget::homol.search()
tomato_Dino$intensity <- '999'

tomato_Dino %>%
  filter(between(mz, 50, 1500)) %>%
  sample_frac(size = 0.1) %>%
  ggplot(aes(x = logp, y = mz)) +
  geom_point()


tomato_Dino$intensity <- as.numeric(tomato_Dino$intensity)
tomato_Dino <- tomato_Dino %>%
  na.exclude()
tomato_Dino$mz <- as.numeric(tomato_Dino$mz)
tomato_Dino$RT <- as.numeric(tomato_Dino$RT)
tomato_Dino$intensity <- as.numeric(tomato_Dino$intensity)
tomato_Dino <- as.data.frame(tomato_Dino)

col_order <- c("mz", "intensity", "RT")
tomato_Dino <- tomato_Dino[, col_order]

# searching for homologue series in the 3-coloumn-df
homol <- homol.search(
  tomato_Dino,
  isotopes,
  elements = c("C", "H", "O"),
  use_C = T,
  #
  minmz = 10,
  maxmz = 35,
  minrt = 0.1,
  maxrt = 1,
  #0
  ppm = T,
  mztol = 1,
  rttol = 0.1,
  # 0.5
  minlength = 3,
  #5
  mzfilter = F,
  spar = .45,
  R2 = .98,
  plotit = F
) # do the thing...


plothomol(homol)
#######################################################

library(plotly)
library(tibble)
library(ggplot2)

#######################################################
homo_tibble <- tibble::as_tibble(homol[[1]])

homo_tibble <- homo_tibble %>%
  mutate(hsid = homo_tibble$`HS IDs`) %>%
  mutate(mzsplit = homo_tibble$'m/z increment') %>%
  tidyr::separate_rows(c('hsid', 'mzsplit'), sep = "/") %>%
  mutate(hsid = as.factor(hsid))

summary(homo_tibble)
homol_filter <- homo_tibble %>%
  filter(homo_tibble$"HS IDs" != "0")
write.csv(
  homol_filter,
  file = "simulate_lcms_tomato_homolseries.txt",
  sep = "\t",
  row.names = T,
  col.names = T,
  dec = ","
)

# grouping mz increments for plot legend/interactive plotly
homol_filter$mzsplit <- as.numeric(homol_filter$mzsplit)
homol_filter$delta_mz = NA
homol_filter$delta_mz[homol_filter$mzsplit < 14] = "x <14"
homol_filter$delta_mz[homol_filter$mzsplit < 15 &
                        homol_filter$mzsplit > 14] = ">14 x <15"
homol_filter$delta_mz[homol_filter$mzsplit < 18 &
                        homol_filter$mzsplit > 15] = ">15 x <18"
homol_filter$delta_mz[homol_filter$mzsplit < 21 &
                        homol_filter$mzsplit > 18] = ">18 x <21"
homol_filter$delta_mz[homol_filter$mzsplit < 30 &
                        homol_filter$mzsplit > 21] = ">21 x <30"
homol_filter$delta_mz[homol_filter$mzsplit < 32 &
                        homol_filter$mzsplit > 30] = ">30 x <32"
homol_filter$delta_mz[homol_filter$mzsplit < 60 &
                        homol_filter$mzsplit > 32] = ">32 x <60"
# homol_filter$delta_mz[homol_filter$mzsplit < 90 & homol_filter$mzsplit > 60] = ">60 x <90"
# homol_filter$delta_mz[homol_filter$mzsplit < 120 & homol_filter$mzsplit > 90] = ">90 x <120"
# homol_filter$delta_mz[homol_filter$mzsplit < 150 & homol_filter$mzsplit > 120] = ">120 x <150"
# homol_filter$delta_mz[homol_filter$mzsplit < 200 & homol_filter$mzsplit > 150] = ">150 x <200"

neworderLegend = c("x <14",
                   ">14 x <15",
                   ">15 x <18",
                   ">18 x <21",
                   ">21 x <30",
                   ">30 x <32",
                   ">32 x <60")
homol_filter$delta_mz <-
  factor(homol_filter$delta_mz, levels = neworderLegend)

homo_plot <-
  ggplot(data = homol_filter, aes(
    x = RT,
    y = mz,
    text = paste(
      "m/z increment ",
      homol_filter$`m/z increment`,
      "\n",
      "HS ID ",
      homol_filter$`HS IDs`,
      "\n",
      sep = ""
    )
  )) +
  labs(title = "homologue series; subset") +
  geom_line(
    data = homol_filter,
    aes(
      x = RT,
      y = mz,
      group = hsid,
      color = delta_mz,
      alpha = 1
    ),
    lwd = 1
  ) + #x=mz, y=rt
  geom_point(alpha = 1,
             color = "darkgrey",
             size = 0.7) +
  #theme(legend.title = element_text())+
  theme_bw()

ggplotly_homo_plot <-
  ggplotly(homo_plot) %>% plotly::add_annotations(
    text = "m/z increment",
    xref = "paper",
    yref = "paper",
    x =
      1.02,
    xanchor = "left",
    y =
      0.8,
    yanchor = "bottom",
    # Same y as legend below
    legendtitle =
      TRUE,
    showarrow = FALSE
  ) %>%
  plotly::layout(legend = list(y = 0.8, yanchor = "top"))
ggplotly_homo_plot
homo_plot

#############################################################

#combining the ouput of homol.search with the annotations from the start
tomato_reunited <-
  merge(homol_filter,
        tomato_adducts_long,
        by.x = c('mz'),
        by.y = "mz")

tomato_annotated <- tomato_reunited %>%
  select(hsid, mzsplit, everything())

tomato_annotated$delta_mz <- as.character(tomato_annotated$delta_mz)

#write.csv(tomato_reunited, file ="tomato_reunited.csv")
# df <- tomato_annotated %>%
#   arrange(`!Notes:Name_neutral`, -`mzsplit`) %>%
#   filter(duplicated(`!Notes:Name_neutral`) == FALSE)
