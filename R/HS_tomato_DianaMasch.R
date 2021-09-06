# load required libraries ------------------------------------------------------
library(tidyverse)
library(rcdk)
library(MetaboCoreUtils)
library(plotly)
library(tibble)
library(ggplot2)

#######################################################
homo_tibble <- tibble::as_tibble(homol[[1]]) 

homo_tibble <- homo_tibble %>%
  mutate(hsid=homo_tibble$`HS IDs`) %>%
  mutate(mzsplit=homo_tibble$'m/z increment') %>%
  tidyr::separate_rows(c('hsid', 'mzsplit'), sep="/") %>%
  mutate(hsid=as.factor(hsid)) 

summary(homo_tibble)
homol_filter <- homo_tibble %>%
  filter(homo_tibble$"HS IDs" != "0")

# grouping mz increments for plot legend/interactive plotly 
homol_filter$mzsplit <- as.numeric(homol_filter$mzsplit)
homol_filter$delta_mz = NA
homol_filter$delta_mz[homol_filter$mzsplit < 14] = "x <14"
homol_filter$delta_mz[homol_filter$mzsplit < 15 & homol_filter$mzsplit > 14] = ">14 x <15"
homol_filter$delta_mz[homol_filter$mzsplit < 18 & homol_filter$mzsplit > 15] = ">15 x <18"
homol_filter$delta_mz[homol_filter$mzsplit < 21 & homol_filter$mzsplit > 18] = ">18 x <21"
homol_filter$delta_mz[homol_filter$mzsplit < 30 & homol_filter$mzsplit > 21] = ">21 x <30"
homol_filter$delta_mz[homol_filter$mzsplit < 32 & homol_filter$mzsplit > 30] = ">30 x <32"
homol_filter$delta_mz[homol_filter$mzsplit < 60 & homol_filter$mzsplit > 32] = ">32 x <60"
# homol_filter$delta_mz[homol_filter$mzsplit < 90 & homol_filter$mzsplit > 60] = ">60 x <90"
# homol_filter$delta_mz[homol_filter$mzsplit < 120 & homol_filter$mzsplit > 90] = ">90 x <120"
# homol_filter$delta_mz[homol_filter$mzsplit < 150 & homol_filter$mzsplit > 120] = ">120 x <150"
# homol_filter$delta_mz[homol_filter$mzsplit < 200 & homol_filter$mzsplit > 150] = ">150 x <200"

neworderLegend = c("x <14", ">14 x <15", ">15 x <18", ">18 x <21", ">21 x <30", ">30 x <32", ">32 x <60")
homol_filter$delta_mz <- factor(homol_filter$delta_mz, levels=neworderLegend)

homo_plot <- ggplot(data=homol_filter, aes(x = RT, y = mz, text = paste(
  "m/z increment ", homol_filter$`m/z increment`, "\n",
  "HS ID ", homol_filter$`HS IDs`, "\n",
  sep = ""))) +
  labs( title = "homologue series; subset" )+
  geom_line(data=homol_filter, aes(x=RT, y=mz,group = hsid, color=delta_mz, alpha = 1), lwd=1)+ #x=mz, y=rt
  geom_point(alpha = 1, color = "darkgrey", size = 0.7)+
  #theme(legend.title = element_text())+
  theme_bw()

ggplotly_homo_plot <- ggplotly(homo_plot) %>% plotly::add_annotations( text="m/z increment", xref="paper", yref="paper",
                                                               x=1.02, xanchor="left",
                                                               y=0.8, yanchor="bottom",    # Same y as legend below
                                                               legendtitle=TRUE, showarrow=FALSE ) %>%
                                              plotly::layout( legend=list(y=0.8, yanchor="top" ) )
ggplotly_homo_plot 
homo_plot

#############################################################
