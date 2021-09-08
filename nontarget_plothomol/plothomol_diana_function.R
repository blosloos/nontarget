###############################################################################

library(plotly)
library(tibble)
library(ggplot2)

###############################################################################
plothomol_Diana <- function(homol) {
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
  
  ggplotly_homo_plot }
plothomol_Diana(homol)