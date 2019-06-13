plot_gene_expression <- function(gene_name, tpm, sample_sheet) {
  gene_tpm <- tpm[which(row.names(tpm) == gene_name), ] %>%
    as.data.frame() %>% rownames_to_column()
  colnames(gene_tpm) <- c("SampleID", "tpm")
  gene_tpm <- merge(gene_tpm, sample_sheet)
  gene_tpm$Timepoint <- as.numeric(as.character(gene_tpm$Timepoint))
  
  
  fdiff <- gene_tpm %>% group_by(Timepoint) %>% 
    summarize(diff = abs(diff(tpm)))
  gene_tpm <- merge(gene_tpm, fdiff)
  
  gene_tpm$xoffset <- 0
  
  gene_tpm[which(!duplicated(gene_tpm[which(gene_tpm$diff < 3),"Timepoint"])), "xoffset"] <- 3
  
  gene_tpm$xTimepoint <- gene_tpm$Timepoint + gene_tpm$xoffset
  
  max_tpm <- max(gene_tpm$tpm)
  
  g <- ggplot(gene_tpm, aes(x = xTimepoint, y = tpm)) + 
    geom_smooth(method = "loess", color = rgb(215/255, 52/255, 88/255), alpha = 0.5,
                fill = rgb(247/255, 214/255, 222/255), fullrange=TRUE, span = 0.35)+
    geom_point(color = "#E8E8E8", size = 4) +
    geom_point(color = rgb(215/255, 52/255, 88/255), size = 3) +
    theme_tufte(base_size = 18, base_family = 'Helvetica') +
    theme(axis.line = element_line(colour = "#222222"),
          axis.title.x = element_text(colour = "#222222", margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(colour = "#222222", margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.y = element_text(color = "#222222", margin = margin(t = 0, r = 40, b = 0, l = 20)),
          axis.text.x = element_text(color = "#222222", margin = margin(t =40, r =40, b = 40, l = 40)),
          axis.ticks.length=unit(.3, "cm"),
          plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ 
    # scale_x_continuous(expand = c(0,0), limits = c(-4,333), breaks = seq(0,330,30)) +
    # scale_y_continuous(expand = c(0,0), limits = c(-1*max_tpm*0.15, max_tpm+max_tpm*0.15)) +
    coord_cartesian(xlim = c(-4,333), ylim=c(-1*max_tpm*0.05, max_tpm+max_tpm*0.15)) +
    coord_capped_cart(left = 'both', bottom = 'both') + 
    xlab("Time from induction (min)") +
    ylab(paste0(gene_name, " expression (fpkm)"))+ 
    ggtitle(paste0(gene_name, " expression"))
  
  show(g)
  return(g)
}