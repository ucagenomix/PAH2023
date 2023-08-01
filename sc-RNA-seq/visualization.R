# Adapted from https://github.com/satijalab/seurat/blob/master/R/visualization.R
Dim_Plot <- function(obj, 
                     group.by, 
                     cell.name = NULL, 
                     #cell.number = F,
                     label.show = F, 
                     colors.use = NULL, 
                     figure.plot = F,
                     legend = T, 
                     legend.title = NULL, 
                     base.size = 8, 
                     pt.size = NULL, 
                     key.size = 3.5, 
                     label.size = 4, 
                     theme.bw = F,...) {
  require(Seurat)
  require(ggplot2)
  
  if ("Seurat" != class(obj)[1]) {
    stop("object should be of class Seurat")
  }
  
  if (!missing(x = group.by)) {
    obj <- SetIdent(object = obj, value = group.by)
  }
  
  plot <- DimPlot(object = obj, group.by = group.by, cols = colors.use, label = label.show, label.size = label.size,
                  shuffle = T, raster = F, pt.size = pt.size, repel = T,...)
  
  plot <- plot + theme(text = element_text(size = base.size, family = "Helvetica"),
                       axis.text = element_text(size = base.size, family = "Helvetica"))
  
  if(theme.bw){ 
    plot <- plot + 
      theme_bw(base_line_size = 0) +
      theme(axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black")
      )
  } 
  
  # if (!cell.number) {
  #   plot
  # }
  
  if (is.null(cell.name)) {
    Idents(obj) <- group.by
    # Calculate number of cells per cluster 
    cell.nb <- table(obj@active.ident)
    # Add cell number per cluster to cluster labels
    ClusterLabels = paste(names(cell.nb), paste0("(n = ", cell.nb, ")"))
    ClusterBreaks = names(cell.nb)
  }else{  
    Idents(obj) <- cell.name
    cell.nb <- table(obj@active.ident)
    ClusterLabels = paste(names(cell.nb), paste0("(n = ", cell.nb, ")"))
    ClusterBreaks = names(cell.nb)
  }
  
  plot <- plot + theme(
    plot.title = element_text(hjust = 0.5, size = base.size, family = "Helvetica", face = "bold"),
    text = element_text(size = base.size, family = "Helvetica"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.direction ="vertical",
    legend.position = "right",
    legend.text = element_text(size = base.size),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm')
  )  
  plot <- plot + guides(colour = guide_legend(override.aes = list(size = key.size),
                                              nrow = length(unique(obj@active.ident)))) 
  plot <- plot + scale_color_discrete(labels = ClusterLabels, 
                                      #breaks = ClusterBreaks, 
                                      type = colors.use) 
  plot <- plot + labs(title = "", colour = paste(legend.title))
  
  if (figure.plot && theme.bw) {
    stop(message = "Cannot return plot if figure.plot = T and theme.bw = T.",
         "i" = "Set theme.bw = F")
  }
  
  if (figure.plot) {
    
    if(!legend){
      plot <- plot & NoLegend()
    }
    
    # pull axis labels
    x.lab.reduc <- plot$labels$x
    y.lab.reduc <- plot$labels$y
    
    plot <- plot & NoAxes()
    
    axis.plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_point() +
      xlim(c(0, 10)) + ylim(c(0, 10)) +
      theme_classic() +
      ylab(y.lab.reduc) + xlab(x.lab.reduc) +
      theme(plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent"),
            text = element_text(size = base.size, family = "Helvetica", colour = "black"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(
              arrow = arrow(angle = 15, length = unit(.5, "cm"),  type = "closed")
            )
      )
    
    figure.layout <- c(
      area(t = 1, l = 1, b = 11, r = 11),
      area(t = 10, l = 1, b = 11, r = 2))
    
    plot.figure <- plot + axis.plot +
      patchwork::plot_layout(design = figure.layout)
    return(plot.figure)
  } else {
    return(plot)
  }
}


# From https://github.com/xmc811/Scillus/blob/master/R/visualization.R

Plot_Stats <- function(sobj, 
                       plot.type, 
                       ident = 'seurat_clusters',
                       group.by = "sample",
                       pal.setup = 'Set1',
                       ncol = 4,
                       plot.ratio = 1,
                       text.size = 10,
                       tilt.text = FALSE) {
  
  
  set_colors <- function(pal, n) {
    
    if (all(pal %in% rownames(RColorBrewer::brewer.pal.info))) {
      num <- c()
      for (i in seq(length(pal))) {
        num[i] <- RColorBrewer::brewer.pal.info[pal[i],][[1]]
      }
      full.pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
    } else if (all(are.colors(pal))) {
      full.pal <- pal
    } else {
      stop('Incorrect custom.colors setup. Please input valid RColorBrewer custom.colors names or color names.')
    }
    
    if (n <= length(full.pal)) {
      return(full.pal[1:n])
    } else {
      warning("Number of colors required exceeds custom.colors capacity. RdYlBu spectrum will be used instead.", 
              immediate. = TRUE)
      #return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
      #return(colorRampPalette(rev(brewer.pal(11, "PiYG")))(n))
      return(colorRampPalette(rev(brewer.pal(14, "RdYlGn")))(n))
      
    }
  }
  
  are.colors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  if (is.data.frame(pal.setup)) {
    pal <- pal.setup[pal.setup[[1]] == group.by,][[2]]
  } else {
    pal <- pal.setup
  }
  
  stat <- tibble::tibble(group = sobj[[group.by]][[1]], 
                         ident = sobj[[ident]][[1]])
  stat %<>%
    group_by(.data$group, 
             .data$ident) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  ncolors <- if (plot.type == 'prop_fill') {
    length(unique(sobj[[ident]][[1]]))
  } else {
    length(unique(sobj[[group.by]][[1]]))
  }
  
  colors <- set_colors(pal, ncolors)
  
  thm <- theme(aspect.ratio = plot.ratio,
               legend.title = element_text(size = text.size, color = "black", family = "Helvetica"),
               legend.text = element_text(size = text.size, color = "black", family = "Helvetica"),
               axis.title = element_text(size = 8, color = "black", family = "Helvetica"),
               axis.text = element_text(size = text.size, color = "black", family = "Helvetica"),
               axis.title.x = element_blank(),
               axis.text.x = element_text(size = text.size,color = "black", family = "Helvetica"),
               axis.text.y = element_text(size = text.size, color = "black", family = "Helvetica")
  ) + theme_classic()
  
  thm2 <- theme(axis.text.x = element_text(size = text.size, color = "black", family = "Helvetica"),
                axis.text.y = element_text(size = text.size, color = "black", family = "Helvetica"),
                axis.title = element_text(size = 8, color = "black", family = "Helvetica"),
                legend.position = "none")
  
  thm3 <- theme(axis.text.x = element_text(angle = 45, size = text.size,
                                           color = "black", family = "Helvetica",hjust = 1,vjust = 1),
                axis.text.y = element_text(size = text.size, color = "black", family = "Helvetica"),
                axis.title = element_text(size = 8, color = "black", family = "Helvetica")
  ) 
  
  switch(plot.type,
         group.count = stat %>%
           group_by(.data$group) %>%
           summarise(sum(n)) %>%
           ggplot(aes(x = .data$group, 
                      y = .data$`sum(n)`)) +
           geom_col(aes(fill = .data$group)) +
           geom_text(aes(label = .data$`sum(n)`), 
                     vjust = -0.5, 
                     size = text.size * 0.35) +
           scale_fill_manual(values = colors) + 
           labs(x = group.by, y = "Number of Cells") + 
           thm + thm2 + if (tilt.text) {thm3},
         
         prop_fill = ggplot(stat) + 
           geom_bar(aes(x = .data$group, 
                        y = .data$freq, 
                        fill = .data$ident), 
                    position = "fill", 
                    stat = "identity") +
           #scale_y_continuous(labels = scales::percent) +
           scale_y_continuous(labels = waiver()) +
           scale_fill_manual(values = colors, 
                             name = "Cluster") +
           labs(x = group.by, y = "Proportions (%)") + 
           thm + if (tilt.text) {thm3},
         
         prop_multi = stat %>%
           mutate(freq = round(.data$freq, 3)) %>%
           ggplot() + 
           geom_bar(aes(x = .data$group,
                        y = .data$freq, 
                        fill = .data$group), 
                    stat = "identity") +
           # geom_text(aes(x = .data$group, 
           #               y = .data$freq, 
           #               label = scales::percent(.data$freq)), 
           #           vjust = -0.5, 
           #           size = text.size * 0.35) +
           scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                              #labels = scales::percent.format()
                              labels = waiver()
           ) +
           facet_wrap(~ ident, 
                      ncol = ncol, 
                      scales = "free") +
           scale_fill_manual(values = colors, 
                             name = "Group") +
           labs(x = NULL, 
                y = "Proportion (%)") + 
           theme(strip.text.x = element_text(size = text.size),
                 strip.background = element_blank(), 
                 strip.placement = "outside") + 
           thm + thm2 + if (tilt.text) {thm3},
         
         stop("Unknown plot type")
  )
}

