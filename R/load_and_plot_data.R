#' @import data.table
#' @import ggplot2
#' @import patchwork
#' @export
load_and_plot_data <- function(file.name, p.val.col, tool.name, interaction_var, p.sig=1E-6, p.sugg=NA, label.qq = T) {

  # Read data file
  result.table <- fread(file.name)

  # Convert all col names to lc to avoid stupid conflicts:
  setnames(result.table, names(result.table), str_to_lower(names(result.table)))
  p.val.col <- str_to_lower(p.val.col)

  # Set gene and p. value column to something standard:
  setnames(result.table, p.val.col, "p.value.selected")
  result.table[,log.p:=-log10(p.value.selected)] # Add -log10 p value
  if ('chr' %in% names(result.table)) {  # Account for BOLT-LMM results...
    setnames(result.table, 'chr', 'chrom')
    result.table[,chrom:=as.character(chrom)]
  }
  result.table[,chrom:=factor(chrom, levels = mean_chr_pos[,chrom])] # Make sure chromosomes are standardised to those in transcripts

  masks <- validate_row(result.table, 'mask', 'ALL')
  mafs <- validate_row(result.table, 'maf', 'ALL')
  tests <- validate_row(result.table, 'test', 'ALL')
  validate_row(result.table,'ac',40)  # This is just to make sure imputed variants aren't removed due to how I filter gene masks. It automatically adds the 'ac' column if it doesn't exist with a value of 40 (to GWAS results)

  # Make standard plots
  # First part just holds a list of possible plots:
  plots <- list()

  # And then loop through all possible combinations
  for (curr_mask in masks) {
    for (curr_maf in mafs) {
      for (curr_test in tests) {

        # Set names that are actually variable in the dataset
        plot_vars <- c(curr_mask,curr_maf,curr_test)
        plot_vars <- plot_vars[str_detect(plot_vars,'ALL',negate=T)]
        if (length(plot_vars) == 0) {
          name <- tool.name
        } else {
          name <- paste(tool.name,paste(plot_vars,collapse='-'),sep="-")
        }

        # For the test of the interaction variant, we DON'T want to label ALL genes that
        # pass significance, because this should be the vast majority of them. We
        # want to label only ones that might be putative 'fails'. Also don't need a QQ plot
        # as it doesn't make sense for this particular test

        # Get the y max value empirically from the data:
        # Also set to a minimum of 10 so we don't get squashed plots
        if (result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test,max(log.p,na.rm = T)] > 10) {
          ymax <- result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test,max(log.p,na.rm = T)]
          ymax <- ymax + (ymax * 0.05)
        } else {
          ymax <- 10
        }

        if (curr_test == paste0('ADD-INT_', interaction_var)) {
          manh.plot <- plot_manh(result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test],
                                 ymax=ymax, p.sig = p.sig, p.sugg = p.sugg, invert_label = T)
          qq.plot <- plot_qq(result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test], ymax=ymax, is.null = T, label.markers = label.qq)
        } else {
          manh.plot <- plot_manh(stats = result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test],
                                 ymax=ymax, p.sig = p.sig, p.sugg = p.sugg)
          qq.plot <- plot_qq(result.table[ac > 30 & !is.infinite(log.p) & mask == curr_mask & maf == curr_maf & test == curr_test], ymax=ymax, label.markers = label.qq)
        }

        comb.plot <- manh.plot[[1]] + qq.plot +
          plot_layout(ncol = 2, nrow = 1, widths = c(3,1.3)) +
          plot_annotation(title = str_replace(name, '-', ' - '), theme = theme(plot.title=element_text(size=18,face="bold",colour="black")))

        plots[[name]] <- list('manh.plot' = manh.plot[[1]],
                              'clusters' = manh.plot[[2]],
                              'qq.plot' = qq.plot,
                              'comb.plot' = comb.plot)
      }
    }
  }

  return(list('gene.table' = result.table,
              'plots' = plots,
              'masks' = masks))

}

#' @import data.table
validate_row <- function(table, row_name, default) {

  if (row_name %in% names(table)) {
    unique_items <- table[get(row_name) != "" & !is.na(get(row_name)), unique(get(row_name))]
  } else {
    table[,eval(row_name):=default]
    unique_items <- c(default)
  }
  return(unique_items)

}

#' @import ggplot2
#' @import data.table
#' @import lemon
plot_manh <- function(stats, ymax, p.sig, p.sugg, label_data = data.table(), invert_label=F) {

  if (!is.na(p.sugg)) {
    sig.lines <- -log10(c(p.sig,p.sugg))
    sig.colours <- c('red','orange')
  } else {
    sig.lines <- -log10(c(p.sig))
    sig.colours <- c('red')
  }

  # text.sig <- min(sig.lines)  # Set minimum value to label genes
  #
  # if (invert_label) {
  #   label_data <- sig_cluster(stats[log.p < text.sig])
  # } else {
  #   label_data <- sig_cluster(stats[log.p > text.sig])
  # }

  manh.plot <- ggplot(stats, aes(manh.pos, log.p, colour = chrom)) +
    geom_point(size=0.25) +
    geom_hline(yintercept = sig.lines, colour = sig.colours, linetype = 2) +
    scale_x_continuous(name = "Chromosome", label = mean_chr_pos[,chrom], breaks = mean_chr_pos[,mean_pos]) +
    scale_y_continuous(name = expression(bold(-log[10](italic(p)))), limits = c(0,ymax)) +
    scale_colour_manual(values = rep(get_hg_palette()[1:2],12)[1:23], breaks = c(as.character(1:22),"X")) +
    coord_capped_cart(bottom="both",left="both") + # Comes from the "lemon" package
    get_hg_theme() + theme(panel.grid.major = element_blank())

  if (nrow(label_data) != 0) {
    manh.plot <- manh.plot + geom_text(inherit.aes = F, data = label_data, aes(manh.pos, log.p, label = symbol),hjust=0,angle=45, size=3)
  }

  return(list(manh.plot,label_data))

}

#' @import ggplot2
#' @import data.table
plot_qq <- function(stats, ymax, label.y = FALSE, is.null = FALSE, label.markers = TRUE) {

  ## QQplot
  qqplot.data <- data.table(observed = stats[,p.value.selected],
                            symbol = stats[,symbol])
  setkey(qqplot.data,observed)
  qqplot.data <- qqplot.data[!is.na(observed)]
  lam <- qqplot.data[,median(qchisq(1 - observed, 1))] / qchisq(0.5, 1)
  qqplot.data[,observed:=-log10(observed)]
  qqplot.data[,expected:=-log10(ppoints(nrow(qqplot.data)))]

  # make sure we set a reasonable y-max based on number of points (e.g. GWAS/Burden
  # will be different)
  expected_ymax <- qqplot.data[,ceiling(max(expected))]

  if (is.null) {

    qq.plot <- ggplot() +
      annotate(geom='text', x=expected_ymax / 2, y = ymax / 2, label = 'NULL QQ PLOT', hjust = 0.5, vjust = 0.5,size=10)

  } else {

    qq.plot <- ggplot(qqplot.data, aes(expected, observed)) +
      geom_abline(slope = 1, intercept = 0,colour="red") +
      geom_point(size=0.25)

    if (label.markers) {
      qq.plot <- qq.plot +
        geom_text(inherit.aes = F, data = qqplot.data[observed > -log10(1.4e-6)], aes(expected, observed, label = symbol), position = position_nudge(-0.04,0),hjust=1, size=3)
    }
  }

  qq.plot <- qq.plot +
    scale_x_continuous(name = expression(bold(Expected~-log[10](italic(p)))), limits = c(0,expected_ymax)) +
    scale_y_continuous(name = ifelse(label.y,expression(bold(Observed~-log[10](italic(p)))),expression('')), limits = c(0,ymax)) +
    # Note to future dev: Make sure to leave list(bquote) otherwise you get confusing warnings
    annotate('text', x = 1, y = expected_ymax * 0.9, hjust=0, vjust=0.5, label=list(bquote(lambda==.(lam))), parse=TRUE) +
    get_hg_theme() + theme(panel.grid.major = element_blank())

  if (!label.y) {
    qq.plot <- qq.plot + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }

  return(qq.plot)

}

