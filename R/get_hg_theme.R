#' @import ggplot2
#' @export
get_hg_theme <- function(add_legend=FALSE) {

  theme <- theme(panel.background=element_rect(fill="white"),
                 line=element_line(size=1,colour="black",lineend="round"),
                 axis.line=element_line(size=1),
                 text=element_text(size=16,face="bold",colour="black"),
                 axis.text=element_text(colour="black"),
                 axis.ticks=element_line(size=1,colour="black"),
                 axis.ticks.length=unit(.1,"cm"),
                 strip.background=element_rect(fill="white"),
                 axis.text.x=element_text(angle=45,hjust=1),
                 legend.position=ifelse(add_legend, "right", "blank"),
                 panel.grid.major=element_line(colour="grey",size=0.5),
                 legend.key=element_blank())

  return(theme)

}
