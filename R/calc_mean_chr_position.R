#' @import data.table
calc_mean_chr_position <- function(transcripts) {

  mean_chr_pos <- transcripts[,mean(manh.pos),by="chrom"]
  setnames(mean_chr_pos, "V1", "mean_pos")
  return(mean_chr_pos)

}
