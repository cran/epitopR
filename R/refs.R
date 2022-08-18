#' @name prep_ref_hu
#' @title prep human reference tables based on IMGT database
#' @param url_imgt, string, github link of IMGT database
#' @export
#' @return data frame, MHC II sequence of officially named alleles
#' @importFrom
#' seqinr read.fasta
#' @importFrom
#' dplyr slice
#' @examples
#' \donttest{
#' out <- prep_ref_hu()
#' }
prep_ref_hu <- function(url_imgt = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/") {
# locations of IPD-IMGT/HLA Database -
# github: https://github.com/ANHIG/IMGTHLA/tree/Latest/fasta

  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop(
      "Package \"seqinr\" must be installed to use this function.",
      call. = FALSE
    )
  }

# list of available locus names
nm_lc <- c("A", "B", "C",
           "DPA1", "DPB1",
           "DQA1", "DQA2",
           "DQB1",
           "DRB1", "DRB345")

for(i in 1:length(nm_lc)) {
  tmp<-  read.fasta(paste0(url_imgt,
                           nm_lc[i], "_prot.fasta"),
                    seqtype = "AA",
                    as.string = TRUE)

  seqs <- data.frame(id	= "", allele	= "", length	= "", seq	= "",
                     allele_short_nm	= "",	allele_query_nm	= "")

for (j in 1:length(tmp)) {
    # 1. create allele info based on annotation string
    # annotation string: ">HLA:HLA00401 C*01:02:01:01 366 bp"
    annot <- attr(tmp[[j]],"Annot")
    allele <- gsub(".* ", "", str_remove(annot, "\\ [0-9]{1,3}+\\ bp"))
    allele_short <- str_extract(allele,"^[A-Za-z]{0,3}+[0-9]{0,3}+\\*+[0-9]{1,3}+\\:+[0-9]{1,3}")

    seqs[j, ]$id <- str_remove(gsub(" .*", "", annot), ">")
    seqs[j, ]$allele <- allele
    seqs[j, ]$length <- gsub(".* ", "", trimws(str_remove(annot, "bp"), which = "right"))

    seqs[j, ]$allele_short_nm <- str_replace_all(allele_short,
                                                 "\\*|\\:", "_")
    seqs[j, ]$allele_query_nm <- paste0("HLA-", allele_short)

    # 2. pull out AA sequence for each allele
    seqs[j, ]$seq <- tmp[[j]][1]

    }
    assign(nm_lc[i], seqs)
}

# 3. prep reference table for IEDB prediction. Only longest sequence of the allele are needed in the final table
  final <- bind_rows(A, B, C, DPA1, DPB1, DQA1, DQA2, DQB1, DRB1, DRB345) %>%
            mutate(allele = allele_short_nm) %>%
            select(allele, seq, length) %>%
            distinct() %>%
            group_by(allele) %>%
            dplyr::slice(which.max(length)) %>%
            select(allele, seq) %>%
            data.frame()

  row.names(final) <- seq(1:dim(final)[1])
  # write.table(final, gzfile("inst/extdata/ref/ref_human.gz"))

  return(final)
}




