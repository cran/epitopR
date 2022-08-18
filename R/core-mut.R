#' @name core_mut
#' @title Determine presence of mutation in core binding sequence
#' @description The core_mut() function appends a new column to the peptide dataframe, identifying those that have a mutation in the core binding sequence. Since MHCII has an open binding pocket, it presents peptides that may be several amino acids longer than the core sequence pattern required to bind a particular MHC. In some cases, the user may want to filter their results, in order to keep only peptides with a mutation in the core binding sequence (when compared to the equivalent self peptide). In order to achieve this, the stimulating and self antigens are aligned using a multiple sequence alignment tool, from the bioconductor package msa. The sequence positions of the core in the stimulating peptide are determined, and sequences are kept only if there is a sequence difference between stimulating and self peptides at that position.
#' @param dat_in
#' dataframe, output of mhcII_hu(). The stimulating peptide, core pattern, start, and end positions will be pulled from this dataframe.
#' @param ag_stim
#' string, amino acid sequence of the stimulating antigen
#' @param ag_self
#' string, amino acid sequence of the self antigen
#' @export
#' @return data frame, peptide with flag of whether or not a mutation is in the core binding sequence
#' @import
#' tidyverse
#' @importFrom
#' Biostrings BString
#' @importFrom
#' Biostrings matchPattern

core_mut <- function(dat_in, ag_stim, ag_self) {
  if (!requireNamespace("Biostrings", quietly = TRUE) | !requireNamespace("msa", quietly = TRUE)) {
    stop(
      "Bioconductor package \"Biostrings\" or \"msa\" not found. Please install.",
      call. = FALSE
    )
  }

  # step 1: validations on aa code, seq aligement, and colnames
  if (str_detect(ag_stim, "B|J|O|U|X|Z") | str_detect(ag_self, "B|J|O|U|X|Z")) {
    stop("Input sequence contains invalid amino acid code. Please check.")
  }

  if (nchar(ag_stim) != nchar(ag_self)) {
    message("It seems input sequences are not aligned. aligning ....")
    algn <- msa(c(ag_stim, ag_self), type = "protein")
    ag_stim <- toupper(toString(unmasked(algn)[[1]]))
    ag_self <- toupper(toString(unmasked(algn)[[2]]))
    Sys.sleep(10)
    message("Alignment is done.")
  }

  if (!all(c("start", "end", "pep_stim", "core") %in% names(dat_in))) {
    stop("Please check if your nput data contains these columnes: start, end, pep_stim, core")
  }

  # step 2: working dataframe prep
  # add ag_stim, ag_self, and flag of whether or not pep_stim contains dashes
  nmer <- unique(nchar(dat_in$pep_stim))
  num_pep <- dim(dat_in)[1]

  # add ag_stim and ag_self to the data
  dat_in <- dat_in %>%
            mutate(ag_stim = ag_stim,
                   ag_self = ag_self)

  # step 3: pull out ag_self
  out0 <- list()
  cnt <- 1

  # for each pep_stim, pull out substr of ag_self, with dash flag
  # start and end are positions are changed to aligned ag_stim after pull_ag_self, are not the same as in iedb prediction tables
  for(i in 1:num_pep) {
    out0[[cnt]] <- pull_ag_self(dat_in[i,])
    cnt <- cnt + 1
  }

  out0 <- bind_rows(out0)

  # step 4: seperate data to with or without dashes
  # note: start and end positions from pull_ag_self call are refer to aligned ag_stim, not the same as in iedb prediction tables anymore
  wdash <- out0 %>% filter(dash == "yes") %>% select(-dash) %>% mutate(core_mut = "yes")
  wodash <- out0 %>% filter(dash == "no") %>% select(-dash)

  # step 5: for those rows have exact match between pep_stim and ag_stim, find mutation to core, then create final output table
  out <- find_core_mut(wodash) %>% rbind(wdash)
  return(out)
}



