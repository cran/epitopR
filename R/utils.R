#' Basic functions of the package
#'
#' preproc_huII() format and validate allele names
#' @param allele_in a vector contains allele name(s)
#'
#' pull_seq_huII() pull out sequence of each allele based on ref table
#' @param alleles_in vector, allele names
#' @param tbl_ref_in dataframe, reference table, default is human_all.csv from github
#'
#' prep_lbl_huII() concatenate allele strings for prediction table names
#' @param alleles_in vector, allele names
#'
#' comb_pred_tbl() combine individual prediction tables by method, exclude none-binders and keep strong and weak binders only
#' @param nm_method string, prediction method used for IEDB prediction
#' @param nm_sht string, short name of alleles
#' @param nm_fd string, folder name which contains predict tables from IEDB
#' @param thold_score list of vectors, binder thresholds by score
#' @param thold_rank vector, binder thresholds by rank
#'
#' find_nonself_huII() find nonself binding peptides
#' @param dat_in dataframe, combined prediction tables from either of netmhciipan or recommended method
#'
#' pull_ag_self() pull out aligned ag_self based on aligned ag_stim position, add a core mutation flag
#' @param dat_in dataframe with pep_stim, core, aligned ag_stim and ag_self columns
#'
#' find_core_mut() find mutation position to core
#' @param dat_in dataframe with pep_stim, core, pep_self selected from pull_ag_self
#'
#' align_seq() align protein sequences
#' @param seq1 path to the fasta or fastq file to read
#' @param seq2 path to the fasta or fastq file to read
#'
#' @name utils
#' @import
#' dplyr
#' fs
#' janitor
#' stringr
#' tidyverse
#' utils
#' readr
#' @importFrom
#' stats setNames
#' @importFrom
#' purrr map_df
#' @importFrom
#' Biostrings readAAStringSet
#' @importFrom
#' Biostrings unmasked
#' @importFrom
#' msa msa
#' @importFrom
#' msa msaClustalW
#' @importFrom
#' seqinr read.fasta
#'
#'
NULL
#> NULL

#' @rdname utils
preproc_huII <- function(allele_in) {
  # 1. format conversion: a*01"01 to A_01_01
  allele_out <- toupper(allele_in) %>%
                str_replace_all(., "\\*|\\:", "_") %>%
                unique()

  # 2. letter code of protein sequences
  prt_seq <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  # 3. validation
  # if input without number but NOT all chars are in prt_seq, then print a warning message; else, pass
  # if input is allele but not in a good format ( c|cc|ccc+0or1digit+"_"+d|dd|ddd+"_"+d|dd|ddd ), then print a message; else, pass
  if(!all(str_detect(allele_out, "[0-9]"))) {
    if (!all(unique(unlist(strsplit(allele_out, split = ""))) %in% prt_seq)) {
        warning("Input sequence contains invalid character. Please check.")
    } else {
      return(allele_out)
    }
  } else {
    if (!all(str_detect(allele_out, "^[A-Za-z]{1,3}+[0-9]{1}+\\_+[0-9]{1,3}+\\_+[0-9]{1,3}$"))) {
      warning(paste0("Input allele is not in the expected format. Please check."))
    } else {
      return(allele_out)
    }
  }
}

#' @rdname utils
pull_seq_huII <- function(alleles_in,
                     tbl_ref_in) {
  tbl_seq <- tbl_ref_in %>% filter(allele %in% c(alleles_in))

  return(tbl_seq)
}

#' @rdname utils
prep_lbl_huII <- function(alleles_in) {
  lbl_out <- alleles_in %>%
    str_to_lower() %>%
    str_c(sep = "_", collapse = T) %>%
    str_replace_all(., "TRUE", "-")

  return(lbl_out)
}

#' @rdname utils
comb_pred_tbl <- function(nm_method, nm_sht, nm_fd, thold_score, thold_rank) {
  # list all of output files, and filter out empty ones (403 Forbidden, file starts with "<!DOCTYPE")
  fl <- dir_ls(nm_fd, glob =  paste0("*",nm_method, ".txt"))
  fl <- fl[sapply(fl, file.size) > 500]

  # map all non-empty predict tables into 1
  if (length(fl) > 1) {
    datout <- fl %>%
      setNames(nm = .) %>%
      map_df(~read_tsv(.x, col_types = cols(), col_names = TRUE), .id = "id")
  } else {
    datout <- data.frame()
  }

  ###*** start of netmhciipan ***###
  if (nm_method == "netmhciipan") {
    # score_val: IC50
    datout <- datout %>%
      mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      filter(ic50 != "-") %>%
      mutate(pep_stim = peptide,
             method = nm_method,
             core = core_peptide,
             score_val = as.numeric(ic50),
             rank_val = as.numeric(rank)) %>%
      mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_netpan) ~ "strong",
                                       score_val <= max(thold_score$cutoff_netpan) ~ "weak",
                              TRUE ~ "no"),
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                       rank_val <= max(thold_rank) ~ "weak",
                                       TRUE ~ "no")) %>%
      filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak"))

    if (dim(datout)[1] > 0) {
      datout <- datout %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
      datout <- data.frame()
    }
  }
  ###*** end of netmhciipan ***###

  ###*** start of recommended ***###
  if  (str_detect(nm_method, "rec")) {
    datout <- datout %>%
      mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      mutate(pep_stim = peptide)

    # start of comblib #
    # score_val: IC50 (column "comblib_score")
    comblib <- datout %>%
      filter(comblib_score != "-") %>%
      mutate(method = "comblib",
             core = comblib_core,
             score_val = as.numeric(comblib_score),
             rank_val = as.numeric(comblib_rank)) %>%
      mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_comblib) ~ "strong",
                                       score_val <= max(thold_score$cutoff_comblib) ~ "weak",
                                       TRUE ~ "no"),
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                       rank_val <= max(thold_rank) ~ "weak",
                                       TRUE ~ "no")) %>%
      filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak"))

    if(dim(comblib)[1] > 0 ) {
      comblib <- comblib %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method)  %>%
        distinct()
    } else {
      comblib <- comblib %>%
        mutate(ag_type = "") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method)
    }
    # end of comblib #

    # start of nn_align #
    # score_val: IC50
    nn_align <- datout %>%
      filter(nn_align_ic50 != "-") %>%
      mutate( method = "nn_align",
              core = nn_align_core,
              score_val = as.numeric(nn_align_ic50),
              rank_val = as.numeric(nn_align_rank)) %>%
       mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_nn_align) ~ "strong",
                                        score_val <= max(thold_score$cutoff_nn_align) ~ "weak",
                                        TRUE ~ "no"),
      strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                rank_val <= max(thold_rank) ~ "weak",
                                TRUE ~ "no")) %>%
      filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak"))

    if (dim(nn_align)[1] > 0) {
      nn_align <- nn_align %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
      nn_align <- nn_align %>%
        mutate(ag_type = "") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method)
    }
    # end of nn_align #

    # start of sturniolo #
    # score_val: sturniolo_score; only label strong binders
    sturniolo <- datout %>%
      filter(sturniolo_score != "-") %>%
      mutate(method = "sturniolo",
             core = sturniolo_core,
             score_val = as.numeric(sturniolo_score),
             rank_val = as.numeric(sturniolo_rank)) %>%
      mutate(strength_ic50 = case_when(score_val >= max(thold_score$cutoff_sturniolo) ~ "strong",
                                       TRUE ~ "no"),
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                       rank_val <= max(thold_rank) ~ "weak",
                                       TRUE ~ "no"))  %>%
      filter(strength_ic50 %in% c("strong") | strength_rank %in% c("strong"))

    if (dim(sturniolo)[1] > 0) {
      sturniolo <- sturniolo %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
      sturniolo <- sturniolo %>%
        mutate(ag_type = "") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method)
    }
    # end of sturniolo #

    if(dim(comblib)[1] > 0 | dim(nn_align)[1] > 0 | dim(sturniolo)[1] > 0) {
      datout <- rbind(comblib, nn_align, sturniolo)
    } else {
      datout <- data.frame()
    }

  }
  ###*** end of recommended ***###

  ###*** start of netpna el ***###
  if (str_detect(nm_method, "el")) {
    # score_val: percentile rank
    datout <- datout %>%
      mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      filter(rank != "-") %>%
      mutate(pep_stim = peptide,
             method = nm_method,
             core = core_peptide,
             score_val = as.numeric(score),
             rank_val = as.numeric(rank)) %>%
      mutate(strength_ic50 = "na",
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                      rank_val <= max(thold_rank) ~ "weak",
                                      TRUE ~ "no")) %>%
      filter(strength_rank %in% c("strong", "weak"))

    if (dim(datout)[1] > 0) {
      datout <- datout %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%

        distinct()
    } else {
      datout <- data.frame()
    }
  }
  ###*** end of netpan el ***###

  return(datout)
}

#' @rdname utils
find_nonself_huII <- function(dat_in) {
    self_peps <- dat_in %>%
                  filter(ag_type == "self")
    allo_peps <- dat_in %>%
                  filter(ag_type == "stim") %>%
                  anti_join(self_peps, by = "pep_stim") %>%
      select(-ag_type)
    return(allo_peps)
}

#' @rdname utils
pull_ag_self <- function(dat_in) {
  nmer <- nchar(dat_in$pep_stim)[1]
  # for each pep_stim, find starting position of best matches to aligned ag_stim, max allowed mismatch is 10
  dat_out <- Biostrings::matchPattern(dat_in$pep_stim[1],
                                      BString(dat_in$ag_stim[1]),
                                      fixed=TRUE,
                                      max.mismatch = 10)@ranges %>%
    data.frame() %>%
    rowwise() %>%
    # for each range of best match result
    mutate(pep_stim = dat_in$pep_stim[1],
           ag_stim = dat_in$ag_stim[1],
           ag_self = dat_in$ag_self[1],
           core = dat_in$core[1]) %>%
    mutate(pep1 = str_sub(ag_stim, start, end), # substr based on starting position
           cnt = str_count(pep1, "-")) %>% # then count number of dashes of each substr
    # remove existing dashes and fill out removed ones with number of none-dashes chars
    # we only do this once to simplify the process
    mutate(pep2 = str_sub(ag_stim, start, end+cnt), "-") %>%
    # find distance between each pair of substr and pep_stim, and keep the closest one
    mutate(cnt_match = adist(pep_stim, pep2)) %>%
    data.frame() %>%
    filter(cnt_match == min(cnt_match)) %>%
    # pull out substr of pep_stim and pep_self, it could be w or w/o dashes
    mutate(pep_stim_alg = str_sub(ag_stim, start, start + nmer - 1),
           pep_self = str_sub(ag_self, start, start + nmer - 1)) %>%
    # flag of yes or no of dash
    mutate(dash = ifelse(str_detect(pep_stim_alg, "-") | str_detect(pep_self, "-"),
                         "yes",
                         "no")) %>%
    #select(-pep_self) %>%
    #mutate(pep_self = pep2) %>%
    select(pep_stim, pep_self, core, start, end, dash)

  return(dat_out)
}

#' @rdname utils
find_core_mut <- function(dat_in) {
  # step 1. constants
  # length and number of peptide
  nmer <- unique(nchar(dat_in$pep_stim))
  num_pep <- length(dat_in$pep_stim)

  # step 2: add start and end positions of core to pep_stim
  dat_in <- dat_in %>%
    rowwise() %>%
    mutate(core_st = str_locate(pep_stim, core)[1],
           core_ed = str_locate(pep_stim, core)[2])

  # step 3: split pep_stim and pep_self into char, and combinde them into one matrix
  tmp_stim <- str_split_fixed(dat_in$pep_stim, "", nmer) %>%
    data.frame() %>%
    setNames(c(paste0("stim", seq(1:nmer))))

  tmp_self <- str_split_fixed(dat_in$pep_self, "", nmer) %>%
    data.frame() %>%
    setNames(c(paste0("self", seq(1:nmer))))

  mutat <- cbind(tmp_stim, tmp_self)

  # step 4: for each pep_stim, find mutation position
  for(i in 1:num_pep){
    for (j in 1:nmer){
      if(mutat[i,j] != mutat[i, j+nmer]){
        mutat[i,j] = j
      }
    }
  }

  mutat <- mutat %>%
    select(names(.)[str_detect(names(.), "stim")]) %>%
    mutate_if(is.character, str_replace_all, pattern = "[A-Za-z]", replacement = '0') %>%
    mutate_all(as.numeric)

  # step 5: add flag of mutation to core position
  dat_out <- cbind(dat_in, mutat) %>%
    rowwise() %>%
    mutate(core_mut = ifelse(any(across(names(.)[str_detect(names(.), "stim")]) %in% c(seq(core_st, core_ed, 1))), "yes", "no")) %>%
    select(pep_stim, pep_self, core, start, end, core_mut)

  return(dat_out)
}

#' @rdname utils
align_seq <- function(seq1, seq2) {
  # file validation
  if (file.exists(seq1) & file.exists(seq2) ){
    nms <- c(seq1, seq2)
  } else {
    stop("invaild file(s), please check.")
  }

  # load sequences to AA string object, and align the sequences
  algn <- msa(readAAStringSet(nms))

  # pull out aligned sequences
  algn.seq <- c(toString(unmasked(algn)[[1]]), toString(unmasked(algn)[[2]]))

  return(c(toupper(algn.seq[1]), toupper(algn.seq[2])))
}

