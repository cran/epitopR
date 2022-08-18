#' @name mhcII_hu
#' @title Prediction of Nonself Peptide Presentation on Human MHC Class II
#' @description It determines which non-self peptides can be presented by a given HLA class II allele. This function takes a sequence for a stimulating antigen and the corresponding self antigen, and given a defined sequence length, queries the IEDB API with the user's choice of peptide binding prediction method. The set of peptides present in the results for the stimulating antigen but not the self antigen are then carried forward as non-self peptides. If desired, the user can adjust the default thresholds (by IC50 binding affinity or percentile rank) used to define "strong" and "weak" binders. The output is a dataframe of non-self peptides that are predicted to bind to the presenting allele.
#' @param ag_present
#' character vector, presenting allele, formatted with either "_", "*", or":" separating loci, antigen, and allele. For example, "DRB1_08_01".
#' @param ag_stim
#' character vector, stimulating antigen, can either be an HLA class II allele entered in the same format as ag_present, or a character vector of the amino acid sequence of the protein
#' @param ag_self
#' character, self antigen, can either be an HLA class II allele entered in the same format as ag_present, or a character vector of the amino acid sequence of the protein
#' @param seq_len
#' string, length of peptides to consider
#' @param fd_out
#' string, output folder name; default output is current working directory
#' @param method
#' string, IEDB prediction method to be used. Options are "netmhciipan', "netmhciipan_el" or "recommended." Default is netmhciipan.
#' @param cutoff_score
#' list of vectors. Defines the thresholds required to be included in results, and to be labeled, "strong" or "weak" binder. Multiple prediction methods are used, each of which provide different raw outputs (i.e. IC50, "strength", "score"). Our justification for the default thresholds is listed in the mhcII_hu vignette, however the user may choose to specify alternate cutoffs if desired.
#' @param cutoff_rank
#' vector, IEDB adjusts all outputs in comparison to a set of random natural peptides in order to determine an normalized adjusted percentile rank. With normalized ranks, the same thresholds can be used across different methods. Default thresholds are <2\% for strong binders and <10\% for weak binders.
#' @param url_iedb
#' string, iedb api url
#'
#' @export
#' @return data frame, MHC II binding prediction result table
#' @import
#' dplyr
#' fs
#' here
#' httr
#' janitor
#' stringr
#' tibble
#' tidyverse
#' @importFrom
#' stats setNames
#' @examples
#' \donttest{
#' mhcII_hu(ag_present=c("DRB1_08_01"),ag_stim=c("DQA1_01_01","DQA1_04_01"),ag_self=c("DQA1_02_01"))
#' }

mhcII_hu <- function(ag_present, ag_stim, ag_self,
                     seq_len = '15',
                     #fd_out = as.character(paste0(here(), "/", "outputs", "/")),
                     fd_out = as.character(paste0(tempdir(), "/", "outputs", "/")),
                     method = "netmhciipan",
                     cutoff_score = list(cutoff_netpan = c(50, 500),
                                   cutoff_comblib = c(50, 500),
                                   cutoff_nn_align = c(50, 500),
                                   cutoff_sturniolo = c(2),
                                   cutoff_el = c(2, 10)),
                     cutoff_rank = c(2, 10),
                     url_iedb = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/') {

  #* pre check: method and folder *#
  # allow ambiguous string of method
  if (str_detect(method, "net|pan|mhcii") & !str_detect(method, "el")) {
    # as of 2002/06, netpan default versions: api is v4.0, website only up to v3.2.
    # v4.0 generates identical result as netpan_ba
    method <- "netmhciipan"
  } else if (str_detect(method, "el")) {
    method <- "netmhciipan_el"
  } else if (str_detect(method, "rec|com")){
    method <- "recommended"
  } else {
    method = "netmhciipan"
  }

  if (!requireNamespace("httr", quietly = TRUE)) {
    stop(
      "Package \"httr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  #* end of pre check *#

  #* 1. setup: ref table and folders *#
  #ref_hu <- read.table(gzfile(here("inst/extdata/ref/ref_human.gz")),row.names=1)
  #ref_hu <- read.table(gzfile("~/projects/epitopR/inst/extdata/ref/ref_human.gz"),row.names=1)
  ref_hu <- read.table(gzfile(system.file("extdata/ref", "ref_human.gz", package = "epitopR")), row.names=1)

  # if output folder doesn't exist, create one; else, clean up the folder
  if(!dir.exists(fd_out)) {
    dir.create(file.path(fd_out))
    #Sys.sleep(2)
  } else {
    file.remove(fs::dir_ls(fd_out, glob = "*.txt"))
    #Sys.sleep(2)
  }
  #* end of setup *#

  #* 2. preprocess *#
  ag_present <- preproc_huII(ag_present)
  ag_self <- preproc_huII(ag_self)
  # keep the elements of stim that are not in self for computational efficiency
  # this is different concept with find_nonself_huII()
  ag_stim <- preproc_huII(ag_stim)
  ag_stim <- setdiff(ag_stim, ag_self)

  if(all(ag_stim %in% ag_self)) {
    warning("stim antigens are identical with self antigens, no prediction will be made.")
  }

  if(all(ag_present %in% ag_stim)) {
    warning("presenting antigens are identical with stim antigens, no prediction will be made.")
  }
  #* end of preprocess *#

  #* 3. prep for api calls *#
  # 3.a pull out sequences
  if (!all(str_detect(ag_stim, "[0-9]"))) {
    ref_self <-  as.data.frame(ag_self) %>%
      mutate(allele = paste0("self_seq_",seq(nrow(.))),
             seq = ag_self,
             allele_query_nm = allele) %>%
      select(allele, seq, allele_query_nm)
  } else {
    ref_self <- pull_seq_huII(alleles_in = ag_self,
                              tbl_ref_in = ref_hu) %>%
      mutate(allele_query_nm = gsub("_", ":", sub("_", "*", allele)))
  }

  if (!all(str_detect(ag_stim, "[0-9]"))) {
    ref_stim <-  as.data.frame(ag_stim) %>%
      mutate(allele = paste0("stim_seq_",seq(nrow(.))),
             seq = ag_stim,
             allele_query_nm = allele) %>%
      select(allele, seq, allele_query_nm)
  } else {
    ref_stim <- pull_seq_huII(alleles_in = ag_stim,
                              tbl_ref_in = ref_hu) %>%
      mutate(allele_query_nm = gsub("_", ":", sub("_", "*", allele)))
  }

  ref_self_stim <- rbind(ref_self, ref_stim)

  if (!all(str_detect(ag_present, "[0-9]"))) {
    stop("please specify MHC allele name of your presenting antigen.")
  } else {
    ref_present <- pull_seq_huII(alleles_in = ag_present,
                                 tbl_ref_in = ref_hu) %>%
      mutate(allele_query_nm = gsub("_", ":", sub("_", "*", allele)))
  }

  # 3.b queries
  seq_names <- ref_self_stim %>%
    rowid_to_column() %>%
    mutate(antigen = allele,
           allele_query_nm = paste0("HLA-", str_replace_all(antigen, "_", "*")),
           iedb_query_seq = seq,
           length = nchar(iedb_query_seq),
           allele = str_replace_all(antigen, "_", "*"),
           ag_type = if_else(str_detect(antigen, "self"), "self",
                             if_else(str_detect(antigen, "stim"), "stim",
                                     if_else(antigen %in% ag_self, "self", "stim"))),
           seq_num = rowid) %>%
    select(seq_num, allele, length, antigen, allele_query_nm, iedb_query_seq, ag_type)

  seq_names_sht <- seq_names %>%
    select(antigen, ag_type)

  #* end of prep for api *#

  #* 4. iedb prediction *#
  for(j in 1:nrow(seq_names)) { # for each antigen
    # set parameter list
    param2httr = list('method' = method,
                      'sequence_text' = as.character(seq_names$iedb_query_seq[j]),
                      'allele' = str_c(ref_present$allele_query_nm, collapse = ","),
                      'length' = seq_len)
    # send request to api
    res_api <- httr::POST(url = url_iedb,
                          body = param2httr)

    # extract content from the request
    cat(httr::content(res_api, "text"),
        file = paste0(fd_out, seq_names$antigen[j], "_", method, ".txt"))
  }
  rm(j)
  #* end of prediction *#

  #* 5. combine all prediction result into one *#
  # call to combine individual prediction report of each method
  tmp <- comb_pred_tbl(nm_method = method,
                       nm_sht = seq_names_sht,
                       nm_fd = fd_out,
                       thold_score = cutoff_score,
                       thold_rank = cutoff_rank)
  # assign method name to the combined report
  assign(method,tmp)
  rm(tmp)

  #* 6. peptide table *#
  # if na report, then empty dataframe; else, return dataframe without self peptides
  if (exists('netmhciipan')) { # start of if method = netpan
    if(dim(netmhciipan)[1] == 0) {
      peps_all <- data.frame()
    } else {
      peps_all <- netmhciipan %>%
        find_nonself_huII()
    } # end of netpan
  } else if (exists('netmhciipan_el')) { # start of if method = el
    if(dim(netmhciipan_el)[1] == 0) {
      peps_all <- data.frame()
    } else {
      peps_all <- netmhciipan_el %>%
        find_nonself_huII()
    }
  } else if (exists('recommended')) { # start of if method = recommended
    if(dim(recommended)[1] == 0) {
      peps_all <- data.frame()
    } else {
      peps_all <- recommended %>%
        find_nonself_huII()
    }
  } else {
    peps_all <- data.frame()
  }

  #* 7. final table *#
  # empty table if no predict peptide; else return table sore by rank
  if (dim(peps_all)[1] == 0) {
    final <- data.frame()
  } else if (dim(peps_all)[1] == 1 & is.na(peps_all[1,]$pep_stim)) {
    final <- peps_all
  } else {
    final <- peps_all %>%
      group_by(pep_stim, core) %>%
      arrange(desc(score_val), .by_group = TRUE) %>%
      ungroup()
  }

  unlink(fd_out, recursive = TRUE)

  return(final)
}
