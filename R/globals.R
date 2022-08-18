# global variables
utils::globalVariables(
  c(# vars for refs.R
    "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQA2", "DQB1", "DRB1", "DRB345",

  # vars for mhcII-hu.R
    ".", "seq_num", "allele", "antigen", "ag_type", "binder", "adjusted_rank", "peptide", "net_core",
    "iedb_core", "net_rank", "dupe_count", "net_antigen", "net_allele", "presenting_allele", "stim_antigen",
    "net_binder", "iedb_binder", "iedb_rank", "iedb_antigen", "iedb_allele", "iedb_query_id", "iedb_query_seq",
    "allele_short_nm", "rowid", "query_fasta", "netmhciipan", "core_peptide", "recommended",
    "core", "method", "stim", "present_allele", "netmhciipan_core", "allele_query_nm", "comblib",
    "sturniolo", "netmhciipan_el",  "thold",

    # vars for comb_pred_tbl()
    "comblib_score",  "nn_align_ic50",  "sturniolo_score",
    "comblib_core", "nn_align_core", "sturniolo_core",
    "sturniolo_rank", "nn_align_rank", "comblib_rank",
    "ic50",  "score", "score_var", "rank_var",
    "strength_ic50", "strength_rank",
    "score_val", "ranl_val",

  # vars for core mutation related
  "start", "end", "unmasked", "core_ed", "core_st", "flag", "self", "start.tmp", "start.x.algn",
  "sub.cnt", "x.algn", "y.algn", "nmer", "num_pep", "out", "out0", "cnt", "wdash", "wodash",
  "ag_self", "ag_stim", "cnt_match", "dash", "map_df", "pep1", "pep2", "pep_self", "pep_stim",
  "pep_stim_alg", "rank_val"
  )
)


