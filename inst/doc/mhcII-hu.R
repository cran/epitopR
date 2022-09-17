## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----libraries, message=FALSE-------------------------------------------------
library(epitopR)
library(tidyverse)
library(here)
library(fs)
library(Biostrings)

## -----------------------------------------------------------------------------

# input antigen name
out_named_ag <- mhcII_hu(ag_present = c("DRB1_08_01"),
                 ag_stim = c("DQA1_01_01","DQA1_04_01"),
                 ag_self = c("DQA1_02_01"),
                 seq_len = 15,
                 method = "net")


## -----------------------------------------------------------------------------
# input protein sequence
dqa0101 <- "MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDQPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIQGLRSVGASRHQGPL"

dqa0202 <- "MILNKALMLGALALTTVMSPCGGEDIVADHVASYGVNLYQSYGPSGQFTHEFDGDEEFYVDLERKETVWKLPLFHRLRFDPQFALTNIAVLKHNLNILIKRSNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGHSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVLIIRGLRSVGASRHQGPL"

dqa0401 <- "MILNKALLLGALALTTVMSPCGGEDIVADHVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPVLRQFRFDPQFALTNIAVTKHNLNILIKRSNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGHSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIRGLRSVGASRHQGPL"

out_stored_ag <- mhcII_hu(ag_present =  c("DRB1_08_01"),
                  ag_stim = c(dqa0101, dqa0401),
                  ag_self = c(dqa0202),
                  seq_len = 15,
                  method = "net")


## -----------------------------------------------------------------------------
# run prediction of dqa0202 and dqa0401
out_single_ag <- mhcII_hu(ag_present =  c("DRB1_08_01"),
                  ag_stim = c(dqa0101),
                  ag_self = c(dqa0202),
                  seq_len = 15,
                  method = "net")

dqa0101_algn <- "MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDQPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIQGLRSVGASRHQGPL"

dqa0202_algn <- "MILNKALMLGALALTTVMSPCGGEDIVADHVASYGVNLYQSYGPSGQFTHEFDGDEEFYVDLERKETVWKLPLFHRLR-FDPQFALTNIAVLKHNLNILIKRSNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGHSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVLIIRGLRSVGASRHQGPL"

core_mut_result <- core_mut(out_single_ag, ag_stim=dqa0101_algn, ag_self=dqa0202_algn)

core_mut_result <- core_mut_result %>%
  filter(core_mut=="yes") %>%
  select(-c(core_mut))

## -----------------------------------------------------------------------------
out_stored_ag %>%
   mutate(antigen=as.factor(antigen)) %>%
   ggplot(aes(antigen, rank_val, color=as.factor(strength_rank))) +
   geom_jitter(width=0.2) +
   xlab("Stimulating Antigen") +
   ylab("Adjusted Rank (Percentile)") +
   labs(color = "Predicted Binding") +
   theme_light()

# Plot of rank x IC50

out_stored_ag %>%
  ggplot(aes(score_val, rank_val, color=as.factor(antigen), group=as.factor(antigen))) +
  geom_point() +
  facet_wrap(~as.factor(antigen)) +
  labs(color="Antigen")

out_stored_ag <- mhcII_hu(ag_present =  c("DRB1_08_01"),
                 ag_stim = c(dqa0101, dqa0401),
                 ag_self = c(dqa0202),
                 seq_len = 15,
                 method = "net",
                 cutoff_score = list(cutoff_netpan = c(50, 200),
                                   cutoff_comblib = c(50, 200),
                                   cutoff_nn_align = c(50, 200),
                                   cutoff_sturniolo = c(2)),
                 cutoff_rank= c(2, 10))

# Plot of rank x IC50

out_stored_ag %>%
  ggplot(aes(score_val, rank_val, color=as.factor(antigen), group=as.factor(antigen))) +
  geom_point() +
  facet_wrap(~as.factor(antigen)) +
  labs(color="Antigen")


## -----------------------------------------------------------------------------
dqa_01_peps <- out_named_ag %>%
  filter(antigen=="DQA1_01_01") %>%
  filter(strength_rank%in% c("weak", "strong"))

# library(ggseqlogo)
#ggseqlogo(dqa_01_peps$pep_stim)


