test_that("test mhcII_hu prediction works", {
  pred_table <- mhcII_hu(ag_present = c("DRB1_08_01"),
           ag_stim = c("DQA1_01_01","DQA1_04_01"),
           ag_self = c("DQA1_02_01"),
           seq_len = 15,
           method = "net")
})
