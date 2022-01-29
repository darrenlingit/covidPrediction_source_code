rm(list=ls())
input_dir <- "model/"
input_files_bkwr <-
  list.files(path = input_dir,
             pattern = "bkwr_*",
             full.names = TRUE)
input_files_forw <-
  list.files(path = input_dir,
             pattern = "forw_*",
             full.names = TRUE)
input_files <- c(input_files_bkwr, input_files_forw)
input_files <- input_files[-9]
cal.chisq.stats <- NULL
for(u in 1:length(input_files)){
  load(input_files[u])
  file.nm <- print(load(input_files[u]))
  if(grepl("i=0", input_files[u]) == TRUE){
    tmp <- get(file.nm)
    datau <- tmp[[1]]$data
  }else{
    tmp <- get(file.nm)
    datau <- tmp[[1]]$data
  }
  
  datau <- datau %>% mutate(newinterval = ifelse(day %in% c(1,2), "1",
                                                 ifelse(day %in% c(3:5), "2",
                                                        ifelse(day %in% c(6:10), "3", "4"))))
  datau  <-  datau %>% group_by(newinterval, event) %>% summarize(sum.obs = sum(obs.tot), 
                                                                  sum.exp = sum(exp.tot))
  Obs <- datau %>% dplyr::select(-sum.exp) %>% spread(event, sum.obs)
  Exp <- datau %>% dplyr::select(-sum.obs) %>% spread(event, sum.exp)
  
  Obs <- as.matrix(Obs[,-1]) 
  Exp <- as.matrix(Exp[,-1]) 
  
  df <- length(unique(datau$newinterval))
  tmps <- data.frame(stats = sum((Obs-Exp)^2/Exp)/df, model = file.nm)
  cal.chisq.stats <- rbind(cal.chisq.stats, tmps)
}

calib_result <- cal.chisq.stats %>% separate(model, c("model", "U"), "\\.save") %>% spread(model, stats)
save(calib_result, file = "paperPrep/finalizedCode/figs/calibration/calib_result_table.RData")



## exclude patients have event on day 0 (add in paper)
