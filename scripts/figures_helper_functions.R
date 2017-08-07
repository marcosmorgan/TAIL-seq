rename_cell <- function(cell){
  type <- if(grepl("bm",   cell)){"BM"
  } else if(grepl("escs",  cell)){"ESCs"
  } else if(grepl("mefs",  cell)){"MEFs"
  } else if(grepl("liver", cell)){"Liver"
  } else if(grepl("Liver", cell)){"liver"
  } else if(grepl("GV",    cell)){"gv"
  } else if(grepl("ESCs",  cell)){"escs"
  } else if(grepl("MEFs",  cell)){"mefs"
  } else if(grepl("BM",    cell)){"bm"
  } else if(grepl("gv",    cell)){"GV"
  } else {"Others"}
  return(type)
}

rename_treatment <- function(cell, treatment){
  change <- if(grepl("ctrl", treatment)){"CTL"
  } else if(grepl("dcko", treatment) & grepl("gv", cell)){"cKO"    
  } else if(grepl("dcko", treatment) & grepl("GV", cell)){"cKO"  
  } else if(grepl("dcko", treatment) & !grepl("gv", cell) & !grepl("GV", cell)){"iKO"
  } else if(grepl("CTL",  treatment)){"ctrl"
  } else if(grepl("iKO",  treatment)){"dcko"
  } else if(grepl("cKO",  treatment)){"dcko"
  } else {"other"}
  return(change)
}  
  
rename_Ts <- function(mod){
  type <-if("T"     == mod){"U"
  } else if("TT"    == mod){"UU"
  } else if("TTT"   == mod){"UUU"
  } else if("TTTT"  == mod){"UUUU"
  } else if("TTTTT" == mod){"UUUUU"
  } else {"other"}
  return(type)
}

correct_recovery <- function(a_length){
  corrected = 2**-((a_length*-0.018828)-0.783679)
  return(corrected)
}

missing_oligo_legend <- function(reads_stats){
        reads_summary_temp <- reads_stats %>% filter(type == "A") %>% select(-type, -trans, -gen)       %>%
                              mutate(monoT  = "monoT", monoC  = "monoC", monoG  = "monoG",
                                          oligoT = "oligoT", oligoC = "oligoC", oligoG = "oligoG")      %>%
                              gather(type, typeDup, c(monoT, monoC, monoG, oligoT, oligoC, oligoG))     %>%
                              select(-typeDup) %>% anti_join(reads_stats)                               %>% 
                              mutate(trans = 0, gen = 0)                                                %>%
                              bind_rows(reads_stats)
        return(reads_summary_temp)
}

missing_oligo <- function(reads_stats){
  reads_summary_temp <- reads_stats %>% filter(type == "A") %>% select(-type, -frequency)       %>%
                        mutate(monoT  = "monoT", monoC  = "monoC", monoG  = "monoG",
                               oligoT = "oligoT", oligoC = "oligoC", oligoG = "oligoG")         %>%
                        gather(type, typeDup, c(monoT, monoC, monoG, oligoT, oligoC, oligoG))   %>%
                        select(-typeDup) %>% anti_join(reads_stats) %>% mutate(frequency = 0)   %>%
                        bind_rows(reads_stats)
  return(reads_summary_temp)
}

missing_T <- function(reads_stats){
  reads_summary_temp <- reads_stats %>% filter(type_short == "A") %>% select(-type_short, -frequency)  %>%
    mutate(T  = "T", C  = "C", G  = "G")  %>% gather(type_short, typeDup, c(T, C, G))                  %>%
    select(-typeDup) %>% anti_join(reads_stats) %>% mutate(frequency = 0)                              %>%
    bind_rows(reads_stats)
  return(reads_summary_temp)
}

pval <- function(dcko, ctrl){
  res <- t.test(dcko, ctrl, var.equal=TRUE, alternative = 'less')$p.value
  return(res)
}

trypval <- function(dcko, ctrl){
  out <- tryCatch(
    {pval(dcko, ctrl)},
    error=function(cond){return(NA)}
  )
  return(out)
}

pval_2 <- function(dcko, ctrl){
  res <- t.test(dcko, ctrl, var.equal=TRUE)$p.value
  return(res)
}

trypval_2 <- function(dcko, ctrl){
  out <- tryCatch(
    {pval_2(dcko, ctrl)},
    error=function(cond){return(NA)}
  )
  return(out)
}

do_chi <- function(a, b, c, d){
  cont <- matrix(c(a, b, c-a, d-b), ncol = 2)
  pval <- chisq.test(cont)$p.value
  return(pval)
}

sig <- function(i){
  res <- "NS"
  if(i < 0.05) res = "*"  
  if(i < 0.01) res = "**" 
  if(i < 0.001) res = "***" 
  return( res)
}
