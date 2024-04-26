calculo_de_blanco <- function (Blancos, metodo) {
  library(dplyr)
  library(readxl)
  library(readr) 
  library(ranger)
  library(plyr)
  switch(metodo,
         "promedio" = {
           Blanco_para_restar =  as.data.frame(t(apply(Blancos[,c(5:28)], 2, mean, na.rm=TRUE)))
           Blanco_para_restar <- Blanco_para_restar  %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
           Blanco_para_restar$K = 0 # esto es porque quiero ver que pasa con una serie con blanco no detectado
           Blanco_para_restar <- bind_rows(replicate(3, Blanco_para_restar, simplify = FALSE))
           Blanco_para_restar$Lote=as.factor(levels(Blancos$Lote))
         },
         "promedioporlote"={
           Blanco_para_restar <- aggregate(Blancos, list(as.numeric(as.character(Blancos$Lote))), FUN=mean, na.rm=TRUE)
           Blanco_para_restar <- Blanco_para_restar[,c(-31,-30)]
           Blanco_para_restar<- dplyr::rename(Blanco_para_restar, Lote=Group.1)
           Blanco_para_restar$Lote <- as.factor(Blanco_para_restar$Lote)
           Blanco_para_restar[,-1] <- Blanco_para_restar[,-1]  %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
           Blanco_para_restar <- Blanco_para_restar[,c(6:29,1)]
         },
         "promediosinmaxnimin"={
           # Blanco_para_restar <- Blancos %>% 
           #   summarise(across(5:28, ~ {
           #     x <- .[!.%in% c(min(., na.rm = TRUE), max(., na.rm = TRUE))]
           #     mean(x, na.rm = TRUE)
           #   }))
           Blanco_para_restar <- Blancos %>% 
             summarise(across(5:28, ~ {
               min_val <- min(., na.rm = TRUE)
               max_val <- max(., na.rm = TRUE)
               x <- .[! (. == min_val | . == max_val) | (duplicated(.) & (. == min_val | . == max_val))]
               mean(x, na.rm = TRUE)
             }))
           Blanco_para_restar <- lapply(Blanco_para_restar, function(x) replace(x, is.nan(x), NA))
           Blanco_para_restar <- as.data.frame(Blanco_para_restar)
           # Blanco_para_restar$K = 0 # esto es porque quiero ver que pasa con una serie con blanco no detectado
           Blanco_para_restar <- bind_rows(replicate(3, Blanco_para_restar, simplify = FALSE))
           Blanco_para_restar$Lote=as.factor(levels(Blancos$Lote))
         })
  }