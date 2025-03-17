library(dplyr)
mergeddata$iscalma=mergeddata$wd==0&mergeddata$ws==0
mergeddata$datestart=as.Date(mergeddata$date-12*3600)
calmas=mergeddata[,c(43:44)]


# Contar la cantidad de calmas por día
calmas_por_dia <- calmas %>%
  group_by(datestart) %>%
  summarise(n_calmas = sum(iscalma, na.rm = TRUE))

# Ver resultado
print(calmas_por_dia)

write_csv(calmas_por_dia, "calmaspordia.csv")
