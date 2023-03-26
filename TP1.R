# TP1 de Linear Time Series
# (TD4 selon la numérotation des TD)

rm(list = ls())
require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
library("ggplot2")


# Question 1
donnees1 = read.csv(file="./Données_TP1/Donnees1.csv")

n = length(donnees1$XM) - 4
xm = zoo(donnees1$XM[1:n])

# Question 2
index = 1:n
data_display = as.data.frame(cbind(index,xm))
p1 = ggplot(data = data_display) + geom_line(aes(x=index,y=xm))
p1

# On observe une saisonnalité très importante