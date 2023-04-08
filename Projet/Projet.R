# Projet - Linear time series

#----------------- Partie I : Les données -------------------

rm(list = ls())

library("ggplot2")
library("zoo")
# library("tseries")

# Préparation des données
source = read.csv(file = "valeurs_mensuelles_industrie_alimentaire.csv", sep=";", dec = ".")
source2 = source[4:length(source$Libellé),1:2]
colnames(source2) = c("Date","Indice")

Date_num = seq(1990,1990+397/12,1/12)
source2$Indice = as.numeric(source2$Indice)
source2 = cbind(Date_num,source2)



# 2. Transformation de la série pour la rendre stationnaire
p = ggplot(data=source2) + geom_line(aes(x=Date_num,y=Indice))
p

# On constate une série avec un tendance haussière. On différencie donc une première fois :
data = as.data.frame(cbind(Date_num[2:length(Date_num)],diff(source2$Indice)))

# Cette transformation de série semble stationnaire. (Cf. représentation graphique dans la prochaine question)

# 3. Représentation graphique
p

colnames(data)= c("Date","FD_indice")
p_diff = ggplot(data=data) + geom_line(aes(x=Date,y=FD_indice))+ylab("Différence première de l'indice")
p_diff

#---------------- Partie II : Modèles ARMA ------------------
# 4. Choix du modèle ARMA
acf(data$FD_indice)
pacf(data$FD_indice)

# ACF => q = 1
# PACF => p = 7


# 5. Modélisation ARIMA

#---------------- Partie III : Prévision --------------------