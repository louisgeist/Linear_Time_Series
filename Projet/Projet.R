# Projet - Linear time series

#----------------- Partie I : Les données -------------------

rm(list = ls())

library("ggplot2")
library("zoo")
library("dplyr")
library("fUnitRoots")
# library("tseries")


# Préparation des données
source = read.csv(file = "valeurs_mensuelles_industrie_alimentaire.csv", sep=";", dec = ".")
source2 = source[4:length(source$Libellé),1:2]
colnames(source2) = c("Date","Indice")

Date_num = seq(1990,1990+397/12,1/12)
source2$Indice = as.numeric(source2$Indice)
source2 = source2[-1] # on supprime la variable "date" qui n'est pas numérique
source2 = cbind(Date_num,source2)
colnames(source2)=c("Date","Indice")
source2 = dplyr::filter(source2,source2$Date<2022) #on arrête la série avant le Covid


### 2. Transformation de la série pour la rendre stationnaire

# Etape 1 : choix de la spécification du test ADF
summary(lm(source2$Indice ~ source2$Date))


# Etape 2 : choix du nombre de lags
Qtests <- function(series, k, fitdf=0) { #réalise le test de Ljung-Box pour les horizons jusqu'à 24 de la série "series" mis en argument
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

exogeneisation_residus = function (series,specification){
  lag_max = 36
  for(lag in 0:lag_max){
    my_adf = fUnitRoots::adfTest(series, lags = lag, type = specification)
    tab_p_val_autocorr = Qtests(my_adf@test$lm$residuals, 24, fitdf = length(my_adf@test$lm$coefficients))
    
    non_rejet = c((tab_p_val_autocorr[,2]> 0.05) | is.na(tab_p_val_autocorr[,2]))
    if(sum(non_rejet)==24){ #test si tous les tests de Ljung-Box sont soit non rejetés, soit n'ont pas été réalisés car le lag était trop faible pour interprétation 
      return(lag)
    }
  }
  return(paste0("Pas d'absence d'autocorrélation trouvé jusqu'à l'ajout du lag = ",lag_max))
}

lag_pour_adf_valide = exogeneisation_residus(source2$Indice,"ct")


# Etape 3 : réalisation du test ADF
test_adf = fUnitRoots::adfTest(source2$Indice, lags = lag_pour_adf_valide, type = "ct")
print(paste0("La pvaleur du test ADF est : ", round(test_adf@test$p.value,digits = 2)))



# Ccl: la série source2$indice admet une racine unitaire

### On réitère les 3 étapes  pour la série différenciée au premier ordre
d_Indice = diff(source2$Indice)
Date = source2$Date[2:(length(d_Indice)+1)]
source3 = as.data.frame(cbind(Date,d_Indice))

# Etape 1bis :
summary(lm(source3$d_Indice ~ source3$Date))

# Etape 2bis :
lag_pour_adf_valide_d = exogeneisation_residus(source3$d_Indice,"nc")

# Etape 3bis :
test_adf_d = fUnitRoots::adfTest(source3$d_Indice, lags = lag_pour_adf_valide, type = "nc")

print(paste0("La pvaleur du test ADF est : ", round(test_adf_d@test$p.value,digits = 2)))
# Ccl : la sérience source3$d_Indice n'admet pas de racine unitaire et donc est stationnaire



### 3. Représentation graphique
p = ggplot(data=source2) + geom_line(aes(x=Date,y=Indice))
p
ggsave("Serie_brute.png",path="./Images_pour_rapport",width = 10, height = 5)


p_diff = ggplot(data=source3) + geom_line(aes(x=Date,y=d_Indice))
p_diff
ggsave("Serie_differenciee.png",path="./Images_pour_rapport",width = 10, height = 5)


#---------------- Partie II : Modèles ARMA ------------------
statio_series = source3$d_Indice
# 4. Choix du modèle ARMA

acf(statio_series)
# ACF => q = 1

pacf(statio_series)
# PACF => p = 7

q_max = 1
p_max = 7


arima701 = arima(statio_series,c(1,0,1),include.mean = F)
p_vals701 = Qtests(arima701$residuals,24,fitdf = length(arima701$coef))
p_vals701
# On choisit donc un modèle ARMA(7,1) pour la série de l'indice différenciée une fois.

# On veut tester l'ajustement de tous les sous-modèles de l'ARIMA(7,0,1).
modele = arima701


#test_ajustement : #renvoie un tableau de Booléens, qui indique quel modèle est correctement ajusté, 
#c'est-à-dire les modèles tels que les coefficients d'ordre p et q soient tous les deux significatifs au seuil de 5% (test de Student)

test_ajustement = function(statio_series, p_max, q_max){ 
  res = matrix(NA, nrow = (p_max+1), ncol = (q_max+1))
  rownames(res) <- paste0("p=",0:p_max)
  colnames(res) <- paste0("q=",0:q_max)
  
  for(p in 0:(p_max)){
    for(q in 0:(q_max)){
      model = arima(x = statio_series,order = c(p,0,q),include.mean = F)
      p_vals = 2*(1-pnorm(abs(model$coef)/diag(model[["var.coef"]])**(1/2))) #test de Student
      
      if(p==0 | q==0){
        if(p!=0 | q!=0){ # on laisse p=0, q=0 en NA.. 
          if(p==0){ # si p ou q est nul, alors on ne peut évidemment que tester la significativité de l'autre coefficient
            if(p_vals[q]<0.05){
              res[1,q+1]=TRUE
            }
            else{
              res[1,q+1]=FALSE
            }
          }
          if(q==0){
            if(p_vals[p]<0.05){
              res[p+1,1]=TRUE
            }
            else{
              res[p+1,1]=FALSE
            }
          }
        }
      }
      else{
        if(p_vals[p]<0.05 & p_vals[p+q]<0.05){
          res[p+1,q+1] = TRUE
        }
        else {res[p+1,q+1] = FALSE}
      }
    }
  }
 return(res)
}

tab_ajustement = test_ajustement(statio_series,p_max,q_max)
print("Liste des modèles correctement ajustés : ")
print(tab_ajustement)



modele_valide = function (series, tab_ajustement, p_max,q_max){
  validation_matrix = tab_ajustement
  
  
  for(p in 0:(p_max)){
    for(q in 0:(q_max)){
      if (is.na(tab_ajustement[p+1,q+1])==FALSE){
        
      
      if (tab_ajustement[p+1,q+1]==TRUE) {# on regarde la valité des modèles qui sont bien ajustés
        model = arima(x = statio_series, order = c(p,0,q), include.mean = F)
        tab_p_val_autocorr = Qtests(model$residuals, 24, fitdf = p+q)
        
        non_rejet = c((tab_p_val_autocorr[,2]> 0.05) | is.na(tab_p_val_autocorr[,2]))
        
        if(sum(non_rejet)!=24){ #test si tous les tests de Ljung-Box sont soit non rejetés, soit n'ont pas été réalisés car le lag était trop faible pour interprétation 
          validation_matrix[p+1,q+1] = FALSE
        }
        # sinon, on laisse TRUE dans le tableau validation_matrix
      }
      }
    }
  }
  return(validation_matrix)
}

tab_valide = modele_valide(statio_series,tab_ajustement,p_max,q_max)
print("Parmi les modèles ajustés, voici le tableau des modèles valides (c'est-à-dire tels que leurs résidus ne sont pas autocorrélés)")
print(tab_valide)


# Conclusion : on garde seulement les modèles ARIMA(7,0,0) et ARIMA(0,0,1)


# 5. Modélisation ARIMA

# Ce qui suit provient du TP1_Données2 corrigé par Saurel et c'est ce qui est attendu de nou
#install.packages("polynom")
library("polynom")

ar2ma1 = arima(dxm-mean(dxm),c(2,0,1))
print(ar2ma1)


model = ar2ma1
phi = polynomial(c(1,-model$coef[1:(model$arma[1])]))
racines = polyroot(phi)

print(Mod(racines))

print(phi)

#votre_pol = polynomial()
#racines = solve()
#mod(racines)
polyroot(c(1,-1.22,0.25))


#---------------- Partie III : Prévision --------------------