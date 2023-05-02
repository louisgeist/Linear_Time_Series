# Projet - Linear time series

#----------------- Partie I : Les données -------------------

rm(list = ls())

library("ggplot2")
library("zoo")
library("dplyr")
library("fUnitRoots")
# library("tseries")
library("polynom") #pour la question 5 (justification ARIMA)


# Préparation des données
source = read.csv(file = "valeurs_mensuelles_industrie_alimentaire.csv", sep=";", dec = ".")
source2 = source[4:length(source$Libellé),1:2]
colnames(source2) = c("Date","Indice")

Date_num = seq(1990,1990+397/12,1/12)
source2$Indice = as.numeric(source2$Indice)
source2 = source2[-1] # on supprime la variable "date" qui n'est pas numérique
source2 = cbind(Date_num,source2)
colnames(source2)=c("Date","Indice")

# On met deux valeurs de côté pour la prévision à la fin
source_end = slice_tail(source2, n=2)
source2 = slice(source2, 1:(length(source2$Indice)-2))

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
val_a_prevoir = as.data.frame(cbind(source_end$Date,c(source_end$Indice[2]-source_end$Indice[1],source_end$Indice[1]-source2$Indice[396])))

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
dev.print(device = png, file = "./Images_pour_rapport/acf.png", width = 600)
# ACF => q = 1

pacf(statio_series)
dev.print(device = png, file = "./Images_pour_rapport/pacf.png", width = 600)
# PACF => p = 7

q_max = 1
p_max = 7


#test_ajustement : #renvoie un tableau de Booléens, qui indique quel modèle est correctement ajusté, 
#c'est-à-dire les modèles tels que les coefficients d'ordre p et q soient tous les deux significatifs au seuil de 5% (test de Student)

test_ajustement = function(statio_series, p_max, q_max){ 
  res = matrix(NA, nrow = (p_max+1), ncol = (q_max+1))
  rownames(res) <- paste0("p=",0:p_max)
  colnames(res) <- paste0("q=",0:q_max)
  
  for(p in 0:(p_max)){
    for(q in 0:(q_max)){
      model = arima(x = statio_series,order = c(p,0,q),include.mean = FALSE)
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

tab_ajustement = test_ajustement(statio_series, p_max, q_max)
print("Liste des modèles correctement ajustés : ")
print(tab_ajustement)

modele_valide = function (series, tab_ajustement, p_max,q_max){
  validation_matrix = tab_ajustement
  
  
  for(p in 0:(p_max)){
    for(q in 0:(q_max)){
      if (is.na(tab_ajustement[p+1,q+1])==FALSE){
        
      
      if (tab_ajustement[p+1,q+1]==TRUE) {# on regarde la valité des modèles qui sont bien ajustés
        model = arima(x = statio_series, order = c(p,0,q), include.mean = FALSE)
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

# On garde seulement les modèles ARIMA(7,0,0) et ARIMA(1,0,1) pour la série différenciée à l'ordre 1

arima700 = arima(x = statio_series, order = c(7,0,0), include.mean = FALSE)
arima101 = arima(x = statio_series, order = c(1,0,1), include.mean = FALSE)

print(paste0("AIC pour ARIMA(7,0,0) : ", AIC(arima700)))
print(paste0("AIC pour ARIMA(1,0,1) : ", AIC(arima101)))

print(paste0("BIC pour ARIMA(7,0,0) : ", BIC(arima700)))
print(paste0("BIC pour ARIMA(1,0,1) : ", BIC(arima101)))

# On conserve donc le modèle ARIMA(1,0,1), car il minimise les deux critères d'information (AIC et BIC).

# 5. Modélisation ARIMA

# Il s'agit de montrer que le modèle ARMA(1,1) qu'on a pour la série différenciée est bien causal.
# Or un ARMA est causal ssi pas de racine dans le disque unité du polynôme phi

arma_causal = function(model){
  if(model$arma[1]==0){# gère le cas trivial d'un modèle MA
    return(TRUE)
  }
  else{
    phi = polynomial(c(1,-model$coef[1:(model$arma[1])]))
    racines = polyroot(phi) #les coefficients polynomiaux sont donnés dans l'ordre CROISSANT
    
    for(i in 1:length(racines)){
      if (abs(racines[i])<=1) {
        return(FALSE)
      }
    }
    return(TRUE) #si toutes les racines sont de module strictement supérieur à 1, alors ARMA causal
  }
}

arma_causal(arima101) # renvoie TRUE

# Remarque : dans ce cas précis, où le polynôme est de degré 1, il était clair que la racine est en dehors
# du disque unité, mais l'écriture de cette fonction avait pour but d'écrire une routine généralisable

# Le modèle ARMA pour la série différenciée est causal.
# Donc, par définition d'un ARIMA, la série non transformée suit un modèle ARIMA(1,1,1).


#---------------- Partie III : Prévision --------------------
# 6. Equation de la région de confiance de niveau alpha
phi = arima101$coef[1]
psi = arima101$coef[2]
sigma = 1 #A DETERMINER

mat_sigma = matrix( data = c(1,phi-psi,phi-psi,1+(phi-psi)**2) *sigma, nrow=2,ncol=2)

mat_sigma %*% solve(mat_sigma)

diag = eigen(solve(mat_sigma))

rac_inverse = (diag$vectors %*% diag(sqrt(c(diag$values)),2) %*% t(diag$vectors))

verif = rac_inverse %*% rac_inverse %*% mat_sigma
verif
