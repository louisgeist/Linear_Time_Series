source = read.csv(file = "valeurs_mensuelles_industrie_alimentaire.csv", sep=";", dec = ".")
source = source[4:length(source$Libellé),]

statio_series = diff(as.numeric(source[,2]))

model = arima(x = statio_series,order = c(1,0,4),include.mean = F)
model[["var.coef"]]

# Les termes diagonaux sont négatifs.

p_vals = 2*(1-pnorm(abs(model$coef)/diag(model[["var.coef"]])**(1/2))) #test de Student
