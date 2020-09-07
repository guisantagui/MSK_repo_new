library(knitr)
library(randomForest)


myurl <- "https://raw.githubusercontent.com/spitakiss/Data607_Pres1/master/FullCoverage.csv"
full.cov = read.csv(myurl,header=TRUE)


full.cov$men <- factor(full.cov$men)
full.cov$urban <- factor(full.cov$urban)
full.cov$private <- factor(full.cov$private)
full.cov$y <- factor(full.cov$y)

kable(head(full.cov))


kable(summary(full.cov))


attach(full.cov)
marital = relevel(marital,ref="S")

# Fit to Logistic regression 
FullcovModel = glm(y~men+urban+private+factor(marital)+age+seniority,family=binomial(link=logit))

# Model Output
summary(FullcovModel)


# Calculate Pseudo R-2
FullcovModelRed <- glm(y~1, family=binomial(link=logit))
1-(logLik(FullcovModel))/(logLik(FullcovModelRed))


rf <- randomForest(y~men+urban+private+factor(marital)+age+seniority,data=full.cov)

# Model output
print(rf)



# Calculate sensitivity and false positive measures for logit model

fity_ypos <- FullcovModel$fitted[y == 1]
fity_yneg <- FullcovModel$fitted[y == 0]

sort_fity <- sort(FullcovModel$fitted.values)

sens <- 0
spec_c <- 0

for (i in length(sort_fity):1){
        sens <- c(sens, mean(fity_ypos >= sort_fity[i]))
        spec_c <- c(spec_c, mean(fity_yneg >= sort_fity[i]))
        
} 

# Calculate sensitivity and false positive measure for random forest model

fity_ypos2 <- as.numeric(rf$pred[y == 1]) - 1
fity_yneg2 <- as.numeric(rf$pred[y == 0]) - 1

sort_fity2 <- as.numeric(sort(rf$pred)) - 1

sens2 <- 0
spec_c2 <- 0

for (i in length(sort_fity2):1){
        sens2 <- (c(sens2, mean(fity_ypos2 >= sort_fity2[i])))
        spec_c2 <- (c(spec_c2, mean(fity_yneg2 >= sort_fity2[i])))
} 

# plot ROC curves

plot(spec_c, sens, xlim = c(0, 1), ylim = c(0, 1), type = "l", 
     xlab = "false positive rate", ylab = "true positive rate", col = 'blue')
abline(0, 1, col= "black")

lines(spec_c2, sens2, col='green')
legend("topleft", legend = c("logit","random forest") , pch = 15, bty = 'n', col = c("blue","green"))