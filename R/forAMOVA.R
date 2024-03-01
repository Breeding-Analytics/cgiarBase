#Procedure amova
forAMOVA<-function(agc.dist,agc.env){
     suppressWarnings(library(vegan))
     agc.adonis <- adonis2(agc.dist ~ Pop, data=agc.env, permutations=9999)
     adonis.data <- data.frame(agc.adonis)
     adonis.data$MS <- adonis.data$SumOfSqs / adonis.data$Df
     pop.sizes <- table(agc.env$Pop)
     n1 <- mean(pop.sizes)
     # how much variance is detected at each stratification.
     #Variations  Within samples
     Sigma.withinSamples <- adonis.data["Residual", "MS"]
     #Variations  Between samples
     Sigma.betweenSamples <- (adonis.data[1, "MS"] - Sigma.withinSamples) / n1
     Sigma.total <- Sigma.withinSamples + Sigma.betweenSamples
     #Phi-pop-total
     PHI.ST <- Sigma.betweenSamples / Sigma.total
     variance=cbind(adonis.data[,c(1,2,6,4,5)],
                    Sigma=c(Sigma.betweenSamples,Sigma.withinSamples,Sigma.total),
                    PercVar=c(100*(Sigma.betweenSamples/Sigma.total),100*(Sigma.withinSamples/Sigma.total),100),
                    Phi=c(PHI.ST,NA,NA)
                    )
   variance=cbind(c("Between Groups","Within Groups","Total"),variance)
	 colnames(variance)[1]=c("source")
	 return(variance)
}