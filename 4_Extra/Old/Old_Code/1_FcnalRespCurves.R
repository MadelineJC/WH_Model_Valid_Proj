# For Bayesian Stats class final project, start by simulating types II and III fcnal responses

library(ggplot2)
library(dplyr)

# Modeling the basic fcnal response curves
a <- 0.1
h <- 0.2

TypeIModel <- function(R1){
  C1 <- a*R1
  return(C1)
}

TypeIIModel <- function(R2){
  C2 <- a*R2/(1 + a*h*R2)
  return(C2)
}

TypeIIIModel <- function(R3){
  C3 <- a*R3^2/(1 + a*h*R3^2)
  return(C3)
}

# Plotting the basic fcnal response curves in base
plot(TypeIModel(1:50),lty=2,xlab="Prey Density",ylab="Prey Items Consumed",main="Type I Functional Response")
plot(TypeIIModel(1:2000),lty=2,xlab="Prey Density",ylab="Prey Items Consumed",main="Type II Functional Response")
plot(TypeIIIModel(1:2000),lty=2,xlab="Prey Density",ylab="Prey Items Consumed",main="Type III Functional Response")

# Type I: Output as DF
thing <- c()
for(R1 in 1:50){
  C1 <- TypeIModel(R1)
  thing <- c(thing,C1)
}
thing=data.frame(thing)
c.vector <- c(1:50)
thing$new.col <- c.vector
colnames(thing) <- c("PreyItemsConsumed","PreyDensity")
dfI <- data.frame(thing)
dfI <- dfI[c("PreyDensity","PreyItemsConsumed")]

# Type II: Output as DF
result <- c()
for(R2 in 1:2000){
  C2 <- TypeIIModel(R2)
  result <- c(result,C2)
}
result=data.frame(result)
a.vector <- c(1:2000)
result$new.col <- a.vector
colnames(result) <- c("PreyItemsConsumed","PreyDensity")
dfII <- data.frame(result)
dfII <- dfII[c("PreyDensity","PreyItemsConsumed")]

# Type III: Output as DF
output <- c()
for(R3 in 1:2000){
  C3 <- TypeIIIModel(R3)
  output <- c(output,C3)
}
output=data.frame(output)
b.vector <- c(1:2000)
output$new.col <- b.vector
colnames(output) <- c("PreyItemsConsumed","PreyDensity")
dfIII <- data.frame(output)
dfIII <- dfIII[c("PreyDensity","PreyItemsConsumed")]

# Combining data frames but in a bird-brain way
result$new.col <- output$PreyItemsConsumed
colnames(result) <- c("PreyItemsConsumed_TII","PreyDensity","PreyItemsConsumed_TIII")
df_TII_TIII <- data.frame(result)
df_TII_TIII <- df_TII_TIII[c("PreyDensity","PreyItemsConsumed_TII","PreyItemsConsumed_TIII")]

df_All <- left_join(df_TII_TIII,thing,by=c("PreyDensity"="PreyDensity"))
colnames(df_All) <- c("PreyDensity","PreyItemsConsumed_TII","PreyItemsConsumed_TIII","PreyItemsConsumed_TI")
df_AllT <- df_AllT[c("PreyDensity","PreyItemsConsumed_TI","PreyItemsConsumed_TII","PreyItemsConsumed_TIII")]

# Plotting w/ ggplot
TI <- ggplot(dfI,aes(dfI$PreyDensity,dfI$PreyItemsConsumed))
TI <- TI + geom_line(colour="green")
TI <- TI + labs(title="Type I Functional Response",
     x="Prey Density",
     y="Prey Items Consumed")
TI

TII <- ggplot(dfII,aes(dfII$PreyDensity,dfII$PreyItemsConsumed))
TII <- TII + geom_line(colour="pink")
TII <- TII + labs(title="Type I Functional Response",
                x="Prey Density",
                y="Prey Items Consumed")
TII

TIII <- ggplot(dfIII,aes(dfIII$PreyDensity,dfIII$PreyItemsConsumed))
TIII <- TIII + geom_line(colour="blue")
TIII <- TIII + labs(title="Type I Functional Response",
                x="Prey Density",
                y="Prey Items Consumed")
TIII

# Combining plots into one ggplot plot
ggplot(df_AllT,aes(x=df_AllT$PreyDensity))+
  geom_line(aes(y=df_All$PreyItemsConsumed_TI,colour="green"))+ # Add Type I
  geom_line(aes(y=df_All$PreyItemsConsumed_TII,colour="pink"))+ # Add Type II
  geom_line(aes(y=df_All$PreyItemsConsumed_TIII,colour="blue"))+ # Add Type III
  labs(title="Type I, II, and II Functional Responses",
       x="Prey Density",
       y="Prey Items Consumed")+
  scale_color_manual(name="Legend",
                     labels=c("Type III",
                              "Type I",
                              "Type II"),
                     values=c("blue","green","pink"))
