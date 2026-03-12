source("http://thegrantlab.org/misc/cdc.R")
View(cdc)
head(cdc$height)
tail(cdc$height, n = 20)

plot(x=cdc$height, y=cdc$weight, col="purple", lwd=3, xlab="height(inches)", ylab="weight(pounds)")
cor(cdc$height, cdc$weight)
hist(cdc$height, xlab = "height(inches)")

weight_kg <- cdc$weight * 0.454
height_m <- cdc$height * 0.0254
BMI <- (weight_kg) / (height_m^2)
plot(x=cdc$height, y=BMI, col="pink", lwd=2, xlab="height(meters)", ylab="BMI(kg/m^2)")
cor(cdc$height, BMI)

head(BMI >= 30, 100)
sum(BMI >= 30)
eg <- c(TRUE, TRUE, FALSE, FALSE)
sum(eg)
sum(BMI >= 30)
sum(BMI >= 30) / length(BMI)

plot(x=cdc[1:100,]$height, y=cdc[1:100,]$weight, col="maroon", lwd=3, xlab="height(inches)", ylab="weight(pounds)")

gender <- cdc[,9]
summary(gender)
df <- data.frame(gender, BMI >= 30)
df <- df[BMI >= 30,]
male_df <- subset(df, gender == "m")
table(male_df)
