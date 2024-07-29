########################################################################
# Machine Learning Methods in R: Boot Camp
# July 30, 2024
########################################################################

# Do this only once:
# install.packages(c("caTools", "dendextend", "ggbeeswarm", "ggplot2", "pheatmap", "pROC", "randomForest", "readr", "rpart", "rpart.plot", "Rtsne"))

library(caTools)
library(dendextend)
library(ggbeeswarm)
library(ggplot2)
library(pheatmap)
library(pROC)
library(randomForest)
library(readr)
library(rpart)
library(rpart.plot)
library(Rtsne)


########################################################################

# Set your own path - different from mine

setwd('~/Desktop/bioinfo-tutorials/full-course-2024-ML-workshop/')

# Explore association between gene expression vs. sex and/or smoking status

library(readr)

df <- read_csv("df_ml.2024.csv")
df = data.frame(df)

head(df)


boxplot( geneA ~ sex  , data = df )
boxplot( geneA ~ smoke  , data = df ) # looks like a strong association
boxplot( geneA ~ sex + smoke  , data = df) # more complexity with 2 variables


# Visualize a linear regression models

library(ggplot2)

ggplot(data = df, aes(x = age, y=geneA)) + geom_point()

ggplot(data = df, aes(x= smoke, y = geneA)) + 
  geom_boxplot() + 
  facet_grid( . ~ sex) 

library(ggbeeswarm)
ggplot(data = df, aes(x= smoke, y = geneA, color = smoke)) +
  geom_beeswarm(priority ='random') +
  facet_grid( . ~ sex) 

ggplot(data = df, aes(x = age, y=geneA)) + geom_point() + 
  geom_smooth(method = "lm", se = 0)

ggplot(data = df, aes(x = age, y=geneA, color = smoke)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = 0)

ggplot(data = df, aes(x = age, y= geneA, color = smoke)) + 
  geom_point(alpha = .8) + 
  facet_grid( ~ sex) + 
  geom_smooth(method = "lm", se = 0)

########################################################################

# Fit a linear regression model for gene expression 

lm( geneA ~ age + smoke, data = df) # lower expression with more smoking  

mod = lm( geneA ~ age + smoke, data = df ) # save the fitted model as an object

summary(mod) # more information about the fitness and significance

coef(mod)

# summary(mod) is an object with its own components e.g. $coefficients
summary(mod) $ coefficients      
coef(summary(mod))

# ... and those components are also objects e.g. $coefficients is a matrix
# so the p-values can be extracted as its last 4th column
summary(mod) $ coefficients [,4] 
summary(mod) $ coefficients [,'Pr(>|t|)'] 


fitted.values(mod)
residuals(mod)


# Looking at Adjusted R-squared and the significance of coefficients:
# Adjusted R-squared tries to balance better fit with fewer variables 

colnames(df)

summary( lm( geneA ~ . , data = df ))

summary( lm( geneA ~ age + smoke + sex + geneB + geneC + class , data = df )) 
summary( lm( geneA ~ age + smoke + sex + geneB +       + class , data = df )) 
summary( lm( geneA ~ age + smoke + sex +               + class , data = df )) 
summary( lm( geneA ~ age + smoke + sex                         , data = df )) 
summary( lm( geneA ~ age + smoke                               , data = df )) 
summary( lm( geneA ~ age                                       , data = df )) 


summary( lm( geneA ~ age + smoke + sex + geneB + geneC + class , data = df )) $ adj.r.squared
summary( lm( geneA ~ age + smoke + sex + geneB +       + class , data = df )) $ adj.r.squared
summary( lm( geneA ~ age + smoke + sex +               + class , data = df )) $ adj.r.squared
summary( lm( geneA ~ age + smoke + sex                         , data = df )) $ adj.r.squared
summary( lm( geneA ~ age + smoke                               , data = df )) $ adj.r.squared
summary( lm( geneA ~ age                                       , data = df )) $ adj.r.squared


########################################################################
# Prediction

library(caTools)

set.seed(555)
split = sample.split(df$geneA, SplitRatio = 0.8)
split

train_set = subset(df, split == TRUE)
test_set  = subset(df, split == FALSE)

# First use THE SAME training set for model building and for prediction

mod = lm(geneA ~ age, data = train_set)
summary(mod)

fitted.values(mod)                # fitted values from model building
predict(mod, newdata = train_set) # more general approach

#Visualize - use THE SAME training set

df_pred = cbind(train_set, y_pred = predict(mod, newdata = train_set))

ggplot(data = df_pred, aes(x = age, y=geneA)) + geom_point(col = 'blue')

ggplot(data = df_pred, aes(x = age, y=y_pred)) + geom_point(col = 'red')

ggplot(data = df_pred, aes(x = age, y=geneA)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = df_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 


#Visualize - test set

df_pred = cbind(test_set, y_pred = predict(mod, newdata = test_set))

ggplot(data = df_pred, aes(x = age, y=geneA)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = df_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 


###############################################################################
# Beware of incompatible training and test sets

train_set = subset(df, smoke == 'No')
test_set  = subset(df, smoke == 'Yes')

mod = lm(geneA ~ age, data = train_set)

df_pred = cbind(test_set, y_pred = predict(mod, newdata = test_set))

ggplot(data = df_pred, aes(x = age, y=geneA)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = df_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 


# check on full data
ggplot(data = df, aes(x = age, y=geneA)) +  geom_point() +  facet_grid( ~ smoke) 

###############################################################################
# Tree-based  models

ggplot(data = df, aes(x = age, y=geneB )) +  geom_point() 

set.seed(555)
split = sample.split(df$geneB, SplitRatio = 0.8)
train_set = subset(df, split == TRUE)
test_set = subset(df, split == FALSE)

library(rpart)
library(rpart.plot)

# regression model
mod = rpart(geneB ~ age, data = train_set)

mod
plot(mod)
rpart.plot(mod)

# on training set

df_pred = cbind(train_set, y_pred = predict(mod, newdata = train_set))

ggplot(data = df_pred, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = df_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 

# on test set

df_pred = cbind(test_set, y_pred = predict(mod, newdata = test_set))

ggplot(data = df_pred, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = df_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 

# on a grid

grid_set = data.frame(age = seq(0, 50, by = .1))
grid_pred = cbind(grid_set, y_pred = predict(mod, newdata = grid_set))

ggplot(data = train_set, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_smooth(method = 'lm', col = 'green', se = 0) +
  geom_point(data = grid_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) 


###############################################################################
# More tree-based  models

mod = rpart(geneA ~ age + smoke + sex, data = train_set)
rpart.plot(mod)

mod = rpart(geneA ~ age + smoke + sex, data = train_set, control = rpart.control(minsplit = 100))
rpart.plot(mod)

mod = rpart(geneB ~ age + smoke + sex, data = train_set, control = rpart.control(maxdepth = 30))
rpart.plot(mod)

# classification tree - predict a categorical variable

mod = rpart(smoke ~ geneA + age + sex, data = train_set)
rpart.plot(mod)

mod = rpart(class ~ geneA + age + sex + smoke, data = train_set, control = rpart.control(minsplit = 20))
rpart.plot(mod)

###############################################################################
# Random forest

library(randomForest)

set.seed(555)
split = sample.split(df$geneB, SplitRatio = 0.8)
train_set = subset(df, split == TRUE)
test_set  = subset(df, split == FALSE)

mod = randomForest(geneB ~ age, data = train_set, ntree = 10, maxnodes = 3)

# on grid

grid_set = data.frame(age = seq(0, 50, by = .1))
grid_pred = data.frame(grid_set, y_pred = predict(mod, newdata = grid_set))

ggplot(data = train_set, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_point(data = grid_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) + 
  geom_smooth(method = 'lm', col = 'green', se = 0)


# try ntree = 1, maxnodes = 300
# try ntree = 300, maxnodes = 3


###############################################################################
# Loess regression

mod = loess(geneB ~ age, data = train_set)

# on grid

grid_set = data.frame(age = seq(0, 50, by = .1))
grid_pred = data.frame(grid_set, y_pred = predict(mod, newdata = grid_set))

ggplot(data = train_set, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_point(data = grid_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) + 
  geom_smooth(method = 'lm', col = 'green', se = 0)

ggplot(data = train_set, aes(x = age, y=geneB)) + 
  geom_point(col = 'blue') +  
  geom_point(data = grid_pred, aes(x = age, y=y_pred), col = 'red', pch = 1) + 
  geom_smooth( col = 'green', se = 0)


#################################################################
# Logistic regression

df$class
df$class = factor(df$class, levels = c('Benign', 'Malignant'))

set.seed(555)
split = sample.split(df$class, SplitRatio = 0.5)
train_set = subset(df, split == TRUE)
test_set = subset(df, split == FALSE)

# tree

mod = rpart(class ~ age + sex + smoke + geneA, data = train_set, control = rpart.control(minsplit = 20))
rpart.plot(mod)

test_malignant = test_set$class == 'Malignant'
pch = ifelse (test_malignant, 2, 1)
col = ifelse (test_malignant, 'red', 'blue')

plot(test_set[,c('age', 'geneA')], pch = pch, col = col)



# logistic regession model

mod = glm(class ~ age + sex + smoke + geneA,  data = train_set, family = binomial)
summary(mod)

y_pred_link = predict(mod, newdata = test_set, type = 'link')
y_pred_prob = predict(mod, newdata = test_set, type = 'response')

plot(y_pred_link, y_pred_prob, pch = pch, col = col)
points(y_pred_link, test_malignant, pch = '|', col = col)
abline(h=.5, col = 'green')
abline(v=0, col = 'green')


###############################################################################
# estimate the performance metrics


thresh = .5

y_pred_class = factor( ifelse(y_pred_prob > thresh, 'Malignant', 'Benign'), 
                       levels = levels (df$class))

cm = table(predicted = y_pred_class, true = test_set$class)  # confusion matrix
cm

sensitivity = round( cm[2,2]/sum(cm[,2]), d=3)
specificity = round( cm[1,1]/sum(cm[,1]), d=3)
precision   = round( cm[2,2]/sum(cm[2,]), d=3)
recall      = round( cm[2,2]/sum(cm[,2]), d=3)
F1          = round( 2 * precision * recall / (precision + recall), d=3)

cat(thresh, sensitivity, specificity, precision, recall, F1, '\n', sep = '\t')

plot(y_pred_link, y_pred_prob, pch = pch, col = col)
points(y_pred_link, test_malignant, pch = '|', col = col)
abline(h=thresh, col = 'green')


# Now try all thresh in  seq(0, 1, by = .02)) 



###############################################################################
# estimate the area under ROC

library(pROC)

result.roc <- roc(test_set$class, y_pred_prob)#, levels=levels(test_set$class)) 
plot(result.roc)
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")

auc(result.roc)

# Now remove age and repeat

########################################################################
# END OF SCRIPT

