#!/share/apps/R/3.3.3/bin/R
# letters_train <- letters[1:16000, ]
# letters_test <- letters[16001:20000,]

# library(kernlab)
# classifier <- ksvm(X231 ~ .,  data=Te,kernel = "vanilladot")
# letter_classifier

classifier = svm(formula = X231 ~ .,
                 data = Te,
                 type = 'C-classification',
                 kernel = 'linear')

# library(tibble)
pred   = predict(classifier, newdata = Te[-231], type = "response")

library(pROC)
test_roc = roc(as.numeric(as.character(Te$X231)) ~ as.numeric(as.character(pred)), plot = FALSE, print.auc = FALSE)
test_roc$auc[1]





# #predictions:
# letter_predictions <- predict(letter_classifier, letters_test)
# #Check the accuracy:
# caret::confusionMatrix(letter_predictions,letters_test$letter)
# agreement <- letter_predictions == letters_test$letter
# prop.table(table(agreement))
# #We get an accuracy of 84% with a simple linear kernel,let's try with RBF kernel:
# letter_classifier_rbf <- ksvm(letter ~ ., data = letters_train,kernel = "rbfdot")

# #predictions:
# letter_predictions <- predict(letter_classifier_rbf,letters_test)
# #Check the accuracy:
# caret::confusionMatrix(letter_predictions,letters_test$letter)
# agreement <- letter_predictions == letters_test$letter
# prop.table(table(agreement))



# cd /mnt/BioScratch/edahi
setwd('/mnt/BioScratch/edahi')
# install.packages("ISLR")
# install.packages("tibble")
# install.packages("lattice")
# install.packages("caret")
# install.packages("pROC")
# install.packages("nnet")
test.x  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_X.txt")
test.y  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_Y.txt")
train.x = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_X.txt")
train.y = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_Y.txt")

Tr = data.frame(cbind(t(train.x),t(train.y)))
Te = data.frame(cbind(t(test.x),t(test.y)))
library(tibble)

model_glm = glm(X231 ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28+X29+X30+X31+X32+X33+X34+X35+X36+X37+X38+X39+X40+X41+X42+X43+X44+X45+X46+X47+X48+X49+X50+X51+X52+X53+X54+X55+X56+X57+X58+X59+X60+X61+X62+X63+X64+X65+X66+X67+X68+X69+X70+X71+X72+X73+X74+X75+X76+X77+X78+X79+X80+X81+X82+X83+X84+X85+X86+X87+X88+X89+X90+X91+X92+X93+X94+X95+X96+X97+X98+X99+X100+X101+X102+X103+X104+X105+X106+X107+X108+X109+X110+X111+X112+X113+X114+X115+X116+X117+X118+X119+X120+X121+X122+X123+X124+X125+X126+X127+X128+X129+X130+X131+X132+X133+X134+X135+X136+X137+X138+X139+X140+X141+X142+X143+X144+X145+X146+X147+X148+X149+X150+X151+X152+X153+X154+X155+X156+X157+X158+X159+X160+X161+X162+X163+X164+X165+X166+X167+X168+X169+X170+X171+X172+X173+X174+X175+X176+X177+X178+X179+X180+X181+X182+X183+X184+X185+X186+X187+X188+X189+X190+X191+X192+X193+X194+X195+X196+X197+X198+X199+X200+X201+X202+X203+X204+X205+X206+X207+X208+X209+X210+X211+X212+X213+X214+X215+X216+X217+X218+X219+X220+X221+X222+X223+X224+X225+X226+X227+X228+X229+X230, data = Tr, family = "binomial")
model_glm_pred = ifelse(predict(model_glm, type = "link") > 0, "Yes", "No")
obs = ifelse(Tr$X231 > 0, "Yes", "No")
train_tab = table(predicted = model_glm_pred, actual = obs)

library(caret)
train_con_mat = confusionMatrix(train_tab, positive = "Yes")
c(train_con_mat$overall["Accuracy"], 
  train_con_mat$byClass["Sensitivity"], 
  train_con_mat$byClass["Specificity"])

library(pROC)
test_prob = predict(model_glm, newdata = Te, type = "response")
test_roc = roc(Te$X231 ~ test_prob, plot = FALSE, print.auc = FALSE)

cat("test2",as.character(round(test_roc$auc[1],5)),"\n",file = "output.txt", append = TRUE)













#
# #
# # #
# # # #
# # # # #
# # # # # #
# randomForest
test.x  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_X.txt")
test.y  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_Y.txt")
train.x = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_X.txt")
train.y = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_Y.txt")

Tr = data.frame(cbind(t(train.x),t(train.y)))
Te = data.frame(cbind(t(test.x),t(test.y)))
library(randomForest)
Tr$X231 = factor(Tr$X231)
rf <- randomForest(
  X231 ~ .,
  data=Tr
)
Te$X231 = factor(Te$X231)
pred = predict(rf, newdata=Te[-231], type = "response")

library(pROC)
test_roc = roc(as.numeric(as.character(Te$X231)) ~ as.numeric(as.character(pred)), plot = FALSE, print.auc = FALSE)
test_roc$auc[1]

#
# #
# # #
# # # #
# # # # #
# # # # # #
# Support Vector Machines
# install.packages('e1071')
test.x  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_X.txt")
test.y  = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Test_Y.txt")
train.x = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_X.txt")
train.y = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_Bcell_Resting_Rep1_Train_Y.txt")

Tr = data.frame(cbind(t(train.x),t(train.y)))
Te = data.frame(cbind(t(test.x),t(test.y)))
Tr$X231 = factor(Tr$X231)
Te$X231 = factor(Te$X231)

library(e1071)
classifier = svm(formula = X231 ~ .,
                 data = Tr,
                 type = 'C-classification',
                 kernel = 'linear')

pred   = predict(classifier, newdata = Te[-231], type = "response")

library(pROC)
test_roc = roc(as.numeric(as.character(Te$X231)) ~ as.numeric(as.character(pred)), plot = FALSE, print.auc = FALSE)
test_roc$auc[1]






prefix="/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets/Classes2_Practical/mm10_WholeBrain_Rep1"
test.x  = read.table(paste0(prefix,"_Test_X.txt") )
test.y  = read.table(paste0(prefix,"_Test_Y.txt") )
train.x = read.table(paste0(prefix,"_Train_X.txt"))
train.y = read.table(paste0(prefix,"_Train_Y.txt"))

Tr = data.frame(cbind(t(train.x),t(train.y)))
Te = data.frame(cbind(t(test.x),t(test.y)))
library(tibble)
model_glm = glm(X231 ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28+X29+X30+X31+X32+X33+X34+X35+X36+X37+X38+X39+X40+X41+X42+X43+X44+X45+X46+X47+X48+X49+X50+X51+X52+X53+X54+X55+X56+X57+X58+X59+X60+X61+X62+X63+X64+X65+X66+X67+X68+X69+X70+X71+X72+X73+X74+X75+X76+X77+X78+X79+X80+X81+X82+X83+X84+X85+X86+X87+X88+X89+X90+X91+X92+X93+X94+X95+X96+X97+X98+X99+X100+X101+X102+X103+X104+X105+X106+X107+X108+X109+X110+X111+X112+X113+X114+X115+X116+X117+X118+X119+X120+X121+X122+X123+X124+X125+X126+X127+X128+X129+X130+X131+X132+X133+X134+X135+X136+X137+X138+X139+X140+X141+X142+X143+X144+X145+X146+X147+X148+X149+X150+X151+X152+X153+X154+X155+X156+X157+X158+X159+X160+X161+X162+X163+X164+X165+X166+X167+X168+X169+X170+X171+X172+X173+X174+X175+X176+X177+X178+X179+X180+X181+X182+X183+X184+X185+X186+X187+X188+X189+X190+X191+X192+X193+X194+X195+X196+X197+X198+X199+X200+X201+X202+X203+X204+X205+X206+X207+X208+X209+X210+X211+X212+X213+X214+X215+X216+X217+X218+X219+X220+X221+X222+X223+X224+X225+X226+X227+X228+X229+X230, data = Tr, family = "binomial")
model_glm_pred = ifelse(predict(model_glm, type = "link") > 0, "Yes", "No")
obs = ifelse(Tr$X231 > 0, "Yes", "No")
train_tab = table(predicted = model_glm_pred, actual = obs)
library(caret)
train_con_mat = confusionMatrix(train_tab, positive = "Yes")
c(train_con_mat$overall["Accuracy"], 
  train_con_mat$byClass["Sensitivity"], 
  train_con_mat$byClass["Specificity"])

library(pROC)
test_prob = predict(model_glm, newdata = Te, type = "response")
test_roc = roc(Te$X231 ~ test_prob, plot = FALSE, print.auc = FALSE)
test_roc$auc[1]

cat("mm10_WholeBrain_Rep1",as.character(round(test_roc$auc[1],5)),"\n",file = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/21.Baselines/mm10_WholeBrain_Rep1.txt", append = TRUE)
