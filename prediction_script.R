#Call and install the dependencies 
library("caret")
library("protr")
library("Boruta")
library("ranger")

#Take user input as Primary protein sequence as character
#for example, https://www.uniprot.org/uniprotkb/P0C5H3/entry
X <- "MKIIIFLIVSSLMLIGVKTDNGYLLNKATGCKVWCVINNASCNSECKLRRGNYGYCYFWKLACYCEGAPKSELWAYATNKCNGKL"
# Calculate the amino acid composition based features
aa_composition <- extractAAC(X)
dpc_composition <-  extractDC(X)
tpc_composition <- extractTC(X)
ctriad_feature <- extractCTriad(X)
compo <- extractCTDC(X)
distri <- extractCTDD(X)
transi <- extractCTDT(X)

aac<-as.data.frame(aa_composition)
dpc<-as.data.frame(dpc_composition)
tpc<-as.data.frame(tpc_composition)
ctriad<-as.data.frame(ctriad_feature)
composition <-as.data.frame(compo)
transition <- as.data.frame(transi)
distribution <- as.data.frame(distri)
colnames(aac)='Column'
colnames(dpc)='Column'
colnames(tpc)='Column'
colnames(ctriad)='Column'
colnames(composition)='Column'
colnames(transition)='Column'
colnames(distribution)='Column'
XX<-rbind(aac, dpc, tpc, ctriad, composition, transition, distribution)
XX<- t(XX)

#Load all three best performing models saved as pickle
#Load model saved as pickle
RF1 <- readRDS("best_perf_cardio_model_RF_w_comp.pickle")
RF2 <- readRDS("best_perf_neuro_model_RF_w_comp+ctd.pickle")
SVM <- readRDS("best_perf_entero_model_SVM_w_comp+cctriad.pickle")

##Result of user input
toxpred1<-predict(RF1, newdata = XX)
toxpred2<-predict(RF2, newdata = XX)
toxpred3<-predict(SVM, newdata = XX)

#Export results in csv/xlsx file
if(toxpred1==1){term1='Cardiotoxic'}else{term1='Non-Cardiotoxic'}
if(toxpred2==1){term2='Neurotoxic'}else{term2='Non-Neurotoxic'}
if(toxpred3==1){term3='Enterotoxic'}else{term3='Non-Enterotoxic'}
XXX=data.frame(Seq=X,Cardio=term1,Neuro=term2,Entero=term3)
write.csv(XXX,file = 'Result.csv')
