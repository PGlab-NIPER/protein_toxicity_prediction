#Call and install the dependencies 
library("caret")
library("protr")
library("Boruta")
library("ranger")

#Take user input as Primary protein sequence as character
#for example, https://www.uniprot.org/uniprotkb/P0C5H3/entry
X <- "MINKIKILFSFLALLLSFTSYAKAEDLHDKSELTDLALANAYGQYNHPFIKENIKSDEISGEKDLIFRNQGDSGNDLRVKFATADLAQKFKNKNVDIYGASFYYKCEKISENISECLYGGTTLNSEKLAQERVIGANVWVDGIQKETELIRTNKKNVTLQELDIKIRKILSDKYKIYYKDSEISKGLIEFDMKTPRDYSFDIYDLKGENDYEIDKIYEDNKTLKSDDISHIDVNLYTKKKV"
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

#Calculate PnGT based features
text=X
library(readxl)
#Grep list of amino acids one-letter codes from pattern.xlsx file to create 11 distinctive PnGT features
Pattern <- read_excel("Pattern.xlsx")
chars <- strsplit(text, "")[[1]]
combinations <- vector("list", length = length(chars) - 1)
for (i in 2:length(chars)) {
  combinations[[i - 1]] <- paste(chars[i - 1], chars[i], sep = "")
}
combinations <- unlist(combinations)


aa=matrix(NA,ncol = 12,nrow = 1)
aa=as.data.frame(aa)
colnames(aa)=c('Seq',c(1:11))
aa[,2:12]<-0

for (j in 1:length(combinations)) {
  text_to_check <- combinations[j]
  
  for (k in 1:nrow(Pattern)) {
    possible_characters= as.character(Pattern[k,2])
    is_present <- unlist(strsplit(text_to_check, '')) %in% unlist(strsplit(possible_characters, ','))
    
    sum1=sum(as.numeric(is_present))
    if(sum1==2){
      aa[1,k+1]=aa[1,k+1]+1
    }
  }
}
aa$Seq=text

aa=aa[,-1]
colnames(aa)=c('Tiny','Small','Aliphatic','Nonpolar','Aromatic','Polar','Charged','Basic','Acidic','Hydrophobic','Hydrophilic')

#Combine all features
XX=cbind(XX,aa)

#Load all three best performing models saved as pickle
#Load model saved as pickle
RF1 <- readRDS("cardio_new.pickle")
RF2 <- readRDS("neuro_new.pickle")
SVM1 <- readRDS("entero_new.pickle")

##Result of user input
which(colnames(XX)=="NA")
colnames(XX)[23]="NA."
toxpred1<-predict(RF1, newdata = XX)
toxpred2<-predict(RF2, newdata = XX)
toxpred3<-predict(SVM1, newdata = XX)

#Export results in csv/xlsx file
if(toxpred1==0){term1='Cardiotoxic'}else{term1='Non-Cardiotoxic'}
if(toxpred2==0){term2='Neurotoxic'}else{term2='Non-Neurotoxic'}
if(toxpred3==0){term3='Enterotoxic'}else{term3='Non-Enterotoxic'}
XXX=data.frame(Seq=X,Cardio=term1,Neuro=term2,Entero=term3)
write.csv(XXX,file = 'Result.csv')
