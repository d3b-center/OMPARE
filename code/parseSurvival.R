############################################
#Purpose: Code to create Survival Data
#Author: Pichai Raman
#Date: 5/2/2019
############################################


parseSurvival <- function()
{
rawSurvData <- read.delim("../data/Reference/CBTTC_PullApril21.txt", stringsAsFactors = F);
formatData <- rawSurvData[,c("CBTTC.Specimen.Group.ID", "Last.Known.Clinical.Status", "Age.at.Diagnosis", "Age.at.Last.Known.Clinical.Status")];

#Remove places with no data
formatData <- unique(formatData[grepl("7316", formatData[,"CBTTC.Specimen.Group.ID"]),]);
formatData <- unique(formatData[!grepl("Not Reported", formatData[,"Age.at.Last.Known.Clinical.Status"]),]);

#Get overall survival time but remove tiems where diagnosis and last known clinical status are same
formatData[,"OS_Time"] <- as.numeric(formatData[,"Age.at.Last.Known.Clinical.Status"])-as.numeric(formatData[,"Age.at.Diagnosis"]);
formatData <- formatData[!formatData[,"OS_Time"]==0,];


#They have comma's to separate CBTTC Specimen ID's, I HATE COMMAS
formatDataNoComma <- formatData[!grepl(",", formatData[,"CBTTC.Specimen.Group.ID"]),];
formatDataComma <- formatData[grepl(",", formatData[,"CBTTC.Specimen.Group.ID"]),];
turnTS <- function(x)
{
  tmpList <- x[[1]]
  out <- strsplit(tmpList, split=",")[[1]]
  output <- data.frame(out, t(x));
  output <- output[-2];
  return(output)
}

formatDataTS <- do.call("rbind", apply(formatDataComma, FUN=turnTS, MARGIN=1));
colnames(formatDataTS) <- colnames(formatDataNoComma)
formatData <- na.omit(unique(rbind(formatDataNoComma, formatDataTS)));

#Final formatting
formatData[,"OS_Event"] <- ifelse(formatData[,"Last.Known.Clinical.Status"]=="Alive", 0, 1)
clinData <- read.delim("../data/Reference/study_view_clinical_data.txt", stringsAsFactors=F);
mapping <- read.delim("../data/Reference/mappingFile.txt", header=F, stringsAsFactors=F);

getID <- function(x)
{
  out <- strsplit(x, "_")[[1]]
  out <- paste(out[1])
  return(out);
}

mapping[,"samps"] <- sapply(mapping[,2], FUN=getID);
formatData <- merge(formatData, mapping, by.x="CBTTC.Specimen.Group.ID", by.y="samps")

formatData <- formatData[,c("V2", "OS_Time", "OS_Event")]
colnames(formatData)[1] <- "samps";
write.table(formatData, "../data/Reference/survData.txt", row.names=F, sep="\t")
}

