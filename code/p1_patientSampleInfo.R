# Patient & Sample Information
patientSampleInfo <- function() {
  df1 <- data.frame(c("Subject ID", "Sex", "Age at Diagnosis (days)", "Ethnicity"), c(getSubjectID(), getSex(), getAge(), getEthnicity()))
  df2 <- data.frame(c("Age at Collection (days)", "Tumor Location", "Tumor Type", ""), c(getCollectionDate(), getTumorLocation(), getTumorType(), ""))
  return(cbind(df1,df2))
}
