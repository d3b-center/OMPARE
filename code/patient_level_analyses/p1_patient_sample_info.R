# Patient & Sample Information
patient_sample_info <- function() {
  df1 <- data.frame(V1 = c("Subject ID", "Sex", "Age at Diagnosis (days)", "Ethnicity"), V2 = c(sampleInfo$subjectID, sampleInfo$sex, sampleInfo$age_diagnosis_days, sampleInfo$ethnicity))
  df2 <- data.frame(V3 = c("Age at Collection (days)", "Tumor Location", "Tumor Type", ""), V4 = c(sampleInfo$age_collection_days, sampleInfo$tumorLocation, sampleInfo$tumorType, ""))
  df <- cbind(df1, df2)
  return(df)
}
