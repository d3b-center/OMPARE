# Patient & Sample Information
patientSampleInfo <- function() {
  df1 <- data.frame(c("Subject ID", "Sex", "Age", "Ethnicity"), c(getSubjectID(), getSex(), getAge(), getEthnicity()))
  df2 <- data.frame(c("Medical Facility", "Primary Physician", "Pathologist", "Lab Director"), c(getMedicalFacility(), getPrimPhysician(), getPathologist(), getLabDirector()))
  df3 <- data.frame(c("Collection Date", "Tumor Location", "Tumor Type", "P/R"), c(getCollectionDate(), getTumorLocation(), getTumorType(), getPrimRelapse()))
  return(cbind(df1,df2, df3))
}