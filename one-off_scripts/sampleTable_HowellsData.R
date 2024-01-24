setwd("/Users/nikolai/data")
d.orig <- read.csv("Howell.csv")
d <- d.orig
women <- d[d$Sex == "F",]
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
d <- cbind(d[,1:4], d[,linMeasNames])
p <- matrix(nrow = 30, ncol = 4, byrow = T, data = c(
  "NORSE", "Northern Europe", "Norway", "Norway - Norse",
  "ZALAVAR", "Central Europe", "Hungary", "Hungary - Zalavar",
  "BERG", "Central Europe", "Austria", "Austria - Berg",
  "EGYPT", "North Africa", "Egypt", "Egypt - Egypt",
  "TEITA", "East Africa", "Kenya", "Kenya - Teita",
  "DOGON", "West Africa", "Mali", "Mali - Dogon",
  "ZULU", "South Africa", "Zulu", "South Africa - Zulu",
  "BUSHMAN", "South Africa", "San", "South Africa - Sān",
  "AUSTRALI", "Australia", "Australia", "Australia - Australia",
  "TASMANIA", "Australia", "Tasmania", "Australia - Tasmania",
  "TOLAI", "Melanesia", "Papua New Guinea", "Papua New Guinea - Tolai",
  "MOKAPU", "Polynesia", "Hawaii", "Hawaii - Mokapu",
  "BURIAT", "North Asia", "Siberia", "Siberia - Buriat",
  "ESKIMO", "Greenland", "Greenland", "Greenland - Eskimo",
  "PERU", "South America", "Peru", "Peru - Peru",
  "ANDAMAN", "Indonesia", "Andaman Islands", "Indonesia - Andaman",
  "EASTERI", "Polynesia", "Easter Island", "Polynesia - Easter Island",
  "ARIKARA", "North America", "USA", "USA - Arikara",
  "AINU", "East Asia", "North Japan", "Japan - Ainu",
  "NJAPAN", "East Asia", "North Japan", "Japan - North Japan",
  "SJAPAN", "East Asia", "South Japan", "Japan - South Japan",
  "HAINAN", "East Asia", "China", "China - Hainan",
  "ANYANG", "East Asia", "China", "China - Anyang",
  "ATAYAL", "East Asia", "Taiwan", "Taiwan - Atayal",
  "PHILLIPI", "Southeast Asia", "Philippines", "Phillipines - Phillipines",
  "GUAM", "Micronesia", "Guam", "Guam - Guam",
  "MORIORI", "Polynesia", "Chatham Islands", "Chatham Islands - Moriori",
  "SMAORI", "New Zealand", "New Zealand", "New Zealand - South Maori",
  "NMAORI",  "New Zealand", "New Zealand", "New Zealand - North Maori",
  "SANTACR", "North America", "USA", "USA - Santa Cruz"
))
p

total <- table(d$Population)
mss <- table(d$Population[d$Sex == "M"])
fss <- table(d$Population[d$Sex == "F"])
total - fss - mss

sample_data <- cbind(Females = fss, Males = mss, Total = total)
rownames(sample_data) <- gsub(pattern = " ", "", rownames(sample_data))
rownames(sample_data) <- p[match(rownames(sample_data), p[,1]),4]

tableprint<-xtable::xtable(sample_data,label="tab:howellsDataComposition",
                           caption=c(paste0("A table detailing the composition of Howells' linear measurement dataset used in 
                                            our empirical analysis. Elements of the table represent numbers of individuals in each 
                                            population corresponding to each 
                                            estimated column sex. Further details regarding the nature of this sample can be found in \\citep{howellsWhoWhoSkulls1995}."), 
                                     "Composition of Linear Measurement Data Used in Empirical Analysis"))
sink("~/dissertation/tables/linear_data.tex")
print(tableprint,tabular.environment="longtable", floating = F) #,width="\\textwidth")
sink()
