# Information ####################
# Nana Kwarteng
# 2024-10-01: All functions moved in re-structuring. Replaced with src statement
# 2024-09-26: standardize treat 1 and treat 2 with netmeta pairwise. TE=1-2 (treat 2 is reference)
# 2024-08-27: Add visualization for arm based
# 2024-05-20: Add contrast Based
# 2024-05-23: Add V matrix for SMDs
# Call Metaregression functions
# Functions for model matrix for network meta regression
# Functions for performing network meta regression using multivariate meta-analysis
# all_funcs<-list.files("functions/",full.names = TRUE,recursive = TRUE)
all_funcs<-list.files("functions/",full.names = TRUE,recursive = TRUE)

for(i in all_funcs){
  source(i)
}

# unmask functions
chknumeric <- netmeta:::chknumeric
chklevel <- netmeta:::chklevel
chkclass <- netmeta:::chkclass
updateversion <- netmeta:::updateversion
chklogical <- netmeta:::chklogical
setchar <- netmeta:::setchar
setref <- netmeta:::setref
catch <- netmeta:::catch
calcV <- netmeta:::calcV
runNN <- netmeta:::runNN
formatN<-netmeta:::formatN
contrmat <- metafor::contrmat
bldiag <- metafor::bldiag
rma.mv <- metafor::rma.mv
formatCI<-netmeta:::formatCI
replaceNA<-netmeta:::replaceNA
formatPT<-netmeta:::formatPT
summary.netmetareg<-summary.netmetareg_old
