###########################################################################
## AGE-DEPENDANT SEIRS MODEL WITH 5-YEAR AGE CLASSES USING UN DEMOG DATA ##
## Modified for Thailand by sompob@tropmedres.ac (8 May 2020)
###########################################################################

library("deSolve")
library("dplyr")


#########  INCIDENCE DATA
incdata_X<-read.csv("Thcovidcases.csv")
incdata_X[,1]<-as.Date(incdata_X[,1],"%d-%m-%y")


########## POP Structure
load('THpopstruct.RData')
load('mort_sever_default.Rda')

popstruc <- THpop %>% 
  select(age_category, pop) %>% 
  rename(agefloor = age_category) %>% 
  as.data.frame()

popbirth <- THpop %>% 
  select(age_category, birth) %>% 
  as.data.frame() # unit should be per person per day

mort <- THpop %>% 
  pull(death) # unit should be per person per day

ihr <- mort_sever_default %>% 
  select(age_category, ihr) %>% 
  as.data.frame()

ifr <- mort_sever_default %>% 
  select(age_category, ifr) %>% 
  as.data.frame()


##########   CONTACT DATA
load('THcontacts.RData')


#########    POP AGEING
# per year ageing matrix
A <- length(age_categories)
dd<-seq(1:A)/seq(1:A)
ageing <- t(diff(diag(dd),lag=1)/(5*365.25))
ageing<-cbind(ageing,0*seq(1:A)) # no ageing from last compartment


#########   INITIALISE SIMULATION/INTERVENTION START TIMES
startdate<-as.Date("2020-02-15") 
# stopdate<-Sys.Date() # today
stopdate<-as.Date("2020-12-31")
# stopdate<-as.Date("2020-03-18")
day_start <- as.numeric(startdate-startdate)
day_stop <- as.numeric(stopdate-startdate)
times <- seq(day_start, day_stop)

tin<-as.numeric(startdate-as.Date("2020-01-01"))/365.25
initP<-sum(popstruc[,2])       # population size 
ageindcase<-20                 # age of index case (years)
aci <- floor((ageindcase/5)+1) # age class of index case


######    THESE ARE JUST VARIABLE DEFINITIONS - PLEASE DO NOT CHANGE   #################################
######     GO TO LINE 652 TO Choose Interventions' start dates         #################################

# date to start the self isolation intervention
date_selfis_on<-as.Date("2021-12-15")
# date to start the social distancing intervention
#date_dist_on<-as.Date("2020-03-17")
date_dist_on<-as.Date("2021-12-15")
# date to start the handwashing intervention
# date_hand_on<-as.Date("2020-02-01")
date_hand_on<-as.Date("2021-12-31")
# date to start the working from home
#date_work_on<-as.Date("2020-03-19")
date_work_on<-as.Date("2021-12-15")
# date to start the school closure
#date_school_on<-as.Date("2020-03-23")
date_school_on<-as.Date("2021-12-15")
# date to start cocooning the elderly
date_cocoon_on<-as.Date("2021-12-14")
# date to start vaccinating
date_vaccine_on<-as.Date("2021-12-31")
# date to start international travel ban
date_travelban_on<-as.Date("2021-12-31")
# date to start screening
date_screen_on<-as.Date("2021-12-21")
# date to start voluntary quarantine
#date_quarantine_on<-as.Date("2020-03-15")
date_quarantine_on<-as.Date("2021-12-19")
# date to start lockdown low 
date_lockdown_low_on<-as.Date("2021-12-31")
# date to start lockdown mid 
# date_lockdown_mid_on<-as.Date("2020-03-23")
date_lockdown_mid_on<-as.Date("2021-12-18")
# date to start lockdown high 
date_lockdown_high_on<-as.Date("2021-12-20")


#############   DEFINE PARAMETERS
parameters <- c(
  # Transmission instrinsic
  p=0.042,           # probabilty of infection given a contact 
  rho = 25,          # relative infectiousness of incubation phase (%) min 0 max 100 step 0.5 
  omega=200,         # average duration of immunity (years) min 0.5 max 100 step 0.5 
  gamma=3.5,         # average incubation period (days) min 1 max 7 step 0.5 
  nui=4.5,             # average duration of symptomatic infection period (days) min 1 max 7 step 0.5
  report=0,          # percentage of all asymptomatic infections that are reported (%) min 0 max 100 step 1
  reportc=0,         # percentage of all symptomatic infections that are reported (%) min 0 max 100 step 1
  reporth=95,        # percentage of all infections requiring hospitalisation that are actually admitted to hospital (%) min 0 max 100 step 1
  beds_available = sum(popstruc[,2])*2.54/1000,#80000, # maximum number of hospital beds - numeric 
  icu_beds_available = sum(popstruc[,2])*6.6/10000,#8000, # maximum number of hospital beds - numeric 
  ventilators_available = 10000, # maximum number of ventilators - numeric
  give = 85 ,        # system capacity stressor  
  pdeath_h = 30,     # probability of dying when hospitalised 
  pdeath_hc = 35,    # probability of dying when denied hospitalisation 
  pdeath_icu = 50,   # probability of dying when admitted to ICU 
  pdeath_icuc = 75,  # probability of dying when admission to ICU denied 
  pdeath_vent = 75,  # probability of dying when ventilated 
  pdeath_ventc = 85, # probability of dying when ventilator denied 
  ihr_scaling = 1,   # scaling factor for infection hospitalisation rate
  nus = 10,          # duration of non-fatal hospitalised infection (days) min 1 max 20 step 0.5
  nusc = 10,         # duration of non-fatal denied hospitalisation infection (days) min 1 max 20 step 0.5
  nu_icu = 10,       # duration of non-fatal icu infection (days) min 1 max 20 step 0.5
  nu_icuc = 10,      # duration of non-fatal denied icu infection (days) min 1 max 20 step 0.5
  nu_vent = 10,      # duration of non-fatal ventilated infection (days) min 1 max 20 step 0.5
  nu_ventc = 10,     # duration of non-fatal denied ventilation infection (days) min 1 max 20 step 0.5
  rhos= 5,           # relative level of contacts from severely ill patients (%) min 0 max 100 step 1
  amp=80,            # relative amplitude of seasonal forcing (%) min 0 max 100 step 1
  phi=12,            # month of peak in seasonal forcing
  pclin=15,          # probability upon infection of developing clinical symptoms
  prob_icu = 70,     # probability upon hospitalisation of requiring icu admission   
  prob_vent = 80,    # probability upon admission to the UCI of requiring a ventilator
  # INTERVENTIONS
  # self isolation
  selfis_on=as.numeric(date_selfis_on-startdate),
  selfis_dur=16,    # duration of self-isolation protocol (weeks) min 1 max 52 step 1
  selfis_cov=50,    # coverage of self isolation (%) min 0 max 100 step 1
  selfis_eff=50,    # adherence to self isolation (%) min 0 max 100 step 1
  # social distancing
  dist_on=as.numeric(date_dist_on-startdate),
  dist_dur=12,      # duration of social distancing protocol (weeks) min 1 max 52 step 1
  dist_cov=50,      # coverage of social distancing (%) min 0 max 100 step 1
  dist_eff=100,     # adherence to social distancing (%) min 0 max 100 step 1
  # hand washing
  hand_on=as.numeric(date_hand_on-startdate),
  hand_dur=12,      # duration of increased hand hygiene protocol (weeks) min 1 max 52 step 1
  hand_eff=5,       # efficacy of hand hygiene  (%) min 0 max 100 step 1
  # working at home
  work_on=as.numeric(date_work_on-startdate),
  work_dur=12,      # duration of working from home protocol (weeks) min 1 max 52 step 1
  work_cov=50,      # coverage of working from home (%) min 0 max 100 step 1
  work_eff=85,      # efficacy of working from home (%) min 0 max 100 step 1
  w2h = 10,         # work contacts that get attibuted to home when working from home (%) min 0 max 100 step 1
  # school closures
  school_on=as.numeric(date_school_on-startdate),
  school_dur=12,    # duration of school closure (weeks) min 1 max 52 step 1
  school_eff=90,    # efficacy of school closure (%) min 0 max 100 step 1
  s2h = 20,         # school contacts that get attibuted to home when school closes (%) min 0 max 100 step 1
  # cocooning the elderly
  cocoon_on = as.numeric(date_cocoon_on-startdate), 
  cocoon_dur=16,    # duration of elderly cocoon protocol (weeks) min 1 max 52 step 1
  cocoon_eff=35,    # efficacy of elderly cocoon (%) min 0 max 100 step 1
  cocoon_cov=75,    # coverage of elderly cocoon (%) min 0 max 100 step 1
  age_cocoon=70,    # minimum age for elderly cocoon min 0 max 100 step 5
  # vaccination
  vaccine_on= as.numeric(date_vaccine_on-startdate),
  vaccine_eff=90,   # vaccine efficacy (%)- min 0 max 100 step 1
  vaccine_cov=0,    # vaccine coverage (%)- min 0 max 100 step 1
  vac_campaign = 4, # Number of weeks it takes to reach maximum coverage - min 1 max 8 step 1
  # imported cases 
  mean_imports = 0,           # user defined - mean number of infectious migrants per day (number) - min 0 max 500 step 1
  travelban_on= as.numeric(date_travelban_on-startdate),
  travelban_dur = 16,         # duration of internation travel restrictions (weeks) - min 1 max 52 step 1
  travelban_eff=50,           # travel restriction efficacy (%) - min 0 max 100 step 1
  # screening - increases the rate of isolation of infectious people in the model
  screen_on = as.numeric(date_screen_on-startdate), 
  screen_dur = 12,            # duration of intensified screening (week) - min 1 max 52 step 1
  screen_cov = 90,            # sensitivity of screening test min 25 max 100 step 1
  screen_overdispersion = 4,  # overdispersion of cases around index case. If  1 likelihood same as general population min 1 max 5 step 0.2 
  screen_contacts = 4,        # number of contacts screened per index case min 1 max 10 step 1
  # quarantine - This is the bi-product of increasing testing of suspected cases with a certain false positivity rate and voluntary home quarantining of people sharing a house with an infectious case
  quarantine_on = as.numeric(date_quarantine_on-startdate),
  quarantine_cov = 70,        # coverage of quarantine (%)- min 0 max 100 step 1
  quarantine_dur = 24,        # duration of quarantine (weeks) - min 1 max 52 step 1
  quarantine_days = 14,       # days in isolation for average person (days)  - min 5 max 21 step 1
  quarantine_effort = 2,      # days to implement maximum quarantine coverage - min 1 max 5
  quarantine_eff_home = 50,   # increase in the number of contacts at home when quarantined (%) - min 0 max 100 step 52
  quarantine_eff_other = 90,  # reduction in the number of other contacts when quarantined (%) - min 0 max 100 step 52
  # lockdown
  lockdown_low_on=as.numeric(date_lockdown_low_on-startdate),
  lockdown_low_dur = 16,
  lockdown_mid_on=as.numeric(date_lockdown_mid_on-startdate),
  lockdown_mid_dur = 16,
  lockdown_high_on=as.numeric(date_lockdown_high_on-startdate),
  lockdown_high_dur = 16,
  # mean household size
  household_size = 2          # mean household size (number) - min 1 max 10 step 1 
)


# Scale parameters to percentages/ rates
parameters["rho"]<-parameters["rho"]/100
parameters["omega"]<-(1/(parameters["omega"]*365))
parameters["gamma"]<-1/parameters["gamma"]
parameters["nui"]<-1/parameters["nui"]
parameters["report"]<-parameters["report"]/100
parameters["reportc"]<-parameters["reportc"]/100
parameters["reporth"]<-parameters["reporth"]/100
parameters["nus"]<-1/parameters["nus"]
parameters["rhos"]<-parameters["rhos"]/100
parameters["amp"]<-parameters["amp"]/100
parameters["selfis_dur"]<-parameters["selfis_dur"]*7
parameters["selfis_cov"]<-parameters["selfis_cov"]/100
parameters["selfis_eff"]<-parameters["selfis_eff"]/100
parameters["dist_dur"]<-parameters["dist_dur"]*7
parameters["dist_cov"]<-parameters["dist_cov"]/100
parameters["dist_eff"]<-parameters["dist_eff"]/100
parameters["hand_dur"]<-parameters["hand_dur"]*7
parameters["hand_eff"]<-parameters["hand_eff"]/100
parameters["work_dur"]<-parameters["work_dur"]*7
parameters["work_cov"]<-parameters["work_cov"]/100
parameters["work_eff"]<-parameters["work_eff"]/100
parameters["w2h"]<-parameters["w2h"]/100
parameters["school_dur"]<-parameters["school_dur"]*7
parameters["schoolcov"]<-parameters["schoolcov"]/100
parameters["school_eff"]<-parameters["school_eff"]/100
parameters["s2h"]<-parameters["s2h"]/100
parameters["cocoon_dur"]<-parameters["cocoon_dur"]*7
parameters["cocoon_cov"]<-parameters["cocoon_cov"]/100
parameters["cocoon_eff"]<-parameters["cocoon_eff"]/100
parameters["age_cocoon"]<-floor((parameters["age_cocoon"]/5)+1)
parameters["travelban_eff"]<-parameters["travelban_eff"]/100
parameters["vaccine_eff"]<-parameters["vaccine_eff"]/100
parameters["vaccine_cov"]<-parameters["vaccine_cov"]/100
parameters["vac_campaign"]<-parameters["vac_campaign"]*7
parameters["travelban_dur"]<-parameters["travelban_dur"]*7
parameters["screen_dur"]<-parameters["screen_dur"]*7
parameters["screen_cov"]<-parameters["screen_cov"]/100
parameters["quarantine_cov"]<-parameters["quarantine_cov"]/100
parameters["quarantine_dur"]<-parameters["quarantine_dur"]*7
parameters["quarantine_days"]<-parameters["quarantine_days"]
parameters["quarantine_effort"]<-1/parameters["quarantine_effort"]
parameters["quarantine_eff_home"]<-parameters["quarantine_eff_home"]/-100
parameters["quarantine_eff_other"]<-parameters["quarantine_eff_other"]/100
parameters["give"]<-parameters["give"]/100
parameters["pdeath_h"]<-parameters["pdeath_h"]/100
parameters["pdeath_hc"]<-parameters["pdeath_hc"]/100
parameters["pdeath_icu"]<-parameters["pdeath_icu"]/100
parameters["pdeath_icuc"]<-parameters["pdeath_icuc"]/100
parameters["pdeath_vent"]<-parameters["pdeath_vent"]/100
parameters["pdeath_ventc"]<-parameters["pdeath_ventc"]/100
parameters["nusc"]<-1/parameters["nusc"]
parameters["nu_icu"]<-1/parameters["nu_icu"]
parameters["nu_icuc"]<-1/parameters["nu_icuc"]
parameters["nu_vent"]<-1/parameters["nu_vent"]
parameters["nu_ventc"]<-1/parameters["nu_ventc"]
parameters["pclin"]<-parameters["pclin"]/100
parameters["prob_icu"]<-parameters["prob_icu"]/100
parameters["prob_vent"]<-parameters["prob_vent"]/100
parameters["lockdown_low_dur"]<-parameters["lockdown_low_dur"]*7
parameters["lockdown_mid_dur"]<-parameters["lockdown_mid_dur"]*7
parameters["lockdown_high_dur"]<-parameters["lockdown_high_dur"]*7

# parameters2<-parameters
# 
# #########    SEVERITY AND MORTALITY
# # age dependent hosp and mort - correction parameters allow to change the shape of the IFR curve
# ifr_correction_young<-2
# ifr_correction_old<-1.75
# ihr <- read.csv("covidagehosp_X.csv",header=TRUE) # hospitalisation rate given infection
# ifr <- read.csv("covidagefrpi_X.csv",header=TRUE) # fatality rate given infection
# ihr<- parameters["ihr_scaling"]*ihr/100   # csv data is in percentages
# ifr_original<-ifr/100   # csv data is in percentages
# for (i in 1:A){
#   ifr[i,2]=ifr[i,2]/max(ifr[,2])    # transform ifr into a normalised age profile (highest value turns into 1)
# }
# ifr[1:14,2]<-ifr[1:14,2]/ifr_correction_young
# ihr$severe[15:21]<-ihr$severe[15:21]*ifr_correction_old

###########################################################################
# Define the indices for each variable
Sindex<-1:A
Eindex<-(A+1):(2*A)
Iindex<-(2*A+1):(3*A)
Rindex<-(3*A+1):(4*A)
Xindex<-(4*A+1):(5*A)
Hindex<-(5*A+1):(6*A)
HCindex<-(6*A+1):(7*A)
Cindex<-(7*A+1):(8*A)
CMindex<-(8*A+1):(9*A)
Vindex<-(9*A+1):(10*A)
QSindex<-(10*A+1):(11*A)
QEindex<-(11*A+1):(12*A)
QIindex<-(12*A+1):(13*A)
QRindex<-(13*A+1):(14*A)
CLindex<-(14*A+1):(15*A)
QCindex<-(15*A+1):(16*A)
ICUindex<-(16*A+1):(17*A)
ICUCindex<-(17*A+1):(18*A)
Ventindex<-(18*A+1):(19*A)
VentCindex<-(19*A+1):(20*A)
CMCindex<-(20*A+1):(21*A)

###########################################################################
# MODEL INITIAL CONDITIONS
initI<-0*popstruc[,2]  # Infected and symptomatic
initE<-0*popstruc[,2]  # Incubating
initE[aci]<-1          # place random index case in E compartment
initR<-0*popstruc[,2]  # Immune
initX<-0*popstruc[,2]  # Isolated 
initV<-0*popstruc[,2]  # Vaccinated 
initQS<-0*popstruc[,2] # quarantined S 
initQE<-0*popstruc[,2] # quarantined E  
initQI<-0*popstruc[,2] # quarantined I  
initQR<-0*popstruc[,2] # quarantined R  
initH<-0*popstruc[,2]  # hospitalised 
initHC<-0*popstruc[,2] # hospital critical 
initC<-0*popstruc[,2]  # Cumulative cases (true)
initCM<-0*popstruc[,2] # Cumulative deaths (true)
initCL<-0*popstruc[,2] # symptomatic cases
initQC<-0*popstruc[,2] # quarantined C 
initICU<-0*popstruc[,2]   # icu
initICUC<-0*popstruc[,2]  # icu critical
initVent<-0*popstruc[,2]  # icu vent
initVentC<-0*popstruc[,2] # icu vent crit
initCMC<-0*popstruc[,2]   # Cumulative deaths (true)
initS<-popstruc[,2]-initE-initI-initR-initX-initV-initH-initHC-initQS-initQE-initQI-initQR-initCL-initQC-initICU-initICUC-initVent-initVentC  # Susceptible (non-immune)


# set up a function to solve the equations
covid<-function(t, Y, parameters) 
{
  with(as.list(c(Y, parameters)),
       {
         S <- Y[Sindex]
         E <- Y[Eindex]
         I <- Y[Iindex]
         R <- Y[Rindex]
         X <- Y[Xindex]
         H <- Y[Hindex]
         HC <- Y[HCindex]
         C <- Y[Cindex]
         CM <- Y[CMindex]
         V <- Y[Vindex]
         QS <- Y[QSindex]
         QE <- Y[QEindex]
         QI <- Y[QIindex]
         QR <- Y[QRindex]
         CL <- Y[CLindex]
         QC <- Y[QCindex]
         ICU <- Y[ICUindex]
         ICUC <- Y[ICUCindex]
         Vent <- Y[Ventindex]
         VentC <- Y[VentCindex]
         CMC <- Y[CMCindex]
         P <- (S+E+I+R+X+V+H+HC+QS+QE+QI+QR+CL+QC+ICU+ICUC+Vent+VentC)
         # print(sum(P))
         
         # health system performance
         f <- c(1,(1+give)/2,(1-give)/2,0)
         KH<-beds_available
         KICU<- icu_beds_available
         Kvent<- ventilators_available
         x.H <- c(0,(1+give)*KH/2,(3-give)*KH/2,2*KH)
         x.ICU <- c(0,(1+give)*KICU/2,(3-give)*KICU/2,2*KICU)
         x.Vent <- c(0,(1+give)*Kvent/2,(3-give)*Kvent/2,2*Kvent)
         fH <- splinefun(x.H, f, method = "hyman")
         fICU <- splinefun(x.ICU, f, method = "hyman")
         fVent<- splinefun(x.Vent, f, method = "hyman")
         critH<-min(1-fH(sum(H)+sum(ICUC))+(1-reporth),1)
         crit<-min(1-fICU(sum(ICU)+sum(Vent)+sum(VentC)),1)
         critV<-min(1-fVent(sum(Vent)),1)
         # print(fH(sum(H)))
         
         # interventions
         isolation<-(t>=selfis_on)*(t<=selfis_on+selfis_dur)
         distancing<-(t>=dist_on)*(t<=(dist_on+dist_dur))
         handwash<-(t>=hand_on)*(t<=(hand_on+hand_dur))
         workhome<-(t>=work_on)*(t<=(work_on+work_dur))
         schoolclose<-(t>=school_on)*(t<=(school_on+school_dur))
         cocoon<-(t>=cocoon_on)*(t<=(cocoon_on+cocoon_dur))*cocoon_cov
         vaccine<-(t>=(vaccine_on))*(t<=vaccine_on+vac_campaign)
         travelban<-(t>=travelban_on)*(t<=(travelban_on+travelban_dur))
         screen<-(t>=screen_on)*(t<=(screen_on+screen_dur))
         quarantine<-(t>=quarantine_on)*(t<=(quarantine_on+quarantine_dur))
         lockdown_low<-(t>=lockdown_low_on)*(t<=(lockdown_low_on+lockdown_low_dur))
         lockdown_mid<-(t>=lockdown_mid_on)*(t<=(lockdown_mid_on+lockdown_mid_dur))
         lockdown_high<-(t>=lockdown_high_on)*(t<=(lockdown_high_on+lockdown_high_dur))

         screen_eff<-0
         selfis<-0
         school<-1
         dist<-1
         hand<-0
         vaccinate<-0
         trvban_eff<-0
         quarantine_rate<-0
         
         if (lockdown_low || lockdown_mid || lockdown_high){
           if(lockdown_low){
             selfis<-0.5
             dist<-0.25
             school<-0
             trvban_eff<-0
             quarantine_rate<-0
             work<-0
             cocoon<-0.95
             hand<-0.05
             vaccinate<-0
           }
           if(lockdown_mid){
             selfis<-0.5
             dist<-0.35
             school<-0.85
             trvban_eff<-0
             quarantine_rate<-0.05
             work<-0.5
             cocoon<-0.95
             hand<-0.05
             vaccinate<-0
           }
           if(lockdown_high){
             selfis<-0.95
             dist<-0.95
             school<-0.85
             trvban_eff<-0.95
             quarantine_rate<-0.9
             work<-0.75
             cocoon<-0.95
             hand<-0.075
             vaccinate<-0
           }
         }
         else{
           if (workhome){
             work<-work_cov*work_eff
           }else{work<-1}
           if (isolation){
             selfis<-selfis_cov
             if(screen){
               screen_eff<-min((report*I+reportc*(CL)+H+ICU+Vent+reporth*HC+ICUC+VentC)*screen_contacts*(screen_overdispersion*I/P)*screen_cov/P,1) 
             }
           }
           if (schoolclose){
             school<-school_eff
           }
           if(distancing){
             dist<-dist_cov*dist_eff
           }
           if(handwash){
             hand<-hand_eff
           }
           if(vaccine){
             vac_rate <- (-log(1-vaccine_cov)/vac_campaign)
             vaccinate <- vac_rate
           }
           if(travelban){
             trvban_eff<-travelban_eff
           }
           if(quarantine){
             quarantine_rate<-min(((I+CL+H+ICU+Vent+HC+ICUC+VentC)*(household_size-1)/P),1)*quarantine_cov*quarantine_effort
           }
         }
         
         
         # cocooning the elderly
         cocoon_mat<-matrix((1-cocoon_eff),nrow = length(popstruc$pop),ncol = length(popstruc$pop))
         cocoon_mat[1:(age_cocoon-1),1:(age_cocoon-1)]<-1
         
         # contact matrices
         cts<-(contact_home+distancing*(1-dist)*contact_other+(1-distancing)*contact_other
               +(1-schoolclose)*contact_school # school on
               +schoolclose*(1-school)*contact_school # school close
               +schoolclose*contact_home*school*s2h # inflating contacts at home when school closes
               +(1-workhome)*contact_work  # normal work
               +workhome*(1-work)*contact_work # people not working from home when homework is active
               +contact_home*workhome*work*w2h # inflating contacts at home when working from home
         )
         
         # Final transmission related parameters
         contacts <- (1-cocoon)*cts+cocoon*cts*cocoon_mat+cocoon*(1+school*(1-school_eff)+work*(1-work_eff))*contact_home*(1-cocoon_mat)
         seas <- 1+amp*cos(2*3.14*(t-(phi*365.25/12))/365.25)
         importation <- mean_imports*(1-trvban_eff)
         HH<-H+ICU+Vent
         HHC<-HC+ICUC+VentC
         lam <- (1-hand)*p*seas*(contacts%*%((rho*E+(I+CL+importation)+(1-selfis_eff)*(X+HHC)+rhos*(HH))/P))
         # contacts under home quarantine
         lamq<-(1-hand)*p*seas*((1-quarantine_eff_home)*contact_home%*%(((1-selfis_eff)*(X+HHC))/P))+(1-hand)*p*seas*(1-quarantine_eff_other)*(contact_other%*%((rho*E+(I+CL+importation)+(1-selfis_eff)*(X+HHC)+rhos*(HH))/P))
         
         # birth/death
         b1<-sum(popbirth[,2]*popstruc[,2])
         birth<-0*popbirth[,2]
         birth[1]<-b1
         
         # ODE system
         dSdt <- -S*lam-S*vaccinate+omega*R+ageing%*%S-mort*S+birth-quarantine_rate*S +(1/quarantine_days)*QS
         dEdt <- S*lam-gamma*E+ageing%*%E-mort*E + (1-vaccine_eff)*lam*V-quarantine_rate*E+(1/quarantine_days)*QE
         dIdt <- gamma*(1-pclin)*(1-screen_eff)*(1-ihr[,2])*E-nui*I+ageing%*%I-mort*I + (1/quarantine_days)*QI - quarantine_rate*I
         dCLdt<- gamma*pclin*(1-selfis)*(1-ihr[,2])*E-nui*CL+ageing%*%CL-mort*CL + (1/quarantine_days)*QC
         dRdt <- nui*I-omega*R+nui*X+nui*CL+ageing%*%R-mort*R + (1/quarantine_days)*QR + nus*(1-pdeath_h*ifr[,2])*H + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC
         dXdt <- gamma*selfis*pclin*(1-ihr[,2])*E+gamma*(1-pclin)*screen_eff*(1-ihr[,2])*E-nui*X+ageing%*%X-mort*X 
         dVdt <- vaccinate*S -(1-vaccine_eff)*lam*V +ageing%*%V - mort*V
         
         dQSdt <- quarantine_rate*S+ ageing%*%QS-mort*QS - (1/quarantine_days)*QS - lamq*QS
         dQEdt <- quarantine_rate*E - gamma*QE + ageing%*%QE-mort*QE - (1/quarantine_days)*QE + lamq*QS 
         dQIdt <- quarantine_rate*I + gamma*(1-ihr[,2])*(1-pclin)*QE-nui*QI+ageing%*%QI-mort*QI - (1/quarantine_days)*QI
         dQCdt <- gamma*(1-ihr[,2])*pclin*QE-nui*QC+ageing%*%QC-mort*QC - (1/quarantine_days)*QC
         dQRdt <- nui*QI+nui*QC+ageing%*%QR-mort*QR - (1/quarantine_days)*QR
         
         dHdt <- gamma*ihr[,2]*(1-prob_icu)*(1-critH)*E + gamma*ihr[,2]*(1-prob_icu)*(1-critH)*QE - nus*H + ageing%*%H-mort*H  # all pdeath have to be lower than
         dHCdt <- gamma*ihr[,2]*(1-prob_icu)*critH*E + gamma*ihr[,2]*(1-prob_icu)*critH*QE - nusc*HC + ageing%*%HC-mort*HC 
         dICUdt <- gamma*ihr[,2]*prob_icu*(1-crit)*(1-prob_vent)*E + gamma*ihr[,2]*prob_icu*(1-crit)*(1-prob_vent)*QE - nu_icu*ICU +ageing%*%ICU - mort*ICU 
         dICUCdt <- gamma*ihr[,2]*prob_icu*crit*(1-prob_vent)*E + gamma*ihr[,2]*prob_icu*crit*(1-prob_vent)*QE - nu_icuc*ICUC +ageing%*%ICUC - mort*ICUC 
         dVentdt <- gamma*ihr[,2]*prob_icu*(1-crit)*(1-critV)*prob_vent*E + gamma*ihr[,2]*prob_icu*(1-crit)*(1-critV)*prob_vent*QE + (1-critV)*VentC*1/2 - nu_vent*Vent +ageing%*%Vent - mort*Vent 
         dVentCdt <- gamma*ihr[,2]*prob_icu*prob_vent*(1-crit)*critV*E +gamma*ihr[,2]*prob_icu*prob_vent*crit*E+
           gamma*ihr[,2]*prob_icu*prob_vent*(1-crit)*critV*QE + gamma*ihr[,2]*prob_icu*prob_vent*crit*QE - 
           (1-critV)*VentC*1/2-nu_ventc*VentC +ageing%*%VentC - mort*VentC 
         
         dCdt <- report*gamma*(1-pclin)*(1-ihr[,2])*(E+QE)+reportc*gamma*pclin*(1-ihr[,2])*(E+QE)+ 
           gamma*ihr[,2]*(1-critH)*(1-prob_icu)*(E+QE)+gamma*ihr[,2]*critH*reporth*(1-prob_icu)*(E+QE)+
           gamma*ihr[,2]*prob_icu*(E+QE)
         dCMdt<- nus*pdeath_h*ifr[,2]*H + nusc*pdeath_hc*ifr[,2]*HC + nu_icu*pdeath_icu*ifr[,2]*ICU + nu_icuc*pdeath_icuc*ifr[,2]*ICUC +  nu_vent*pdeath_vent*ifr[,2]*Vent + nu_ventc*pdeath_ventc*ifr[,2]*VentC + 
           mort*H + mort*HC + mort*ICU + mort*ICUC + mort*Vent + mort*VentC 
         dCMCdt <- nusc*pdeath_hc*ifr[,2]*HC+nu_icuc*pdeath_icuc*ifr[,2]*ICUC + nu_ventc*pdeath_ventc*ifr[,2]*VentC + 
           mort*HC + mort*ICUC + mort*VentC

         # return the rate of change
         list(c(dSdt,dEdt,dIdt,dRdt,dXdt,dHdt,dHCdt,dCdt,dCMdt,dVdt,dQSdt,dQEdt,dQIdt,dQRdt,dCLdt,dQCdt,dICUdt,dICUCdt,dVentdt,dVentCdt,dCMCdt))
       }
  ) 
}


process_ode_outcome <- function(out){
  # Start Bridge ----
  critH<-c()
  crit<-c()
  critV<-c()
  # End Bridge ----
  
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  f <- c(1,(1+parameters["give"])/2,(1-parameters["give"])/2,0) 
  KH<-parameters["beds_available"]
  KICU<- parameters["icu_beds_available"]
  Kvent<- parameters["ventilators_available"]
  x.H <- c(0,(1+parameters["give"])*KH/2,(3-parameters["give"])*KH/2,2*KH) 
  x.ICU <- c(0,(1+parameters["give"])*KICU/2,(3-parameters["give"])*KICU/2,2*KICU) 
  x.Vent <- c(0,(1+parameters["give"])*Kvent/2,(3-parameters["give"])*Kvent/2,2*Kvent) 
  fH <- splinefun(x.H, f, method = "hyman") 
  fICU <- splinefun(x.ICU, f, method = "hyman") 
  fVent<- splinefun(x.Vent, f, method = "hyman") 
  for (i in 1:length(times)){
    critH[i]<-min(1-fH(sum(out[i,(Hindex+1)]))+(1-parameters["reporth"]),1)
    crit[i]<-min(1-fICU((sum(out[i,(ICUindex+1)]))+(sum(out[i,(Ventindex+1)]))+(sum(out[i,(VentCindex+1)]))))
    critV[i]<-min(1-fVent((sum(out[i,(Ventindex+1)]))),1)
  }
  
  # total population
  pop1<-out[,(Sindex+1)]+out[,(Eindex+1)]+out[,(Iindex+1)]+out[,(CLindex+1)]+out[,(Rindex+1)]+out[,(Xindex+1)]+out[,(Vindex+1)]+
    out[,(QSindex+1)]+out[,(QEindex+1)]+out[,(QIindex+1)]+out[,(QCindex+1)]+out[,(QRindex+1)]+
    out[,(Hindex+1)]+out[,(HCindex+1)]+out[,(ICUindex+1)]+out[,(ICUCindex+1)]+out[,(Ventindex+1)]+out[,(VentCindex+1)] 
  tpop1<-rowSums(pop1)
  time<-as.Date(out[,1]+startdate)
  # daily incidence
  inc1 <- parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(QEindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(QEindex+1)]%*%(1-ihr[,2])
  
  inc1h<- parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*critH*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*critH*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["prob_icu"]+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["prob_icu"]
  
  dailyinc1<-rowSums(inc1)+rowSums(inc1h)      # daily incidence
  cuminc1<-colSums(inc1)+colSums(inc1h)        # cumulative incidence
  previcureq1<-rowSums(out[,(Hindex+1)])+ rowSums(out[,(ICUCindex+1)])      # requirement for beds
  previcureq21<-rowSums(out[,(ICUindex+1)])+rowSums(out[,(VentCindex+1)])    # requirement for icu
  previcureq31<-rowSums(out[,(Ventindex+1)])   # requirement for icu
  cmortality1<-rowSums(out[,(CMindex+1)])      # cumulative mortality
  overloadH1<-rowSums(out[,(HCindex+1)])       # requirement for beds
  overloadICU1<-rowSums(out[,(ICUCindex+1)])   # requirement for beds
  overloadVent1<-rowSums(out[,(VentCindex+1)]) # requirement for beds
  ccases1<-rowSums(out[,(Cindex+1)])           # cumulative cases
  
  inc_overloadH1<-((parameters["gamma"]*(1-parameters["prob_icu"])*out[,(Eindex+1)]))
  inc_overloadICU1<-((parameters["gamma"]*parameters["prob_icu"]*(1-parameters["prob_vent"])*out[,(Eindex+1)]))
  for (i in 1:length(times)) {
    inc_overloadH1[i,]<-inc_overloadH1[i,]*critH[i]*ihr[,2]
    inc_overloadICU1[i,]<-inc_overloadICU1[i,]*crit[i]*ihr[,2]
  }
  inc_overloadH1<-cumsum(rowSums(inc_overloadH1))
  inc_overloadICU1<-cumsum(rowSums(inc_overloadICU1))
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  ##########################    CALCULATE MORTALITY 
  pdeath_hc<-parameters["pdeath_hc"]
  prob_icu<-parameters["prob_icu"]
  prob_vent<-parameters["prob_vent"]
  pdeath_icuc<-parameters["pdeath_icuc"]
  pdeath_ventc<-parameters["pdeath_ventc"]
  
  
  cinc_mort_H1 <- cumsum(rowSums(parameters["nus"]*parameters["pdeath_h"]*(out[,(Hindex+1)]%*%ifr[,2])+ out[,(Hindex+1)]%*%mort))
  cinc_mort_HC1 <- cumsum(rowSums(parameters["nusc"]*parameters["pdeath_hc"]*(out[,(HCindex+1)]%*%ifr[,2]) + out[,(HCindex+1)]%*%mort))
  cinc_mort_ICU1 <- cumsum(rowSums(parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUindex+1)]%*%ifr[,2] + out[,(ICUindex+1)]%*%mort))
  cinc_mort_ICUC1 <- cumsum(rowSums(parameters["nu_icuc"]*parameters["pdeath_icuc"]*out[,(ICUCindex+1)]%*%ifr[,2] + out[,(ICUCindex+1)]%*%mort))
  cinc_mort_Vent1 <- cumsum(rowSums(parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(Ventindex+1)]%*%ifr[,2] + out[,(Ventindex+1)]%*%mort))
  cinc_mort_VentC1 <- cumsum(rowSums(parameters["nu_ventc"]*parameters["pdeath_ventc"]*out[,(VentCindex+1)]%*%ifr[,2] + out[,(VentCindex+1)]%*%mort))
  base_mort_H1 <- cumsum(rowSums(out[,(Hindex+1)]%*%mort))
  base_mort_HC1 <- cumsum(rowSums(out[,(HCindex+1)]%*%mort))
  base_mort_ICU1 <- cumsum(rowSums(out[,(ICUindex+1)]%*%mort))
  base_mort_ICUC1 <- cumsum(rowSums(out[,(ICUCindex+1)]%*%mort))
  base_mort_Vent1 <- cumsum(rowSums(out[,(Ventindex+1)]%*%mort))
  base_mort_VentC1 <- cumsum(rowSums(out[,(VentCindex+1)]%*%mort))
  base_mort_S1 <- cumsum(rowSums(out[,(Sindex+1)]%*%mort))
  base_mort_E1 <- cumsum(rowSums(out[,(Eindex+1)]%*%mort))
  base_mort_I1 <- cumsum(rowSums(out[,(Iindex+1)]%*%mort))
  base_mort_CL1 <- cumsum(rowSums(out[,(CLindex+1)]%*%mort))
  base_mort_X1 <- cumsum(rowSums(out[,(Xindex+1)]%*%mort))
  base_mort_QS1 <- cumsum(rowSums(out[,(QSindex+1)]%*%mort))
  base_mort_QE1 <- cumsum(rowSums(out[,(QEindex+1)]%*%mort))
  base_mort_QI1 <- cumsum(rowSums(out[,(QIindex+1)]%*%mort))
  base_mort_QC1 <- cumsum(rowSums(out[,(QCindex+1)]%*%mort))
  base_mort_QR1 <- cumsum(rowSums(out[,(QRindex+1)]%*%mort))
  base_mort_R1 <- cumsum(rowSums(out[,(Rindex+1)]%*%mort))
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  
  Rt <- NULL
  for (i in (ceiling(1/parameters["nui"])+1):length(times)){
    Rt[i] <- cumsum(sum(parameters["gamma"]*out[i,(Eindex+1)]))/cumsum(sum(parameters["gamma"]*out[(i-1/parameters["nui"]),(Eindex+1)]))
    if(Rt[i] >= 7) Rt[i] <- NA
  }
  
  # Export in a cohesive format ----
  results <- list()
  results$time <- startdate + times  # dates
  results$Rt <- Rt
  results$cum_mortality <- round(cmortality1)  # cumulative mortality
  results$pct_total_pop_infected <- round(100 * tail(cumsum(rowSums(parameters["gamma"]*out[,(Eindex+1)])),1)/sum(popstruc[,2]), 1)  # proportion of the  population that has been infected at the end of the simulation
  results$doubling_time <- round(log(2)*7 / (log(dailyinc1[2+7] / dailyinc1[2])), 2)  # (Baseline only) to double the number of infections at inception
  results$required_beds <- round(previcureq1)  # required beds
  results$saturation <- parameters["beds_available"]  # saturation
  results$daily_incidence <- round(dailyinc1)  # daily incidence (Reported)
  results$daily_total_cases <- round(rowSums(parameters["gamma"]*out[,(Eindex+1)]+parameters["gamma"]*out[,(QEindex+1)])) # daily incidence (Reported + Unreported)  # daily incidence (Reported + Unreported)
  results$hospital_surge_beds <- round(previcureq1)
  results$icu_beds <- round(previcureq21)
  results$ventilators <- round(previcureq31)
  
  results$death_natural_non_exposed <- round(base_mort_S1)
  results$death_natural_exposed <- round(base_mort_E1 + base_mort_I1 + base_mort_CL1 + base_mort_X1 + base_mort_QS1 + 
                                           base_mort_QE1 + base_mort_QI1 + base_mort_QC1 + base_mort_QR1 + base_mort_R1)
  results$death_treated_hospital <- round(cinc_mort_H1)
  results$death_treated_icu <- round(cinc_mort_ICU1)
  results$death_treated_ventilator <- round(cinc_mort_Vent1)
  results$death_untreated_hospital <- round(cinc_mort_HC1)
  results$death_untreated_icu <- round(cinc_mort_ICUC1)
  results$death_untreated_ventilator <- round(cinc_mort_VentC1)
  results$total_deaths <- results$death_treated_hospital + results$death_treated_icu + results$death_treated_ventilator +
    results$death_untreated_hospital + results$death_untreated_icu + results$death_untreated_ventilator
  results$total_deaths_end <- last(results$total_deaths)
  
  
  # !!!! code re-using of variable names but with different str() !!! - request some cleaning
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  ## AGE DEPENDENT MORTALITY
  cinc_mort_H1 <- parameters["nus"]*parameters["pdeath_h"]*(out[,(Hindex+1)])
  cinc_mort_HC1 <- parameters["nusc"]*parameters["pdeath_hc"]*(out[,(HCindex+1)])
  cinc_mort_ICU1 <- parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUindex+1)]
  cinc_mort_ICUC1 <- parameters["nu_icuc"]*parameters["pdeath_icuc"]*out[,(ICUCindex+1)] 
  cinc_mort_Vent1 <- parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(Ventindex+1)] 
  cinc_mort_VentC1 <- parameters["nu_ventc"]*parameters["pdeath_ventc"]*out[,(VentCindex+1)] 
  totage1<-as.data.frame(cinc_mort_H1+cinc_mort_HC1+cinc_mort_ICU1+cinc_mort_ICUC1+cinc_mort_Vent1+cinc_mort_VentC1)
  basemort_H1<-(out[,(Hindex+1)])
  basemort_HC1<-(out[,(HCindex+1)])
  basemort_ICU1<-(out[,(ICUindex+1)])
  basemort_ICUC1<-(out[,(ICUCindex+1)])
  basemort_Vent1<-(out[,(Ventindex+1)])
  basemort_VentC1<-(out[,(VentCindex+1)])
  totbase1<-as.data.frame(basemort_H1+basemort_HC1+basemort_ICU1+basemort_ICUC1+basemort_Vent1+basemort_VentC1)
  tc<-c()
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  for (i in 1:dim(cinc_mort_H1)[1]) {
    for (j in 1:dim(cinc_mort_H1)[2]) {
      tc<-rbind(tc,c(i, j, totage1[i,j]*ifr[j,2]+totbase1[i,j]*mort[j])) 
    }
  }
  tc<-as.data.frame(tc)
  colnames(tc)<-c("Day","Age","value")
  
  results$tc <- tc %>%
    mutate(Date = startdate + Day,
           age_cat = case_when(
             Age >=  1 & Age <= 6   ~ "≤ 30 y.o.",
             Age >  6 & Age <= 8    ~ "30-40 y.o.",
             Age >  8 & Age <= 10    ~ "40-50 y.o.",
             Age >  10 & Age <= 12    ~ "50-60 y.o.",
             Age >  12 & Age <= 14    ~ "60-70 y.o.",
             Age >=  15  ~ "≥ 70 y.o.")) %>%
    mutate(age_cat = factor(age_cat, levels = rev(c("≤ 30 y.o.", "30-40 y.o.",
                                                    "40-50 y.o.", "50-60 y.o.", "60-70 y.o.", "≥ 70 y.o."))))
  
  mortality_lag <- data.frame(Age = popstruc$agefloor)
  if(nrow(out) >= 30)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day30 = out[30,CMindex+1]/out[30,Cindex+1]) %>%
                                                    mutate(day30 = ifelse(is.infinite(day30), 0, day30)) %>%
                                                    rename(`Day 30` = day30))
  if(nrow(out) >= 60)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day60 = out[60,CMindex+1]/out[60,Cindex+1]) %>%
                                                    mutate(day60 = ifelse(is.infinite(day60), 0, day60)) %>%
                                                    rename(`Day 60` = day60))
  if(nrow(out) >= 90)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day90 = out[90,CMindex+1]/out[90,Cindex+1]) %>%
                                                    mutate(day90 = ifelse(is.infinite(day90), 0, day90)) %>%
                                                    rename(`Day 90` = day90))
  if(nrow(out) >= 120)  mortality_lag <- bind_cols(mortality_lag, 
                                                   data.frame(day120 = out[120,CMindex+1]/out[120,Cindex+1]) %>%
                                                     mutate(day120 = ifelse(is.infinite(day120), 0, day120)) %>%
                                                     rename(`Day 120` = day120))
  
  results$mortality_lag <- mortality_lag
  
  return(results)
}


###########    
Y<-c(initS,initE,initI,initR,initX,initH,initHC,initC,initCM,initV, initQS, initQE, initQI, initQR, initCL, initQC, initICU, initICUC, initVent, initVentC, initCMC) # initial conditions for the main solution vector
out0 <- ode(y = Y, times = times, method = "euler", hini = 0.05, func = covid, parms = parameters)


### using Rcpp
Rcpp::sourceCpp('comoTH.cpp')
out1 <- ode(y = Y, times = times, func = covidTH, parms = parameters)


results <- process_ode_outcome(out0)

##### list of variables that can be visualised
names(results)

### daily incidence cases
plot(x=results$time,y=results$daily_incidence,ylim=c(0,200),col="red",type='l',xlab="date",ylab="daily incidence cases")
lines(incdata_X,type='l')

### R(t)
plot(x=results$time, y=results$Rt, type='l',xlab = "date",ylab = "predicted R(t)")

