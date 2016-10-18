#In the simplest form of my simulations I set the prevalence at 10% 
#I multiply the estimates of beta by a factor such that 10% is the equilibrium prevalence. 
#The initial time is year 1995. The initial average viral load is vinit = 4.748596 (the value in the data in 1995). 
#The phenotypic variance is 1. 
#Three values of heritability (0.06, 0.36, 0.66) are tested. 
#I adjust the mutational variance to 0.08, 0.23, 0.38 for these three values of heritability respectively, 
#such that phenotypic variance and heritability stay about stable over time. 
#Total population size is 20x10^6 and I tune the susceptible birth rate such that the population size stays about constant.


N <- 200000
Infected <- 28000
tfinal <- 365*40		# number of DAYS;  5 years = 1825 days; 20 years = 7300 days; 50 years = 18250 days
xgraphlimit = round(round(10*((tfinal+365)/365), digits = 0)/10)
tprintHerit <- 365*10
AverageDegree = 0.6			# Average number of partners per person
MinDuration = 0.25*365.0    	# Minimum duration of any relationship
MaxDuration = 7*365.0     	# Maximum duration of any relationship
ProbSex = 1             	# Probability that any given couple will have sex per day
MaxLinks = 1             	# Maximum number of concurrent partners (1 = monogamy) 

BirthRate = 0.002
#NaturalDeathRate = 0.00025 				# OK value
#NaturalDeathRate = 0.0003 				# virulence_source.updated.7May.Oct2015.Mittler.14Nov.prev8.spvl4.pdf3
#NaturalDeathRate = 0.0001 				# Test value
#NaturalDeathRate = 0.00005
NaturalDeathRate = 0.000025

Start_Treatment_Campaign = 700*365  	# Day of epidemic when public health authorities start treatment campaign
Percent_Under_Care = 0.0 		    	# Percentage of people who will get treated if VL exceeds a threshold
Increment_Prob_Under_Care_High_spVL = 0  # Increase in probabity of going under care for folks with high setpoint VLs
DefinitionHighSetPoint = 1e5            # VL Threshold above which someone is more likely to come under care
TreatmentTimeInPatients = 4*365			# Time (within host) following infection when people get treated (assumed the same for all patients)
Treatment_Threshold = 1			        # VL threshold for triggering treatment
VL_After_Treatment = 50	     			# Viral load after treatment
CD4_After_Treatment = 1         		# CD4 load after treatment
TreatmentThreshold_CD4 = 1    			# threshold CD4 category for treatment
CD4_Determines_Treatment = 1  			# treatment based on CD4 (CD4_flag=1) or VL (CD4_flag=0)
CD4_Exp_Flag = 1              			# given cd4 dynamics, are times per cd4 category exp. distributed
Increment_Prob_Under_Care_High_CD4=0
DefinitionHighSetPoint_CD4=0
Start_SexReduction_Campaign = 700*365   # Day of epidemic when public health authorities start campaign to reduce # unprotected partners
Reduction_Mean_Degree = 1.3        		# Reduction in mean degree resulting from this campaign
Start_Condom_Campaign = 700*365   		# Day of epidemic when public health authorities start campaign to reduce probability of unprotected sex
Percent_using_condums_after_condom_campaign = 0.30   # Reduction in the probablity of having unprotected sex resulting from this campaign
Start_Faithfulness_Campaign = 700*365 	# Day of epidemic when public health authorities start campaign to reduce partnership turnover
Increased_Duration = 200               	# Increase in partnership durations resulting from this campaign
Inf_Rate_Function <- 0					# Use "1" for exponential function, Lingappa 2010, "0" for Fraser asymptotic function
InfRateBaseline <- 0.0000003			# 0.0001178/365		# Lingappa 2010 
InfRateExponent <- 3.52					# Lingappa 2010 
MaxInfRate <- 0.0053					# 0.001044; Asymptotic function from Fraser 2007, Assuming P_inf(1 year) = 1 - (1 _ P_inf(1 day))^365
VHalfMaxInfRate <- 13938				# Fraser 2007 
HillCoeffInfRate <- 1.02  				# Fraser 2007 

CD4_Determined_Death <- 1     			# Death rate determined by unpublished CD4 progression model considered by Fraser and colleagues (over-rides Gamma_Death)
Gamma_Death <- 1
Donors_RandomTimeToAIDS_Flag <- 1 		# if 1, GetGammaDelay(shape_parameter, theta);
                                 		# else  Donors_RandomTimeToAIDS[i] =  ExpectedDelayTime;
shape_parameter <- 3.46					# 3.46 is from Fraser

death_rate_constant <- 0.000000003 		# 0.000003 is from CASCADE
death_rate_exponent <- 6.45           	# 6.45 is from CASCADE

Dmax <- 9271 							# from Fraser, maximum time in days of asymptomatic state
D50 <- 3058 							# from Fraser, spVL at which duration is half maximum
Dk <- 0.41 								# 0.41 is from Fraser, Hill coefficient

V0 <- 100 
Progression_Model <- 2        			# 1 = model we have been using all along, 2 = VL jumps to X at time of AIDS
V_peak <- 1e7 
t_peak <- 21 
t_acute <- 90 
AverageLogSP0 <- 4.8
VarianceLogSP0 <- 0.91
MutationVariance <- 0.15				
prog_rate <- 0.20 
Heritability <- 0.36

#Three values of heritability (0.06, 0.36, 0.66) are tested. 
#I adjust the mutational variance to 0.08, 0.23, 0.38 for these three values of heritability respectively, 

idum <- trunc(10000*runif(1))
#idum = 5             # Set to a fixed value for debugging purposes

printOutput=1  #turn on(=1)/off(=0) printStats() and other printf commands - if debugging, want them off

# Specify times and bin sizes for printing viral loads for histograms
VL_print_time = 365  # Print VL distribution very x days
VL_print_lower = -0.25   # Lower limit of VL distribution (Log scale)
VL_print_upper = 8.25    # Upper limit of VL distribution (Log scale)
VL_print_interval = 0.5  # Bin size for printing out viral load distribution (Log scale)

### V_peak at AIDS is currently set at 0.5*V_peak (defined in the .c file)


# Recreate the input file with the parameters specified above
cat(N, "N\n", file="Parameters.txt")
cat(Infected, "Infected\n", file="Parameters.txt", append=TRUE)
cat(AverageDegree, "AverageDegree\n", file="Parameters.txt", append=TRUE)
cat(MinDuration, "MinDuration\n", file="Parameters.txt", append=TRUE)
cat(MaxDuration, "MaxDuration\n", file="Parameters.txt", append=TRUE)
cat(ProbSex, "ProbSex\n", file="Parameters.txt", append=TRUE)
cat(MaxLinks, "MaxLinksPerPerson\n", file="Parameters.txt", append=TRUE)
cat(Start_Treatment_Campaign, "Time_Treatment\n", file="Parameters.txt", append=TRUE)
cat(Percent_Under_Care, "Percent_Under_Care\n", file="Parameters.txt", append=TRUE)
cat(Increment_Prob_Under_Care_High_spVL, "Increment_Prob_Under_Care_High_spVL\n", file="Parameters.txt", append=TRUE)
cat(DefinitionHighSetPoint, "DefinitionHighSetPoint\n", file="Parameters.txt", append=TRUE)
cat(TreatmentTimeInPatients, "Start_Treatments\n", file="Parameters.txt", append=TRUE)
cat(Treatment_Threshold, "Treatment_Threshold\n", file="Parameters.txt", append=TRUE)
cat(VL_After_Treatment, "VL_Treatment\n", file="Parameters.txt", append=TRUE)
cat(CD4_Determines_Treatment,"CD4_Death\n", file="Parameters.txt", append=TRUE)
cat(TreatmentThreshold_CD4, "TreatmentThreshold_CD4\n", file="Parameters.txt", append=TRUE)
cat(CD4_After_Treatment, "CD4_After_Treatment\n", file="Parameters.txt", append=TRUE)
cat(Start_SexReduction_Campaign, "Start_SexReduction_Campaign\n", file="Parameters.txt", append=TRUE)
cat(Reduction_Mean_Degree, "Reduction_Mean_Degree\n", file="Parameters.txt", append=TRUE)
cat(Start_Condom_Campaign, "Start_Condom_Campaign\n", file="Parameters.txt", append=TRUE)
cat(Percent_using_condums_after_condom_campaign, "Percent_using_condums_after_condom_campaign\n", file="Parameters.txt", append=TRUE)
cat(Start_Faithfulness_Campaign, "Start_Faithfulness_Campaign\n", file="Parameters.txt", append=TRUE)
cat(Increased_Duration, "Increased_Duration\n", file="Parameters.txt", append=TRUE)
cat(Inf_Rate_Function, "Linear_Inf_Rate_Function\n", file="Parameters.txt", append=TRUE)
cat(InfRateBaseline, "LinearInfRateBaseline\n", file="Parameters.txt", append=TRUE)
cat(InfRateExponent, "LinearInfRateIncrease\n", file="Parameters.txt", append=TRUE)
cat(MaxInfRate, "MaxInfRate\n", file="Parameters.txt", append=TRUE)
cat(VHalfMaxInfRate, "VHalfMaxInfRate\n", file="Parameters.txt", append=TRUE)
cat(HillCoeffInfRate, "HillCoeffInfRate\n", file="Parameters.txt", append=TRUE)
cat(BirthRate, "BirthRate\n", file="Parameters.txt", append=TRUE)
cat(NaturalDeathRate, "NaturalDeathRate\n", file="Parameters.txt", append=TRUE)
cat(Progression_Model, "Progression_Model\n", file="Parameters.txt", append=TRUE)
cat(CD4_Determined_Death,"CD4_Death\n", file="Parameters.txt", append=TRUE)
cat(Gamma_Death, "Gamma_Death\n", file="Parameters.txt", append=TRUE)
cat(shape_parameter, "shape_parameter\n", file="Parameters.txt", append=TRUE)
cat(death_rate_constant, "death_rate_constant\n", file="Parameters.txt", append=TRUE)
cat(death_rate_exponent,"death_rate_exponent\n", file="Parameters.txt", append=TRUE)
cat(Dmax,  "Dmax\n", file="Parameters.txt", append=TRUE)
cat(D50, "D50\n", file="Parameters.txt", append=TRUE)
cat(Dk, "Dk\n", file="Parameters.txt", append=TRUE)
cat(V0, "V0\n", file="Parameters.txt", append=TRUE)
cat(V_peak, "V_peak\n", file="Parameters.txt", append=TRUE)
cat(t_peak, "t_peak\n", file="Parameters.txt", append=TRUE)
cat(t_acute, "t_acute\n", file="Parameters.txt", append=TRUE)
cat(AverageLogSP0, "AverageLogSP0\n", file="Parameters.txt", append=TRUE)
cat(VarianceLogSP0, "VarianceLogSP0\n", file="Parameters.txt", append=TRUE)
cat(MutationVariance,"MutationVariance\n",file="Parameters.txt",append=TRUE)
cat(tfinal, "tfinal\n", file="Parameters.txt", append=TRUE)
cat(prog_rate, "Prog_Rate\n", file="Parameters.txt", append=TRUE)
cat(Heritability, "Heritability\n", file="Parameters.txt", append=TRUE)
cat(tprintHerit, "tprintHerit\n", file="Parameters.txt", append=TRUE)
cat(VL_print_time, "VL_Print_time\n", file="Parameters.txt", append=TRUE)
cat(VL_print_lower, "VL_Print_time_lower\n", file="Parameters.txt", append=TRUE)
cat(VL_print_upper, "VL_Print_time_upper\n", file="Parameters.txt", append=TRUE)
cat(VL_print_interval, "VL_Print_time_interval\n", file="Parameters.txt", append=TRUE)
cat(idum, "idum\n", file="Parameters.txt", append=TRUE)
cat(printOutput, "printOutput\n", file="Parameters.txt", append=TRUE)
cat(CD4_Exp_Flag, "CD4_Exp_Flag\n", file="Parameters.txt", append=TRUE)
cat(Donors_RandomTimeToAIDS_Flag, "Donors_RandomTimeToAIDS_Flag\n", file="Parameters.txt", append=TRUE)
cat(Increment_Prob_Under_Care_High_CD4, "Increment_Prob_Under_Care_High_CD4\n", file="Parameters.txt", append=TRUE)
cat(DefinitionHighSetPoint_CD4, "DefinitionHighSetPoint_CD4\n", file="Parameters.txt", append=TRUE)

#####################################################################

cat("**** About to compile the C code ****\n")
system("rm a.out spVL")
system("gcc -O3 virulence_source.updated.Oct2016.alpha+0.057.c")
system("cp a.out spVL")
cat("**** About to run the C program ****")
system("./spVL")

spVL.run <- read.table("ViralLoadOutput.txt", header=TRUE)
patient.run <- read.table("PatientRegistryOutput.txt", header=TRUE)

