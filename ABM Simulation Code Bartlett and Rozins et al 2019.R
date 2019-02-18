#### Simulation Code for Bartlett, Rozins et al. 2019 ' Industrial bees: the impact of agricultural intensification on local disease prevalence'

# Agent based model  built by Dr Lewis J. Bartlett. Contact: ljbartl@emory.edu, or see online for up to date comtact details

# Set file writeout location
OUTDIR <- ''

## Set parameters and parameter ranges

# Number of repeats for each unique parameter set
RepeatCount <- 
# Length of simulations
Days <- 
# Apiary Sizes
HiveCount <- c()
# Define spatial apiary arrangement - can be a line, a circle, or a square grid; mixed for panmitic
LayoutTypes <- c('Linear','Ring','Grid','Mixed')
# Define birthrate (recruitment) - this will be a mean - interpret as 'new bees per day'
MBirth <- 1600
# Calculate deathrate (likelihood per bee per day) based on a mean max hive size
MaxHive <- 58200
DeathRate <- MBirth / MaxHive
# Define drift rate (likelihood per bee per day)
DriftRates <- c()

# Define disease mortality (likelihood per day) and maximum contagion (probability of infecting a suspectible in a day)
BetaRange <- c()
AlphaRange <- c()


### Create functions which do the 'hard-lifting' of the repeating ecology steps

#Now We get into making our simulation functions

# Create function to fill / update cols TotalNS and TotalNI
# Two versions: one just to calculate number of neighbours at initialisation
# The second is the workhorse function

InitialiseCalcN <- function(frame, arrangement, drift){
  
  # Check if there are appropriate ways to calculate neighbours for the supplied
  # spatial arrangement
  if(arrangement != 'Linear' & arrangement != 'Ring' & arrangement != 'Grid' & arrangement != 'Mixed'){
    
    stop('Spatial Arrangement not Recognised: must be one of Mixed, Grid, Ring or Linear')
    
  }
  
  # For spatially structured:
  if(arrangement != 'Mixed'){
    
    # Extract number of hives
    NumHives <- NROW(frame)
    
    #Neighbours for linear are simple (non-existent hives e.g. 0 are filtered later)
    for(H in 1:NumHives){
      if(arrangement=='Linear'){
        Nbrs <- c(H-1, H+1)
      }
      #Neighbours for a ring require special conditions for first and last
      if(arrangement=='Ring'){
        
        if(H != 1 & H != NumHives){
          Nbrs <- c(H-1, H+1)
        }
        if(H == 1){
          Nbrs <- c(NumHives, 2)
        }
        if(H == NumHives){
          Nbrs <- c(H-1, 1)
        }
        
      }
      
      # Neighbours in a grid extracted using an expanded grid buffered by an 
      # outer row and column of 0 values and matching using the wonders of 
      # arr.ind=T in the which() function
      # again the 0 values get filtered out after
      if(arrangement=='Grid'){
        ExpGrid <- matrix(0, ncol=(sqrt(NumHives)+2), nrow=(sqrt(NumHives)+2))
        
        IDGrid <- matrix(1:NumHives, ncol=sqrt(NumHives), byrow=T)
        
        ExpGrid[2:(sqrt(NumHives)+1), 2:(sqrt(NumHives)+1)] <- IDGrid
        
        Nbrs <- c(ExpGrid[which(ExpGrid==H, arr.ind=T)+c(1,0)],
                  ExpGrid[which(ExpGrid==H, arr.ind=T)+c(-1,0)],
                  ExpGrid[which(ExpGrid==H, arr.ind=T)+c(0,1)],
                  ExpGrid[which(ExpGrid==H, arr.ind=T)+c(0,-1)])
      }
      
      #Filter out hypothetical neighbours who dont actually exist (e.g. hive -4, hive 0)
      NbrHives <-(match(Nbrs, Ref$Hive, nomatch=NA))[
        (!(is.na(match(Nbrs, Ref$Hive, nomatch=NA))))]
      
      frame$Nbrs[[H]] <- NbrHives
      
      #Sum up total number neighbours
      frame$CountN[H] <- length(NbrHives)
      
    }
    
  }
  
  # Mixed case
  if(arrangement == 'Mixed'){
    
    Nbrs <- list(frame$Hive)
    
    frame$CountN <- NROW(frame)
    frame$Nbrs <- rep.int(Nbrs, times = NROW(frame))
    
  }
  
  #Return the frame to replace the supplied one
  return(frame)
  
}

# Derivative version for repetitive updating (speed things up, above version for intialisation only)

CalcN <- function(frame, arrangement, drift){
  
  
  # For spatially structured:
  if(arrangement != 'Mixed'){
    
    # Extract number of hives
    NumHives <- NROW(frame)
    
    for(H in 1: NumHives){
      
      #Pull out neighbours
      NbrHives <- frame$Nbrs[[H]]
      
      #Sum up totals in neighbours
      frame$TotalNS[H] <- sum(frame$TotalS[NbrHives])
      frame$TotalNI[H] <- sum(frame$TotalI[NbrHives])
      
      #Sum up numbers of neighbouring drifting bees drifting into this give, scaled for number of neighbours those hives have (to coarsely correct for bees moving into 'two hives at once', though note this is approximate)
      frame$DriftNS[H] <- sum(rbinom(n = length (NbrHives), size = rbinom(1,frame$TotalS[NbrHives],drift), prob = 1/(frame$CountN[NbrHives])))
      frame$DriftNI[H] <- sum(rbinom(n = length (NbrHives), size = rbinom(1,frame$TotalI[NbrHives],drift), prob = 1/(frame$CountN[NbrHives])))
    }
    
  }
  
  
  # Mixed case
  if(arrangement == 'Mixed'){
    
    #Totals are populations totals minus own colony
    frame$TotalNS <- sum(frame$TotalS) - frame$TotalS
    frame$TotalNI <- sum(frame$TotalI) - frame$TotalI
    
    # Draw total number of drifting bees then randomly assign across all hives
    
    frame$DriftNS <- rbinom(NROW(frame), rbinom(1, sum(frame$TotalS), drift), 1/NROW(frame))
    frame$DriftNI <- rbinom(NROW(frame), rbinom(1, sum(frame$TotalI), drift), 1/NROW(frame))
    
    
  }
  
  #Return the frame to replace the supplied one
  return(frame)
  
}


# Create function which calculates and modifies TotalS and TotalI cols

StepMod <- function(frame, arrangement, death, drift, d.mort, d.cont){
  
  # Extract number of hives
  NumHives <- NROW(frame)
  
  # Create the frame which we will update for export
  # The old frame should not change as we go along
  NewFrame <- frame
  
  PTNI <- 0
  
  for( H in 1:NumHives){
    # Get change in suscs due to demography
    # Births, then the proportion of bees which die
    DeltaS <- frame$BR[H] - rbinom(1, frame$TotalS[H], death)
    
    
    # Get change in infecteds due to death
    DeltaI <- -1 * rbinom(1, frame$TotalI[H], (death + d.mort))
    
    
    # Now we can calculate infection
    
    #Drifters out of the colony
    SDrift <- rbinom(1, frame$TotalS[H], drift)
    IDrift <- rbinom(1, frame$TotalI[H], drift)
    
    # internal infection
    # and infection from susc bees drifting into other hives (take a probability averaging liberty with this bit)
    
    IntS <- frame$TotalS[H] - SDrift
    IntI <- (frame$TotalI[H] - IDrift) + frame$DriftNI[H]
    
    NewI <- (
      sum((replicate(n = IntS, expr = rbinom(n = 1, size = IntI, prob = d.cont))) > 0) 
      +
      sum((replicate(n = SDrift, expr = rbinom(n = 1, size = sum(frame$TotalI[(frame$Nbrs[[H]])]), prob = d.cont/length(frame$Nbrs[[H]])))) > 0)
    )
    

    #Use these deltas to caculate the new number of bees in each class
    
    NewFrame$TotalS[H] <- frame$TotalS[H] + (DeltaS - NewI)
    NewFrame$TotalI[H] <- frame$TotalI[H] + (DeltaI + NewI)
    
    PTNI <- PTNI + NewI 
    
  }
  # Return the modified frame with new totals of each class
  # HOWEVER THE NEIGHBOUR TOTALS ARE OUTDATED AT THIS STEP
  
  RetList <- list(NewFrame, PTNI)
  
  return(RetList)
  
}


#### Model runs (nested for loops)

## Run tracker
RunNumber <- 0

for(Repeats in RepeatCount){
  
  for(Day in Days){
    
    for(Cont in BetaRange){
      
      for(Mort in AlphaRange){
        
        for(Layout in LayoutTypes){
          
          for(ASize in HiveCount){
            
            for(DriftRate in DriftRates){
              
              ### Individual Simulation run starts here
              
              #Tracks total infection across the apiary at each time step
              ITrack <- as.data.frame(matrix(NA, ncol = 4, nrow = Day))
              colnames(ITrack) <- c('Time', 'TI', 'TS', 'NI')
              
              SimLength <- Day
              
              Spat <- Layout
              
              # For a Grid spatial design, only accept a square number
              
              if( (sqrt(ASize)%%1) != 0 & Spat=='Grid'){
                print('Run Aborted - Apiary Size not Square')
                break
              }
              
              
              ######## Intialise the model #######
              
              # Create reference table of hives (ID, brithrate, suscs, infecteds, neighbours)
              Ref <- as.data.frame(matrix(NA, ncol= 10, nrow=ASize))
              colnames(Ref) <- c('Hive','BR','TotalS','TotalI','TotalNS','TotalNI','DriftNS','DriftNI','CountN', 'Nbrs')
              
              Ref$Nbrs <- as.list(Ref$Nbrs)
              
              #Give each hive a number (functionally equivalent to row number)
              Ref$Hive <- 1:ASize
              
              #Populate with hives of laying rates forming a normal distribition about the mean rate
              Ref$BR <- round(rnorm(n = ASize, mean = MBirth, sd = (MBirth/10)))
              
              # Initialise with single individual in one hive infected
              # Initial colony size related to mean birth rate with some variation around that
              Ref$TotalS <- round(rnorm(n = ASize, mean = 9*MBirth, sd = (9/8)*MBirth))
              Ref$TotalI <- 0
              Ref$TotalI[sample(x=(1:NROW(Ref)), size=1, prob = NULL)] <- 1
              
              Ref <-  InitialiseCalcN(frame = Ref, arrangement = Spat, drift = DriftRate)
              
              # Run 1st function to calc number of neighbouring S and I
              Ref <- CalcN(frame = Ref, arrangement = Spat, drift = DriftRate)
              
              
              # Set the simulations running!
              
              DayTrack <- 0
              
              while(DayTrack < SimLength){
                
                DayTrack <- DayTrack +1
                
                Temp <- StepMod(frame = Ref, arrangement = Spat, death = DeathRate, drift = DriftRate,
                                d.mort = Mort, d.cont = Cont)
                Ref <- Temp[[1]]
                
                ITrack$Time[DayTrack] <- DayTrack
                ITrack$TI[DayTrack] <- sum(Ref$TotalI)
                ITrack$TS[DayTrack] <- sum(Ref$TotalS)
                ITrack$NI[DayTrack] <- Temp[[2]]
                
                rm(Temp)
                
                Ref <- CalcN(frame = Ref, arrangement = Spat, drift = DriftRate)
                
                
                #Reset if pathogen goes extinct (see inititalisation code above)
                if(sum(Ref$TotalI) == 0){
                  Ref <- as.data.frame(matrix(NA, ncol= 9, nrow=ASize))
                  colnames(Ref) <- c('Hive','BR','TotalS','TotalI','TotalNS','TotalNI','DriftNS','DriftNI','CountN')
                  
                  Ref$Hive <- 1:ASize
                  
                  Ref$BR <- round(rnorm(n = ASize, mean = MBirth, sd = (MBirth/10)))
                  
                  Ref$TotalS <- round(rnorm(n = ASize, mean = 9*MBirth, sd = (9/8)*MBirth))
                  Ref$TotalI <- 0
                  Ref$TotalI[sample(x=(1:NROW(Ref)), size=1, prob = NULL)] <- 1
                  
                  Ref <-  InitialiseCalcN(frame = Ref, arrangement = Spat, drift = DriftRate)
                  
                  Ref <- CalcN(frame = Ref, arrangement = Spat, drift = DriftRate)
                  
                  DayTrack <- 0
                  
                }
                # Print daytrack for troubleshooting purposes, not recommended for many repeats
                #print(DayTrack)
                  
              }
              
              RunNumber <- RunNumber + 1
              
              save.image(file = paste(OUTDIR,'Output',RunNumber,'.RData', sep = ''))
              
              print(RunNumber)
              
            }
          }
        }
      }
    }
  }
}

