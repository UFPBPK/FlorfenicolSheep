## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(ggprism)
library(ggplot2)
library(tibble)
library(rlang)
library(patchwork)

## Set working direction to the data files




Tol_liver=3
Tol_kidney=0.3
Tol_muscle=0.2



Tol_liver_fsis=3.7
Tol_kidney_fsis=0.3
Tol_muscle_fsis=0.3



Tol_liver_tol=0.15
Tol_kidney_tol=0.15
Tol_muscle_tol=0.15



## Physiological global parameter

BW.mean        = 58.73                 ## kg; Bodyweight of market age sheep, Li 2021

## Dosing parameters global parameters

PDOSEiv = 0                       ## IV dose:: (mg/kg)
PDOSEsc = 0                       ## SC dose:: (mg/kg)
PDOSEim = 20                      ## IM dose:: (mg/kg)
PDOSEoral = 0                     ## Oral gavage:: (mg/kg)

## Time parameters for simulation 
N         <- 1000                 ## Number of iteration

end_time=150

TDoses=2

tinterval=48

dt = 2

route="im"


tsamp=tgrid(0,tinterval*(TDoses-1)+end_time*24,dt)




library(mrgsolve)  

pbpk.code <- '
$PARAM @annotated

 //IM absorption rate constants
 
  Kim       : 0.3971582  :                       Absorption rate Yang et al. 2019
  Kdissim   : 0.1183221  :                        Rate of chemical becoming available for absorption from the injection 
  
 //SC absorption rate constants
  
  Ksc       : 0.208137913     :                    Absorption rate Yang et al. 2019
  Kdisssc   : 0.009089343     :                    Rate of chemical becoming available for absorption from the injection site 
  
  
  
  // Cardiac Output and Blood flow
  
  QCC       :  7.15      :                       Cardiac output (L/h/kg), Li 2021
  QLCa       : 0.3054    :                       Fraction of flow to the liver, Li 2021
  QKCa       : 12.98/100 :                       Fraction of flow to the kidneys, Li 2021
  QFCa       : 6.18/100  :                       Fraction of flow to the fat, Li 2021
  QMCa       : 10.09/100 :                       Fraction of flow to the muscle, Li 2021
  QLuC      : 1         :                       Fraction of blood flow to the lung (Achenbach, 1995), Li 2021
  QRestCa    : 0.4021    :                       Fraction of blood flow to the rest of te body
  Htc       :34.69/100  :                        Percentage of Hematocrit in the blood, Li 2021
  PCV       :36.15/100  :                       Serum as a percentage of whole blood, Li 2021
  
  // Tissue Volume
  
  VBloodC   : 0.0486    :                       Unitless, Fraction vol. of blood; Li,2021               
  VLCa       : 1.27/100  :                       Unitless, Fraction vol. of liver; Li 2021
  VKCa       : 0.23/100  :                       Unitless, Fraction vol. of kidney; Li 2021
  VFCa       : 20.98/100 :                       Unitless, Fractional fat tissue, Li 2021
  VMCa       : 24.78/100 :                       Unitless, Fractional muscle tissue, Li 2021
  VLuCa      : 1.07/100  :                       Unitless, Fractional lung tissue, Li 2021
  VabCa      : 0.011664  :                       Unitless, Fractional arterial blood, Li 2017
  VvbCa      : 0.036936  :                       Unitless, Fractional venous blood, Li 2017
  VRestCa    : 0.4681    :                       Unitless, Fractional rest of the body, Li 2017
  
  // Mass Transfer Parameters (Chemical-specific parameters) 
  // Partition coefficients(PC, tissue:plasma)

  
  PLu       : 0.9        :                     Lung: plasma PC 
  PL        : 1.3        :                     Liver plasma PC Yang et al. 2019
  PK        : 0.9        :                     Kidney plasma PC Yang et al. 2019
  PM        : 0.9        :                     Muscle plasma PC Yang et al. 2019
  PF        : 0.1137     :                     Fatp plasma PC Yang et al. 2019
  PRest     : 1.13       :                     Rest of the body plasma PC
  
  PLum      : 12.741375 :                   Lung: plasma PC 
  PLm       : 12.741375 :                   Liver: plasma PC Yang et al. 2019
  PKm       : 2.902587 :                   Kidney:plasma PC Yang et al. 2019
  PMm       : 0.4778312:                   Muscle:plasma PC Yang et al. 2019
  PFm       : 0.1      :                   Fat:plasma PC Yang et al. 2019
  PRestm    : 0.9      :                   Rest of the body:plasma PC
  
  // Absorption and elimination parameters
  
  GEC       : 0.182   :               1/(h); Gastric emptying time; 
  K0C       : 1.991   :               1/(h); Rate of uptake from the stomach into the liver
  KabsC     : 0       :               1/(h); Rate of absorption of chemical from small intestine to liver Lin 2015
  KunabsC   : 0       :               1/(h); Rate of unabsorbed dose to appear in feces from small intestine
  
  KbileC    : 0       :               1/(h*kg); Biliary elimination rate from liver to feces   
  KurineC   : 1.391167:               L/(h*kg); Rate of urine elimination from urine storage
  KurineCm  : 0.0042  :               L/(h*kg); Rate of urine elimination from urine storage
  
  Kffa      : 1.52    :                Rate of transfer to florfenicol amine
  Kmo       : 0.1014  :                Rate of transfer to other matabolite
  
  KbileCm   :   0     :               Biliary excretion rate
  
  FracFFA   :   0.5   :               Fractions of FFA 
  
  PB        :  0.2    :               Plasma protein binding coefficient for florfenicol
  PBm       :  0.2    :               Plasma protein binding coefficient for florfenicol amine

  
  BW        : 58.73      :                    kg, Bodyweight of market age goats, Li 2021
  PDOSEsc   : 0          :                    SC dose, (mg/kg)
  PDOSEim   : 0          :                    IM dose, (mg/kg)


  Fracim    : 0.45        :                    Fractions readily available for absorption for IM dose Yang et al 2019
  Fracsc    : 0.26        :                    Fractions readily available for absorption for SC doseYang et al 2019

  

$MAIN


// Cardiac output and blood flows
      
      double Freem =1-PBm;
      double Free  = 1-PB;                     // Unitless; fraction of chemical free in blood after plasma protein binding


 double sumQ = QLCa + QKCa + QMCa + QFCa + QRestCa; // sum up cardiac output fraction
  double QLC = QLCa/sumQ;                     // adjusted blood flow rate fraction to liver
  double QKC = QKCa/sumQ;                     // adjusted blood flow rate fraction to kidney
  double QMC = QMCa/sumQ;                     // adjusted blood flow rate fraction to muscle
  double QFC = QFCa/sumQ;                     // adjusted blood flow rate fraction to fat
  double QRestC = QRestCa/sumQ;               // adjusted blood flow rate fraction to rest of body
  
  double sumV = VLCa + VKCa + VMCa + VFCa + VLuCa + VRestCa + VvbCa + VabCa; //sum up the tissue volumes
  
  double VLC = VLCa/sumV;                      //adjusted fraction of tissue volume of liver
  double VKC = VKCa/sumV;                      //adjusted fraction of tissue volume of kidney
  double VMC = VMCa/sumV;                      //adjusted fraction of tissue volume of muscle
  double VFC = VFCa/sumV;                      //adjusted fraction of tissue volume of fat
  double VLuC = VLuCa/sumV;                    //adjusted fraction of tissue volume of lung
  double VabC = VabCa/sumV;
  double VvbC = VvbCa/sumV;
  double VRestC = VRestCa/sumV;                //adjusted fraction of tissue volume of rest of body


  
  // Cardiac output and blood flows

double QC = QCC*(BW)*(1-Htc);            //  L/h Cardiac output adjusted for plasma
double QK = (QKC*QC);                    //  L/h Blood flow to kidney            
double QL = (QLC*QC);                    //  L/h Blood flow to liver
double QM = (QMC*QC);                    //  L/h Blood flow to muscle
double QF = (QFC*QC);                    //  L/h Blood flow to fat
double QLu = (QLuC*QC);                  //  L/h Blood flow to lung
double QRest = QC-(QK+QL+QM+QF);         //  L/h Blood flow to rest of the body                                     
double QBal = QC - (QK+QL+QM+QF+QRest);  //  L/h Balance check of blood flows should equal zero

double VL = VLC*BW;                                  // L,   Volume of liver 
double VK = VKC*BW;                                  // L,   Volume of kidney 
double VM = VMC*BW;                                  // L,   Volume of muscle 
double VF = VFC*BW;                                  // L,   Volume of fat
double VLu = VLuC*BW;                                // L,   Volume of lung 
double VRest = VRestC*BW;  // L,   Volume of Rest
double VBal = BW- Vab- Vvb - VK - VL-VM- VF -VRest - VLu;

double Vvb=VvbC*BW*((1-Htc)) ;                       // L,   Volume of plasma 
double Vab=VabC*BW*(1-Htc) ;                         // L,   Volume of plasma



      
      // Elimination and absorption parameters
      
      double Kbile        = KbileC*BW;     // Kbile:       1/h; Billiary elimination; liver to feces storage                         
      double Kurine       = KurineC*BW;   // Kurine:      1/h; Urinary elimination;
      
      double Kbilem        = KbileCm*BW;     // Kbile:       1/h; Billiary elimination; liver to feces storage                         
      double Kurinem       = KurineCm*BW;    // Kurine:      1/h; Urinary elimination;   
      
      double GE           = GEC;           // GE:          1/h; Gastric emptying time                            
      double K0           = K0C;           // K0:          1/h; Rate of uptake from the stomach into the liver                                
      double Kabs         = KabsC;         // Kabs:        1/h; Rate of absorption from small intestines to liver
      double Kunabs       = KunabsC;       // Kunabs:      1/h; Rate of elimination from small intestines




$CMT ADOSE Aab Avb Absorbim Absorbsc Amtsite Asctsite 
                     ARest ALu AK AF AM Aurine AST ASI AabsST AabsSI  
                     Afeces AL AUCK AUCCR AUCCL AUCCF AUCCM AUCCLu Aabm Avbm 
                     ARestm ALum AKm AFm AMm Aurinem ALm AUCKm AUCCRm AUCCLm AUCCFm AUCCMm AUCCLum Affa AotherMO ADOSEim ADOSEsc AUCP            

             
             
$ODE

dxdt_ADOSE=0;


// Dosing for IM and SC intramuscular 
// The chaning rate of the amount of dose via IM

double RIM          = Kim*Amtsite;                              // Rim, drug absorption rate of intramuscular route (mmol/h)
double Rsite        = -RIM + Kdissim*ADOSEim;                   // Rsiteim, changing rate of drug at the intramuscular injection site (mmol/h)
dxdt_Absorbim       = RIM;                                      // Absorbim, amount of drug absorbed by intramuscular injection (mmol)
dxdt_Amtsite        = Rsite;                                    // Amtsiteim, amount of the drug at the intramuscular injection site (mmol)
dxdt_ADOSEim        = -Kdissim*ADOSEim;

// The chaning rate of the amount of dose via IM  
double RSC          = Ksc*Asctsite;                              // Rim, drug absorption rate of intramuscular route (mmol/h)
double Rscsite      = -RSC + Kdisssc*ADOSEsc;                    // Rsiteim, changing rate of drug at the intramuscular injection site (mmol/h)
dxdt_Absorbsc       = RSC;                                       // Absorbim, amount of drug absorbed by intramuscular injection (mmol)
dxdt_Asctsite       = Rscsite;                                   // Amtsiteim, amount of the drug at the intramuscular injection site (mmol)
dxdt_ADOSEsc        = -Kdisssc*ADOSEsc;


// Concentration of the chemical in blood compartment

      
      double Cv=Avb/Vvb;                 // Concentration of Chemical in venous plasma
      double Ca=Aab/Vab;                 // Concentration of Chemical in arterial plasma
      
      double Cvm=Avbm/Vvb;               // Concentration of Chemical in venous plasma
      double Cam=Aabm/Vab;               // Concentration of Chemical in arterial plasma
      
      // Concentration of the chemical in vein and tissue compartments
      
      double CL = AL/(VL);              // Concentration of Chemical in liver 
      double CVL = CL/(PL);             // Concentration of Chemical in venous plasma leaving liver 
      double CK = AK/(VK);              // Concentration of Chemical in Kidney  
      double CVK = CK/(PK);             // Concentration of Chemical in venous plasma leaving Kidney 
      double CF = AF/(VF);              // Concentration of Chemical in fat 
      double CVF = CF/(PF);             // Concentration of Chemical in venous plasma leaving fat 
      double  CRest = ARest/(VRest);    // Concentration of Chemical in Rest 
      double CVRest = CRest/(PRest);    // Concentration of Chemical in venous plasma leaving Rest 
      double CLu = ALu/(VLu);           // Concentration of Chemical in Lung
      double CVLu = CLu/(PLu);          // Concentration of Chemical in venous plasma leaving Lung
      double CM=AM/(VM);                // Concentration of Chemical in muscle
      double CVM=CM/(PM);               // Concentration of Chemical venous plasma leaving in muscle
      
      
      double CLm = ALm/(VL);              // Concentration of FFA in liver 
      double CVLm = CLm/(PLm) ;           // Concentration of FFA in venous plasma leaving liver 
      double CKm = AKm/(VK);              // Concentration of FFA in Kidney  
      double CVKm = CKm/(PKm);            // Concentration of FFA in venous plasma leaving Kidney 
      double CFm = AFm/(VF);              // Concentration of FFA in fat 
      double CVFm = CFm/(PFm);            // Concentration of FFA in venous plasma leaving fat 
      double CRestm = ARestm/(VRest);     // Concentration of FFA in Rest 
      double CVRestm = CRestm/(PRestm);   // Concentration of FFA in venous plasma leaving Rest 
      double CLum = ALum/(VLu );          // Concentration of FFA in Lung
      double CVLum = CLum/(PLum);         // Concentration of FFA in venous plasma leaving Lung
      double CMm=AMm/(VM);                // Concentration of FFA in muscle
      double CVMm=CMm/(PMm);              // Concentration of FFA venous plasma leaving in muscle      

      double Cv1=Cv*1000*358.21;               // Concentration of FFC in venous plasma in mg/L
      double Ca1=Ca*1000*358.21;               // Concentration of FFC in arterial plasma in mg/L
      
      double Cvm1=Cvm*1000*247.28;           // Concentration of FFA in venous plasma in mg/L
      double Cam1=Cam*1000*247.28;           // Concentration of FFA in arterial plasma in mg/L
      
      // Concentration of the chemical in vein and tissue compartments
      
      double CL1 = CL*1000*358.21;              // Concentration of FFC in liver in mg/L
      double CVL1 = CVL*1000*358.21;            // Concentration of FFC in venous plasma leaving liver in mg/L
      double CK1 = CK*1000*358.21;              // Concentration of FFC in Kidney in mg/L 
      double CVK1 = CVK*1000*358.21 ;           // Concentration of FFC in venous plasma leaving Kidney in mg/L
     double  CF1 = CF*1000*358.21;              // Concentration of FFC in fat in mg/L
      double CVF1 = CVF*1000*358.21 ;           // Concentration of FFC in venous plasma leaving fat in mg/L
      double CRest1 = CRest*1000*358.21;        // Concentration of FFC in Rest in mg/L
      double CVRest1 = CVRest*1000*358.21;      // Concentration of FFC in venous plasma leaving Rest in mg/L
      double CLu1 = CLu*1000*358.21;            // Concentration of FFC in Lung in mg/L
      double CVLu1 = CVLu*1000*358.21;          // Concentration of FFC in venous plasma leaving Lung in mg/L
      double CM1=CM*1000*358.21;                // Concentration of FFC in muscle in mg/L
      double CVM1=CVM*1000*358.21;              // Concentration of FFC venous plasma leaving in muscle in mg/L
      
      
      
      double CLm1 = CLm*1000*247.28;              // Concentration of FFA in liver in mg/L
      double CVLm1 = CVLm*1000*247.28;            // Concentration of FFA in venous plasma leaving liver in mg/L
      double CKm1 = CKm*1000*247.28;              // Concentration of FFA in Kidney in mg/L 
      double CVKm1 = CVKm*1000*247.28 ;           // Concentration of FFA in venous plasma leaving Kidney in mg/L
      double CFm1 = CFm*1000*247.28;              // Concentration of FFA in fat in mg/L
      double CVFm1 = CVFm*1000*247.28;            // Concentration of FFA in venous plasma leaving fat in mg/L
      double CRestm1 = CRestm*1000*247.28;        // Concentration of FFA in Rest in mg/L
      double CVRestm1 = CVRestm*1000*247.28;      // Concentration of FFA in venous plasma leaving Rest in mg/L
      double CLum1 = CLum*1000*247.28 ;           // Concentration of FFA in Lung in mg/L
      double CVLum1 = CVLum*1000*247.28;          // Concentration of FFA in venous plasma leaving Lung in mg/L
      double CMm1 =CMm*1000*247.28;               // Concentration of FFA in muscle in mg/L
      double CVMm1 =CVMm*1000*247.28 ;            // Concentration of FFA venous plasma leaving in muscle in mg/L        

             
             
     // Rate of FFC concentration changes in different compartments
      double Rurine     = Kurine*CVK;                            // Rurine:     mg/h; Rate of change in urine elimination
      double RK         = QK*(Ca-CVK)*Free - Rurine;             // RK:         mg/h; Rate of change in kidney 
      double RST        = - K0*AST - GE*AST ;           // RabsST:     mg/h; Rate of absorption in stomach          
      double RabsST     = K0*AST;                                // RabsSI:     mg/h; Rate of absorption in small intesine       
      double RSI        = GE*AST - Kabs*ASI - Kunabs*ASI;        // RSI:        mg/h; Rate of chnage in small intestine                 
      double  RabsSI     = Kabs*ASI;                              // RST:        mg/h; Rate of chnage in stomach      
      double RL         = QL*(Ca-CVL)*Free - Kbile*AL + 
        Kabs*ASI + K0*AST- FracFFA*Kffa*AL -(1-FracFFA)*Kmo*AL;  // RL:  mg/h; Rate of change in liver compartment
      
      
      
      double RM         = QM*(Ca-CVM)*Free;                      // RM:         mg/h; Rate of change in muscle compartment
      double RF         = QF*(Ca-CVF)*Free;                      // RF:         mg/h; Rate of change in fat compartment
      double RLu        = QLu*(Cv-CVLu)*Free;                    // RLu:        mg/h; Rate of change in lung compartment
      double RRest      = QRest*(Ca - CVRest)*Free;              // RRest:      mg/h; Rate of change in rest of body                                
      double Rvb = (QRest*CVRest*Free) + (QK*CVK*Free) + 
        (QL*CVL*Free) + (QF*CVF*Free) +
        (QM*CVM*Free) + RIM + RSC - (QC*Cv*Free);          // Rvb:         mg/h; Rate of change in venous plasma
      double Rab=QC*(CVLu-Ca)*Free;                              // Rab:         mg/h; Rate of change in arterial plasma
      
      double Rfeces     = Kbile*AL + Kunabs*ASI;                 // Rfeces:      mg/h; Rate of change in feces elimination
      
      
      
      
      // Rate of FFA concentration changes in different compartments
      
      double Rffa        = FracFFA*Kffa*AL;                          // Rffa:        mg/h; Rate of FFA input to the liver 
      double RotherMO    = (1-FracFFA)*Kmo*AL;
      double Rurinem     = Kurinem*CVKm;                             // Rurinem:     mg/h; Rate of change in urine elimination
      double RKm         = QK*(Cam-CVKm)*Freem - Rurinem;            // RKm:         mg/h; Rate of change in kidney 
      double RLm         = QL*(Cam-CVLm)*Freem + Rffa - Kbilem*ALm;  // RLm:         mg/h; Rate of change in liver compartment
      double RMm         = QM*(Cam-CVMm)*Freem;                      // RMm:         mg/h; Rate of change in muscle compartment
      double RFm         = QF*(Cam-CVFm)*Freem;                      // RFm:         mg/h; Rate of change in fat compartment
      double RLum        = QLu*(Cvm-CVLum)*Freem;                    // RLum:        mg/h; Rate of change in lung compartment
      double RRestm      = QRest*(Cam - CVRestm)*Freem;              // RRestm:      mg/h; Rate of change in rest of body                                
      double Rvbm = (QRest*CVRestm*Freem) + (QK*CVKm*Freem) + 
        (QL*CVLm*Freem) + (QF*CVFm*Freem) +
        (QM*CVMm*Freem) - (QC*Cvm*Freem);                            // Rvbm:         mg/h; Rate of change in venous plasma
      double Rabm=QC*(CVLum-Cam)*Freem;                              // Rabm:         mg/h; Rate of change in arterial plasma
      
      
     // ODE Equation for the amount of FFC in tissues
            
            
      dxdt_ARest     = RRest;              // amount of FFC in rest of the body                                       
      dxdt_AL        = RL;                 // amount of FFC in liver
      dxdt_AST       = RST;                // amount of FFC in stomach                                           
      dxdt_ASI       = RSI;                // amount of FFC in small intestine                                            
      dxdt_AabsSI    = RabsSI;             // amount of FFC absorbed in small intestine
      dxdt_AabsST    = RabsST;             // amount of FFC absorbed in small intestine
      dxdt_AK        = RK;                 // amount of FFC in kidney                                   
      dxdt_Aurine    = Rurine;             // amount of FFC in urine                              
      dxdt_Afeces    = Rfeces;             // amount of FFC in feces     
      dxdt_AF        = RF;                 // amount of FFC in fat     
      dxdt_AM        = RM;                 // amount of FFC in muscle     
      dxdt_ALu       = RLu;                // amount of FFC in lung      
      dxdt_Aab       = Rab;                // amount of FFC in arterial plasma      
      dxdt_Avb       = Rvb;                // amount of FFC in venous plasma
      
      
      
      // ODE Equation for the amount of FFA in tissues
      
      
      dxdt_Affa       = Rffa;
      dxdt_AotherMO   = RotherMO;
      dxdt_ARestm     = RRestm;              // amount of FFA in rest of the body                                       
      dxdt_ALm        = RLm;                 // amount of FFA in liver
      dxdt_AKm        = RKm;                 // amount of FFA in kidney                                   
      dxdt_Aurinem    = Rurinem;             // amount of FFA in urine                              
      dxdt_AFm        = RFm;                 // amount of FFA in fat     
      dxdt_AMm        = RMm;                 // amount of FFA in muscle     
      dxdt_ALum       = RLum;                // amount of FFA in lung      
      dxdt_Aabm       = Rabm;                // amount of FFA in arterial plasma      
      dxdt_Avbm       = Rvbm;                // amount of FFA in venous plasma        

      
    
      // Equation for the AUC of FFC in tissue compartment
      
      dxdt_AUCK      = CK;                 // AUC of FFC in kidney compartment                         
      dxdt_AUCCR     = CRest;              // AUC of FFC in rest of the body compartment                              
      dxdt_AUCCL     = CL;                 // AUC of FFC in liver compartment  
      dxdt_AUCCF     = CF;                 // AUC of FFC in fat compartment
      dxdt_AUCCM     = CM;                 // AUC of FFC in muscle compartment
      dxdt_AUCCLu    = CLu;                // AUC of FFC in lung compartment
      
      
      // Equation for the AUC of FFA in tissue compartment
      dxdt_AUCP       = Cv1;                 // AUC of FFA in Plasma compartment 
      dxdt_AUCKm      = CKm;                 // AUC of FFA in kidney compartment                         
      dxdt_AUCCRm     = CRestm;              // AUC of FFA in rest of the body compartment                              
      dxdt_AUCCLm     = CLm;                 // AUC of FFA in liver compartment  
      dxdt_AUCCFm     = CFm;                 // AUC of FFA in fat compartment
      dxdt_AUCCMm     = CMm;                 // AUC of FFA in muscle compartment
      dxdt_AUCCLum    = CLum;                // AUC of FFA in lung compartment

      // Mass Balance Check
      
      double Atissue =  ARest + AK + AL + AST + ASI + AF + AM + ALu + Aab + Avb;   // Amount of chemical in different tissue
      double Aloss   = Aurine + Afeces+ Affa+ AotherMO;                            // Amount of lost chemical via feces and Urine
      double Atotal  = Atissue + Aloss;                                           // Amount of total chemical in the body
      double BAL     = Absorbim + Absorbsc - Atotal;                   // Amount balance
      
      //  Mass Balance Check
      
      double Atissuem =  ARestm + AKm + ALm + AFm + AMm + ALum + Aabm + Avbm;       // Amount of chemical in different tissue
      double Alossm   = Aurinem;                                                    // Amount of lost chemical via feces and Urine
      double Atotalm  = Atissuem + Alossm;                                          // Amount of total chemical in the body
      double BALm     = Affa - Atotalm;                                             // Amount balance


$TABLE
capture Plasma = Cv1;
capture Liver  = CL1;
capture Kidney = CK1;
capture Lung  = CLu1;
capture Muscle = CM1;
capture Fat  = CF1;
capture Rest = CRest1;
capture Balance = BAL;

capture Plasma_FFA = Cvm1;
capture Liver_FFA  = CLm1;
capture Kidney_FFA = CKm1;
capture Lung_FFA  = CLum1;
capture Muscle_FFA = CMm1;
capture Fat_FFA  = CFm1;
capture Rest_FFA = CRestm1;
capture Balance_FFA = BALm;

capture AUC_P       = AUCP;
capture AUCC_R     = AUCCR;                 // AUC of FFA in rest of the body compartment                              
capture AUCC_L     = AUCCL;                 // AUC of FFA in liver compartment  
capture AUCC_F     = AUCCF;                 // AUC of FFA in fat compartment
capture AUCC_M     = AUCCM;                 // AUC of FFA in muscle compartment
capture AUCC_Lu    = AUCCLu;                // AUC of FFA in lung compartment

'


## Build mrgsolve-based PBPK Model

mod <- mcode_cache("pbpk", pbpk.code)





## Dosing fraction constant for slow and fast absorption
Frac.mean = 0.45                     ## 0.95 for conventional formulation and 0.5 for long-acting formulation Lin 2015
Fracsc.mean = 0.26                     ## 0.95 for conventional formulation and 0.5 for long-acting formulation; same as Lin 2017

#IM absorption rate constats
Kim.mean = 0.397                     # Absorption coefficient for IM dose fitted
Kdiss.mean = 0.118                # Dissolution coefficient for IM dose fitted

#SC absorption rate constats
Ksc.mean = 0.208                      # 0.15 for conventional formulation# 0.3 for long-acting formulation Lin 2016
Kdisssc.mean = 0.009                  # 1/h Rate of chemical becoming available for absorption from the injection site (need to be fitted)
KmC.mean=0.002                        # Rate of metabolism 

## Cardiac Output and Blood flow
QCC.mean = 7.15                       # SD(2.47%)# Cardiac output (L/h/kg), Li 2021
QLC.mean = 0.3054                      # SD(0.445) # Fraction of flow to the liver, Li 2021
QKC.mean = 12.98/100                  # SD(5.67%) Fraction of flow to the kidneys, Li 2021
QFC.mean = 6.18/100                   # SD(4.75%) Fraction of flow to the fat, Li 2021
QMC.mean = 10.09/100                  # SD(4.206%) , Fraction of flow to the muscle, Li 2021
QLuC.mean = 1                         # Fraction of blood flow to the lung (Achenbach, 1995), Li 2021
Htc=34.69/100                    # SD (3.09%) #  Percentage of Hematocrit in the blood, Li 2021
PCV=36.15/100                    # Serum as a percentage of whole blood, Li 2021
QrestC.mean=1-(QLC.mean+QKC.mean+QMC.mean+QFC.mean)       # Blood flow in the rest of the body


####################### Tissue Volume ######################################################
############################################################################################

VBloodC.mean   = 0.0486               ## L/kg BW; Fraction vol. of plasma; Li,2021               
VLC.mean       = 1.27/100             ## L/kg BW; Fraction vol. of liver; LI 2021
VKC.mean       = 0.23/100             ## L/kg BW; Fraction vol. of kidney; LI2021
VFC.mean       = 20.98/100            ##  Fractional fat tissue, Li 2021
VMC.mean       = 24.78/100            ## Fractional muscle tissue, Li 2021
VLuC.mean      = 1.07/100             ## Fractional lung tissue, Li 2021
VabC.mean=VBloodC.mean*0.24                ## Fractional arterial blood, Yang 2019
VvbC.mean=VBloodC.mean*0.76                ## Fractional venous blood, Yang 2019
VRestC.mean=1-(VLC.mean+VKC.mean+VMC.mean+VFC.mean+VLuC.mean+VabC.mean+VvbC.mean)

############# Mass Transfer Parameters (Chemical-specific parameters) ########################
############### Partition coefficients(PC, tissue:plasma) ####################################


PLu.mean = 0.9                   ## Lung: plasma PC
PL.mean = 1.3                    ## Liver: plasma PC
PK.mean = 0.9                    ## Kidney:plasma PC
PM.mean = 0.9                   ## Muscle:plasma PC
PF.mean = 0.1137                    ## Fat:plasma PC
Prest.mean = 1.13                ## Rest of the body:plasma PC Yang et al. 2019


PLum.mean    = 12.741375             # Lung: plasma PC 
PLm.mean     = 12.741375             # Liver: plasma PC Yang et al. 2019
PKm.mean     = 2.902587              # Kidney:plasma PC Yang et al. 2019
PMm.mean     = 0.4778312             # Muscle:plasma PC Yang et al. 2019
PFm.mean     = 0.1                   # Fat:plasma PC Yang et al. 2019
Prestm.mean  = 0.9                   # Rest of the body:plasma PC



GEC.mean       = 0.034              ## 1/h/kg; Gastric emptying time; Lin 2016
K0C.mean       = 0.0017             ## 1/h/kg; Rate of uptake from the stomach into the liver Lin 2016
KabsC.mean     = 0.0002             ## 1/h/kg; Rate of absorption of chemical from small intestine to liver Lin 2016
KunabsC.mean   = 0.003              ## 1/h/kg; Rate of unabsorbed dose to appear in feces from small intestine
KbileC.mean    = 0.026              ## 1/h/kg; Biliary elimination rate from liver to feces   
KurineC.mean   = 1.391167           ## 1/h/kg; Rate of urine elimination from urine storage Lin 2016
KurineCm.mean   = 0.0042            ## 1/h/kg; Rate of urine elimination from urine storage Lin 2016
Free.mean      = 1                  ## Unitless; fraction of chemical free in blood after plasma protein binding


Kffa.mean      = 1.52               # Rate of transfer to florfenicol amine
Kmo.mean       = 0.1014             # Rate of transfer to other matabolite

KbileCm.mean   =   0                # Biliary excretion rate

FracFFA.mean   =   0.5              # Fractions of FFA 

PB.mean        =  0.2                    # Plasma protein binding coefficient for florfenicol
PBm.mean       =  0.2                    # Plasma protein binding coefficient for florfenicol amine




######################## Standard deviations for each parameters ######################################
#######################################################################################################

PC.CV=0.2                      ## Coefficient of variation for partition coefficients sd/mean Li 2017
K.CV=0.3                       ## Coefficient of variation for absorption/excretion parameters sd/mean Li 2017
P.CV=0.3                       ## Coefficient of variation for physiological parameters sd/mean Li 2017
Frac.CV=0.1                    ## Coefficient of variation for fractional parameters sd/mean Li 2017

QCC.sd = 2.47                  ## Standard deviation for cardiac output Li 2021              
QLC.sd = 0.445/100             ## Standard deviation for fractions of blood flow to liver Li 2021
QKC.sd = 5.6/100               ## Standard deviation for fractions of blood flow to Kidney Li 2021
QMC.sd = 4.206/100             ## Standard deviation for fractions of blood flow to Muscle Li 2021
QFC.sd = P.CV*QFC.mean              ## Standard deviation for fractions of blood flow to fat
QrestC.sd = P.CV*QrestC.mean        ## Standard deviation for fractions of blood flow to rest of the body


BW.sd = 5.18                   ## Standard deviation for bodyweight Li 2021
VLC.sd = 0.24/100              ## Standard deviation for fractions of volume of liver Li 2021
VKC.sd = 0.03/100              ## Standard deviation for fractions of volume of kidney Li 2021
VMC.sd = 2.2/100               ## Standard deviation for fractions of volume of muscle Li 2021
VFC.sd = 2.68/100              ## Standard deviation for fractions of volume of fat Li 2021
VLuC.sd = 0.24/100             ## Standard deviation for fractions of volume of lung Li 2021
VRestC.sd = P.CV*VRestC.mean       ## Standard deviation for fractions of volume of rest of the body
VvbC.sd = P.CV*VvbC.mean           ## Standard deviation for fractions of volume of arterial blood
VabC.sd = P.CV*VabC.mean           ## Standard deviation for fractions of volume of venous blood


PL.sd = PC.CV*PL.mean               ## Standard deviation for liver:Plasma Partition coefficient
PK.sd = PC.CV*PK.mean               ## Standard deviation for kidney:Plasma Partition coefficient
PM.sd =PC.CV*PM.mean                ## Standard deviation for muscle:Plasma Partition coefficient
PF.sd = PC.CV*PF.mean               ## Standard deviation for fat:Plasma Partition coefficient
PLu.sd = PC.CV*PLu.mean             ## Standard deviation for lung:Plasma Partition coefficient
Prest.sd = PC.CV*Prest.mean         ## Standard deviation for rest of the body:Plasma Partition coefficient



PLm.sd = PC.CV*PLm.mean               ## Standard deviation for liver:Plasma Partition coefficient
PKm.sd = PC.CV*PKm.mean               ## Standard deviation for kidney:Plasma Partition coefficient
PMm.sd =PC.CV*PMm.mean                ## Standard deviation for muscle:Plasma Partition coefficient
PFm.sd = PC.CV*PFm.mean               ## Standard deviation for fat:Plasma Partition coefficient
PLum.sd = PC.CV*PLum.mean             ## Standard deviation for lung:Plasma Partition coefficient
Prestm.sd = PC.CV*Prestm.mean         ## Standard deviation for rest of the body:Plasma Partition coefficient


Kim.sd = K.CV*Kim.mean              ## Standard deviation for absorption coefficient for IM dose
Frac.sd = Frac.CV*Frac.mean         ## Standard deviation for fractions available for immediate absorption during IM dose
Kdiss.sd = K.CV*Kdiss.mean          ## Standard deviation for dissolution coefficient for IM dose
Ksc.sd = K.CV*Ksc.mean              ## Standard deviation for absorption coefficient for SC dose
Fracsc.sd = Frac.CV*Fracsc.mean     ## Standard deviation for fractions available for immediate absorption during SC dose
Kdisssc.sd = K.CV*Kdisssc.mean      ## Standard deviation for dissolution coefficient for SC dose
KmC.sd = K.CV*KmC.mean              ## Standard deviation for metabolism in stomach for Oral dose
KurineC.sd = K.CV*KurineC.mean      ## Standard deviation for urinary excretion rate
KurineCm.sd = K.CV*KurineCm.mean    ## Standard deviation for urinary excretion rate


Kffa.sd      = Kffa.mean*K.CV       # Rate of transfer to florfenicol amine
Kmo.sd       = Kmo.mean*K.CV    # Rate of transfer to other matabolite

KbileCm.mean   =   0                # Biliary excretion rate
KbileCm.mean   =   0                # Biliary excretion rate

FracFFA.sd= FracFFA.mean* Frac.CV   # Fractions of FFA 

PB.sd= PB.mean* Frac.CV             # Plasma protein binding coefficient for florfenicol
PBm.sd= PBm.mean* Frac.CV           # Plasma protein binding coefficient for florfenicol amine




m.log.PL= log(PL.mean^2/(PL.sd^2+PL.mean^2)^0.5) 

sd.log.PL = (log(1+PL.sd^2/PL.mean^2))^0.5 

m.log.PK = log(PK.mean^2/(PK.sd^2+PK.mean^2)^0.5) 

sd.log.PK = (log(1+PK.sd^2/PK.mean^2))^0.5 

m.log.PM = log(PM.mean^2/(PM.sd^2+PM.mean^2)^0.5)

sd.log.PM = (log(1+PM.sd^2/PM.mean^2))^0.5

m.log.PF= log(PF.mean^2/(PF.sd^2+PF.mean^2)^0.5)

sd.log.PF = (log(1+PF.sd^2/PF.mean^2))^0.5

m.log.PLu= log(PLu.mean^2/(PLu.sd^2+PLu.mean^2)^0.5)

sd.log.PLu = (log(1+PLu.sd^2/PLu.mean^2))^0.5

m.log.Prest= log(Prest.mean^2/(Prest.sd^2+Prest.mean^2)^0.5)

sd.log.Prest = (log(1+Prest.sd^2/Prest.mean^2))^0.5




m.log.PLm= log(PLm.mean^2/(PLm.sd^2+PLm.mean^2)^0.5) 

sd.log.PLm = (log(1+PLm.sd^2/PLm.mean^2))^0.5 

m.log.PKm = log(PKm.mean^2/(PKm.sd^2+PKm.mean^2)^0.5) 

sd.log.PKm = (log(1+PKm.sd^2/PKm.mean^2))^0.5 

m.log.PMm = log(PMm.mean^2/(PMm.sd^2+PMm.mean^2)^0.5)

sd.log.PMm = (log(1+PMm.sd^2/PMm.mean^2))^0.5

m.log.PFm= log(PFm.mean^2/(PFm.sd^2+PFm.mean^2)^0.5)

sd.log.PFm = (log(1+PFm.sd^2/PFm.mean^2))^0.5

m.log.PLum= log(PLum.mean^2/(PLum.sd^2+PLum.mean^2)^0.5)

sd.log.PLum = (log(1+PLum.sd^2/PLum.mean^2))^0.5

m.log.Prestm= log(Prestm.mean^2/(Prestm.sd^2+Prestm.mean^2)^0.5)

sd.log.Prestm = (log(1+Prestm.sd^2/Prestm.mean^2))^0.5





m.log.Kim= log(Kim.mean^2/(Kim.sd^2+Kim.mean^2)^0.5)

sd.log.Kim = (log(1+Kim.sd^2/Kim.mean^2))^0.5

m.log.Fracim= log(Frac.mean^2/(Frac.sd^2+Frac.mean^2)^0.5)

sd.log.Fracim = (log(1+Frac.sd^2/Frac.mean^2))^0.5

m.log.Kdiss= log(Kdiss.mean^2/(Kdiss.sd^2+Kdiss.mean^2)^0.5) 

sd.log.Kdiss = (log(1+Kdiss.sd^2/Kdiss.mean^2))^0.5

m.log.Ksc= log(Ksc.mean^2/(Ksc.sd^2+Ksc.mean^2)^0.5)

sd.log.Ksc = (log(1+Ksc.sd^2/Ksc.mean^2))^0.5

m.log.Fracsc= log(Fracsc.mean^2/(Fracsc.sd^2+Fracsc.mean^2)^0.5) 

sd.log.Fracsc = (log(1+Fracsc.sd^2/Fracsc.mean^2))^0.5 

m.log.Kdisssc= log(Kdisssc.mean^2/(Kdisssc.sd^2+Kdisssc.mean^2)^0.5)

sd.log.Kdisssc = (log(1+Kdisssc.sd^2/Kdisssc.mean^2))^0.5

m.log.KmC= log(KmC.mean^2/(KmC.sd^2+KmC.mean^2)^0.5)

sd.log.KmC = (log(1+KmC.sd^2/KmC.mean^2))^0.5

m.log.KurineC= log(KurineC.mean^2/(KurineC.sd^2+KurineC.mean^2)^0.5)

sd.log.KurineC = (log(1+KurineC.sd^2/KurineC.mean^2))^0.5 

m.log.KurineCm= log(KurineCm.mean^2/(KurineCm.sd^2+KurineCm.mean^2)^0.5)

sd.log.KurineCm = (log(1+KurineCm.sd^2/KurineCm.mean^2))^0.5 



m.log.PB= log(PB.mean^2/(PB.sd^2+PB.mean^2)^0.5)

sd.log.PB = (log(1+PB.sd^2/PB.mean^2))^0.5 


m.log.PBm= log(PBm.mean^2/(PBm.sd^2+PBm.mean^2)^0.5)

sd.log.PBm = (log(1+PBm.sd^2/PBm.mean^2))^0.5 


m.log.Kmo= log(Kmo.mean^2/(Kmo.sd^2+Kmo.mean^2)^0.5)

sd.log.Kmo = (log(1+Kmo.sd^2/Kmo.mean^2))^0.5 


m.log.Kffa= log(Kffa.mean^2/(Kffa.sd^2+Kffa.mean^2)^0.5)

sd.log.Kffa = (log(1+Kffa.sd^2/Kffa.mean^2))^0.5 


m.log.FracFFA= log(FracFFA.mean^2/(FracFFA.sd^2+FracFFA.mean^2)^0.5)

sd.log.FracFFA = (log(1+FracFFA.sd^2/FracFFA.mean^2))^0.5 


sd.log.KmC <-
  sqrt(log(1 + KmC.sd ^ 2 / KmC.mean ^ 2)) # standard deviation of lognormal distribution
m.log.KmC <-
  log(KmC.mean) - 0.5 * sd.log.KmC ^ 2 # mean of lognormal distribution


set.seed(324)#+i) # set random seed so that the simulation result is reproducible, because randomly generated data is same if you set same random seed.
idata <- 
  data_frame(ID=1:N) %>% 
  mutate(
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QCC.mean, sd = QCC.sd),
      b = qnorm(0.975, mean = QCC.mean, sd = QCC.sd),
      mean = QCC.mean,
      sd = QCC.sd
    ),  # Cardiac output index (L/h/kg)
    # Fracion of blood flow to organs (unitless) 
    QLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QLC.mean, sd = QLC.sd),
      b = qnorm(0.975, mean = QLC.mean, sd = QLC.sd),
      mean = QLC.mean,
      sd = QLC.sd
    ),
    # Fraction of blood flow to the kidneys (2016 Lin)
    QKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QKC.mean, sd =  QKC.sd),
      b = qnorm(0.975, mean = QKC.mean, sd =  QKC.sd),
      mean = QKC.mean,
      sd =  QKC.sd
    ),
    #QMC = 0.180			# Fraction of blood flow to the muscle (2016 Lin)
    QMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QMC.mean, sd =  QMC.sd),
      b = qnorm(0.975, mean = QMC.mean, sd =  QMC.sd),
      mean = QMC.mean,
      sd =  QMC.sd
    ),
    # Fraction of blood flow to the fat (2016 Lin)
    QFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QFC.mean, sd =  QFC.sd),
      b = qnorm(0.975, mean = QFC.mean, sd =  QFC.sd),
      mean = QFC.mean,
      sd =  QFC.sd
    ),
    # Fraction of blood flow to the rest of body (total sum equals to 1)
    QRestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QrestC.mean, sd = QrestC.sd),
      b = qnorm(0.975, mean = QrestC.mean, sd = QrestC.sd),
      mean = QrestC.mean,
      sd = QrestC.sd
    ),
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = BW.mean, sd = BW.sd), 
                    b = qnorm(0.975, mean = BW.mean, sd = BW.sd), 
                    mean = BW.mean, sd = BW.sd), 
    # Fractional organ tissue volumes (unitless)
    # Fractional liver tissue (1933 Swett)
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VLC.mean, sd =  VLC.sd),
      b = qnorm(0.975, mean = VLC.mean, sd =  VLC.sd) ,
      mean = VLC.mean,
      sd =  VLC.sd
    ),
    # Fractional kidney tissue (1933 Swett)
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VKC.mean, sd =  VKC.sd),
      b = qnorm(0.975, mean = VKC.mean, sd =  VKC.sd),
      mean = VKC.mean,
      sd =  VKC.sd
    ),
    # Fractional fat tissue (2016 Lin, 2014 Leavens)
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VFC.mean, sd =  VFC.sd),
      b = qnorm(0.975, mean = VFC.mean, sd =  VFC.sd),
      mean = VFC.mean,
      sd =  VFC.sd
    ),
    # Fractional muscle tissue (2016 Lin, 2014 Leavens)
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VMC.mean, sd =  VMC.sd),
      b = qnorm(0.975, mean = VMC.mean, sd =  VMC.sd),
      mean = VMC.mean,
      sd =  VMC.sd
    ),
    # Fractional lung tissue (2016 Lin, 2014 Leavens)
    VLuCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VLuC.mean, sd =  VLuC.sd),
      b = qnorm(0.975, mean = VLuC.mean, sd =  VLuC.sd),
      mean = VLuC.mean,
      sd =  VLuC.sd
    ),
    # Venous blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VvbCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VvbC.mean, sd =  VvbC.sd),
      b = qnorm(0.975, mean = VvbC.mean, sd =  VvbC.sd),
      mean = VvbC.mean,
      sd =  VvbC.sd
    ),
    # Arterial blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VabCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VabC.mean, sd =  VabC.sd),
      b = qnorm(0.975, mean = VabC.mean, sd =  VabC.sd),
      mean = VabC.mean,
      sd =  VabC.sd
    ),
    # Fractional rest of body (total sum equals to 1)
    VRestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VRestC.mean, sd =  VRestC.sd),
      b = qnorm(0.975, mean = VRestC.mean, sd =  VRestC.sd),
      mean = VRestC.mean,
      sd =  VRestC.sd
    ),
    PL = rlnormTrunc(
      N,
      meanlog = m.log.PL,
      sdlog = sd.log.PL,
      min = qlnorm(0.025, meanlog = m.log.PL, sdlog = sd.log.PL),
      max = qlnorm(0.975, meanlog = m.log.PL, sdlog = sd.log.PL)
    ),
    PK = rlnormTrunc(
      N,
      meanlog = m.log.PK,
      sdlog = sd.log.PK,
      min = qlnorm(0.025, meanlog = m.log.PK, sdlog = sd.log.PK),
      max = qlnorm(0.975, meanlog = m.log.PK, sdlog = sd.log.PK)
    ),
    PM = rlnormTrunc(
      N,
      meanlog = m.log.PM,
      sdlog = sd.log.PM,
      min = qlnorm(0.025, meanlog = m.log.PM, sdlog = sd.log.PM),
      max = qlnorm(0.975, meanlog = m.log.PM, sdlog = sd.log.PM)
    ),
    PF = rlnormTrunc(
      N,
      meanlog = m.log.PF,
      sdlog = sd.log.PF,
      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF),
      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)
    ),
    PLu = rlnormTrunc(
      N,
      meanlog = m.log.PLu,
      sdlog = sd.log.PLu,
      min = qlnorm(0.025, meanlog = m.log.PLu, sdlog = sd.log.PLu),
      max = qlnorm(0.975, meanlog = m.log.PLu, sdlog = sd.log.PLu)
    ),
    Prest = rlnormTrunc(
      N,
      meanlog = m.log.Prest,
      sdlog = sd.log.Prest,
      min = qlnorm(0.025, meanlog = m.log.Prest, sdlog = sd.log.Prest),
      max = qlnorm(0.975, meanlog = m.log.Prest, sdlog = sd.log.Prest)
    ),
    
    PLm = rlnormTrunc(
      N,
      meanlog = m.log.PLm,
      sdlog = sd.log.PLm,
      min = qlnorm(0.025, meanlog = m.log.PLm, sdlog = sd.log.PLm),
      max = qlnorm(0.975, meanlog = m.log.PLm, sdlog = sd.log.PLm)
    ),
    PKm = rlnormTrunc(
      N,
      meanlog = m.log.PKm,
      sdlog = sd.log.PKm,
      min = qlnorm(0.025, meanlog = m.log.PKm, sdlog = sd.log.PKm),
      max = qlnorm(0.975, meanlog = m.log.PKm, sdlog = sd.log.PKm)
    ),
    PMm = rlnormTrunc(
      N,
      meanlog = m.log.PMm,
      sdlog = sd.log.PMm,
      min = qlnorm(0.025, meanlog = m.log.PMm, sdlog = sd.log.PMm),
      max = qlnorm(0.975, meanlog = m.log.PMm, sdlog = sd.log.PMm)
    ),
    PFm = rlnormTrunc(
      N,
      meanlog = m.log.PFm,
      sdlog = sd.log.PFm,
      min = qlnorm(0.025, meanlog = m.log.PFm, sdlog = sd.log.PFm),
      max = qlnorm(0.975, meanlog = m.log.PFm, sdlog = sd.log.PFm)
    ),
    PLum = rlnormTrunc(
      N,
      meanlog = m.log.PLum,
      sdlog = sd.log.PLum,
      min = qlnorm(0.025, meanlog = m.log.PLum, sdlog = sd.log.PLum),
      max = qlnorm(0.975, meanlog = m.log.PLum, sdlog = sd.log.PLum)
    ),
    Prestm = rlnormTrunc(
      N,
      meanlog = m.log.Prestm,
      sdlog = sd.log.Prestm,
      min = qlnorm(0.025, meanlog = m.log.Prestm, sdlog = sd.log.Prestm),
      max = qlnorm(0.975, meanlog = m.log.Prestm, sdlog = sd.log.Prestm)
    ),
    
    
    #{Kinetic Constants}
    # IM Absorption Rate Constants
    # /h, IM absorption rate constant
    Kim = rlnormTrunc(
      N,
      meanlog = m.log.Kim,
      sdlog = sd.log.Kim,
      min = qlnorm(0.025, meanlog = m.log.Kim, sdlog = sd.log.Kim),
      max = qlnorm(0.975, meanlog = m.log.Kim, sdlog = sd.log.Kim)
    ),
    Fracim = rlnormTrunc(
      N,
      meanlog = m.log.Fracim,
      sdlog = sd.log.Fracim,
      min = qlnorm(0.025, meanlog = m.log.Fracim, sdlog = sd.log.Fracim),
      max = qlnorm(0.975, meanlog = m.log.Fracim, sdlog = sd.log.Fracim)
    ),
    Kdiss = rlnormTrunc(
      N,
      meanlog = m.log.Kdiss,
      sdlog = sd.log.Kdiss,
      min = qlnorm(0.025, meanlog = m.log.Kdiss, sdlog = sd.log.Kdiss),
      max = qlnorm(0.975, meanlog = m.log.Kdiss, sdlog = sd.log.Kdiss)
    ),
    # SC Absorption Rate Constants
    # /h, SC absorption rate constant
    Ksc = rlnormTrunc(
      N,
      meanlog = m.log.Ksc,
      sdlog = sd.log.Ksc,
      min = qlnorm(0.025, meanlog = m.log.Ksc, sdlog = sd.log.Ksc),
      max = qlnorm(0.975, meanlog = m.log.Ksc, sdlog = sd.log.Ksc)
    ),
    Fracsc = rlnormTrunc(
      N,
      meanlog = m.log.Fracsc,
      sdlog = sd.log.Fracsc,
      min = qlnorm(0.025, meanlog = m.log.Fracsc, sdlog = sd.log.Fracsc),
      max = qlnorm(0.975, meanlog = m.log.Fracsc, sdlog = sd.log.Fracsc)
    ),
    Kdisssc = rlnormTrunc(
      N,
      meanlog = m.log.Kdisssc,
      sdlog = sd.log.Kdisssc,
      min = qlnorm(0.025, meanlog = m.log.Kdisssc, sdlog = sd.log.Kdisssc),
      max = qlnorm(0.975, meanlog = m.log.Kdisssc, sdlog = sd.log.Kdisssc)
    ),
    # Percentage Plasma Protein Binding unitless
    PB = rlnormTrunc(
      N,
      meanlog = m.log.PB,
      sdlog = sd.log.PB,
      min = qlnorm(0.025, meanlog = m.log.PB, sdlog = sd.log.PB),
      max = qlnorm(0.975, meanlog = m.log.PB, sdlog = sd.log.PB)
    ),
    
    PBm = rlnormTrunc(
      N,
      meanlog = m.log.PBm,
      sdlog = sd.log.PBm,
      min = qlnorm(0.025, meanlog = m.log.PBm, sdlog = sd.log.PBm),
      max = qlnorm(0.975, meanlog = m.log.PBm, sdlog = sd.log.PBm)
    ),
    
    #{Metabolic Rate Constant}
    KmC = rlnormTrunc(
      N,
      meanlog = m.log.KmC,
      sdlog = sd.log.KmC,
      min = qlnorm(0.025, meanlog = m.log.KmC, sdlog = sd.log.KmC),
      max = qlnorm(0.975, meanlog = m.log.KmC, sdlog = sd.log.KmC)
    ),
    
    FracFFA = rlnormTrunc(
      N,
      meanlog = m.log.FracFFA,
      sdlog = sd.log.FracFFA,
      min = qlnorm(0.025, meanlog = m.log.FracFFA, sdlog = sd.log.FracFFA),
      max = qlnorm(0.975, meanlog = m.log.FracFFA, sdlog = sd.log.FracFFA)),
    
    GEC       = 0,   #               1/(h); Gastric emptying time; 
    K0C       = 0,   #               1/(h); Rate of uptake from the stomach into the liver
    KabsC     = 0,       #               1/(h); Rate of absorption of chemical from small intestine to liver Lin 2015
    KunabsC   = 0,       #               1/(h); Rate of unabsorbed dose to appear in feces from small intestine
    
    KbileC    = 0,       #               1/(h*kg); Biliary elimination rate from liver to feces   
    
    Kffa = rlnormTrunc(
      N,
      meanlog = m.log.Kffa,
      sdlog = sd.log.Kffa,
      min = qlnorm(0.025, meanlog = m.log.Kffa, sdlog = sd.log.Kffa),
      max = qlnorm(0.975, meanlog = m.log.Kffa, sdlog = sd.log.Kffa)),
    
    Kmo = rlnormTrunc(
      N,
      meanlog = m.log.Kmo,
      sdlog = sd.log.Kmo,
      min = qlnorm(0.025, meanlog = m.log.Kmo, sdlog = sd.log.Kmo),
      max = qlnorm(0.975, meanlog = m.log.Kmo, sdlog = sd.log.Kmo)),
    
    
    
    
    KurineCm = rlnormTrunc(
      N,
      meanlog = m.log.KurineCm,
      sdlog = sd.log.KurineCm,
      min = qlnorm(0.025, meanlog = m.log.KurineCm, sdlog = sd.log.KurineCm),
      max = qlnorm(0.975, meanlog = m.log.KurineCm, sdlog = sd.log.KurineCm)),
    
    
    # Urinary Elimination Rate Constants
    KurineC = rlnormTrunc(
      N,
      meanlog = m.log.KurineC,
      sdlog = sd.log.KurineC,
      min = qlnorm(0.025, meanlog = m.log.KurineC, sdlog = sd.log.KurineC),
      max = qlnorm(0.975, meanlog = m.log.KurineC, sdlog = sd.log.KurineC)
    ), DOSEim = PDOSEim*BW*Fracim/1000/358.21, DOSEimslow = PDOSEim*BW*(1-Fracim)/1000/358.21, 
    DOSEsc = PDOSEsc*BW*Fracsc/1000/358.21, DOSEscslow = PDOSEsc*BW*(1-Fracsc)/1000/358.21,
    DOSEiv=PDOSEiv*BW/1000/358.21)



if (route == "iv") {
  
  ev1 <- ev (ID   = 1:N, amt  = idata$DOSEiv, ii = tinterval, tinf = 0.01, addl = TDoses-1, cmt  = "Avb", replicate = FALSE)
  ex =  ev1
  
}




if (route == "im") {
  ev1 <- ev(ID=1:N, amt= idata$DOSEim, ii=tinterval, addl=TDoses-1, tinf = 0.1, cmt="Amtsite", replicate = FALSE)
  ev2 <- ev(ID=1:N, amt= idata$DOSEimslow, ii=tinterval, addl=TDoses-1, tinf = 0.1, cmt="ADOSEim", replicate = FALSE)
  ex =  ev1 + ev2
}

if (route == "sc") {
  ev1 <- ev(ID=1:N, amt= idata$DOSEsc, ii=tinterval, addl=TDoses-1, tinf = 0.1, cmt="Asctsite", replicate = FALSE)
  ev2 <- ev(ID=1:N, amt= idata$DOSEscslow, ii=tinterval, addl=TDoses-1, tinf = 0.1, cmt="ADOSEsc", replicate = FALSE)
  ex =  ev1 + ev2
}



set.seed(11009)
{
  out <- 
    mod %>% 
    data_set(ex) %>%
    update(atol = 1E-15,maxsteps=50000) %>%
    idata_set(idata) %>%
    mrgsim(obsonly=TRUE, tgrid=tsamp)
}



out.sum = out %>% as_data_frame %>%
  mutate(Time1 = time/24 - (TDoses-1)*tinterval/24) %>%
  dplyr::select(ID, Time1, Plasma:Fat) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma, 0.01), CV50 = median(Plasma), CV99 = quantile(Plasma, 0.99), 
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL99 = quantile(Liver, 0.99), 
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK99 = quantile(Kidney, 0.99), 
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM99 = quantile(Muscle, 0.99), 
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF99 = quantile(Fat, 0.99))



out.sum.FFA = out %>% as_data_frame %>%
  mutate(Time1 = time/24 - (TDoses-1)*tinterval/24) %>%
  dplyr::select(ID, Time1, Plasma_FFA:Fat_FFA) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma_FFA, 0.01), CV50 = median(Plasma_FFA), CV99 = quantile(Plasma_FFA, 0.99), 
                   CL1 = quantile(Liver_FFA, 0.01), CL50 = median(Liver_FFA), CL99 = quantile(Liver_FFA, 0.99), 
                   CK1 = quantile(Kidney_FFA, 0.01), CK50 = median(Kidney_FFA), CK99 = quantile(Kidney_FFA, 0.99), 
                   CM1 = quantile(Muscle_FFA, 0.01), CM50 = median(Muscle_FFA), CM99 = quantile(Muscle_FFA, 0.99), 
                   CF1 = quantile(Fat_FFA, 0.01), CF50 = median(Fat_FFA), CF99 = quantile(Fat_FFA, 0.99))



output_plot=as.data.frame(out.sum)

output_plot_FFA=as.data.frame(out.sum.FFA)






p<-ggplot(data=output_plot, aes(x=Time1, y=CL50))  + geom_line(size=1)
p1<-p+geom_ribbon(aes(ymin=CL1, ymax=CL99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_liver),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),10)+ # add tolerance line
  ggtitle("(A) FF in Liver") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + geom_vline(aes(xintercept = min(subset(output_plot, CL99<= Tol_liver & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                       size = 1, color = "red", linetype = 2,
                                                                                       show.legend = F)+
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10(limits = c(0.000000001, 1000), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")



p<-ggplot(data=output_plot, aes(x=Time1, y=CK50))  + geom_line(size=1)
p2<-p+geom_ribbon(aes(ymin=CK1, ymax=CK99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_kidney),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),10)+ geom_vline(aes(xintercept = min(subset(output_plot, CK99<= Tol_kidney & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                 size = 1, color = "red", linetype = 2,
                                                                                                                                                 show.legend = F)+ # add tolerance line
  ggtitle("(B) FF in Kidney") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")




p<-ggplot(data=output_plot, aes(x=Time1, y=CM50))  + geom_line(size=1)
p3<-p+geom_ribbon(aes(ymin=CM1, ymax=CM99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_muscle),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),10)+ geom_vline(aes(xintercept = min(subset(output_plot, CM99<= Tol_muscle & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                 size = 1, color = "red", linetype = 2,
                                                                                                                                                 show.legend = F)+ # add tolerance line# add tolerance line
  ggtitle("(C) FF in Muscle") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")



















######################################################################        FFA                  ##################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################







p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CL50))  + geom_line(size=1)
p4<-p+geom_ribbon(aes(ymin=CL1, ymax=CL99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_liver),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CL99<= Tol_liver & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                      size = 1, color = "red", linetype = 2,
                                                                                                                                                      show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(A) FFA in Sheep Liver") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")



p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CK50))  + geom_line(size=1)
p5<-p+geom_ribbon(aes(ymin=CK1, ymax=CK99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_kidney),color = 'black',size = 1, linetype = 'twodash', show.legend = T, legend = "EU_Tol") + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CK99<= Tol_kidney & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                                          size = 1, color = "red", linetype = 2,
                                                                                                                                                                          show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(B) FFA in Sheep Kidney") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")




p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CM50))  + geom_line(size=1)
p6<-p+geom_ribbon(aes(ymin=CM1, ymax=CM99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_muscle),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CM99<= Tol_muscle & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                       size = 1, color = "red", linetype = 2,
                                                                                                                                                       show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(C) FFA in Sheep Muscle") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")









p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CL50))  + geom_line(size=1)
p4_fsis<-p+geom_ribbon(aes(ymin=CL1, ymax=CL99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_liver_fsis),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CL99<= Tol_liver_fsis & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                           size = 1, color = "red", linetype = 2,
                                                                                                                                                           show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(A) FFA in Sheep Liver 2X20 IM") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10(limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                     axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                     legend.position = "none")



p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CK50))  + geom_line(size=1)
p5_fsis<-p+geom_ribbon(aes(ymin=CK1, ymax=CK99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_kidney_fsis),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CK99<= Tol_kidney_fsis & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                            size = 1, color = "red", linetype = 2,
                                                                                                                                                            show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(B) FFA in Sheep Kidney 2X20 IM") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")




p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CM50))  + geom_line(size=1)
p6_fsis<-p+geom_ribbon(aes(ymin=CM1, ymax=CM99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_muscle_fsis),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CM99<= Tol_muscle_fsis & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                            size = 1, color = "red", linetype = 2,
                                                                                                                                                            show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(B) FFA in Sheep Muscle 2X20 IM") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")



p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CL50))  + geom_line(size=1)
p4_tol<-p+geom_ribbon(aes(ymin=CL1, ymax=CL99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_liver_tol),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CL99<= Tol_liver_tol & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                          size = 1, color = "red", linetype = 2,
                                                                                                                                                          show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(G) FFA in Sheep Liver") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")



p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CK50))  + geom_line(size=1)
p5_tol<-p+geom_ribbon(aes(ymin=CK1, ymax=CK99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_kidney_tol),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CK99<= Tol_kidney_tol & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                           size = 1, color = "red", linetype = 2,
                                                                                                                                                           show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(H) FFA in Sheep Kidney") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                      axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                      axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                      strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                      legend.position = "none")




p<-ggplot(data=output_plot_FFA, aes(x=Time1, y=CM50))  + geom_line(size=1)
p6_tol<-p+geom_ribbon(aes(ymin=CM1, ymax=CM99), linetype=2, alpha=0.01, colour="blue", size=1)+ 
  geom_line(aes(y = Tol_muscle_tol),color = 'black',size = 1, linetype = 'twodash', show.legend = T) + xlim(-(TDoses*tinterval/24-1),end_time)+ geom_vline(aes(xintercept = min(subset(output_plot_FFA, CM99<= Tol_muscle_tol & Time1 > 1) %>% dplyr::select(Time1))), 
                                                                                                                                                           size = 1, color = "red", linetype = 2,
                                                                                                                                                           show.legend = F)+ # add tolerance line # add tolerance line
  ggtitle("(I) FFA in Sheep Muscle") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +scale_shape_prism() + 
  theme_prism() + 
  theme() +scale_y_log10( limits = c(0.000000001, 100), expand = c(0, 0))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                             axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                             axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                             strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                             strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                             legend.position = c(0.8, 0.2), legend.title = element_text(colour="blue", size=10, 
                                                                                                                                        face="bold"), legend.text = element_text(colour="blue", size=10, 
                                                                                                                                                                                 face="bold"))
windows()

figure <- ggarrange( #p4, p5, p6,
                    p4_fsis,  p6_fsis,nrow=1)
                    #p4_tol, p5_tol, p6_tol)

figure


ggsave(
  "Figure 41_Sheep2X20IM.TIFF",
  plot = figure,
  width = 10,
  height = 3,
  dpi = 300
)

WDI_Liver_EU=output_plot_FFA$Time1[min(which(output_plot_FFA$CL99<= Tol_liver & output_plot_FFA$Time1>1))]
WDI_Liver_FSIS=output_plot_FFA$Time1[min(which(output_plot_FFA$CL99<= Tol_liver_fsis & output_plot_FFA$Time1>1))]
WDI_Liver_tol=output_plot_FFA$Time1[min(which(output_plot_FFA$CL99<= Tol_liver_tol & output_plot_FFA$Time1>1))]


WDI_kidney_EU=output_plot_FFA$Time1[min(which(output_plot_FFA$CK99<= Tol_kidney & output_plot_FFA$Time1>1))]
WDI_kidney_FSIS=output_plot_FFA$Time1[min(which(output_plot_FFA$CK99<= Tol_kidney_fsis & output_plot_FFA$Time1>1))]
WDI_kidney_tol=output_plot_FFA$Time1[min(which(output_plot_FFA$CK99<= Tol_kidney_tol & output_plot_FFA$Time1>1))]


WDI_muscle_EU=output_plot_FFA$Time1[min(which(output_plot_FFA$CM99<= Tol_muscle & output_plot_FFA$Time1>1))]
WDI_muscle_FSIS=output_plot_FFA$Time1[min(which(output_plot_FFA$CM99<= Tol_muscle_fsis & output_plot_FFA$Time1>1))]
WDI_muscle_tol=output_plot_FFA$Time1[min(which(output_plot_FFA$CM99<= Tol_muscle_tol & output_plot_FFA$Time1>1))]

WDIs<- cbind.data.frame(WDI_Liver_EU, WDI_Liver_FSIS, WDI_Liver_tol, WDI_kidney_EU, WDI_kidney_FSIS, WDI_kidney_tol, WDI_muscle_EU, WDI_muscle_FSIS, WDI_muscle_tol)


File="WDIs_IM2X20_Sheep.csv"

write.csv(ceiling(WDIs), File)
