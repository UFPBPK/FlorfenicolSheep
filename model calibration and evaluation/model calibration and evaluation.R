## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
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

SolvePBPKFF <- '
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
  
  PLum      : 1.3026    :                      Lung: plasma PC 
  PLm       : 8.7414    :                      Liver: plasma PC Yang et al. 2019
  PKm       : 1.3026    :                      Kidney:plasma PC Yang et al. 2019
  PMm       : 0.2778    :                      Muscle:plasma PC Yang et al. 2019
  PFm       : 0.1       :                      Fat:plasma PC Yang et al. 2019
  PRestm    : 0.0121    :                      Rest of the body:plasma PC
  
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



library(mrgsolve)  

## Build mrgsolve-based PBPK Model
mod <- mcode ("MrgSolve", SolvePBPKFF)



data_Palma_2011_IV_20_FFC<-as.data.frame(read.csv(file="Palma_2011_IV_20_FFC_FFA.csv"))       ### Time in hours 
data_Palma_2011_IV_20_FFC<- cbind.data.frame(data_Palma_2011_IV_20_FFC$Time, data_Palma_2011_IV_20_FFC$Plasma_FF, data_Palma_2011_IV_20_FFC$upper)
names(data_Palma_2011_IV_20_FFC) = c("Time", "Plasma", "Plasma_sd")


data_Balcomb_SC_40_FFC<-read.csv(file="Balcomb_SC_40_FFC.csv")                     ### Time in hours, single SC dose 40 mg/kg
data_Balcomb_SC_40_FFC<- cbind.data.frame(data_Balcomb_SC_40_FFC$Time, data_Balcomb_SC_40_FFC$Plasma_FormA, data_Balcomb_SC_40_FFC$Plasma_FormAsd)
names(data_Balcomb_SC_40_FFC) = c("Time", "Plasma", "Plasma_sd")

data_Balcomb_IM_20_FFC<- read.csv(file="Balcomb_IM_20_FFC.csv")
data_Balcomb_IM_20_FFC<- cbind.data.frame(data_Balcomb_IM_20_FFC$Time, data_Balcomb_IM_20_FFC$Plasma, data_Balcomb_IM_20_FFC$upper, data_Balcomb_IM_20_FFC$lower)
names(data_Balcomb_IM_20_FFC) = c("Time", "Plasma","upper","lower")

data_El_Sheikh_IV_IM_SC_20_FFC<-read.csv(file="El_Sheikh_IV_IM_SC_20_FFC.csv")     ### Time in hours, single dose 20 mg/kg
data_El_Sheikh_IV_IM_SC_20_FFC_IV<- read.csv(file="El_Sheikh_IV_20_FFC.csv") 
names(data_El_Sheikh_IV_IM_SC_20_FFC_IV) = c("Time", "Plasma", "Plasma_sd")

data_El_Sheikh_IV_IM_SC_20_FFC_IM<- read.csv(file="El_Sheikh_IM_20_FFC.csv") 
names(data_El_Sheikh_IV_IM_SC_20_FFC_IM) = c("Time", "Plasma", "Plasma_sd")


data_El_Sheikh_IV_IM_SC_20_FFC_SC<- read.csv(file="El_Sheikh_SC_20_FFC.csv") 
names(data_El_Sheikh_IV_IM_SC_20_FFC_SC) = c("Time", "Plasma", "Plasma_sd")


data_Jianzhong_IV_IM_20_30_FFC<-read.csv(file="Jianzhong_IV_IM_20_30_FFC.csv")     ### Time in minutes (Change it to hours, single IV and IM dose of 20 and 30 mg/kg)

data_Jianzhong_IV_20_FFC<- cbind.data.frame(data_Jianzhong_IV_IM_20_30_FFC$Time, data_Jianzhong_IV_IM_20_30_FFC$Plasma_iv20, data_Jianzhong_IV_IM_20_30_FFC$Plasma_iv20sd)
names(data_Jianzhong_IV_20_FFC) = c("Time", "Plasma", "Plasma_sd")

data_Jianzhong_IV_30_FFC<- cbind.data.frame(data_Jianzhong_IV_IM_20_30_FFC$Time, data_Jianzhong_IV_IM_20_30_FFC$Plasma_iv30, data_Jianzhong_IV_IM_20_30_FFC$Plasma_iv30sd)
names(data_Jianzhong_IV_30_FFC) = c("Time", "Plasma", "Plasma_sd")

data_Jianzhong_IM_30_FFC<- cbind.data.frame(data_Jianzhong_IV_IM_20_30_FFC$Time, data_Jianzhong_IV_IM_20_30_FFC$Plasma_im30, data_Jianzhong_IV_IM_20_30_FFC$Plasma_im30sd)
names(data_Jianzhong_IM_30_FFC) = c("Time", "Plasma", "Plasma_sd")
data_Jianzhong_IM_30_FFC<-data_Jianzhong_IM_30_FFC[-1,]


data_Jianzhong_IM_20_FFC<- cbind.data.frame(data_Jianzhong_IV_IM_20_30_FFC$Time, data_Jianzhong_IV_IM_20_30_FFC$Plasma_im20, data_Jianzhong_IV_IM_20_30_FFC$Plasma_im20sd)
names(data_Jianzhong_IM_20_FFC) = c("Time", "Plasma", "Plasma_sd")
data_Jianzhong_IM_20_FFC<-data_Jianzhong_IM_20_FFC[-1,]



################## Main Calibration data ########################################################################################
data_Wetzlich_2006_SC_3X40_FFC_Serum<-read.csv(file="Wetzlich_2006_SC_3X40_FFC_Serum.csv")        ### Time in hours


data_Wetzlich_2006_SC_3X40_FFA<-as.data.frame(read.csv(file="Wetzlich_2006_SC_3X40_FFA.csv"))     ### Time in days



data_Wetzlich_2006_SC_3X40_FFA_Liver<- cbind.data.frame(data_Wetzlich_2006_SC_3X40_FFA$Time+2, data_Wetzlich_2006_SC_3X40_FFA$Liver, data_Wetzlich_2006_SC_3X40_FFA$Liver_sd)
names(data_Wetzlich_2006_SC_3X40_FFA_Liver) = c("Time", "Liver", "Liver_sd")


data_Wetzlich_2006_SC_3X40_FFA_Kidney<- cbind.data.frame(data_Wetzlich_2006_SC_3X40_FFA$Time+2, data_Wetzlich_2006_SC_3X40_FFA$Kidney, data_Wetzlich_2006_SC_3X40_FFA$Kidney_sd)
names(data_Wetzlich_2006_SC_3X40_FFA_Kidney) = c("Time", "Kidney", "Kidney_sd")


data_Wetzlich_2006_SC_3X40_FFA_Muscle<- cbind.data.frame(data_Wetzlich_2006_SC_3X40_FFA$Time+2, data_Wetzlich_2006_SC_3X40_FFA$Muscle, data_Wetzlich_2006_SC_3X40_FFA$Muscle_sd)
names(data_Wetzlich_2006_SC_3X40_FFA_Muscle) = c("Time", "Muscle", "Muscle_sd")


data_Wetzlich_2006_SC_3X40_FFA_Fat<- cbind.data.frame(data_Wetzlich_2006_SC_3X40_FFA$Time+2, data_Wetzlich_2006_SC_3X40_FFA$Fat, data_Wetzlich_2006_SC_3X40_FFA$Fat_sd)
names(data_Wetzlich_2006_SC_3X40_FFA_Fat) = c("Time", "Fat", "Fat_sd")
data_Wetzlich_2006_SC_3X40_FFA_Fat<- data_Wetzlich_2006_SC_3X40_FFA_Fat




#################################################################################################################################
#################################################################################################################################



## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                  ## return a list of exp (parameters) from log domain
  
  ## Define the three exposure scenario
  ## Exposure scenario for oral exposue to 4.2 mg/kg-d
  
  PDose_IV.A  = 20 
  BW          = 25.56
  tinterval   = 24
  TDOSE       = 1
  dt          = 0.01
  DOSE_IV.A   =  PDose_IV.A*BW/1000/358.21
  
  
  
  ev_1.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Avb", replicate = FALSE)
  
  ev_2.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
  
  ex.iv.A <- ev_1.A+ev_2.A
  
  ## set up the exposure time
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)
  
  
  
  
  ## Get a prediction
  ## Out A: oral exposure to 4.2 mg/kg-d, matrix: Plasma, urine and feces
  out.iv <- 
    mod %>%                                               # model object
    param(pars) %>%                       # select model output
    update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.iv.A, tgrid=tsamp)               # Set up the simulation run
  outdf.iv.plasmaFF<-cbind.data.frame(Time=out.iv$time, 
                                      Plasma=out.iv$Plasma)
  
  cost<- modCost  (model=outdf.iv.plasmaFF,obs=data_El_Sheikh_IV_IM_SC_20_FFC_IV, err="Plasma_sd",x="Time")
  
  outdf.iv.plasmaFFA<-cbind.data.frame(Time=out.iv$time, 
                                       Plasma=out.iv$Plasma_FFA)
  
  cost<- modCost  (model=outdf.iv.plasmaFFA,obs=data_Palma_2011_IV_20_FFA, err="Plasma_sd", x="Time", cost=cost)
  
  Frac          =  0.45
  PDose_IM.B    =  20 
  DOSE_IM.B     =  PDose_IM.B*BW/1000/358.21
  Fracim.B      =  0.45
  DOSEfast.B    =  DOSE_IM.B*Frac             
  DOSEslow.B    =  DOSE_IM.B*(1-Frac)
  TDOSE=1
  
  
  ev_1.B <- ev (ID   = 1, amt  = DOSEfast.B, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Amtsite", replicate = FALSE)
  
  ev_2.B <- ev (ID   = 1, amt  = DOSEslow.B, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSEim", replicate = FALSE)
  
  ev_3.B <- ev (ID   = 1, amt  = DOSEfast.B+DOSEslow.B, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
  
  ex.im.B <- ev_1.B + ev_2.B + ev_3.B
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)
  
  
  
  out.im <- 
    mod %>%                                               # model object
    param(pars) %>%                       # select model output
    update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.im.B, tgrid=tsamp) 
  
  
  outdf.im.Plasma<-cbind.data.frame(Time=out.im$time, 
                                    Plasma=out.im$Plasma)
  
  cost<- modCost  (model=outdf.im.Plasma,obs=data_El_Sheikh_IV_IM_SC_20_FFC_IM, err="Plasma_sd",x="Time", cost=cost)
  cost<- modCost  (model=outdf.im.Plasma,obs=data_El_Sheikh_IV_IM_SC_20_FFC_IM, err="Plasma_sd",x="Time")
  
  
  Frac          = 0.26
  PDose_SC.C    = 20 
  DOSE_SC.C     =  PDose_SC.C*BW/1000/358.21
  DOSEfast.C    =  DOSE_SC.C*Frac             
  DOSEslow.C    =  DOSE_SC.C*(1-Frac) 
  TDOSE         =  1
  tinterval     = 24
  dt=0.01
  
  
  
  
  ev_1.C <- ev (ID   = 1, amt  = DOSEfast.C, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "Asctsite", replicate = FALSE)
  
  ev_2.C <- ev (ID   = 1, amt  = DOSEslow.C, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)
  
  
  
  ev_3.C <- ev (ID   = 1, amt  = DOSEfast.C+DOSEslow.C, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
  
  
  
  ex.sc.C <- ev_1.C + ev_2.C + ev_3.C
  
  ## set up the exposure time
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)
  
  out.sc <- 
    mod %>%                                               # model object
    param(pars) %>%                       # select model output
    update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.sc.C, tgrid=tsamp) 
  
  
  outdf.Plasma<-cbind.data.frame(Time=out.sc$time, 
                                 Plasma=out.sc$Plasma)
  
  
  
  
  
  cost<- modCost (model=outdf.Plasma,obs=data_El_Sheikh_IV_IM_SC_20_FFC_SC, err="Plasma_sd",x="Time")
  
  
  
  Frac          = 0.26
  PDose_SC.C    = 40 
  DOSE_SC.C     = PDose_SC.C*BW/1000/358.21
  DOSEfast.C    = DOSE_SC.C*Frac             
  DOSEslow.C    = DOSE_SC.C*(1-Frac) 
  TDOSE         = 3
  tinterval     = 24
  dt            = 0.01
  
  
  
  
  ev_1.C <- ev (ID   = 1, amt  = DOSEfast.C, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Asctsite", replicate = FALSE)
  
  ev_2.C <- ev (ID   = 1, amt  = DOSEslow.C, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)
  
  
  
  ev_3.C <- ev (ID   = 1, amt  = DOSEfast.C+DOSEslow.C, ii = tinterval,tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
  
  
  
  ex.sc.C <- ev_1.C + ev_2.C + ev_3.C
  
  ## set up the exposure time
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*45, dt)
  
  out.sc <- 
    mod %>%                                               # model object
    param(pars) %>%                       # select model output
    update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.sc.C, tgrid=tsamp) 
  
  
  outdf.Plasma<-cbind.data.frame(Time=out.sc$time, 
                                 Plasma=out.sc$Plasma)
  
  # Set up the simulation run
  outdf.Liver<-cbind.data.frame(Time=out.sc$time/24, 
                                Liver=out.sc$Liver_FFA)
  
  outdf.Kidney<-cbind.data.frame(Time=out.sc$time/24, 
                                 Kidney=out.sc$Kidney_FFA)
  
  
  outdf.Muscle<-cbind.data.frame(Time=out.sc$time/24, 
                                 Muscle=out.sc$Muscle_FFA)
  
  outdf.Fat<-cbind.data.frame(Time=out.sc$time/24, 
                              Fat=out.sc$Fat_FFA)
  
  
  
  
  cost<- modCost  (model=outdf.Liver,obs=data_Wetzlich_2006_SC_3X40_FFA_Liver, err="Liver_sd",x="Time")
  cost<- modCost  (model=outdf.Kidney,obs=data_Wetzlich_2006_SC_3X40_FFA_Kidney, err="Kidney_sd",x="Time", cost=cost)
  cost<- modCost  (model=outdf.Muscle,obs=data_Wetzlich_2006_SC_3X40_FFA_Muscle, err="Muscle_sd",x="Time", cost=cost)
  cost<- modCost  (model=outdf.Fat,obs=data_Wetzlich_2006_SC_3X40_FFA_Fat, err="Fat_sd",x="Time", cost=cost)
  
  
  return(cost)
  
}


## Cost fuction (FME) 
## Estimate the model residual by modCost function

MCcost<-function (pars){
  cost<- pred (pars)
  return(cost)
}


theta <- log(c(
  
  #PLu = 0.9,                        ## Lung: plasma PC 
  #PL = 1.3 ,                        ## Liver: plasma PC Yang et al. 2019
  #PK = 0.9,                         ## Kidney:plasma PC Yang et al. 2019
  #PM = 0.9,                         ## Muscle:plasma PC Yang et al. 2019
  #PF = 0.1137,                      ## Fat:plasma PC Yang et al. 2019
  PRest =  1.13,                  ## Rest of the body:plasma PC
  
  PLm = 12.741 ,                  ## Liver: plasma PC Yang et al. 2019
  PKm = 2.902,                   # Kidney:plasma PC Yang et al. 2019
  PMm = 0.478,                 # Muscle:plasma PC Yang et al. 2019
  PFm = 0.1,                       # Fat:plasma PC Yang et al. 2019
  PRestm = 0.9                    # Rest of the body:plasma PC
  
  ## Absorption and elimination parameters
  #GEC       = 0.182,                ## 1/(h); Gastric emptying time; 
  #K0C      = 1.991,                 ## 1/(h); Rate of uptake from the stomach into the liver
  #KabsC     = 0,                    ## 1/(h); Rate of absorption of chemical from small intestine to liver Lin 2015
  #KunabsC   = 0,                    ## 1/(h); Rate of unabsorbed dose to appear in feces from small intestine
  
  #KbileC    = 0,                    ## 1/(h*kg); Biliary elimination rate from liver to feces   
  #KurineC   = 1.45                  ## L/(h*kg); Rate of urine elimination from urine storage
  #KurineCm   = 0.0042,                ## L/(h*kg); Rate of urine elimination from urine storage
  #Kffa=1.52
  #Kmo=0.1014,
  
  #KbileCm=0,
  
  #FracFFA=0.5,
  
  #PB   =  0.2,
  #PBm  =  0.2
  
))

## PBPK model fitting 


#Fit <- modFit(f = MCcost, 
#p = theta,
#lower =theta*0.5, upper = theta*3,
#method ="Marq",control = nls.lm.control(nprint = 1))


#summary(Fit)
#theta=log(theta)


BW          = 58.58



PDose_IV.A  = 20 
tinterval   = 24
TDOSE      = 1
dt          = 0.05
DOSE_IV.A        =  PDose_IV.A*BW/1000/358.21



ev_1.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Avb", replicate = FALSE)

ev_2.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)

ex.iv.A <- ev_1.A+ev_2.A


## set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*40, dt)




## Get a prediction

out.A <- 
  mod %>%                                               # model object
  param((exp(theta))) %>%                       # select model output
  update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.iv.A, tgrid=tsamp)               # Set up the simulation run
outdf.A<-cbind.data.frame(Time=out.A$time[3: length(out.A$time)], 
                          Plasma=out.A$Plasma[3: length(out.A$time)],
                          Plasma_FFA=out.A$Plasma_FFA[3: length(out.A$time)])



p_IV_ElSheikh<- ggplot(data_El_Sheikh_IV_IM_SC_20_FFC_IV, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.A, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(A) FF in Sheep Plasma, 20 mg/kg IV dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  theme_prism() + 
  theme(legend.position = "none") + scale_y_log10(  limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), axis.text    
                                                                                                                                     = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                                     legend.position = "none")


p_IV_ElSheikh
cost_Calibration_sheep<-modCost(obs=data_El_Sheikh_IV_IM_SC_20_FFC_IV, model=outdf.A, x="Time", err="Plasma_sd")



p_Jianzhong_20_IV<- ggplot(data_Jianzhong_IV_20_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.A, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(J) FF in Sheep Plasma, 20 mg/kg IV dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                  axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                  axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                  strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  legend.position = "none")


p_Jianzhong_20_IV
cost_evaluation_sheep<-modCost(obs=data_Jianzhong_IV_20_FFC, model=outdf.A, x="Time", err="Plasma_sd")




p_Palma_20_IV_FFC<- ggplot(data_Palma_2011_IV_20_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma, ymax=Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.A, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(M) FF in Sheep Plasma, 20 mg/kg IV dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (Days)") +scale_shape_prism() + 
  theme_prism() + scale_y_log10(  limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                   axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                   axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                   strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   legend.position = "none")

p_Palma_20_IV_FFC
#cost_evaluation_sheep<-modCost(obs=data_Palma_2011_IV_20_FFC, model=outdf.A, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)


obs_mean <- data.frame(
  Time   = data_Palma_2011_IV_20_FFC$Time,
  Plasma = data_Palma_2011_IV_20_FFC$Plasma
)

cost_evaluation_sheep <- modCost(
  obs   = obs_mean,
  model = outdf.A,
  x     = "Time"
)
################################## 30 mg/kg IV dose #################################################

PDose_IV.A  = 30 
tinterval   = 24
TDOSE      = 1
dt          = 0.05
DOSE_IV.A        =  PDose_IV.A*BW/1000/358.21



ev_1.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Avb", replicate = FALSE)

ev_2.A <- ev (ID   = 1, amt  = DOSE_IV.A, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)

ex.iv.A <- ev_1.A+ev_2.A


## set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*40, dt)




## Get a prediction
## Out A: oral exposure to 4.2 mg/kg-d, matrix: Plasma, urine and feces
out.A <- 
  mod %>%                                               # model object
  param((exp(theta))) %>%                       # select model output
  update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.iv.A, tgrid=tsamp)               # Set up the simulation run
outdf.A<-cbind.data.frame(Time=out.A$time[3: length(out.A$time)], 
                          Plasma=out.A$Plasma[3: length(out.A$time)])



p_Jianzhong_30_IV<- ggplot(data_Jianzhong_IV_30_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.A, aes(x=Time, y=Plasma), size=0.8)+ ggtitle(" (I) FF in Sheep Plasma, 30 mg/kg IV dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                  axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                  axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                  strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  legend.position = "none")

p_Jianzhong_30_IV
cost_evaluation_sheep<-modCost(obs=data_Jianzhong_IV_30_FFC, model=outdf.A, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)





Frac             = 0.45
PDose_IM.B       = 20 
DOSE_IM.B        =  PDose_IM.B*BW/1000/358.21
Fracim.B      =  0.45
DOSEfast.B    =  DOSE_IM.B*Frac             
DOSEslow.B    =  DOSE_IM.B*(1-Frac)
TDOSE=1


ev_1.B <- ev (ID   = 1, amt  = DOSEfast.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Amtsite", replicate = FALSE)

ev_2.B <- ev (ID   = 1, amt  = DOSEslow.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSEim", replicate = FALSE)

ev_3.B <- ev (ID   = 1, amt  = DOSEfast.B+DOSEslow.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)

ex.im.B <- ev_1.B + ev_2.B + ev_3.B

tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)






out.im <- 
  mod %>%                                               # model object
  param(exp(theta)) %>%                       # select model output
  update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.im.B, tgrid=tsamp) 


outdf.im.Plasma<-cbind.data.frame(Time=out.im$time, 
                                  Plasma=out.im$Plasma)




p_El_Sheikh_IM<- ggplot(data_El_Sheikh_IV_IM_SC_20_FFC_IM, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.im.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(B) FF in Sheep Plasma, 20 mg/kg IM dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10(  limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                     axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     legend.position = "none")



p_El_Sheikh_IM
cost_Calibration_sheep= modCost(obs=data_El_Sheikh_IV_IM_SC_20_FFC_IM, model=outdf.im.Plasma, x="Time", err="Plasma_sd", cost=cost_Calibration_sheep)



p_Jianzhong_20_IM<- ggplot(data_Jianzhong_IM_20_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.im.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(L) FF in Sheep plasma, 20 mg/kg IM dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                    axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                    axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                    strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    legend.position = "none")


p_Jianzhong_20_IM
cost_evaluation_sheep<-modCost(obs=data_Jianzhong_IM_20_FFC, model=outdf.im.Plasma, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)


p_Balcomb_20_IM<- ggplot(data_Balcomb_IM_20_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=1, color="black", size=0.5)+
  xlim(c(0,24))+
  geom_line(data=outdf.im.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(N) FF in Sheep plasma, 20 mg/kg IM dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                    axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                    axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                    strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    legend.position = "none")


p_Balcomb_20_IM
cost_evaluation_sheep<-modCost(obs=data_Jianzhong_IM_20_FFC, model=outdf.im.Plasma, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)



Frac             = 0.45
PDose_IM.B       = 30 
DOSE_IM.B        =  PDose_IM.B*BW/1000/358.21
Fracim.B      =  0.45
DOSEfast.B    =  DOSE_IM.B*Frac             
DOSEslow.B    =  DOSE_IM.B*(1-Frac)
TDOSE=1


ev_1.B <- ev (ID   = 1, amt  = DOSEfast.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Amtsite", replicate = FALSE)

ev_2.B <- ev (ID   = 1, amt  = DOSEslow.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSEim", replicate = FALSE)

ev_3.B <- ev (ID   = 1, amt  = DOSEfast.B+DOSEslow.B, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)

ex.im.B <- ev_1.B + ev_2.B + ev_3.B

tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)



out.im <- 
  mod %>%                                               # model object
  param(exp(theta)) %>%                       # select model output
  update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.im.B, tgrid=tsamp) 


outdf.im.Plasma<-cbind.data.frame(Time=out.im$time, 
                                  Plasma=out.im$Plasma)


p_Jianzhong_30_IM<- ggplot(data_Jianzhong_IM_30_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.im.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(K) FF in Sheep plasma, 30 mg/kg IM dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10(  limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                     axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     legend.position = "none")


cost_evaluation_sheep<-modCost(obs=data_Jianzhong_IM_30_FFC, model=outdf.im.Plasma, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)



Frac             = 0.26
PDose_SC.C       = 20 
DOSE_SC.C        =  PDose_SC.C*BW/1000/358.21
DOSEfast.C    =  DOSE_SC.C*Frac             
DOSEslow.C    =  DOSE_SC.C*(1-Frac) 
TDOSE         =  1
tinterval     = 24
dt=0.01




ev_1.C <- ev (ID   = 1, amt  = DOSEfast.C, ii = tinterval,
              
              addl = TDOSE - 1, cmt  = "Asctsite", replicate = FALSE)

ev_2.C <- ev (ID   = 1, amt  = DOSEslow.C, ii = tinterval,
              
              addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)



ev_3.C <- ev (ID   = 1, amt  = DOSEfast.C+DOSEslow.C, ii = tinterval,
              
              addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)



ex.sc.C <- ev_1.C + ev_2.C + ev_3.C

## set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)

out.sc <- 
  mod %>%                                               # model object
  param((exp(theta))) %>%                       # select model output
  update(atol = 1E-5, maxsteps=1000) %>%              # solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.sc.C, tgrid=tsamp) 


outdf.sc.Plasma<-cbind.data.frame(Time=out.sc$time, 
                                  Plasma=out.sc$Plasma)


p_El_Sheikh_SC<- ggplot(data_El_Sheikh_IV_IM_SC_20_FFC_SC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.sc.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(C) FF in Sheep Plasma, 20 mg/kg SC dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10(  limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                     axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                     legend.position = "none")



p_El_Sheikh_SC
cost_Calibration_sheep= modCost(obs=data_El_Sheikh_IV_IM_SC_20_FFC_SC, model=outdf.sc.Plasma, x="Time", err="Plasma_sd", cost=cost_Calibration_sheep)


Frac          = 0.26
PDose_SC.C    = 40 
DOSE_SC.C     =  PDose_SC.C*BW/1000/358.21
DOSEfast.C    =  DOSE_SC.C*Frac             
DOSEslow.C    =  DOSE_SC.C*(1-Frac) 
TDOSE         =  3
tinterval     = 24
dt=0.01




ev_1.C <- ev (ID   = 1, amt  = DOSEfast.C, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Asctsite", replicate = FALSE)

ev_2.C <- ev (ID   = 1, amt  = DOSEslow.C, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)



ev_3.C <- ev (ID   = 1, amt  = DOSEfast.C+DOSEslow.C, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)



ex.sc.C <- ev_1.C + ev_2.C + ev_3.C

## set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*50, dt)

out.sc <- 
  mod %>%                                               ## Model object
  param(exp(theta)) %>%                                 ## Select model output
  update(atol = 1E-5, maxsteps=1000) %>%              ## Solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.sc.C, tgrid=tsamp) 

# Set up the simulation run
outdf.Plasma<-cbind.data.frame(Time=out.sc$time/24, 
                              Plasma=out.sc$Plasma)

outdf.Liver<-cbind.data.frame(Time=out.sc$time/24, 
                              Liver=out.sc$Liver_FFA)

outdf.Kidney<-cbind.data.frame(Time=out.sc$time/24, 
                               Kidney=out.sc$Kidney_FFA)


outdf.Muscle<-cbind.data.frame(Time=out.sc$time/24, 
                               Muscle=out.sc$Muscle_FFA)

outdf.Fat<-cbind.data.frame(Time=out.sc$time/24, 
                            Fat=out.sc$Fat_FFA)



cost_Calibration_sheep= modCost(obs=data_Wetzlich_2006_SC_3X40_FFA_Liver, model=outdf.Liver, x="Time", err="Liver_sd", cost=cost_Calibration_sheep)
cost_Calibration_sheep= modCost(obs=data_Wetzlich_2006_SC_3X40_FFA_Kidney, model=outdf.Kidney, x="Time", err="Kidney_sd", cost=cost_Calibration_sheep)
cost_Calibration_sheep= modCost(obs=data_Wetzlich_2006_SC_3X40_FFA_Muscle, model=outdf.Muscle, x="Time", err="Muscle_sd", cost=cost_Calibration_sheep)
cost_Calibration_sheep= modCost(obs=data_Wetzlich_2006_SC_3X40_FFA_Fat, model=outdf.Fat, x="Time", err="Fat_sd", cost=cost_Calibration_sheep)







p_sc_plasma_wetzlich<- ggplot(data_Wetzlich_2006_SC_3X40_FFC_Serum, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+xlim(c(0,24))+
  geom_line(data=outdf.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("40 mg/kg SC dose in Sheep Plasma, Wetzlich et al 2006") + 
  ylab("Concentration (mg/L)")+ xlab("Time (Days)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                    axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                    axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                    strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    legend.position = "none")


p_sc_plasma_wetzlich
p2<- ggplot(data_Wetzlich_2006_SC_3X40_FFA_Liver, aes(x=Time, y=Liver), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Liver-Liver_sd, ymax=Liver+Liver_sd), width=1, color="black", size=0.5)+xlim(c(0,45))+
  geom_line(data=outdf.Liver, aes(x=Time, y=Liver), size=0.8)+ ggtitle("(D) FFA in Sheep Liver, 3X40 mg/kg SC dose") + 
  ylab("Concentration (mg/kg)")+ xlab("Time (Days)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                  axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                  axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                  strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  legend.position = "none")



p2

p3<- ggplot(data_Wetzlich_2006_SC_3X40_FFA_Kidney, aes(x=Time, y=Kidney), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Kidney-Kidney_sd, ymax=Kidney+Kidney_sd), width=1, color="black", size=0.5)+xlim(c(0,45))+
  geom_line(data=outdf.Kidney, aes(x=Time, y=Kidney), size=0.8)+ ggtitle("(E) FFA in Sheep Kidney, 3X40 mg/kg SC dose") + 
  ylab("Concentration (mg/kg)")+ xlab("Time (Days)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10(  limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                   axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                   axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                   strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   legend.position = "none")


p3

p4<- ggplot(data_Wetzlich_2006_SC_3X40_FFA_Muscle, aes(x=Time, y=Muscle), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Muscle-Muscle_sd, ymax=Muscle+Muscle_sd), width=1, color="black", size=0.5)+xlim(c(0,45))+
  geom_line(data=outdf.Muscle, aes(x=Time, y=Muscle), size=0.8)+ ggtitle("(F) FFA in Sheep Muscle, 3X40 mg/kg SC dose") + 
  ylab("Concentration (mg/kg)")+ xlab("Time (Days)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                  axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                  axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                  strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                  legend.position = "none")

p4

p5<- ggplot(data_Wetzlich_2006_SC_3X40_FFA_Fat, aes(x=Time, y=Fat), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Fat-Fat_sd, ymax=Fat+Fat_sd), width=1, color="black", size=0.5)+xlim(c(0,45))+
  geom_line(data=outdf.Fat, aes(x=Time, y=Fat), size=0.8)+ ggtitle("(G) FFA in Sheep Fat, 3X40 mg/kg SC dose") + 
  ylab("Concentration (mg/kg)")+ xlab("Time (Days)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1,10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                   axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                   axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                   strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                   legend.position = "none")

p5

Frac          = 0.26
PDose_SC.E    = 40 
DOSE_SC.E     =  PDose_SC.E*BW/1000/358.21
DOSEfast.E   =  DOSE_SC.E*Frac             
DOSEslow.E    =  DOSE_SC.E*(1-Frac) 
TDOSE         =  1
tinterval     = 24
dt=0.01


ev_1.E <- ev (ID   = 1, amt  = DOSEfast.E, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Asctsite", replicate = FALSE)

ev_2.E <- ev (ID   = 1, amt  = DOSEslow.E, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)



ev_3.E <- ev (ID   = 1, amt  = DOSEfast.E+DOSEslow.E, ii = tinterval,tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)



ex.sc.E <- ev_1.E + ev_2.E + ev_3.E

## set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + tinterval*10, dt)


out.sc.E <- 
  mod %>%                                               ## Model object
  param(exp(theta)) %>%                                 ## Select model output
  update(atol = 1E-5, maxsteps=1000) %>%              ## Solver setting, atol: Absolute tolerance parameter
  mrgsim_d(data = ex.sc.E, tgrid=tsamp) 

# Set up the simulation run

outdf.sc.plasma<-cbind.data.frame(Time=out.sc.E$time, 
                                  Plasma=out.sc.E$Plasma)

p_Balcomb_SC_40<- ggplot(data_Balcomb_SC_40_FFC, aes(x=Time, y=Plasma), scale="Free", color="red") + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1) +
  geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=1, color="black", size=0.5)+
  xlim(c(0,24))+
  geom_line(data=outdf.im.Plasma, aes(x=Time, y=Plasma), size=0.8)+ ggtitle("(O) FF in Sheep plasma, 40 mg/kg SC dose") + 
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +scale_shape_prism() + 
  
  theme_prism() + scale_y_log10( limits = c(0.001, 10), expand = c(0, 0), breaks= c(0.001, 0.01,0.1,1, 10)) + theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                                    axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                                    axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                                    strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                                    legend.position = "none")


p_Balcomb_SC_40
cost_evaluation_sheep<-modCost(obs=data_Balcomb_SC_40_FFC, model=outdf.im.Plasma, x="Time", err="Plasma_sd", cost=cost_evaluation_sheep)



Fit_calibration_sheep<- lm(obs ~ mod, data = cost_Calibration_sheep$residuals)
summary(Fit_calibration_sheep)


Fit_evaluation_sheep<- lm(obs ~ mod, data = cost_evaluation_sheep$residuals)
summary(Fit_evaluation_sheep)


p6 <- ggplot (cost_Calibration_sheep$residuals, aes(obs, mod)) + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1)  + 
  geom_abline  (intercept = 0, 
                slope     = 1,
                color     ="black",size = 0.8) + ggtitle(" (H) Regression Analysis for Sheep Model Calibration")  + theme_prism() +
  ylab("Prediction (mg/L or mg/kg)")+ xlab("Observation (mg/L or mg/kg)") + 
  theme(legend.position = "none") +scale_x_log10(  limits = c(0.001, 100), expand = c(0, 0), breaks= c( 0.001,0.01, 0.1,1,10,100)) +
  scale_y_log10(limits = c(0.001, 100), expand = c(0, 0), breaks= c(0.001,0.01,0.1,1,10,100))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                                     axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                                     axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                                     strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                     strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                                     legend.position = "none")




p6

p7 <- ggplot (cost_evaluation_sheep$residuals, aes(obs, mod)) + 
  geom_point(shape = 21, colour = "red", fill = "red", size = 2, stroke = 1)  + 
  geom_abline  (intercept = 0, 
                slope     = 1,
                color     ="black",size = 0.8) + ggtitle(" (P) Regression Analysis for Sheep Model Evaluation") + 
  ylab("Prediction (mg/L or mg/kg)")+ xlab("Observation (mg/L or mg/kg)") +    theme_prism() + 
  theme(legend.position = "none") +scale_x_log10(limits = c(0.01, 100), expand = c(0, 0), breaks= c(0.01, 0.1,1,10, 100)) +
  scale_y_log10(limits = c(0.01, 100), expand = c(0, 0), breaks= c( 0.01,0.1,1,10,100))+ theme(plot.title = element_text(size=12, face="bold", family="A"), 
                                                                                               axis.text  = element_text(size = 8,  color = "black", face="bold", family="A"),
                                                                                               axis.title   = element_text(size = 10, face = "bold", color = "black", family="A"),
                                                                                               strip.text.x = element_text(size=11, face="bold", color="black", family="A"),
                                                                                               strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                                               legend.position = "none")



p7
























windows()

figure <- ggarrange( p_IV_ElSheikh, p_El_Sheikh_IM, p_El_Sheikh_SC, 
                     p2,p3,p4, p5, p6, 
                     p_Jianzhong_30_IV,p_Jianzhong_20_IV, 
                     p_Jianzhong_30_IM, p_Jianzhong_20_IM,
                     p_Palma_20_IV_FFC,
                     p_Balcomb_20_IM, p_Balcomb_SC_40,
                     p7)

figure

ggsave("Figure 2.TIFF", width = 20, height=11, dpi = 600)

summary(Fit_calibration_sheep)
summary(Fit_evaluation_sheep)

calc_mape <- function(obs, mod) {
  100 * mean(abs((obs - mod) / obs))
}

MAPE_cal <- calc_mape(cost_Calibration_sheep$residuals$obs, cost_Calibration_sheep$residuals$mod)
MAPE_eval <- calc_mape(cost_evaluation_sheep$residuals$obs, cost_evaluation_sheep$residuals$mod)

print(MAPE_cal) 
print(MAPE_eval) 

