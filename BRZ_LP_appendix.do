*** ONLINE APPENDIX REPLICATION FILE FOR BEN ZEEV, RAMEY, ZUBAIRY, AEA P & P, MAY 2023

**** BRZ_lp_appendix.do 

****   Models asymmetry using a state-dependent model
***      Produces Figures 1-4 of online appendix

**  Version:  March 2, 2023

*** Requires:
***     brzdat.xlsx (Subset of Ramey-Zubairy JPE 2018 data, plus new potential estimate)
********************************************************************************

drop _all
clear all

set more 1

capture log close
log using appendix_results.log, replace

set scheme s1color

*******************************************************************************
** I. CHOOSE MODEL SPECIFICATION - # OF LAGS AND CONTROL VARIABLE SET
*******************************************************************************

/* BASELINE SETTINGS: p=4, horiz=20, shock newsy, negdef=0, potential rgdp_pott6, option A

   HOW TO PRODUCE ROBUSTNESS FIGURES AND TABLE 1
   
	Figure 1. Use baseline settings
	Figure 2. Use "local negdef = 1" plus other baseline settings
	Figure 3. Use "local shock newsy_resid" plus other baseline settings
	           IMPORTANT - to make a readable graph, some values must be truncated
			   Search "truncation" below to see the lines that must be uncommented
			   only for Figure 3.			    
	Figure 4. Use p=8, option D, plus other baseline settings
	Figure 5. Use BRZ_hbr_robustness.do to produce results with Hall-Barro-Redlick transformation
	
    Table 1 uses the output of the summarize command (search "Table 1" below)
	   (use settings for Fig 1, Fig 2).  For FAIR, see BDM_fair_multipliers.xlsx.)
	
*/
	 
  local p = 4 /* p = 4 is baseline lags in regressions */

  local horiz = 20  /* horizon of IRFs, 20 is baseline */

  local shock newsy  /* baseline is newsy itself; alternative is newsy_resid (innovation to newsy) */
				 
  local negdef = 0 /* baseline = 0, 0s included in negative state; if = 1, 0s included in positive state */
  
  local controls = 0 /* 0 is the baseline; 1 adds Barnichon-Debortoli-Matthes tax and trend controls */
    					  
  local potential rgdp_pott6 /* rgdp_pott6 or rgdp_potma */

  ******************************************************************************
  * SPECIFY CONTROL VARIABLES BASED ON CHOSEN OPTION
  
  if `controls' == 0 {  
  	
  	global linxlist L(1/`p').(newsy y g)  
    global nlxlist negnewsy? posnewsy? negy? posy? negg? posg? 	
	local taxes no taxes
	local trends no trends
  }
  
  else {
  	
  	global linxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy t t2 t3 t4 
    global nlxlist negnewsy? posnewsy? negy? posy? negg? posg? negtaxy? postaxy? post negt post? negt?	
	local taxes with taxes
	local trends with trends	
  }
*******************************************************************************

*******************************************************************************
** II. RAW DATA IMPORTATION FROM RAMEY-ZUBAIRY 2018 AND VARIABLE DEFINITONS
*******************************************************************************

  import excel brzdat.xlsx, sheet("brzdat") firstrow

  drop if quarter<1889

  gen qdate = q(1889q1) + _n-1
  tsset qdate, q

  * BASIC VARIABLES

  *trends
  gen t = _n
  gen t2 = t^2
  gen t3 = t^3
  gen t4 = t^4

  label var rgdp_pott6 "(real) potential GDP estimated with a 6th-degree polynomial trend"
  label var rgdp_potma "(real) potential GDP, estimated with MA on Y/hrs and 2 population series"
  label var news "Ramey-Zubairy narrative military news"
  label var ngov "nominal government purchases"
  label var ngdp "nominal GDP"
  
  label var rgdp "real GDP"
  label var pgdp "GDP deflator"
  label var nfedcurrreceipts_nipa "nominal federal current receipts, NIPA"

  gen newsy = news/(L.`potential'*L.pgdp)
  gen taxy = nfedcurrreceipts_nipa/ngdp

  * normalize variables by potential GDP and shorten names

  gen y = rgdp/`potential'
  gen g = ngov/(pgdp*`potential')

  label var newsy "Nominal military news divided by GDP deflator and potential GDP"
  label var y "Nominal GDP divided by GDP deflator and potential GDP"
  label var g "Nominal govt purchases divided by GDP deflator and potential GDP"
  label var taxy "Nominal federal tax receipts divided by nominal GDP"

  
*******************************************************************************
** III.  ESTIMATE INNOVATIONS TO NEWSY - REQUIRED ONLY WITH LOCAL NEWSY_RESID ALTERNATIVE
*******************************************************************************

* Note that in neither specification are the control variables significant
* Thus, one could justify using newsy itself rather than the residual

ivreg2 newsy $linxlist, robust bw(auto)
predict newsy_resid, resid


*******************************************************************************
** IV.  INTERACTION OF SHOCKS AND CONTROLS WITH STATE
*******************************************************************************

  * NEGATIVE SHOCK STATE EITHER INCLUDING OR EXCLUDING 0'S

  if `negdef' == 0 {
	  gen negshockstate = `shock' <=0  /* 0's included with negative news */
	  local definition 0s in negative state
  }

   else {
	  gen negshockstate =`shock'<0 /* 0's included with positive news */
	  local definition 0s in positive state
  }


  local state negshockstate

  foreach var in `shock' { 
     gen neg0`var' = `var'*`state'
     gen pos0`var' = `var'*(1-`state') 
  }
 
  forvalues i = 1/`p' { 

    foreach var in newsy y g taxy {	
	  gen neg`var'`i' = L`i'.`var'*`state'
      gen pos`var'`i' = L`i'.`var'*(1-`state') 
   }
}

  foreach var in t t2 t3 t4 {
    gen neg`var' = `var'*`state'
    gen pos`var' = `var'*(1-`state') 
}

*******************************************************************************
** V. INITIALIZATIONS FOR CREATING IRFS, MULTIPLIERS, CUMULATIVE ENDOGENOUS
**      VARIABLES - FOR 1-STEP METHOD
*******************************************************************************

  gen h = _n - 1 /* h indexes the horizon h */
  
  * Initialize some sums and parameters to 0 or .

  foreach var in liny ling posy posg negy negg {
  	gen sum`var' = 0
  }


  foreach var in bylin bypos byneg bglin bgpos bgneg up95bylin up95bypos up95byneg  ///
    up95bglin up95bgpos up95bgneg  lo95bylin lo95bypos lo95byneg lo95bglin lo95bgpos ///
    lo95bgneg seylin seypos seyneg seglin segpos segneg  multlin multpos multneg { 
       quietly gen `var' = .
  } 

  label var bylin "GDP IRF, linear model"
  label var bypos "GDP IRF to positive news shock"
  label var byneg "GDP IRF to negative news shock"
  label var bglin "Govt IRF, linear model"
  label var bgpos "Govt IRF to positive news shock"
  label var bgneg "Govt IRF to negative news shock"
  label var seylin "standard error of GDP IRF estimate"
  label var seglin "standard error of govt IRF estimate"

  * Create cumulative variables used for 1-step multiplier estimation

  gen cumuly = 0
  gen cumulg = 0
 
  forvalues i = 0/`horiz' {

     gen f`i'cumuly = F`i'.y + cumuly
     gen f`i'cumulg = F`i'.g + cumulg
   
     gen negf`i'cumulg = f`i'cumulg*`state'
     gen posf`i'cumulg = f`i'cumulg*(1-`state')
   
     replace cumuly = f`i'cumuly
     replace cumulg = f`i'cumulg
   
  }

*******************************************************************************
** VI. ESTIMATION OF IRFS
*******************************************************************************

  forvalues i = 0/`horiz' {

   ivreg2 F`i'.y `shock' $linxlist, robust bw(auto)

       gen bylinh`i' = _b[`shock']  
       gen seylinh`i' = _se[`shock']

   ivreg2 F`i'.g `shock' $linxlist, robust bw(auto)

       gen bglinh`i' = _b[`shock']  
       gen seglinh`i' = _se[`shock'] 

   ivreg2 F`i'.y pos0`shock' neg0`shock' `state' $nlxlist, robust bw(auto) 

       gen byposh`i' = _b[pos0`shock']
       gen bynegh`i' = _b[neg0`shock']
  
       gen seyposh`i' = _se[pos0`shock']
       gen seynegh`i' = _se[neg0`shock'] 

   ivreg2 F`i'.g pos0`shock' neg0`shock' `state' $nlxlist, robust bw(auto) 

       gen bgposh`i' = _b[pos0`shock']
       gen bgnegh`i' = _b[neg0`shock']
  
       gen segposh`i' = _se[pos0`shock']
       gen segnegh`i' = _se[neg0`shock']   
 
   * Create the 3-step multipliers from the IRFs
   
   replace sumliny = bylinh`i' + sumliny
   replace sumling = bglinh`i' + sumling
  
   replace sumposy = byposh`i' + sumposy
   replace sumposg = bgposh`i' + sumposg
  
   replace sumnegy = bynegh`i' + sumnegy
   replace sumnegg = bgnegh`i' + sumnegg
  
   gen multlinh`i' = sumliny/sumling
   gen multposh`i' = sumposy/sumposg
   gen multnegh`i' = sumnegy/sumnegg
  
  * Create series using observation h = horizon
  
    foreach var in bylin bypos byneg bglin bgpos bgneg multlin multpos multneg ///
      seylin seglin seypos seyneg segpos segneg { 
       quietly replace `var' = `var'h`i' if h==`i'	
	   drop `var'h`i'
   }
  
   foreach var in ylin glin ypos gpos yneg gneg { 
      quietly replace up95b`var' = b`var' + 1.96*se`var'
	  quietly replace lo95b`var' = b`var' - 1.96*se`var'	
   }

  }

  * rename multipliers to indicate they are 3-step estimates
  rename multlin multlin3
  rename multpos multpos3
  rename multneg multneg3

  label var multlin3 "cumulative multiplier, 3-step method, linear model"
  label var multpos3 "cumulative multiplier, 3-step method, positive shock"
  label var multneg3 "cumulative multiplier, 3-step method, negative shock"


*******************************************************************************
** VII. PRINTOUTS AND GRAPHS OF IMPULSE RESPONSE FUNCTIONS (IRF) AND F-STATISTICS
*******************************************************************************

  display as text "IRFS AND STANDARD ERRORS AND 3-STEP METHOD MULTIPLIERS"
  list h bgpos segpos bgneg segneg bypos seypos byneg seyneg multlin3 multpos3 ///
       multneg3 if h<=`horiz'
  
  label var bglin "Gov, linear model"
  label var bylin "GDP, linear model"
  label var bgpos "Gov, positive"
  label var bypos "GDP, positive"
  label var bgneg "Gov, negative"
  label var byneg "GDP, negative"
  label var up95bglin "95% CI"
  label var up95bylin "95% CI"
  label var lo95bglin "95% CI"
  label var lo95bylin "95% CI"

  label var h "horizon"

  tw (rarea up95bglin lo95bglin h, bcolor(gs12) clw(medthin medthin)) ///
   (scatter bglin h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick)) ///
   if h<=`horiz', legend(on order(1 4)) title("Government Spending") name(irfglin) 

  tw (rarea up95bylin lo95bylin h, bcolor(gs12) clw(medthin medthin)) ///
   (scatter bylin h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick)) ///
    if h<=`horiz', legend(on order(1 4)) title("GDP") name(irfylin)  

  tw (rarea up95bgneg lo95bgneg h, bcolor(red%10) clw(medthin medthin)) ///
    (scatter up95bgpos bgpos lo95bgpos bgneg h, clw(medthin medthick medthin medthick) ///
    c(l l l l l) clp(- l - l) clc(blue blue blue red) ms(i o i i i i) ///
    mc(blue blue blue red)) if h<=`horiz', title("Government Spending") name(irfgnl) ///
    legend(on order(3 5))

  tw (rarea up95byneg lo95byneg h, bcolor(red%10) clw(medthin medthin)) ///
    (scatter up95bypos bypos lo95bypos byneg h,  clw(medthin medthick medthin medthick) ///
    c(l l l l l) clp(- l - l) clc(blue blue blue red) ms(i o i i i i) ///
    mc(blue blue blue red)) if h<=`horiz', title("GDP") name(irfynl) legend(on order(3 5))
  
  graph combine irfglin irfylin irfgnl irfynl, col(2) iscale(0.5) name(irfcombo) ///
    title("IRFs, `shock' shock, `p' lags, `definition', `taxes', `trends'", size(medsmall))


*******************************************************************************
** VIII. ESTIMATION OF 1-STEP MULTIPLIERS USING LP-IV
*******************************************************************************
  
  * preliminaries - initializing variables to missing

  foreach var in multlin1 multpos1 multneg1 multposj1 multnegj1 Fkplin Fkppos Fkpneg ///
    semlin1 sempos1 semneg1 semposj1 semnegj1 ptestdiff Fdifflin Fdiffpos Fdiffneg ///
    up95multlin lo95multlin up95multneg lo95multneg up95multpos lo95multpos ///
    up95multnegj lo95multnegj up95multposj lo95multposj { 
  
        quietly gen `var' = .
  
   } 

   forvalues i = 0/`horiz' { 

     ivreg2 f`i'cumuly (f`i'cumulg = `shock') $linxlist, robust bw(auto) 
        gen Fkplinh`i'= e(widstat) /* Kleibergen-Paap rk Wald F statistic*/
        gen Fdifflinh`i'= Fkplinh`i'- 23.1085 
        gen multlinh`i' = _b[f`i'cumulg]
        gen semlinh`i' = _se[f`i'cumulg] /* HAC robust standard error*/

     ivreg2 f`i'cumuly (posf`i'cumulg = pos0`shock') `state' $nlxlist, robust bw(auto) 
        gen Fkpposh`i'= e(widstat)
        gen Fdiffposh`i'= Fkpposh`i'- 23.1085 
        gen multposh`i' = _b[posf`i'cumulg]
        gen semposh`i' = _se[posf`i'cumulg]
 
     ivreg2 f`i'cumuly (negf`i'cumulg = neg0`shock') `state' $nlxlist, robust bw(auto)
        gen Fkpnegh`i'= e(widstat)  
        gen Fdiffnegh`i'= Fkpnegh`i'- 23.1085 
        gen multnegh`i' = _b[negf`i'cumulg]
        gen semnegh`i' = _se[negf`i'cumulg]
  
  * Estimate both multipliers jointly, "j" suffix indicates "joint"
  
    ivreg2 f`i'cumuly (posf`i'cumulg negf`i'cumulg = pos0`shock' neg0`shock') ///
     `state' $nlxlist, robust bw(auto) 
        gen multposjh`i' = _b[posf`i'cumulg]
        gen semposjh`i' = _se[posf`i'cumulg]
        gen multnegjh`i' = _b[negf`i'cumulg]
        gen semnegjh`i' = _se[negf`i'cumulg]
  
  * Test for equality of positive and negative multipliers
  
  test posf`i'cumulg=negf`i'cumulg	
  gen ptestdiffh`i' = r(p)
  
  * create multiplier, std. error, CI, and F-stat series with h indexing horizon
  	
 foreach var in multlin multpos multneg multposj multnegj ///
        semlin sempos semneg semposj semnegj {  
    quietly replace `var'1 = `var'h`i' if h==`i'
	drop `var'h`i'
  }
  
  foreach var in  ptestdiff Fkplin Fkppos Fkpneg  { 
    quietly replace `var' = `var'h`i' if h==`i'
	drop `var'h`i'
  }
  
      foreach var in lin pos neg posj negj {
    quietly replace up95mult`var' = mult`var'1 + 1.96*sem`var'1 
	quietly replace lo95mult`var' = mult`var'1 - 1.96*sem`var'1 	
  }
 
  
 foreach var in  Fdifflin Fdiffpos Fdiffneg { 
    quietly replace `var' = `var'h`i' if h==`i'
	quietly replace `var' = 30 if `var'>30
	drop `var'h`i'
  }
  
  
 }
 

*******************************************************************************
** IX. GRAPHS OF ONE-STEP MULTIPLIERS AND AVERAGE STANDARD ERRORS OF THE ESTIMATES 
*******************************************************************************

  label var multposj1 "multiplier, positive shock, 1-step joint estimation"
  label var multnegj1 "multiplier, negative shock, 1-step joint estimation"
  label var semposj1 "standard error of mult estimate, pos shock, 1-step joint estimation" 
  label var semnegj1 "standard error of mult estimate, neg shock, 1-step joint estimation" 
  label var ptestdiff "p-value on H0: pos multiplier = neg multiplier"

  label var multlin1 "linear"
  label var multpos1 "Positive shock"
  label var multneg1 "Negative shock"

  label var up95multlin "95% CI"
  label var lo95multlin "95% CI"
  label var up95multpos "95% CI"
  label var lo95multpos "95% CI"
  label var up95multneg "95% CI"
  label var lo95multneg "95% CI"
  
  label var Fkppos "First-stage F-stat for positive shocks"
  label var Fkpneg "First-stage F-stat for negative shocks"
  label var Fdiffpos "F-stat - threshold for weak instrument, positive shocks"
  label var Fdiffneg "F-stat - threshold for weak instrument, negative shocks"
  
  label var h "horizon"

  display as text "MULTIPLIERS FROM 3-STEP AND 1-STEP PLUS P-VALUE FOR TEST"

  list h multneg3 multnegj1 semnegj1 multpos3 multposj1 semposj1  ptestdiff if h<21, abbrev(10)
  
  display as text "FOR TABLE 1: MEANS OF MULTIPLIER ESTIMATE STANDARD ERRORS"
  
  summ semnegj1 semposj1 if h<21
 
  display as text "FIRST STAGE F-STATISTICS "

  list h Fkplin Fkppos Fkpneg Fdifflin Fdiffpos Fdiffneg if h<=`horiz'

  
********************************************************************************
  * NECESSARY TRUNCATION FOR A READABLE FIGURE 3
  
  *  Anticipations effects lead the first few confidence bands to be huge, so the graph
  *     is not readable without truncation
  *  Uncomment the following lines ONLY FOR FIGURE 3 (which uses newsy_resid)
 /* 
  foreach var in multpos1 multneg1  {
    replace `var' = . if h==0 & `var'>3
  }
  
  foreach var in lo95multneg up95multneg lo95multpos up95multpos {
  	 replace `var' = . if `var'>4.5 
	 replace `var' = . if `var'<-2 
  }
  
  * End of special truncation lines for Figure 3 - 
  */
  
********************************************************************************
  
  label var multpos1 "positive shock"
  label var multneg1 "negative shock"
  label var lo95multneg "95% CI"
  label var up95multneg "95% CI"


  tw (rarea up95multlin lo95multlin h, bcolor(gs12) clw(medthin medthin)) ///
    (scatter multlin1 h, clw(medthick) ///
    c(l l l l l) clp(l) clc(black ) ms(i) mc(black)) if h>0 & h<=`horiz', ///
    title("LP-IV 1-Step Multipliers, Linear") name(multlin_1step) 
  
  tw (rarea up95multneg lo95multneg h, bcolor(red%10) clw(medthin medthin)) ///
    (scatter up95multpos multpos1 lo95multpos multneg1 h, clw(medthin medthick ///
	medthin medthick)   c(l l l l l) clp(- l - l) clc(blue blue blue red ) ///
	ms(i o i i i i) mc(blue blue blue red)) if h<=`horiz', ///
	title("1-Step Multipliers,`shock' shock, `p' lags, `definition', `taxes', `trends'", size(medsmall)) ///
	name(multnl_1step)
  
  tw (rarea up95multneg lo95multneg h, bcolor(red%10) clw(medthin medthin)) ///
    (scatter up95multpos multpos1 lo95multpos multneg1 h, clw(medthin medthick ///
	medthin medthick)   c(l l l l l) clp(- l - l) clc(blue blue blue red ) ///
	ms(i o i i i i) mc(blue blue blue red)) if h<=`horiz', ///
	title("Cumulative Multipliers") legend(on order(3 5) cols(1) bmargin(b=22 r=20)) name(multnl)

label var Fdiffpos "positive shock"
label var Fdiffneg "negative shock"

tw scatter Fdiffpos Fdiffneg h if h<=`horiz', c(l l) clp(- l) clc(blue red) ///
   clw(medthick medthick) ms( o i) mc(blue red) title("F-Stats Relative to Threshold") ///
   legend(cols(2)) ylabel(-30 -20 -10 0 10 20 30) yline(0, lp(_) lc(black)) name(fstat)

grc1leg irfgnl irfynl multnl fstat, cols(2) ysize(4) xsize(8) iscale(0.7) legendfrom(fstat) ///
    name(combo4)

*******************************************************************************

capture log close

 
