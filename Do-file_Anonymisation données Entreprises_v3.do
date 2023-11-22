**Calcul des variables transfomées

/* Chargement des données */
global wd =  "C:\Users\" + c(username)+ "\OneDrive - GOUVCI\CAE_INS - Fichiers de Cellule d'Analyses Economiques ( CAE)\CENTRE DE CALCUL\Anonymisation_BDF"

use "$wd\Base test.dta", clear
/***********************************************/
*		Redimensionnement des variables	 		*
/***********************************************/
 
foreach var in ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288 {
	qui replace `var' = `var'/100
}

/*******************************************************************/
* Estimation de la matrice de variance-covariance des pertubations *
/*******************************************************************/
 
mat V = J(33, 33, 0)
local k = 1
foreach var in ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288 {
	 cap drop `var'_t 
	 qui sum `var'
     qui gen `var'_t = r(mean) + (`var' - r(mean)) / r(sd)
	mat V[`k',`k']=r(Var)-1
	local l=1
	foreach var1 in ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288 {
	if `l' > `k' {
		qui correlate `var' `var1', covariance
		mat V[`k',`l']=r(cov_12)-r(cov_12)/sqrt(r(Var_1)*r(Var_2))
		mat V[`l',`k']=r(cov_12)-r(cov_12)/sqrt(r(Var_1)*r(Var_2))
	}
	local l = `l'+1
 }
 local k=`k'+1
}
qui mat list V
 
/***********************************************************************************************/
*	Génération du vecteur aléatoire des bruits et génération des variables pertubées			*
/***********************************************************************************************/
 

foreach var in ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288 {
	gen t_`var'=0
}
forvalues k=1/1000  {
qui drawnorm ZZ2_a ZZ3_a XB_a RG_a RK_a XD_a XI_a Eff_YE_Total_a TA_a RA_a TB_a TC_a RC_a XC_a XF_a XH_a RS_a Sal_YE_a Code_6052_a Code_6053_a Code_6051_a Code_6042_a Code_6221_a Code_6222_a Code_6252_a Code_6253_a Code_6254_a Code_6255_a Code_6257_a Code_6281_a Code_6284_a Code_6282_a Code_6288_a, cov(V)
foreach var in ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288 {
	qui replace t_`var' = t_`var' + (`var'_t + `var'_a)/1000
}
drop *_a
}

/********************************************************************************************/
*	Estimation des divergences de KulBack-Leibler entre les données initiales et les données *
/********************************************************************************************/

/*Pairwise Correlation Difference */
pwcorr ZZ2 ZZ3 XB RG RK XD XI Eff_YE_Total TA RA TB TC RC XC XF XH RS Sal_YE Code_6052 Code_6053 Code_6051 Code_6042 Code_6221 Code_6222 Code_6252 Code_6253 Code_6254 Code_6255 Code_6257 Code_6281 Code_6284 Code_6282 Code_6288
mat define mat_real = r(C)

pwcorr t_ZZ2 t_ZZ3 t_XB t_RG t_RK t_XD t_XI t_Eff_YE_Total t_TA t_RA t_TB t_TC t_RC t_XC t_XF t_XH t_RS t_Sal_YE t_Code_6052 t_Code_6053 t_Code_6051 t_Code_6042 t_Code_6221 t_Code_6222 t_Code_6252 t_Code_6253 t_Code_6254 t_Code_6255 t_Code_6257 t_Code_6281 t_Code_6284 t_Code_6282 t_Code_6288
mat define mat_synth = r(C)

mat diff = mat_real - mat_synth

mata:
mX = st_matrix("diff")
psd = norm(mX,2)
st_matrix("sX",psd)
end

mat list sX

/*LogCluster */
Generation and evaluation of synthetic patient data (Goncalves et al. 2020)
/*Cross-Classification */
Generation and evaluation of synthetic patient data (Goncalves et al. 2020) 
 