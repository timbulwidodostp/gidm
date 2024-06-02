program gidm, properties(svyb svyj svyr swml)
	version 11

	if replay() {
		if (`"`e(cmd)'"' != "gidm") error 301
		Replay `0'	
	}
	else  Estimate `0'
end 


// Syntax:
// ParseEqns (y1 xvars1) (y2 xvars2) ...
//
// Build the -ml- equivalent equations, including the equations for the
// variance matrix.
//
// Saved results:
// s(p) number of dependent variables
// s(k_eq) number of -ml- equations : unconcentrated model
// s(yvars) the dependent variables
// s(xvars) the independent variables from all the equations
// s(xvars‘i’) the independent variables from equation ‘i’
// s(xvars‘i’_cns) "noconstant" or ""
// s(eqns) -ml- equations : unconcentrated model
// s(ceqns) -ml- equations : concentrated model
// s(eqns0) -ml- equations : constant-only concentrated model
program ParseEqns_gidp, sclass
	gettoken eq eqns : 0 , match(paren) parse("()") bind
	local p 0
	while "`eq'" != "" {
		local ++p
		local 0 `eq'
		syntax varlist(fv) [, noCONStant]
		gettoken yi xvarsi : varlist
		local xvars`p' `xvarsi'
		local xvars`p'_cns `constant'
		if `p'==1 {
			local nocon=nocon
			if nocon==1 {
				local meqns  (`yi' : `yi' = `xvarsi', `constant' noconstant) 
				local meqns0 (`yi' : `yi' = ,`constant' noconstant)		 
			}
			else if nocon==2 {
					qui: levelsof `yi', local (lev)
					qui: capture matrix drop levs
					foreach el of local lev {
						matrix levs = nullmat(levs) \ `el'
						}
					capture scalar drop K
					scalar K =`=rowsof(levs)'
					local end=K-1
					forval i=1/`end' {
						local conspit=levs[`i',1]
						if `i'==1 {
						local meqns  (cons_`conspit' : `yi' = `xvarsi', `constant' ) 
						local meqns0 (cons_`conspit' : `yi' = ,`constant' )		 
						}
						else {
						local meqns  `meqns'  (cons_`conspit' : `yi' = `xvarsi', `constant' ) 
						local meqns0 `meqns0' (cons_`conspit' : `yi' = ,`constant' )		 
						}
					}
			}
			else {
			local meqns  (`yi' : `yi' = `xvarsi', `constant' ) 			
			local meqns0 (`yi' : `yi' = ,`constant')		   			 
			}
			
		}
		if `p'>1 {
			local p_=`p'-1
			local point=A[`p_',1]
			local meqns  `meqns'   (inf_at_`point' : `yi' `xvarsi', `constant') 
			local meqns0 `meqns0'  (inf_at_`point' :   )
		}
		local yvars `yvars' `yi'
		local xvars `xvars' `xvarsi'
		local var_eq`p' `yi' `xvarsi'
		gettoken eq eqns : eqns , match(paren) parse("()") bind
	}
	

	forval i = 1/`p' {
		forval j = `i'/`p' {
		local seqns `seqns' /sigma`i'_`j'
		}
	}
	sreturn clear
	sreturn local p `p'
	sreturn local k_eq = `p'*(`p'+3)/2
	sreturn local yvars `yvars'
	sreturn local xvars `: list uniq xvars'
	sreturn local eqns `meqns' `seqns'
	sreturn local ceqns `meqns'
	sreturn local myeqns `myeqns'
	sreturn local eqns0 `meqns0'
	forval i = 1/`p' {
	sreturn local xvars`i' `xvars`i''
	sreturn local xvars`i'_cns `xvars`i'_cns'
	sreturn local var_eq`i' `var_eq`i''
	}
end



program Estimate, eclass sortpreserve
	syntax anything(id="equation(s)" equalok) [if] 					///
	[in][fweight pweight iweight][, 								///
	vce(passthru) Level(cilevel)  INFlation(numlist)					///
	noLOg noINITial	LINk(string) *                 /// 
	]	
	
	//mark the stimation sample
	marksample touse
	markout `touse' 						
	_vce_parse `touse' , opt(Robust oim opg) argopt(CLuster): `wgt' , `vce'  
	
	
	//generate the matrix A to store the inflation
	local nlist `inflation'
	capture matrix drop A 
	foreach el of local nlist {
		matrix A = nullmat(A) \ `el'
		}
	

	//nocon to check whether constant is displayed 
	capture scalar drop nocon
	if "`link'"=="lgcl" | "`link'"=="pbcl" | "`link'"=="lgcp" | "`link'"=="pbcp"  {
		scalar nocon=1
	} 
	else if "`link'"=="lgml" |"`link'"=="pbml" {
		scalar nocon=2
	}
	else {
		scalar nocon=0
	}
	
	//generate constance 
	capture drop _con
	if _rc { 
            qui gen _con=1
    }
	
	//Parse equations to genderate to be used data
	ParseEqns_gidp `anything'
	local p 	`s(p)'
	local k_eq 	`s(k_eq)'
	local yvars `s(yvars)'
	local xvars `s(xvars)'
	local xvars1 `s(xvars1)'
	local eqns  `s(eqns)'
	local eqns0 `s(eqns0)'
	local ceqns `s(ceqns)'
	local depvar: word 1 of `yvars'
	
	forvalue i=1/`p' {
	local var_eq`i' `s(var_eq`i')'
	}
	capture scalar drop keq
	scalar keq=`s(p)'

	//matrix lev to store the level of DV
	//and scalar K to store the number of levels of DV
	qui{

	//scalar K sore the level of DV
	levelsof `depvar', local (lev)
	capture matrix drop levs
	foreach el of local lev {
		matrix levs = nullmat(levs) \ `el'
		}
	capture scalar drop K
	scalar K =`=rowsof(levs)'
	}
	
	

	//validations
	//1. check whether the link function is specified
	 if (mi("`link'")) {
        display as err "Please specify the model type"
		error
    }
	
	//2.weight
	if "`weight'" != "" { 
	local wgt "[`weight'`exp']"
	}
	
	//3.check inflation function number
	local p_1=`p'-1
	if (`=rowsof(A)'!=`p_1' ) {
        display as err "The number of inflated values does not equal to the number of specified inflation functions"
		error
    }
	
	//4.check wether link function is correctly specified
	if "`link'"!="lgcl" & "`link'"!="lgcp"  & "`link'"!="pbcl"  & "`link'"!="pbcp"  ///
	   & "`link'"!="lgml" & "`link'"!="pbml" ///
	   & "`link'"!="lgpo" & "`link'"!="lgnb" & "`link'"!="pbpo" & "`link'"!="pbnb"  {
        display as err "Please specify the correct key words;see help file for a list of the supported models."
		error 
	} 

	

	//fit the Logit-Poisson
	if "`link'"=="lgpo" {
		
		//initial value
		qui: poisson `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
				
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'
		}
				
		*matrix list b0
		local contin init(b0) search(off) 
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// fit the model 
		ml model lf gidm_lgpo_lf `ceqns' `wgt' if `touse', `vce' ///
		missing	`contin'	                	         		///
		nopreserve												///		
		`log'													///
		`mlopts'												///
		maximize			
	
	}
	

	//fit the probit-Poisson
	else if "`link'"=="pbpo" {
		
		//initial value
		qui: poisson `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
				
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: probit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'
		}
				
		*matrix list b0
		local contin init(b0) search(off) 
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// fit the model 
		ml model lf gidm_pbpo_lf `ceqns' `wgt' if `touse', `vce' ///
		missing	`contin'	                 	                ///
		nopreserve												///		
		`log'													///
		`mlopts'												///
		maximize			
	
	}
	
	
	//fit the logist-NB	
	else if "`link'"=="lgnb" {
		
		// initial values
		qui: nbreg `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// fit the NBmodel
		ml model lf gidm_lgnb_lf `ceqns' (lnalpha:) `wgt' if `touse', `vce' 	///
		missing	`contin' 	waldtest(0) 							 					///
		nopreserve															///		
		`log'																///
		`mlopts'															///
		maximize			
	
	}
	
	//fit the probit-NB	
	else if "`link'"=="pbnb" {
		
		// initial values
		qui: nbreg `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// fit the NBmodel
		ml model lf gidm_pbnb_lf `ceqns' (lnalpha:)  `wgt' if `touse', `vce' 	///
		missing	`contin' waldtest(0) 								 					///
		nopreserve															///		
		`log'																///
		`mlopts'															///
		maximize			
	}	
	
	
	//fit the logit-ordered logistic
	else if "`link'"=="lgcl" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		// continue generate constants
		local kcons=K-1
		forval i=1/`kcons' {
			local neweq `neweq' (cut`i':)
		} 
		
		// initial values
		qui: ologit `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=-e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		*matrix list b0
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}

		//esimate the full model
		
		ml model d1 gidm_lgcl_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	 diff `contin' 		     		 				///
		nopreserve													///		
		`log'														///
		`mlopts'													///
		maximize													///
	
	}
	
	//fit the logit-ordered probit
	else if "`link'"=="lgcp" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		// continue generate constants
		local kcons=K-1
		forval i=1/`kcons' {
			local neweq `neweq' (cut`i':)
		} 
		
		// initial values
		qui: oprobit `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=-e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		*matrix list b0
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}

		//esimate the full model
		
		ml model d1 gidm_lgcp_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	 diff `contin' 		     		 				///
		nopreserve													///		
		`log'														///
		`mlopts'													///
		maximize													///
	
	}	
	
	//fit the probit-ordinal logit 
	else if "`link'"=="pbcl" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		// continue generate constants
		local kcons=K-1
		forval i=1/`kcons' {
			local neweq `neweq' (cut`i':)
		} 
		
		// initial values
		qui: ologit `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: probit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=-e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}

		//esimate the full model
		
		ml model d1 gidm_pbcl_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	 diff `contin' 	    		 				            ///
		nopreserve													///		
		`log'														///
		`mlopts'													///
		maximize													///
	
	}	
	
	//fit the probit-ordinal probit 
	else if "`link'"=="pbcp" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		// continue generate constants
		local kcons=K-1
		forval i=1/`kcons' {
			local neweq `neweq' (cut`i':)
		} 
		
		// initial values
		qui: oprobit `depvar' `s(xvars1)'
		capture matrix drop b0
		mat b0=e(b)
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: probit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=-e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}
		
		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}

		//esimate the full model
		
		ml model d1 gidm_pbcp_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	diff `contin' 	    		 				            ///
		nopreserve													///		
		`log'														///
		`mlopts'													///
		maximize													///
	
	}
	
	//fit the logit multninomial 
	else if "`link'"=="lgml" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		
		
		// initial values
		local base=levs[K,1]
		qui: mlogit `depvar' `s(xvars1)',baseoutcome(`base') 
		capture matrix drop b0
		mat b0=e(b)
		local nlev=K-1
		local nvar: word count `var_eq1'
		forvalue i=1/`nlev' {
			forvalue j=1/`nvar' {
			local conspit=levs[`i',1]
			local eqname `eqname' cons_`conspit'
		    }
		}
		
		matrix coleq b0 = `eqname'
	
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: logit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}

		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// continue generate constants
		ml model d1 gidm_lgml_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	diff	 `contin' 			     		 						///
		nopreserve																///		
		`log'																	///
		`mlopts'																///
		maximize			
	
	}
	
		//fit the probit mulinomial
	else if "`link'"=="pbml" {
		//gen new equations by adding nonconstant
		local neweq `ceqns'
		
		
		// initial values
		local base=levs[K,1]
		qui: mprobit `depvar' `s(xvars1)',baseoutcome(`base') 
		capture matrix drop b0
		mat b0=e(b)
		local nlev=K-1
		local nvar: word count `var_eq1'
		forvalue i=1/`nlev' {
			forvalue j=1/`nvar' {
			local conspit=levs[`i',1]
			local eqname `eqname' cons_`conspit'
		    }
		}
		
		matrix coleq b0 = `eqname'
	
		
		forvalue i=1/`=rowsof(A)' {
		local point=A[`i',1]
		capture drop inf_at_`point'
		gen inf_at_`point'=(`depvar'==`point')
		local ip=`i'+1
		qui: probit inf_at_`point' `var_eq`ip''
		capture matrix drop btemp
		matrix btemp=e(b)
		matrix coleq btemp = inf_at_`point'
		matrix b0 = b0,btemp
		capture drop inf_at_`point'

		}

		local contin init(b0) search(off)
		
		if "`initial'" !="`noinitial'" {
		local contin search(on) 
		}
		
		// continue generate constants
		ml model d1 gidm_pbml_d1 `neweq' `wgt' if `touse', `vce' 	    ///
		missing	diff	 `contin' 			     		 						///
		nopreserve																///		
		`log'																	///
		`mlopts'																///
		maximize			
	
	}
	
	
	ereturn local cmd gidm
	Replay, level(`level')

	capture drop _con

end 

program Replay
	syntax[, Level(cilevel)]
	version 9: ml di,  level(`level') 
end

