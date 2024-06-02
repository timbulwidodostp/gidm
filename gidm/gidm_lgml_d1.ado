program define gidm_lgml_d1
        version 9.1
		args todo  b lnf g
qui {

		//gen xb
		tempvar den 
		gen double `den' = 0
		local end=K-1
		forval i = 1/`end'  {
			tempvar xb`i' 
			mleval `xb`i'' = `b',eq(`i')
			replace `den'=`den'+exp(`xb`i'')
		}
		replace `den'=`den'+1

		
		
		// loop to assign inflation parameters to each pai
		tempvar paiall 
		gen double `paiall' = 0
		local j=1
		local sta=K
		local end=K+`=rowsof(A)'-1
		forval i = `sta'/`end'  {
			tempvar pai`j' xb`i'
			mleval `xb`i'' = `b',eq(`i')
			gen double `pai`j'' = invlogit(`xb`i'') 
			replace `paiall'=`paiall'+`pai`j''
			local j=`j'+1
		}
		
		
		// loop to assign p
		local end=K-1
		forval i = 1/`end'  {
		tempvar p`i'
		gen double `p`i'' = exp(`xb`i'')/`den'
		}
		local las=K
		tempvar p`las'
		gen double `p`las'' = 1/`den'

		
		//calculate clval
		tempvar clval 
		gen double `clval'=0
		
		
		*assign mlpart
		local end=K
		forval i = 1/`end' {
		replace `clval'=(1-`paiall')*(`p`i'') ///
				if $ML_y1 ==levs[`i',1] 
		} 
		
		*assign inflation
		forval j = 1/`=rowsof(A)' {
		replace `clval'=`pai`j''+ `clval' if $ML_y1 ==A[`j',1] 
		}
		
		//calculate clval
		mlsum  `lnf'=log(`clval')
		
		
		if (`todo'==0 | `lnf' >= .) exit	
		
		*************************************gradiant
		*gen ds
		forval i = 1/`=rowsof(A)'  {
		
		 local end=K
		 forval j=1/`end' {
		 if (levs[`j',1]==A[`i',1]) {
		  local num=`j'
		  } 
		 }
		 tempvar d`i'
		 gen double `d`i''=(1-`paiall')/(`pai`i''+(1-`paiall')*`p`num''  )
		}	
		
		
		****grad for mlpart 
		local end=K-1
		forval i = 1/`end'  { 
		 tempvar gt`i' g`i'
		 
		 *assign normal
		 gen double `gt`i''=0
		 local end=K
		 forval j=1/`end' {
		 replace `gt`i''=(`i'==`j') -`p`i'' if $ML_y1 ==levs[`j',1] 
		    }
		 
		 *assign inflate
		 forval j = 1/`=rowsof(A)' {
		 
			 local end=K
			 forval x=1/`end' {
			 if (levs[`x',1]==A[`j',1]) {
			  local num=`x'
			  } 
			 }
		 replace `gt`i''=`gt`i''*`p`num''*`d`j''    if $ML_y1 ==A[`j',1] 
		 }
		
		mlvecsum `lnf' `g`i'' = `gt`i'' ,eq(`i')
	    }
	
		****grad for twopai
 
		local w=1
		local sta=K
		local end=K+`=rowsof(A)'-1
		forval i = `sta'/`end'  { 
		 tempvar gt`i' g`i'
		 
		 *assign normal 
		 gen double `gt`i''=0
		 local end=K
		 forval k=1/`end' {
		 replace `gt`i''=`pai`w''*(1-`pai`w'')/(1-`paiall')*-1    if $ML_y1 ==levs[`k',1] 
		 }
		 
		 *assign inflate
		 forval j = 1/`=rowsof(A)' {
		 
			 local end=K
			 forval x=1/`end' {
			 if (levs[`x',1]==A[`j',1]) {
			  local num=`x'
			  } 
			 }
		 replace `gt`i''=`gt`i''*`d`j''*(`p`num''-(`w'==`j'))    if $ML_y1 ==A[`j',1] 
		 }
		 
		 local w=`w'+1
		 mlvecsum `lnf' `g`i'' = `gt`i'' ,eq(`i')
	     }

		****grad for pai 
		local end=K+`=rowsof(A)'-1
		forval i=1/`end' {
		if `i'==1 {
		mat `g' = `g`i''
		} 
		else {
		mat `g' = (`g',`g`i'')
		}
		}
		
		if (`todo'==1) exit
}
	
	end
