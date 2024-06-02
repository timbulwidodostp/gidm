	program define gidm_pbcp_d1
		args todo  b lnf  g
        version 9.1
		
quietly {
		tempvar xb 
		mleval `xb'   = `b',eq(1)
		
		// loop to assign inflation parameters to each pai
		tempvar paiall
		gen double `paiall' = 0
		local end=2+`=rowsof(A)'-1
		local j=1
		forval i = 2/`end'  {
		tempvar xb`i' pai`j' dpaidxb`j'
		mleval `xb`i''   = `b',eq(`i')
        gen double `pai`j''     = normal(-`xb`i'')
		gen double `dpaidxb`j'' = normalden(`xb`i'')*-1
		replace `paiall'=`paiall'+`pai`j''
		local j=`j'+1
		}

		// the first keq(scalar) contains coef
		// the following K contains constants
		// loop to gen constant for each i
		local j=1
		local sta=2+`=rowsof(A)'
		local end=2+`=rowsof(A)'+K-2
		forval i = `sta'/`end'  {
			tempvar con`j'
			mleval `con`j''   = `b',eq(`i')
			local j=`j'+1
		}
		
		
		*gen qs dqdxb dqdcon
		local end=K-1
		forval i=1/`end' {
		 tempvar q`i' dqdxb`i'
		 gen double `q`i''=normal(`con`i''-`xb')
		 gen double `dqdxb`i''=normalden(`con`i''-`xb')*-1
		}	

        // LL
		tempvar clval
		gen double `clval'=0
		
		*assign non-inflate
		local end=K
		forval i = 1/`end' {
	    tempvar p`i' dpdxb`i'  
			if `i'>1 {
				local i_=`i'-1
			} 
			else {
				local i_=`i'
			}
			
		if `i'==1 {
			replace `clval'= (1-`paiall')*(normal(`con1'-`xb')) if $ML_y1 ==levs[1,1] 
			gen double `p`i''=normal(`con1'-`xb') 
			gen double `dpdxb`i'' =`dqdxb`i''
			
			local end=K-1
			forval k=1/`end' {
			tempvar dpdcon`i'`k'
			gen double `dpdcon`i'`k''=0			
			replace `dpdcon`i'`k''=normalden(`con1'-`xb')  if `k'==`i'
			}
			
		 }
			
		if `i' >1 & `i' <K {
			replace `clval'=(1-`paiall')*(normal(`con`i''-`xb')-normal(`con`i_''-`xb')) if $ML_y1 ==levs[`i',1]
			gen double `p`i''=(normal(`con`i''-`xb')-normal(`con`i_''-`xb'))
			gen double `dpdxb`i'' =`dqdxb`i''-`dqdxb`i_''
			
			local end=K-1
			forval k=1/`end' {
			tempvar dpdcon`i'`k'
			gen double `dpdcon`i'`k''=0
			replace `dpdcon`i'`k''=-normalden(`con`i_''-`xb')     if (`i'-`k')==1
			replace `dpdcon`i'`k''=normalden(`con`i''-`xb')        if `k'==`i'
			}
			
			
		 }
		
		if `i'==K {
			replace `clval'=(1-`paiall')*(1-normal(`con`i_''-`xb')) if $ML_y1 ==levs[K,1]
			gen double `p`i''=(1-normal(`con`i_''-`xb')) 
			gen double `dpdxb`i'' =-`dqdxb`i_''
			
			local end=K-1
			forval k=1/`end' {
			tempvar dpdcon`i'`k'
			gen double `dpdcon`i'`k''=0
			replace `dpdcon`i'`k''=-normalden(`con`i_''-`xb')       if (`i'-`k')==1
			replace `dpdcon`i'`k''=0                 if `k'==`i'
			}
			
			
		 }
		}
		*assign inflate
		forval j = 1/`=rowsof(A)' {
		replace `clval'=`pai`j''+ `clval' if $ML_y1 ==A[`j',1] 
		}
		
		replace `clval'=log(`clval')
		
		mlsum  `lnf'=`clval'
		if (`todo'==0 | `lnf' >= .) exit
		
		//*******************************gradiants
		
		
		*gen ds
		forval i = 1/`=rowsof(A)'  {
		
		 local end=K
		 forval j=1/`end' {
		 if (levs[`j',1]==A[`i',1]) {
		  local num=`j'
		  } 
		 }
		 tempvar d`i' v`i'
		 gen double `d`i''=(1-`paiall')/(`pai`i''+(1-`paiall')*`p`num'')
		 gen double `v`i''=1/(`pai`i''+(1-`paiall')*`p`num'')
		}	
		


		
		
		//g1
		tempvar gt1 g1
		gen double `gt1'=0
		
		local end=K
		local j=1
		forval i = 1/`end'  { 
		replace `gt1'=1/`p`i''*`dpdxb`i''   if $ML_y1 ==levs[`i',1]
		local inf=A[`j',1]
		if levs[`i',1]==`inf' {
		replace `gt1'=(`d`j'')*`dpdxb`i''   if $ML_y1 ==`inf'
		 local j=`j'+1
		 }
		}
		
		mlvecsum `lnf' `g1' = `gt1' ,eq(1)

	
		//a...
		
		local end=2+`=rowsof(A)'-1
		forval i = 2/`end'  {
		tempvar gt`i' g`i'  
		local cnt=`i'-1
		gen double `gt`i''=0
			local end=K
			local j=1 //inflation count
			forval k = 1/`end'  { 
			replace `gt`i''=-1/(1-`paiall')*`dpaidxb`cnt''   if $ML_y1 ==levs[`k',1]
			local inf=A[`j',1]
			if levs[`k',1]==`inf' {
			 replace `gt`i''=`v`j''*((`cnt'==`k')-`p`k'')*`dpaidxb`cnt''   if $ML_y1 ==`inf'
			 local j=`j'+1
			 }
			}
		mlvecsum `lnf' `g`i'' = `gt`i'' ,eq(`i')
		}
		
		
		
		
		//con...
		local sta=2+`=rowsof(A)'
		local end=2+`=rowsof(A)'+K-2
		forval i = `sta'/`end'  {
			tempvar gt`i' g`i'  
			gen double `gt`i''=0
			
			local cnt=`i'-1-`=rowsof(A)'
			local end=K
			local j=1 //inflation count
				forval k = 1/`end'  { 
				replace `gt`i''=1/`p`k''*`dpdcon`k'`cnt''   if $ML_y1 ==levs[`k',1]
				local inf=A[`j',1]
				if levs[`k',1]==`inf' {
				 replace `gt`i''=(`d`j'')*`dpdcon`k'`cnt''   if $ML_y1 ==`inf'
				 local j=`j'+1
				 }
				}
		mlvecsum `lnf' `g`i'' = `gt`i'' ,eq(`i')
		}
		
		
		****g
		local end=K+`=rowsof(A)'
		forval i=1/`end' {
			if `i'==1 {
			mat `g' = `g`i''
			} 
			else {
			mat `g' = (`g',`g`i'')
			}
		}
	}
		if (`todo'==1) exit
end	
	
