	program define gidm_lgpo_lf
        version 9.1
		args lnf mu p1 p2 p3 p4 p5 p6 p7 p8 p9 p10
		tempvar emu paiall
		quietly {
		gen double `emu'  	= exp(`mu')
		// loop to assign inflation parameters to each pai
		gen double `paiall' = 0
		local allvar=0
		forval i = 1/`=rowsof(A)'  {
			tempvar pai`i'
			gen double `pai`i'' = invlogit(`p`i'') 
			replace `paiall'=`paiall'+`pai`i''
			local allvar `allvar',A[`i',1]
		}
		
		//sum of pai no more than 1
		*replace `paiall'=min(`paiall',0.999999)
		
		//the poisson value 
		*local poi_val=$ML_y1 *ln(`emu')-`emu'-lngamma($ML_y1 +1)
		
		//calculate LL function
		forval i = 1/`=rowsof(A)' {
			local inflation=A[`i',1]
			replace `lnf'=ln(`pai`i''+(1-`paiall')*exp($ML_y1 *ln(`emu')-`emu'-lngamma($ML_y1 +1))) ///
			if $ML_y1 ==`inflation' & `paiall'<1
		}
		replace `lnf'=ln(1-`paiall')+$ML_y1 *ln(`emu')-`emu'-lngamma($ML_y1 +1)		///
		if !inlist($ML_y1 ,`allvar') & `paiall'<1
}
	end

	
