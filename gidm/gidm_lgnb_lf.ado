	program define gidm_lgnb_lf
        version 9.1
		args lnf mu p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 
		tempvar emu paiall alpha m p
		quietly{
		gen double `emu'  	= exp(`mu')
		// loop to assign inflation parameters to each pai
		gen double `paiall' = 0
		local allvar=0
		forval i = 1/`=rowsof(A)' {
			tempvar pai`i'
			gen double `pai`i'' = invlogit(`p`i'') 
			replace `paiall'=`paiall'+`pai`i''
			local allvar `allvar',A[`i',1]
		}
		
		local lastp=`=rowsof(A)'+1
		gen double `alpha' = `p`lastp''

		//sum of pai no more than 1
		*replace `paiall'=min(`paiall',0.999999)
		
		//the nagative binomial value 
		if `alpha' < -18 { scalar `alpha' = -18 
		}
		if `alpha' >  18 { scalar `alpha' =  18 
		}
		replace `alpha' = exp(`alpha')
		gen double `m' = 1/`alpha'
		gen double `p' = 1/(1+`alpha'*`emu') 
		*local nb_val=lngamma(`m'+$ML_y1 ) -lngamma($ML_y1 +1) - lngamma(`m') + `m'*ln(`p') +$ML_y1 *ln(1-`p') )
		
		//calculate LL function
		forval i = 1/`=rowsof(A)'  {
			local inflation=A[`i',1]
			replace `lnf'=ln(`pai`i''+(1-`paiall' )*exp(lngamma(`m'+$ML_y1 ) -lngamma($ML_y1 +1) - lngamma(`m') + `m'*ln(`p') +$ML_y1 *ln(1-`p') )) ///
			if $ML_y1 ==`inflation' & `paiall'<1
		}
		replace `lnf'=ln(1-`paiall' )+lngamma(`m'+$ML_y1 ) -lngamma($ML_y1 +1) - lngamma(`m') + `m'*ln(`p') +$ML_y1 *ln(1-`p') 		///
			if !inlist($ML_y1 ,`allvar' ) & `paiall'<1
	}

	end

		
	
