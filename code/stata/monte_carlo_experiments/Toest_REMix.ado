capture program drop Toest_REMix
program def Toest_REMix, rclass

/*Program to estimate a Finite discrete mixture diynamic logit model with only a lag of the endogenous variable. */
/*The program has to be run when we have in memory the data.*/

syntax [, ntipos( integer 2) nI( integer 100) nT(integer 4) pry0(real 0.5) ]

  mata: Toest_REMixmata("y", `ntipos', `nI', `nT', `pry0')

return scalar estbeta_REMix=beta_est
return scalar estAME_REMix=AME_est

end
  mata:
  
  void C_lkMix(todo, alfas_pts, lk, g, H) {
	
	external Y, ntipos
	
	nT=cols(Y)
	nI=rows(Y)
	
	beta=alfas_pts[1]
	alfas=alfas_pts[2..ntipos+1]
	exppt_0=exp(alfas_pts[(ntipos+2)..(2*ntipos)])
	pt_0=exppt_0:/(1+(exppt_0*J(ntipos-1,1,1)))
	exppt_1=exp(alfas_pts[(2*ntipos+1)..(3*ntipos-1)])
	pt_1=exppt_1:/(1+(exppt_1*J(ntipos-1,1,1)))
	pt_0=pt_0, 1-pt_0*J(cols(pt_0),1,1)
	pt_0=pt_0'
	pt_1=pt_1, 1-pt_1*J(cols(pt_1),1,1)
	pt_1=pt_1'
	
	proby_t=J(nI,cols(alfas),-33)
	for (g=1; g<=ntipos; g++) {
	  proby_t[.,g]=(Y[.,2]:*logistic(beta:*Y[.,1]:+alfas[g])):+((1:-Y[.,2]):*(1:-logistic(beta:*Y[.,1]:+alfas[g])))
	}
    for (i=3; i<=cols(Y); i++) {
      for (g=1; g<=cols(alfas); g++) {
        proby_t[.,g]=proby_t[.,g]:*(Y[.,i]:*logistic(beta:*Y[.,i-1]:+alfas[.,g])+(1:-Y[.,i]):*(1:-logistic(beta*Y[.,i-1]:+alfas[.,g])))
      }
	}

  piest= (proby_t*pt_0):*(1:-Y[.,1])+(proby_t*pt_1):*Y[.,1]
  lnpiest=ln(piest)
  lk=lnpiest'*J(nI,1,1)
  }

void Toest_REMixmata(string scalar yes, real scalar ntiposh, real scalar nI, real scalar nT, real scalar pry0){
	
	external Y, ntipos
	ntipos=ntiposh
   	Y=st_data(.,yes)
	Y=colshape(Y, nT)

	betas=st_numscalar("beta")
	
	inipara=(betas, -1, 0.5, 0.6, 0.1)

	S=optimize_init()
	optimize_init_evaluator(S,&C_lkMix())
	optimize_init_conv_maxiter(S, 100)
	optimize_init_evaluatortype(S,"d0")
	optimize_init_params(S,inipara)
    optimize_init_tracelevel(S,"params")
	p=optimize(S)
	estb_CMLE=p[1]

	alfas=p[2..ntipos+1]
	exppt_0=exp(p[(ntipos+2)..(2*ntipos)])
	pt_0=exppt_0:/(1+(exppt_0*J(ntipos-1,1,1)))
	exppt_1=exp(p[(2*ntipos+1)..(3*ntipos-1)])
	pt_1=exppt_1:/(1+(exppt_1*J(ntipos-1,1,1)))
	pt_0=pt_0, 1-pt_0*J(cols(pt_0),1,1)
	pt_0=pt_0'
	pt_1=pt_1, 1-pt_1*J(cols(pt_1),1,1)
	pt_1=pt_1'

	AME_est=0
	for (g=1; g<=ntipos; g++) {
	  AME_est=AME_est+(logistic(p[1]+alfas[g])-logistic(alfas[g]))*pt_0[g]*(1-pry0)
	}
	for (g=1; g<=ntipos; g++) {
	  AME_est=AME_est+(logistic(p[1]+alfas[g])-logistic(alfas[g]))*pt_1[g]*pry0
	}
	
	st_numscalar("beta_est", estb_CMLE)
	st_numscalar("AME_est", AME_est)
  }
  end
