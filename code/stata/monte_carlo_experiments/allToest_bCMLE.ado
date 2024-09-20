capture program drop allToest_bCMLE
program def allToest_bCMLE, rclass

syntax [, nI( integer 100) nT(integer 4)]

  mata: allToest_bCMLEmata("y", `nI', `nT')

return scalar all_estb_CMLE=all_estb_CMLE

end
  mata:
  
  void C_lk(todo, beta_pts, lk, g, H) {
	
	external Y
	
	nT=cols(Y)
	nI=rows(Y)
	
	alfas=0 
	
	tocompus3=J(nT,1,1)
	tocompus3[1,1]=0
	S3=Y*tocompus3
	S2=Y[.,nT]
	S1=Y[.,1]
	
	nseq=2^nT
	Yseq=J(nseq,nT,0)
	for (j=1; j<=nseq; j++) {
	  for (t=1; t<=nT; t++) {
	    m=2^(nT - t+1)
		iniv=(m/2)+1
		i=iniv
		do {
		  endv=i+(m/2)-1
		  Yseq[i..endv,t]=J(endv-i+1,1,1)
		  i=i+m
		} while (i<=nseq)
	  }
	}
	S3seq=Yseq*tocompus3
	S2seq=Yseq[.,nT]
	S1seq=Yseq[.,1]

  Posicion=2:*S3seq:+S2seq:+S1seq:*(2*(nT-1)+1)
  Posicion[1]=1
  Posicion[nseq/2]=2*(nT-1)
  Posicion[(nseq/2)+1]=(2*(nT-1))+1
  Posicion[nseq]=2*(nT-1)+(2*(nT-1))+1
  probyseq=J(nseq,1,1)
  probyseq[.]=(Yseq[.,2]:*logistic(beta_pts:*Yseq[.,1]:+alfas)+(1:-Yseq[.,2]):*(1:-logistic(beta_pts:*Yseq[.,1]:+alfas)))
  for (t=3; t<=nT; t++) {
    probyseq[.]=probyseq[.]:*(Yseq[.,t]:*logistic(beta_pts:*Yseq[.,t-1]:+alfas)+(1:-Yseq[.,t]):*(1:-logistic(beta_pts*Yseq[.,t-1]:+alfas)))
  }
  ProbS=J(2*(nT-1)+(2*(nT-1)+1),1,0)
  for (i=1; i<=nseq; i++) {
    ProbS[Posicion[i]]=ProbS[Posicion[i]]+probyseq[i]
  }
    ProbY=J(nI,1,1)
	ProbY[.]=(Y[.,2]:*logistic(beta_pts:*Y[.,1]:+alfas)+(1:-Y[.,2]):*(1:-logistic(beta_pts:*Y[.,1]:+alfas)))
    for (t=3; t<=nT; t++) {
      ProbY[.]=ProbY[.]:*(Y[.,t]:*logistic(beta_pts:*Y[.,t-1]:+alfas)+(1:-Y[.,t]):*(1:-logistic(beta_pts*Y[.,t-1]:+alfas)))
    }

    ProY_S=J(nI,1,0)
    for (i=1; i<=nI; i++) {
	  Posi_i=2:*S3[i]:+S2[i]:+S1[i]:*(2*(nT-1)+1)
	  if (Posi_i==0) {
	    Posi_i=1
	  }
	  if (S3[i]==nT-1 & S1[i]==0 ) {
	    Posi_i=(2*(nT-1))
	  }
	  if (S3[i]==0 & S1[i]==1 ) {
	    Posi_i=(2*(nT-1))+1
	  }
	  if (S3[i]==nT-1 & S1[i]==1 ) {
	    Posi_i=2*(nT-1)+(2*(nT-1))+1
	  }
      ProY_S[i]=ProbY[i]/ProbS[Posi_i]
    }
	
	Suman=(S3-S2:<nT-1 :& S3-S2:>0)
	lk=log(ProY_S')*Suman
  }

void allToest_bCMLEmata(string scalar yes, real scalar nI, real scalar nT){
	
	external Y
  
   	Y=st_data(.,yes)
	Y=colshape(Y, nT)

	betaest=0
	S=optimize_init()
	optimize_init_evaluator(S,&C_lk())
	optimize_init_conv_maxiter(S, 100)
	optimize_init_evaluatortype(S,"d0")
	optimize_init_params(S, (betaest))
    optimize_init_tracelevel(S,"value")
	p=optimize(S)
	all_estb_CMLE=p[1]

	st_numscalar("all_estb_CMLE", all_estb_CMLE)
  }
  end
