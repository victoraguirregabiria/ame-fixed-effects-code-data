capture program drop Toest_AMEfeT4
program def Toest_AMEfeT4, rclass

syntax [, betae(real 1.0) nI( integer 100) nT(integer 4)]

  mata: Toest_AMEfemata("y",`betae', `nI', `nT')

return scalar estAME_fe=EestAME_fe

end
  mata:
  
  void Toest_AMEfemata(string scalar yes, real scalar betae, real scalar nI, real scalar nT){
//Valild only for T=4 (counting the initial observation as t=1).   

    Y=st_data(.,yes)
	Y=colshape(Y, nT)

	tocompus3=J(nT,1,1)
	tocompus3[1,1]=0
	S3=Y*tocompus3
	S2=Y[.,nT]
	S1=Y[.,1]
	ESE=S1,S2,S3

	tmp1=(ESE:==(0,0,1))*J(3,1,1)
	tmp2=(tmp1:==3)
	S001=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(0,0,2))*J(3,1,1)
	tmp2=(tmp1:==3)
	S002=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(0,1,2))*J(3,1,1)
	tmp2=(tmp1:==3)
	S012=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(0,1,3))*J(3,1,1)
	tmp2=(tmp1:==3)
	S013=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(1,0,1))*J(3,1,1)
	tmp2=(tmp1:==3)
	S101=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(1,0,2))*J(3,1,1)
	tmp2=(tmp1:==3)
	S102=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(1,1,2))*J(3,1,1)
	tmp2=(tmp1:==3)
	S112=tmp2'*J(nI,1,1)

	tmp1=(ESE:==(1,1,3))*J(3,1,1)
	tmp2=(tmp1:==3)
	S113=tmp2'*J(nI,1,1)

	estAME_fe=(1/2)*(exp(betae)-1)*((S001+S112)/nI) + ((exp(betae)-1)/(exp(betae)+1)) *((S012+S101)/nI)

	st_numscalar("EestAME_fe", estAME_fe)
  }
  end
