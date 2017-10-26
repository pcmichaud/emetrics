
clear all
capture log close
set more off

cd ~/emetrics/tauchen

global T = 10
global N = 10000
global rho = 0.97
global sig = 0.1

set obs $N

gen id = _n
gen e0 = sqrt($sig/(1-$rho^2))*rnormal()
forvalues t = 1/$T {
	local t1 = `t' - 1
	gen e`t' = $rho*e`t1' + sqrt($sig)*rnormal() 
}

drop e0
	
mata 
  // objective function for minimum distance
  void mdfunc(todo, b, m, W, func, g, H)
  {
  	// map parameters
  	rho  = b[1,1]
  	sige = b[1,2]
	sigv = b[1,3]
	T = (-1 + sqrt(1 + 8*rows(m)))/2; 
  	// Get "True" Covariance Matrix
  	mCovTrue = J(T,T,.)
  	for (r = 1; r<=T; ++r) {
  		for (c=1; c<=T; ++c) {
  			mCovTrue[r,c] = (r==c)*sigv + (rho^abs(r-c))*sige/(1-rho^2)
  		}
  	} 	
	d = m :- vech(mCovTrue)  	  	
  	func = d'*W*d
  }
  
  // distance 
  void gfunc(b, m, d)
  {
  	// map parameters
  	rho= b[1,1]
  	sige = b[1,2]
	sigv = b[1,3]
	T = (-1 + sqrt(1 + 8*rows(m)))/2; 
  	// Get "True" Covariance Matrix
  	mCovTrue = J(T,T,.)
  	for (r = 1; r<=T; ++r) {
  		for (c=1; c<=T; ++c) {
  			mCovTrue[r,c] = (r==c)*sigv + (rho^abs(r-c))*sige/(1-rho^2)
  		}
  	} 	
	d = m :- vech(mCovTrue)  	  	
  }
    
end

* compute covariance matrix and boostrap
capture program drop cov
program cov, eclass
	tempname b
	corr e*, cov
	matamatrix b = vech(r(C))'
	ereturn post b
end program
qui bootstrap _b, reps(100) : cov
matrix V = e(V)
matrix m = e(b)' 

mata {
	// covariance matrix from data
	m = st_matrix("m")
	V = st_matrix("V")
	
	// weight matrix (use identity for equally weighted)
	W = cholinv(V) 
	//W = I(rows(V))
	
	D = optimize_init()
	optimize_init_which(D , "min")
	optimize_init_evaluator(D, &mdfunc())
	optimize_init_evaluatortype(D, "d0")
	optimize_init_technique(D, "bfgs")
	optimize_init_argument(D, 1, m)		
	optimize_init_argument(D, 2, W)	
			
	// Starting values
	b0 = (0.97,0.02,0.0)
	optimize_init_params(D,b0)
	// Optimize
	b = optimize(D)
	Vb = optimize_result_V(D)
	
	// --- covariance matrix ---
	// compute jacobian
	n = rows(m);
	k = 3;
	G = J(n,k,.)
	eps = 1.0e-6;
	
	// get value of moments at max
	d0 = J(n,1,.);
	gfunc(b, m, d0);
	
	// now with step
	dk = J(n,1,.);
	for (k = 1; k<=3; ++k) {
		bk = b;
		bk[k] = bk[k] + eps; 
		gfunc(bk, m, dk);
		G[.,k] = (dk :- d0) :/ eps;
	}
	
	// compute std	(see Chamberlain)
	Vb = cholinv(G'*W*G)*G'*W*V*W*G*cholinv(G'*W*G);
	
	//send pars to Stata parameters
	st_matrix("bcov",b)
	st_matrix("Vcov",Vb)
	
	
}
* prepare for posting
matrix colnames bcov = rho sige sigv	
matrix colnames Vcov = rho sige sigv	
matrix rownames Vcov = rho sige sigv	


* save to file	
preserve	
matrix bcovs = bcov'
svmat bcovs
keep if bcovs1!=.
keep bcovs1
outsheet bcovs1 using params.csv, comma nonames nolabel noquote replace
ereturn post bcov Vcov	
ereturn display
restore	  
