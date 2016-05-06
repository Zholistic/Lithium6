#include "usrdef.h"
#include "gauss.h"
#include "mex.h"

char WORKSHOP[80] = "~/Documents/active/phd/2d_eos/virial/";

// ==> to control the gaussian integrals;
const int NCUT = 4;
const int _MAX_DIM_ = NCUT * _MAX_GAUSS_;

double calc_db2(double be)
{
	double x = 0.5*log(2.0*PI * be);

	int    ct, it, jt;
	double t, tinc;
	double inte, dblp, dblm, tant;

	tinc = PI_2 / NCUT;
	inte = 0.0;

	for(ct = 0 ; ct < _MAX_DIM_ ; ct++){
		it = ct / _MAX_GAUSS_;
		jt = ct % _MAX_GAUSS_;

		t    = (it + (gm_dGaussPoint[jt] + 1.0) * P5) * tinc;
		tant = tan(t);

		if(tant < 100.0){
			dblp  = exp(+ 2.0*tant) / (2.0*PI);
			dblp  = exp(- dblp);
			dblp *= 2.0;
			dblp /= PI*PI + 4.0 * (tant - x) * (tant - x);
		}else{
			dblp  = 0.0;
		}

		dblm  = exp(- 2.0*tant) / (2.0*PI);
		dblm  = exp(- dblp);
		dblm *= 2.0;
		dblm /= PI*PI + 4.0 * (tant + x) * (tant + x);

		dblp /= cos(t)*cos(t);
		dblm /= cos(t)*cos(t);

		inte += gm_dGaussWeight[jt] * (dblp + dblm);
	}

	inte *= (P5*tinc);

	return exp(be) - inte;
}

double calc_db3(double be)
{
	double res;
	double x = 0.5*log(2.0*PI * be);

	// ==> note that the expression is valid for - 2.0 < x < 1.0, or 0.002915 < be < 1.176;
	res  = - 0.45938;
	res += - 0.40400 * pow(x, 1.0);
	res += - 0.31103 * pow(x, 2.0);
	res += - 0.16998 * pow(x, 3.0);
	res += - 0.17801 * pow(x, 4.0);
	res += - 0.23461 * pow(x, 5.0);
	res += - 0.13623 * pow(x, 6.0);
	res += - 0.02685 * pow(x, 7.0);

	return res;
}

double calc_xb2(double be)
{
	double inc = 0.002;
	double xb2;

	xb2  = calc_db2(be + inc) - calc_db2(be - inc);
	xb2 /= 2.0*inc;
	xb2 *= (- be);

	return xb2;
}

double calc_xb3(double be)
{
	double inc = 0.002;
	double xb3;

	xb3  = calc_db3(be + inc) - calc_db3(be - inc);
	xb3 /= 2.0*inc;
	xb3 *= (- be);

	return xb3;
}

double integral_pressure1(double z)
{
	int    ct, it, jt;
	double t, tinc;
	double inte, dbl, tant;

	tinc = PI_2 / NCUT;
	inte = 0.0;

	for(ct = 0 ; ct < _MAX_DIM_ ; ct++){
		it = ct / _MAX_GAUSS_;
		jt = ct % _MAX_GAUSS_;

		t    = (it + (gm_dGaussPoint[jt] + 1.0) * P5) * tinc;
		tant = tan(t);

		dbl  = log(1.0 + z*exp(- tant));
		dbl /= cos(t)*cos(t);

		inte+= gm_dGaussWeight[jt] * dbl;
	}

	inte *= (P5*tinc);

	return inte;
}

double integral_number1(double z)
{
	int    ct, it, jt;
	double t, tinc;
	double inte, dbl, tant;

	tinc = PI_2 / NCUT;
	inte = 0.0;

	for(ct = 0 ; ct < _MAX_DIM_ ; ct++){
		it = ct / _MAX_GAUSS_;
		jt = ct % _MAX_GAUSS_;

		t    = (it + (gm_dGaussPoint[jt] + 1.0) * P5) * tinc;
		tant = tan(t);

		dbl  = z*exp(- tant) / (1.0 + z*exp(- tant));
		dbl /= cos(t)*cos(t);

		inte+= gm_dGaussWeight[jt] * dbl;
	}

	inte *= (P5*tinc);

	return inte;
}

double integral_kappa1(double z)
{
	int    ct, it, jt;
	double t, tinc;
	double inte, dbl, tant;

	tinc = PI_2 / NCUT;
	inte = 0.0;

	for(ct = 0 ; ct < _MAX_DIM_ ; ct++){
		it = ct / _MAX_GAUSS_;
		jt = ct % _MAX_GAUSS_;

		t    = (it + (gm_dGaussPoint[jt] + 1.0) * P5) * tinc;
		tant = tan(t);

		dbl  = z*exp(- tant) / pow(1.0 + z*exp(- tant), 2.0);
		dbl /= cos(t)*cos(t);

		inte+= gm_dGaussWeight[jt] * dbl;
	}

	inte *= (P5*tinc);

	return inte;
}

//int main(int argc, char** argv)

                  
void huicoeffs(double y[], double x1[], double x2[], double x3[])        
{
	double be = x1[0]; //atof converts string to double
    int coeff = (int)x2[0]; //atoi converts string to integer
    int all_coeff = (int)x3[0]; //atoi converts string to integer
    double db2 = 0;
    double db3 = 0;
    double xb2 = 0;
    double xb3 = 0;
    double beb = 0; 
    if (all_coeff > 0)
    {
    	char szDB2[128];
	    char szDB3[128];
        strcpy(szDB2, "db2_out");
        strcat(szDB2, ".dat");
 
        strcpy(szDB3, "db3_out");
        strcat(szDB3, ".dat");
	    
        ofstream coeff_out2(szDB2);
	    ofstream coeff_out3(szDB3);
	    assert(coeff_out2.is_open());
        assert(coeff_out3.is_open());
        for(beb = 0.003; beb < 3.0; beb +=0.001)
        {
            coeff_out2 << setprecision(16) << setw(24);
            coeff_out2 << setiosflags(ios::scientific);
            coeff_out2 << calc_db2(beb)<< "    ";

            coeff_out3 << setprecision(16) << setw(24);
            coeff_out3 << setiosflags(ios::scientific);
            coeff_out3 << calc_db3(beb) << "    ";
        }     
    }

    if (be > 0)
    {
        db2 = calc_db2(be);
        db3 = calc_db3(be);

        if (coeff == 1)
        {
            cout << db2 << endl;
            y[0] = (double)db2;
        }
        else if (coeff == 2)
        {
            cout << db3 << endl; 
            y[0] = (double)db3;
        }
        xb2 = calc_xb2(be);
        xb3 = calc_xb3(be);
    } else
    {
     y[0] = -1;   
        
    }
  
//---------------------------------------------------------------
	char szEOS[128];
    char c[10];
    sprintf(c , "%lf" , x2[0]);
*
	//strcpy(szEOS, WORKSHOP);
                

    strcpy(szEOS, "pk2d_be_");
    strcat(szEOS, c);
    strcat(szEOS, ".dat");

	ofstream eos2d(szEOS);
	assert(eos2d.is_open());
//---------------------------------------------------------------

	double z, teff;
	double num, pressure, entropy, energy, contact, kappa;
	double bmu;

	for(bmu = - 3.0 ; bmu <= 3.0 + 1.0e-012 ; bmu += 0.01){
		z = pow(10.0, bmu);

		// ==> determine the particle density n;
		num = integral_number1(z) + 2.0 * db2 * pow(z, 2.0) + 3.0 * db3 * pow(z, 3.0);
		teff = 1.0/num;

		// ==> pressure, P/P_0, where P_0 is the pressure at T=0;
		pressure  = integral_pressure1(z) + db2 * pow(z, 2.0) + db3 * pow(z, 3.0);
		pressure *= 2.0*teff*teff;

		// ==> entropy, S/(nk_B);
		entropy  = integral_pressure1(z) + db2 * pow(z, 2.0) + db3 * pow(z, 3.0);
		entropy += 0.5*xb2 * pow(z, 2.0) + 0.5*xb3 * pow(z, 3.0);
		entropy *= 2.0*teff;
		entropy -= log(z);

		// ==> energy, E/E_0, where E_0 is the energy at T=0;
		energy  = xb2 * pow(z, 2.0) + xb3 * pow(z, 3.0);
		energy *= 2.0*teff*teff;
		energy += pressure;
        
        contact = energy - pressure;
		// ==> \kappa / \kappa_0;
		kappa = integral_kappa1(z) + 4.0 * db2 * pow(z, 2.0) + 9.0 * db3 * pow(z, 3.0);

		// ==> negative kappa indicates the breakdown of virial expansion;
		//if(kappa < 0.0) continue;

		// ==> output the pk relation;
		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << pressure << "    ";

		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << entropy << "    ";
		
		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << energy << "    ";		
		
		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << kappa << "    ";

		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << z << "    ";
    
        eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << contact  << "    ";

		eos2d << setprecision(16) << setw(24);
		eos2d << setiosflags(ios::scientific);
		eos2d << teff << "    ";
		
		eos2d << endl;
	}

	eos2d.close();
    //y[0] = (double)1;
    
}

//nlhs = number of expected output mxArrays
//nrhs = number of input mxArrays
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x1,*x2,*x3,*y;
  size_t mrows,ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgIdAndTxt( "MATLAB:huicoeffs:invalidNumInputs",
            "Three input required.");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:huicoeffs:maxlhs",
            "Too many output arguments.");
  }
  
  /* The input must be a noncomplex scalar double.*/
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:huicoeffs:inputNotRealScalarDouble",
            "Input must be a noncomplex scalar double.");
  }
  
  //printf("mrows = %d, ncols = %d \n", mrows, ncols);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
  
  /* Assign pointers to each input and output. */
  x1 = mxGetPr(prhs[0]);
  x2 = mxGetPr(prhs[1]);
  x3 = mxGetPr(prhs[2]);
  y = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */
  huicoeffs(y,x1,x2,x3);
}
