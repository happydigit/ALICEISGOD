#include "Header.h"

using namespace std ;

double simpson(double f(double*,double*), double *par, double a, double b, int n)   
{
    double h = (b-a) / n;    // definition of the integration step
    double z = 0;           // initialize the variable
    // we define different pointers that we will use as input parameters for our f function
    double *x = new double;  
    double *x1 = new double;
    double *x2 = new double;
    for(int i = 0; i <= n-1  ; i++)
    {
        x[0] = a + i*h;   // x_i
        x1[0] = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2
        x2[0] = a + (i+1)*h;   // x_{i+1}
        z = z + f(x, par) + 4*f(x1,par) + f(x2 , par);
        //cout << "Simpson iteration " << i << endl;
    }
    return h*z/6;
}


double expo( double x , double *par)
{
	double m_T = sqrt(x*x + par[1]*par[1]);
	double res = par[0]*TMath::Power( par[2]*(par[1]+par[2]), -1)*x*exp(- (m_T - par[1]) / par[2]);
	return res;
}

		
// Exponential law
// par[0] = C
// par[1] = m_0
// par[2] = T
double expo_law(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*TMath::Power( par[2]*(par[1]+par[2]), -1)*x[0]*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Boltzmann distribution
// par[0] = C
// par[1] = m_0
// par[2] = T
double boltzmann(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*x[0]*m_T*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Levy-Tsallis distribution
// par[0] = C
// par[1] = m
// par[2] = T
// par[3] = n

double levyPT(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = x[0]*par[0]*(par[3]-1)*(par[3]-2)/(par[3]*par[2]*(par[3]*par[2]+par[1]*(par[3]-2))) *x[0]*TMath::Power(1+(m_T - par[1]) / (par[3]*par[2]) , - par[3]);
	return res;
}



double levy(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*(par[3]-1)*(par[3]-2)/(par[3]*par[2]*(par[3]*par[2]+par[1]*(par[3]-2))) *x[0]*TMath::Power(1+(m_T - par[1]) / (par[3]*par[2]) , - par[3]);
	return res;
}


// Power law from [5] in biblio
// par[0] = C = dN / dy
// par[1] = p_0
// par[2] = n
double power_law_Five(double *x , double *par)
{
	double res = par[0]*4*(par[2]-1)*(par[2]-2)*TMath::Power(par[2]-3,-2)*TMath::Power(par[1],-2)*x[0]*TMath::Power(1+ x[0] / (par[1]* (par[2]-3)/2 ) , - par[2]);
	return res;
}


// "Power law" distribution  from [7] in biblio
// par[0] = C
// par[1] = p0
// par[2] = n
double power_law_Seven(double *x, double*par) 
{
	double res = par[0]*x[0]*TMath::Power( 1 + TMath::Power(x[0] / par[1] , 2) , - par[2] );
	return res ;
}



// Blast-wave
// par[0] = A, normalization constant
// par[1] = m , the mass of the studied particle (which will be fixed)
// par[2] = T , the kinetic freeze-out temperature
// par[3] = n , the velocity profile
// par[4] = beta_s
double blast_wave(double *x , double *par)
{
    // x[0] = p_T  here
// we define the variables and parameters
	double a = 0.0;
	//double R = 12*TMath::Power(10, -15);
	double R = 12 ;
	double b = 1.0;
	int n = 500;	
	double m_T = TMath::Sqrt(par[1]*par[1] + x[0]*x[0]);
// we integrate over r
	double h = (b-a) / n;
	double z = 0;
	double r_0 , r_1 , r_2;  
	for(int i = 0; i <= n-1  ; i++)
	{
	    r_0 = a + i*h;   // x_i
	    r_1 = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2               
	    r_2 = a + (i+1)*h;   // x_{i+1}
	    z +=  x[0]*m_T*par[0]*R*R*2*TMath::Pi()*(   r_0*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_0 , par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_0 , par[3] ) * par[4] ) ) / par[2])  
	    + 4*r_1*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_1 , par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_1 , par[3] ) * par[4] ) ) / par[2])  
	    + r_2*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_2 , par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_2 , par[3] ) * par[4] ) ) / par[2])  );
	}
	z =  h*z/6;
	return z;
}

// function that returns the number of the bin corresponding to the upper limit of the fitting range (b)
int UpperBin_Fit(double Upper_limit)
{
	double a = hist->GetBinLowEdge(1); // lower edge of fitting range that is fixed	
	double p_T = a;        // p_T is going to go from a to b to explore [a,b]             
	int counter = 1;		// we begin by studying the first bin	
	while(p_T < Upper_limit) 
	{
		p_T += hist->GetBinWidth(counter);
		counter+=1;
		//cout << "counter = " << counter << endl;
		//cout << "p_T = " << p_T << endl;
		//cout << " " << endl;

	}
	//cout << "UpperBin = " << counter << endl;
	return counter-1;
}


void fcn_blast_wave(int &npar, double *gin, double &f, double *par, int iflag)
{
	f = 0.;
	// we can modify the fitting range [a,b]
	double a = hist->GetBinLowEdge(1); // lower limit of fitting range that is fixed	
	double p_T = a;        // p_T is going to go from a to b to explore [a,b]             
	int i = 1;		// we begin by studying the first bin
	// when p_T > b, the minimization stops
	for(int i=1;i<=UpperBin_Fit(upper_limit);i++)
	{
		double *x = new double ;
		x[0] = hist->GetBinCenter(i);
		double measure = hist->GetBinContent(i);
		double error = hist->GetBinError(i);
		double func = blast_wave(x,par);
		double delta = (func - measure)/error;
		f += delta*delta;	
	}
}

void fcn_levy(int &npar, double *gin, double &f, double *par, int iflag)
{
	f = 0.;
	double range = 3.3 ;
	//hist->GetXaxis()->SetRangeUser(0.3,3);
	for(int i=1;i<=hist->GetNbinsX();i++) {
	double *x = new double ;
	x[0] = hist->GetBinCenter(i);
	double measure = hist->GetBinContent(i);
	double error = hist->GetBinError(i);
	//double func = levy(x,par);
	double func = levy(x,par);
	double delta = (func - measure)/error;
	f += delta*delta;
		if (x[0] >= range)
		{
		i = (hist->GetNbinsX()+1 ) ;
		}
	}
}

double mean_p_T( double func_num(double*,double*), double func_denom(double*,double*), double *par , double b){
double numerator = simpson(func_num, par , 0 , b , 10000);
double denominator = simpson(func_denom, par , 0 , b , 10000);
double res = numerator / denominator;
return res ;
}


// Data / Model fit 
// We compute the ratio data/model for each bin of a given histogram hist
// and store the values in tables
void Data_Model(double X[], double Y[], double Xerr[], double Yerr[], double func(double *,double *), double *par, int par_size, double *covar)
{
	double data_value , model_value , model_error;
	double integral, low_edge, up_edge, bin_width ;
	int Upper_Bin = UpperBin_Fit(upper_limit);

// we define a TF1 function and inject the fitted parameters
// we use a TF1 function to have access to Integral() and IntegralError()
	TF1 *f = new TF1("fit function",func,0,upper_limit,par_size);
		f->SetParameters(par);
	
// we go through each bin of the fitting range	
	for(int i=1;i<=Upper_Bin;i++) 
	{
	// for each bin, we get its content
		data_value = hist->GetBinContent(i);
	// we compute the mean value of the bin via the model
	// i.e. we comute the integral of the model on the bin's width and divide by the bin's width
		low_edge = hist->GetBinLowEdge(i);
		up_edge = hist->GetBinLowEdge(i);
		bin_width = hist->GetBinWidth(i);
		up_edge += bin_width;
		//integral = simpson(func,par,low_edge,up_edge,n);
		integral = f->Integral(low_edge,up_edge);
		model_value = integral / bin_width;
		model_error = f->IntegralError(low_edge,up_edge,par,covar) / bin_width;
	// then we compute the ratio data/model fit
		Y[i-1] = data_value / model_value;
	// we compute the uncertainty on this Y value via propagation of uncertainty
	// considering the data and model uncertainties are uncorrelated 
		Yerr[i-1] = TMath::Sqrt(  TMath::Power(hist->GetBinError(i) / model_value,2) + TMath::Power( data_value* model_error /(model_value*model_value) ,2) ) ; 
	// we also get the center value of the bin for the future plot
		X[i-1] = hist->GetBinCenter(i);
		Xerr[i-1] = 0.5*bin_width;  // Xerr will represent the width of the bin (it is not an error strictly speaking)
		cout << "iteration i = " << i << endl;
		cout << "model_value = " << model_value << endl;
		cout << "model error = " << f->IntegralError(low_edge,up_edge,par,covar) << endl;
		cout << "bin error/model = " << hist->GetBinError(i) / model_value << endl;
		cout << "Y " << Y[i-1] << endl;
		cout << "Yerr " << Yerr[i-1] << endl;
		cout << " " << endl;  
	}	
	
}



double Integral(double func(double *,double *), double *par, int par_size, double *covar,double a, double b)
{

	TF1 *f = new TF1("fit function",func,a,b,par_size);
		f->SetParameters(par);	
	double integral = f->Integral(a,b);		
	return integral ;
}

double Integral_E(double func(double *,double *), double *par, int par_size, double *covar,double a, double b)
{

	TF1 *f = new TF1("fit function",func,a,b,par_size);
		f->SetParameters(par);	

	double integral_error = f->IntegralError(a,b,par,covar) ;		
	return integral_error;
}
