#pragma once
#include <TH2.h>
#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <stdbool.h>
#include <vector>
#include <iostream>
#include <vector>
#include <math.h> 
#include <TRandom.h>
#include <TColor.h>
#include <TPaveStats.h>
#include <TList.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TArrow.h>
#include <TMath.h>
#include <TFile.h> 
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <TMultiGraph.h>
#include <fstream>

using namespace std ;

TH1F *hist =0 ;
TH1F *H1_pp_e1 = 0;
TH1F *H1_pp_E1 = 0;

// we define the upper limit of the fitting range
double const upper_limit = 3.5;

// we defin the number of parameters of each function
int const par_size_expo = 3;
int const par_size_boltzmann = 3;
int const par_size_levy = 4;
int const par_size_blast_wave = 5;



double simpson(double f(double*,double*), double *par, double a, double b, int n)   ;

double expo( double x , double *par);

double expo_law(double *x , double *par);

double boltzmann(double *x , double *par);

double levyPT(double *x , double *par);

double levy(double *x , double *par);

double power_law_Five(double *x , double *par);

double power_law_Seven(double *x, double*par) ;

double blast_wave(double *x , double *par);

int UpperBin_Fit(double Upper_limit);

void fcn_blast_wave(int &npar, double *gin, double &f, double *par, int iflag);

void fcn_levy(int &npar, double *gin, double &f, double *par, int iflag);

double mean_p_T( double func_num(double*,double*), double func_denom(double*,double*), double *par , double b);


void Data_Model(double X[], double Y[], double Xerr[], double Yerr[], double func(double *,double *), double *par, int par_size, double *covar);

double Integral(double func(double *,double *), double *par, int par_size, double *covar,double a, double b);

double Integral_E(double func(double *,double *), double *par, int par_size, double *covar,double a, double b);


