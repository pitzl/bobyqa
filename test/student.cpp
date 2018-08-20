
// Daniel Pitzl, 2018
// Student

#include <random>
#include <iostream> // for cout
#include <time.h> // gettimeofday

#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include <TH1I.h> // ROOT
#include <TFile.h> // ROOT
#include <TMath.h> // ROOT

#include "/home/pitzl/ROOT/bobyqa/bobyqa.h"

using namespace std;

//------------------------------------------------------------------------------
double fitStudent( double x, const double *par )
{
  double t = ( x - par[0] ) / par[1]; // mean and sigma

  // exponent:

  double rn = par[2];
  double xn = 0.5 * ( rn + 1.0 );

  // Normalization needs Gamma function:

  double pk = 0.0;

  if( rn > 0.0 && fabs( xn * log( 1.0 + t*t/rn ) ) < 333 ) {

    double pi = 3.14159265358979323846;

    double aa = par[3] / par[1] / sqrt(rn*pi) *
      TMath::Gamma(xn) / TMath::Gamma(0.5*rn); // norm, par[3] = A*dx

    pk = aa * exp( -xn * log( 1.0 + t*t/rn ) );

    // lim n->inf (1+a/n)^n = e^a

  }

  return pk + par[4]; // BG
}

//----------------------------------------------------------------------------
struct Data {
  int n; // data points = length of arrays x,y,e
  double A; // area
  double *x;
  double *y;
  double *e;
};

//----------------------------------------------------------------------------
REAL chisq( const INTEGER m, const REAL* par, void * data ) // bobyqa_objfun
{
  cout << "par";
  for( int i = 0; i < m; ++i )
    cout<< "  " << par[i];

  Data * d { (Data*)data };

  double parf[m];

  parf[0] = par[0]; // mean
  parf[1] = par[1]; // sigma
  parf[2] = par[2]; // nu

  double dx = d->x[1] - d->x[0]; // bin width
  parf[3] = par[3]*dx; // bin area

  double xrange = d->x[d->n-1] - d->x[0] + dx;
  parf[4] = par[4]*dx/xrange; // normed BG

  double c2 = 0;
  for( int i = 0; i < d->n; ++i ) {
    double x = d->x[i];
    double f = d->A * fitStudent( x, parf );
    double y = d->y[i];
    double e = d->e[i];
    if( e > 1e-14 )
      c2 += pow( (y-f)/e, 2 ); // chisq
  }
  cout << "  chisq " << c2 << endl;

  return c2;
}

//----------------------------------------------------------------------------
int main()
{
  TFile rootFile( "student.root", "recreate" );

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.4, "y" );

  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetLabelOffset( 0.022, "xyz" );

  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(111110);
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  gStyle->SetStatBorderSize(1); // no 'shadow'
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.88);

  gStyle->SetPalette(55); // sunset

  TApplication app( "app", 0, 0 );
  TCanvas c1( "c1", "data and fit", 1100, 0, 900, 900 ); // square
  c1.SetTopMargin( 0.12 );
  c1.SetLeftMargin( 0.16 );
  c1.SetRightMargin( 0.05 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // randoms:
  /*
    default_random_engine
    minstd_rand // Minimal Standard linear_congruential
    minstd_rand0 // Minimal Standard
    mt19937 // Mersenne Twister
    mt19937_64 // Mersene Twister 64 bit
    ranlux24_base
    ranlux48_base
    ranlux24
    ranlux48 // factor 3 slower
    knuth_b
  */
  default_random_engine gen;

  gen.seed( time(NULL) ); // seconds since 1.1.1970

  double nu = 3;
  student_t_distribution<double> stu( nu ); // mean zero, sig 1

  // background:

  double xmin = -8;
  double xmax =  8;
  uniform_real_distribution uni( xmin, xmax );
  double fbg = 0.25; // background fraction

  TH1I hd( "hd", Form( "Student(0,1,%f) + %f*BG;x;n/bin", nu, fbg ),
	   100, xmin, xmax );

  // generate:

  int N = 1000*1000;

  for( int ievt = 0; ievt < N; ++ievt )
    hd.Fill( stu(gen) );

  for( int ievt = 0; ievt < fbg*N; ++ievt )
    hd.Fill( uni(gen) );

  hd.Draw("e");
  c1.Update();

  // get data points:

  int n = hd.GetNbinsX();
  double x[n];
  double y[n];
  double e[n];
  for( int i = 0; i < n; ++i ) {
    x[i] = hd.GetBinCenter(i+1); // bins 0 is underflow
    y[i] = hd.GetBinContent(i+1); // counts
    e[i] = hd.GetBinError(i+1); // Poisson error = sqrt(count)
  }
  double A = hd.GetSumOfWeights(); // inside

  Data data;
  data.n = n; // bins
  data.A = A; // area
  data.x = x;
  data.y = y; // measured
  data.e = e; // error

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // prepare fit:

  INTEGER npar = 5;
  REAL par[npar];

  double xm = hd.GetMean();
  double sm = 0.5*hd.GetRMS();
  double bg = hd.GetBinContent(1); // BG
  par[0] = xm;
  par[1] = sm;
  par[2] = 3.3; // nu

  // parameter values should have similar magnitudes (for trust region rho)
  // don't use (huge) A as parameter, but rather fractions for signal and BG

  par[3] = 1-n*bg/A; // signal fraction
  par[4] = n*bg/A; // BG fraction

  cout << "starting at " << chisq( npar, par, &data )
       << " / " << n-npar << endl;

  // plot:

  const int npnt = 500; // for plotting
  vector <double> vx(npnt);
  vector <double> vy(npnt);
  double parf[npar];
  parf[0] = par[0]; // mean
  parf[1] = par[1]; // sig
  parf[2] = par[2]; // nu
  double dx = x[1] - x[0]; // bin width
  parf[3] = par[3]*dx; // signal
  double xrange = x[n-1] - x[0] + dx;
  parf[4] = par[4]*dx/xrange; // normed BG
  for( int ii = 0; ii < npnt; ++ii ) {
    vx[ii] = xmin + ii*(xmax-xmin)/(npnt-1);
    vy[ii] = A*fitStudent( vx[ii], parf );
  }

  TGraph * gp = new TGraph( vx.size(), &vx[0], &vy[0] );
  gp->SetLineColor(6);
  gp->SetLineWidth(3);
  gp->Draw("l");
  c1.Update();

  cout << "enter any key" << endl;
  string any;
  cin >> any;

  REAL parl[npar]; // lower bounds on fit parameters
  parl[0] = x[0]; // mean
  parl[1] = x[1]-x[0]; // min sigma = bin width
  parl[2] = 1; // min n u
  parl[3] = 0; // signal
  parl[4] = 0; // BG

  REAL parh[npar];
  parh[0] = x[n-1]; // mean
  parh[1] = (x[n-1]-x[0])/3; // max sigma
  parh[2] = 333; // max nu
  parh[3] = 1; // signal
  parh[4] = 1; // BG

  INTEGER npt = 9;
  REAL rhobeg = 0.1;
  REAL rhoend = 0.00001;
  INTEGER iprint = 1;
  INTEGER maxfun = 999;

  int nw = (npt+5)*(npt+npar) + 3*npar*(npar+5)/2;
  REAL work[nw];

  int good = bobyqa( npar, npt, chisq, (void*)&data, par, parl, parh,
		     rhobeg, rhoend,
		     iprint, maxfun, work );
  cout << "BOBYQA " << good << endl;
  cout << "final chisq " << work[0] << " / " << n-npar << endl;

  chisq( npar, par, &data );

  parf[0] = par[0]; // xm
  parf[1] = par[1]; // s
  parf[2] = par[2]; // nu
  parf[3] = par[3]*dx; // signal
  parf[4] = par[4]*dx/xrange; // normed BG
  for( int ii = 0; ii < npnt; ++ii ) {
    vx[ii] = xmin + ii*(xmax-xmin)/(npnt-1);
    vy[ii] = A*fitStudent( vx[ii], parf );
  }

  hd.Draw("e");
  TGraph * gf = new TGraph( vx.size(), &vx[0], &vy[0] );
  gf->SetLineColor(6);
  gf->SetLineWidth(3);
  gf->Draw("l");
  c1.Update();

  cout << "enter any key" << endl;
  cin >> any;

  rootFile.Write();
  rootFile.Close();
  cout << "rlq " << rootFile.GetName() << endl;
  cout << ".x fittp0.C+(\"hd\")" << endl;
  return 0;

} // main
