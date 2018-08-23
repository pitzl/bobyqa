
// Daniel Pitzl, 2018
// Student

#include <random>
#include <iostream> // for cout
#include <iomanip> // setw
#include <time.h> // gettimeofday

#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include <TH1I.h> // ROOT
#include <TFile.h> // ROOT
#include <TMath.h> // ROOT
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector.h"

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
REAL ChiSq( const INTEGER m, const REAL* par, void * data ) // bobyqa_objfun
{
  cout << "par";
  for( int i = 0; i < m; ++i )
    cout<< "  " << par[i];

  Data * d { (Data*)data };

  double parf[m];

  parf[0] = par[0]; // mean
  parf[1] = par[1]; // sigma
  parf[2] = par[2]; // nu

  // user parameters are signal and BG fractions
  // integrate over one bin:

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

  cout << "starting at " << ChiSq( npar, par, &data )
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

  int good = bobyqa( npar, npt, ChiSq, (void*)&data, par, parl, parh,
		     rhobeg, rhoend,
		     iprint, maxfun, work );
  cout << "BOBYQA " << good << endl;
  cout << "final chisq " << work[0] << " / " << n-npar << endl;

  double chisq0 = ChiSq( npar, par, &data );

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

  // scan FCN around minimum:
  // chisq = chisq0 + ( (p-p0)/s )^2
  // => dchi = (dp/s)^2
  // first iter: dchi1 for dp1
  // we want dchi2 = 1
  // sqrt(dchi2/dchi1) = dp2/dp1
  // => dp2 = dp1 / sqrt(dchi1)

  double xar[npar];
  for( int j = 0; j < npar; ++j )
    xar[j] = par[j];

  double dp1[npar];

  bool ldb = 0; // debug flag

  for( int j = 0; j < npar; ++j ) {

    double dp = 0.1*par[j]; // initial step
    int iter = 0;
    bool again = 1;
    double dchi = 0;

    do {
      iter++;
      xar[j] = par[j] + dp;
      double chisq = ChiSq( npar, xar, &data );
      dchi = chisq - chisq0;
      if( ldb )
	cout << "   par " << j << " iter " << iter << "  " << xar[j]
	     << " dchi2 " << dchi << endl;

      if( dchi > 2 )
	dp = dp/sqrt(dchi);
      else if( dchi < 0.01 )
	dp = 2*dp;
      else if( dchi < 0.5 )
	dp = dp/sqrt(dchi);
      else
	again = 0;
      if( iter > 9 ) again = 0;
    }
    while( again );

    cout << "par " << j
	 << "  " << par[j]
	 << ", iter " << iter
	 << ", step " << dp
	 << ", dchi " << dchi
	 << endl;

    xar[j] = par[j]; // back to minimum
    dp1[j] = dp;

  } // j par

  // Hessian for 2nd order derivatives:
  // f''(x) = ( f(x + h) − 2f(x) + f(x − h) ) / h^2
  // d2f/dxdy =
  // ( f(a+h1, b+h2) - f(a+h1, b-h2) - f(a-h1, b+h2) + f(a-h1, b-h2) / (4 h1 h2)

  double H[npar][npar];

  for( int j = 0; j < npar; ++j ) {

    double dpj = dp1[j];
    xar[j] = par[j] + dpj;
    double chisq = ChiSq( npar, xar, &data );

    double dchiup = chisq - chisq0;
    xar[j] = par[j] - dpj;
    chisq = ChiSq( npar, xar, &data );

    double dchidn = chisq - chisq0;
    double f2nd = ( dchiup + dchidn ) / (dpj*dpj);
    H[j][j] = f2nd;

    for( int k = j+1; k < npar; ++k ){

      double dpk = dp1[k];

      xar[j] = par[j] + dpj;
      xar[k] = par[k] + dpk;
      double lupup = ChiSq( npar, xar, &data );

      xar[k] = par[k] - dpk;
      double lupdn = ChiSq( npar, xar, &data );

      xar[j] = par[j] - dpj;
      double ldndn = ChiSq( npar, xar, &data );

      xar[k] = par[k] + dpk;
      double ldnup = ChiSq( npar, xar, &data );

      double df2didj = ( lupup - lupdn - ldnup + ldndn ) / ( 4*dpj*dpk);
      H[j][k] = df2didj;
      H[k][j] = df2didj;

      xar[k] = par[k]; // back to minimum

    } // k

    xar[j] = par[j]; // back to minimum

  } // j

  cout << endl;

  // invert H

  TMatrixD hesse( npar, npar );

  for( int j = 0; j < npar; ++j )
    for( int k = 0; k < npar; ++k )
      TMatrixDRow( hesse, j )(k) = 0.5*H[j][k]; // justify factor 1/2

  TDecompSVD svd( npar, npar );
  svd.SetMatrix( hesse );
  svd.Decompose();

  TVectorD eigen = svd.GetSig();
  cout << "Hesse Eigenvalues";
  for( int j = 0; j < npar; ++j )
    cout << "  " << eigen(j);
  cout << endl;
  cout << "Hesse condition number  " << svd.Condition() << endl;

  Double_t d1, d2;
  svd.Det( d1, d2 );
  cout << "Hesse Det  " << d1 * pow( 2, d2 ) << endl;

  TMatrixD covar( npar, npar );
  covar = svd.Invert();
  //covar = hesse.Invert();

  cout << endl;

  double sigma[npar];
  for( int j = 0; j < npar; ++j ) {
    sigma[j] = sqrt( fabs( TMatrixDRow( covar, j )(j) ) );
    cout << "par " << j
	 << ":  " << par[j]
	 << " +- " << sigma[j] // agress with Minuit
	 << endl;
  }

  // print correlations:

  cout << endl;

  for( int k = 0; k < npar; ++k )
    cout << setw(14) << k;
  cout << endl;

  for( int j = 0; j < npar; ++j ) {
    cout << setw(9) << j << "  ";
    cout << setw(j*14) << "";
    for( int k = j; k < npar; ++k ) {
      int iw = 14;
      if( k == j ) iw = 3;
      cout << setw(iw)
	   << TMatrixDRow( covar, j )(k) / sigma[j] / sigma[k];
    }
    cout << endl;
  }

  cout << endl;

  rootFile.Write();
  rootFile.Close();
  cout << rootFile.GetName() << endl;

  return 0;

} // main
