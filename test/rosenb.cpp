
// Daniel Pitzl, 2019
// bobyqa for Rosenbrock

#include <iostream> // for cout
#include <iomanip> // for setw
#include <cmath> // exp
#include <sys/ioctl.h>

#include "TStyle.h" // ROOT
#include "TApplication.h" // ROOT
#include "TCanvas.h"
#include "TFile.h" // ROOT
#include "TProfile2D.h"
#include "TGraph.h"
#include "TSystem.h"

#include "/home/pitzl/ROOT/bobyqa/bobyqa.h"

using namespace std;

vector <double> x0;
vector <double> x1;

//------------------------------------------------------------------------------
double rosenbrock( const INTEGER n, const REAL * x, void * data )
{
  x0.push_back( x[0] );
  x1.push_back( x[1] );
  return
    1 * pow( x[1] - pow( x[0], 2 ), 2 ) + pow( 1-x[0], 2 );
}

//------------------------------------------------------------------------------
bool kbhit()
{
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
}

//----------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
  // set ROOT styles:

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.4, "y" );
  gStyle->SetTitleOffset( 1.9, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 111111 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.99 );
  gStyle->SetStatY( 0.60 );

  gStyle->SetPalette( 55 ); // sunset colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate(  );

  TApplication app( "app", 0, 0 );
  TCanvas c1( "c1", "plane", 0, 0, 936, 837 ); // square, global

  c1.SetBottomMargin( 0.15 );
  c1.SetLeftMargin( 0.15 );
  c1.SetRightMargin( 0.20 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TFile rootFile( "rosenb.root", "recreate" );

  TProfile2D * p2 = new
    TProfile2D( "rosenbrock", "BObyQA for Rosenbrock;x;y;Rosenbrock",
		81, -5.05, 3.05, 81, -2.05, 6.05 );

  for( double xx = -5; xx < 3.05; xx += 0.1 ) {

    double x2 = pow( xx, 2 );
    double a2 = pow( 1-xx, 2 );

    for( double yy = -2; yy < 6.05; yy += 0.1 ) {
      double rb = 1 * pow( yy - x2, 2 ) + a2;
      p2->Fill( xx, yy, rb + 0.01 );
    }

  }
  p2->SetStats(0);
  gStyle->SetOptStat(10);
  gPad->SetLogz(1);
  p2->Draw("colz");
  c1.Update();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // min:

  const INTEGER npar = 2; // fit parameters

  REAL par[npar];
  par[0] = -4; // traditional: -1, 2
  par[1] =  4;
  //par[0] = -1; // traditional: -1, 2
  //par[1] =  2; // 43 callsx

  REAL parl[npar]; // lower bounds on fit parameters
  parl[0] = -99;
  parl[1] = -99;

  REAL parh[npar]; // lower bounds on fit parameters
  parh[0] =  99;
  parh[1] =  99;

  //INTEGER npt = 2*npar; // 170 calls
  //INTEGER npt = 5; // 60 calls
  INTEGER npt = 6; // 57 calls
  REAL rhobeg = 0.1;
  REAL rhoend = 0.001;
  INTEGER iprint = 1;
  INTEGER maxfun = 999;

  int nw = (npt+5)*(npt+npar) + 3*npar*(npar+5)/2;
  REAL work[nw];
  int data = 0; // dummy

  int good = bobyqa( npar, npt, rosenbrock, (void*)&data, par, parl, parh,
		     rhobeg, rhoend,
		     iprint, maxfun, work );

  cout << "BOBYQA " << good << endl;
  cout << "Rosenbrock " << work[0] << endl;
  cout << "at " << par[0] << ", " << par[1] << endl;

  cout << x0.size() << " points" << endl;

  TGraph gp( x0.size(), &x0[0], &x1[0] );
  gp.SetMarkerColor(1);
  gp.SetMarkerStyle(20);
  gp.SetMarkerSize(1.0);
  gp.Draw("p");

  //for( unsigned i = 0; i < x0.size(); ++i ) cout << i << "  " << x0[i] << "  " << x1[i] << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "enter any key" << endl;

  while( !kbhit() )
    gSystem->ProcessEvents(); // ROOT

  string any;
  cin >> any;

  rootFile.Write();
  cout << rootFile.GetName() << endl;

  return 0;

} // main
