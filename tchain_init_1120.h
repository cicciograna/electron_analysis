// General variables
double CX_SR_LOW( 3. ), CX_SR_HIGH( 59.5 );
double CX_NODE( 0.004 ); // 4100 V - TEMPORARY FIX
double CX_TOF( 0.17 ); // GEANT4 + COMSOL BFIELD - ONLY FORWARD ELECTRONS - TEMPORARY FIX
// double CX_NODE( 0.006 ); // 4100 V
// double CX_TOF( 0.166 ); // GEANT4 + COMSOL BFIELD - ONLY FORWARD ELECTRONS
// double CX_TOF( 0.25 ); // TEMPORARY GEANT4 + COMSOL BFIELD - ONLY FORWARD ELECTRONS
// double CX_TOF( 0.178 ); // TEMPORARY GEANT4 + COMSOL BFIELD - ONLY FORWARD ELECTRONS
double CX_TOF_HI_LOW( 0.1675 ); // TEMPORARY 4100 V - ONLY FORWARD ELECTRONS
double CX_TOF_HI_HIGH( 0.1835 ); // TEMPORARY 4100 V - ONLY FORWARD ELECTRONS
double CX_TOF_LO( 0.2055 ); // 4100 V - SOME BACKWARD ELECTRONS

TF1 *pp3 = new TF1( "pp3","pol3",0.18,CX_TOF_LO );

// POL3
double CX_TOF_p0 = 54468.141523626706;
double CX_TOF_p1 = -847942.726966361;
double CX_TOF_p2 = 4.385758174956211e6;
double CX_TOF_p3 = -7.532956079857603e6;


// PARTICLE STRUCTURE

double CF = 931494.095; // Conversion factor amu -> keV

double c = 299792458; // m/s

double Bfield = 8.5;   // Gauss
double auger_mass = 510.998910;   // keV/c²
double omega =  8.98755e3 * Bfield / auger_mass;    // us^-1

double DEsplat_radius_MCP = 0.1;   // mm = 100 um
// double DEsplat_radius_MCP = 0.05;   // mm = 50 um
// double DEtof_MCP = 0.00035;   // us = 0.35 ns
double DEtof_MCP = 0.00085;   // us = 0.85 ns
// double DEtof_MCP = 0.001;   // us = 1 ns
double DEBfield = 8*0./100.;
double DEEfield = 0.4100*0./100.;
double DELelectron = 1000*0./100.;


// Origin position ( mm )
double origin_x = 0.;
double origin_y = 0.;
double origin_z = 1120.;

double vessel_inner_radius = 300.;

double MCP_radius = 60.;

double ion_MCP_position = 159.;

double src_inner_radius = 70.;
double src_outer_radius = src_inner_radius + 25.;

double ele_inner_radius = 171.;
double ele_outer_radius = ele_inner_radius + 75.;
double ele_MCP_position = origin_z + 932.;

// // Default
// // No PMT Shield
// TFile *isos = new TFile( "isos_GEANT4phi1_theta1_KE1-150_1eV_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut_MuFe_WithSleeves_NoPMTSh_pa20210111.root"   ); double Efield = 0.424; double nodes[] = { 0.1255, 0.1663, 0.2032, 0.2475 }; // Manually inserted nodes from COMSOL generated Bfieldmap

// // With PMT Shield
TFile *isos = new TFile( "isos_GEANT4phi1_theta1_KE1-150_1eV_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut_MuFe_WithSleeves_WithPMTSh_pa20210111.root" ); double Efield = 0.424; double nodes[] = { 0.1255, 0.1665, 0.2100, 0.2475 }; // Manually inserted nodes from COMSOL generated Bfieldmap

// // With PMT Shield @ 450mm
// TFile *isos = new TFile( "isos_GEANT4phi1_theta1_KE1-150_1eV_WithPMT450mm_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root" ); double Efield = 0.424; double nodes[] = { 0.1255, 0.1665, 0.2100, 0.2475 }; // Manually inserted nodes from COMSOL generated Bfieldmap

// // dθ = 0.5, 5eV around standard energies, No PMT Shield
// TFile *isos = new TFile( "isos_GEANT4phi1_theta0.5_StandardErange_5eV_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root" ); double Efield = 0.424; double nodes[] = { 0.1255, 0.1665, 0.2100, 0.2475 }; // Manually inserted nodes from COMSOL generated Bfieldmap

// // dθ = 0.5, 5eV around standard energies, dE = 0.25eV, No PMT Shield
// TFile *isos = new TFile( "isos_GEANT4phi1_theta0.5_StandardErange5eV_dE025eV_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root" ); double Efield = 0.424; double nodes[] = { 0.1255, 0.1665, 0.2100, 0.2475 }; // Manually inserted nodes from COMSOL generated Bfieldmap

TGraph2D *iso_init_prad       = ( TGraph2D* )isos->Get( "iso_init_prad" );
TGraph2D *iso_init_plong      = ( TGraph2D* )isos->Get( "iso_init_plong" );
TGraph2D *iso_init_azimuth    = ( TGraph2D* )isos->Get( "iso_init_azimuth" );
TGraph2D *iso_init_elevation  = ( TGraph2D* )isos->Get( "iso_init_elevation" );
TNamed *iso_comment           = ( TNamed* )  isos->Get( "iso_comment" );

// TTree mode
double tfmin(0.086), tfmax(0.206);
double srmin(0.000), srmax(60.00);
double samin(0.000), samax(360.0);
TH3D *htree = new TH3D( "htree","htree", 240, tfmin, tfmax, 60, srmin, srmax, 360, samin, samax );


// ISO_CURVE_THISTO - COMMENTA DA QUI...
TTree *tthist = (TTree*)isos->Get("tthist");
TTree *rpb = (TTree*)isos->Get("run_per_bin");
double ltf, lsr, lsa, lpr;
int tbi;
vector<int> *runids;
/// ...A QUI, VEDI DOPO


// double prad_cov     = iso_init_prad   ->GetHistogram()->GetCovariance();
// double plong_cov    = iso_init_plong  ->GetHistogram()->GetCovariance();
// double azimuth_cov  = iso_init_azimuth->GetHistogram()->GetCovariance();

TRandom3 *ra = new TRandom3( 0 );

// Definition of the "particle" class
class particle{

public:
  double  ionn;           // ion number
  double  tof;            // tof [us]
  double  mass;           // mass of the particle [amu]
  double  charge;         // charge of the particle [e]
  double  init_x;         // initial x position [mm]
  double  init_y;         // initial y position [mm]
  double  init_z;         // initial z position [mm]
  double  init_azm;       // initial azimuth [deg]
  double  init_elev;      // initial angle with the X axis [deg] [NOTE: IT'S DIFFERENT FROM ELEVATION! GOES FROM 0° TO 90°!]
  double  init_dircos;    // initial director cosine
  double  init_vx;        // initial vx ( X-component of the velocity ) [mm/us]
  double  init_vy;        // initial vy ( Y-component of the velocity ) [mm/us]
  double  init_vz;        // initial vz ( Z-component of the velocity ) [mm/us]
  double  init_vtot;      // initial total velocity [mm/us]
  double  init_vrad;      // initial vrad ( radial component of the velocity = √( vy²+vz² ) ) [mm/us]
  double  init_vlong;     // initial vlong ( longitudinal component of the velocity = vx ) [mm/us]
  double  init_px;        // initial px ( X-component of the momentum ) [keV/c]
  double  init_py;        // initial py ( Y-component of the momentum ) [keV/c]
  double  init_pz;        // initial pz ( Z-component of the momentum ) [keV/c]
  double  init_ptot;      // initial total momentum [keV/x]
  double  init_prad;      // initial prad ( radial component of the momentum = √( py²+pz² ) ) [keV/c]
  double  init_plong;     // initial plong ( longitudinal component of the momentum = px ) [keV/c]
  double  init_ke;        // initial kinetic energy ( keV )
  double  splat_x;        // splat x position [mm]
  double  splat_y;        // splat y position [mm]
  double  splat_z;        // splat z position [mm]
  double  splat_azm;      // splat azimuth [deg]
  double  splat_elev;     // splat angle with the X axis [deg] [NOTE: IT'S DIFFERENT FROM ELEVATION! GOES FROM 0° TO 90°!]
  double  splat_dircos;   // splat director cosine
  double  splat_vx;       // splat vx ( X-component of the velocity ) [mm/us]
  double  splat_vy;       // splat vy ( Y-component of the velocity ) [mm/us]
  double  splat_vz;       // splat vz ( Z-component of the velocity ) [mm/us]
  double  splat_vtot;     // splat total velocity [mm/us]
  double  splat_vrad;     // splat vrad ( radial component of the velocity = √( vy²+vz² ) ) [mm/us]
  double  splat_vlong;    // splat vlong ( longitudinal component of the velocity = vx ) [mm/us]
  double  splat_px;       // splat px ( X-component of the momentum ) [keV/c]
  double  splat_py;       // splat py ( Y-component of the momentum ) [keV/c]
  double  splat_pz;       // splat pz ( Z-component of the momentum ) [keV/c]
  double  splat_ptot;     // splat total momentum [keV/x]
  double  splat_prad;     // splat prad ( radial component of the momentum = √( py²+pz² ) ) [keV/c]
  double  splat_plong;    // splat plong ( longitudinal component of the momentum = px ) [keV/c]
  double  splat_ke;       // splat kinetic energy ( keV )
  double  splat_radius;   // splat radius ( = √( y²+z² ) ) [mm]
  double  run_id;         // univocal identifier for the run
  int     ontarget;       // flag to mark that the particle hit the electron MCP
  int     onmcpplane;     // flag to mark that the particle hit the electron MCP plane



  // These vectors contain the represented quantity at each time step in the particle flight
  // Only useful for the fully recorded data
  vector<double> *vec_tof     = new vector<double>;   // tof [us]
  vector<double> *vec_x       = new vector<double>;   // x position [mm]
  vector<double> *vec_y       = new vector<double>;   // y position [mm]
  vector<double> *vec_z       = new vector<double>;   // z position [mm]
  vector<double> *vec_azm     = new vector<double>;   // azimuth [deg]
  vector<double> *vec_elev    = new vector<double>;   // elevation [deg]
  vector<double> *vec_dircos  = new vector<double>;   // director cosine
  vector<double> *vec_vx      = new vector<double>;   // vx ( X-component of the velocity ) [mm/us]
  vector<double> *vec_vy      = new vector<double>;   // vy ( Y-component of the velocity ) [mm/us]
  vector<double> *vec_vz      = new vector<double>;   // vz ( Z-component of the velocity ) [mm/us]
  vector<double> *vec_vtot    = new vector<double>;   // total velocity [mm/us]
  vector<double> *vec_vrad    = new vector<double>;   // vrad ( radial component of the velocity = √( vy²+vz² ) ) [mm/us]
  vector<double> *vec_vlong   = new vector<double>;   // vlong ( longitudinal component of the velocity = vx ) [mm/us]
  vector<double> *vec_px      = new vector<double>;   // px ( X-component of the momentum ) [keV/c]
  vector<double> *vec_py      = new vector<double>;   // py ( Y-component of the momentum ) [keV/c]
  vector<double> *vec_pz      = new vector<double>;   // pz ( Z-component of the momentum ) [keV/c]
  vector<double> *vec_ptot    = new vector<double>;   // total momentum [keV/x]
  vector<double> *vec_prad    = new vector<double>;   // prad ( radial component of the momentum = √( py²+pz² ) ) [keV/c]
  vector<double> *vec_plong   = new vector<double>;   // plong ( longitudinal component of the momentum = px ) [keV/c]
  vector<double> *vec_ke      = new vector<double>;   // kinetic energy [keV]
  vector<double> *vec_radius  = new vector<double>;   // radius [mm]]
  vector<double> *vec_volt    = new vector<double>;   // voltage [V]
  vector<double> *vec_e       = new vector<double>;   // total electric field [V/mm]
  vector<double> *vec_ex      = new vector<double>;   // ex ( X-component of the electric field ) [V/mm]
  vector<double> *vec_ey      = new vector<double>;   // ey ( Y-component of the electric field ) [V/mm]
  vector<double> *vec_ez      = new vector<double>;   // ez ( Z-component of the electric field ) [V/mm]
  vector<double> *vec_b       = new vector<double>;   // total magnetic field [Gauss]
  vector<double> *vec_bx      = new vector<double>;   // ex ( X-component of the magnetic field ) [Gauss]
  vector<double> *vec_by      = new vector<double>;   // ey ( Y-component of the magnetic field ) [Gauss]
  vector<double> *vec_bz      = new vector<double>;   // ez ( Z-component of the magnetic field ) [Gauss]

  ////////////////////////////////////////////////////////
  //    _____  _    _ _____  ______  __   ____     __   //
  //   |  __ \| |  | |  __ \|  ____| \ \ / /\ \   / /   //
  //   | |__) | |  | | |__) | |__     \ V /  \ \_/ /    //
  //   |  ___/| |  | |  _  /|  __|     > <    \   /     //
  //   | |    | |__| | | \ \| |____   / . \    | |      //
  //   |_|     \____/|_|  \_\______| /_/ \_\   |_|      //
  //                                                    //
  ////////////////////////////////////////////////////////
  
  /////////////////////////////////
  // ADDITIONAL "WORK" VARIABLES //
  /////////////////////////////////
  
  // init_elevation() [deg] - returns the initial elevation from 0 to 180
  double init_elevation() { return TMath::ACos( init_vz/init_vtot ) * TMath::RadToDeg(); }


  // init_azimuth() [deg] - returns the initial azimuth from 0 to 360
  double init_azimuth()       {
    double val = TMath::ATan2( init_vy, init_vx ) * TMath::RadToDeg();
    if( val < 0. ){ val += 360.; }  // comment for -pi to +pi
    return val;
  }

  // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
  double splat_azimuth()      {
    double val = TMath::ATan2( splat_y, splat_x ) * TMath::RadToDeg();
    if( val < 0. ){ val += 360.; }  // comment for -pi to +pi
    return val;
  }
  
  double wtfloor() { return ( omega*tof - TMath::Floor( omega*tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2.; }

  
  //////////////////////////////////////
  // THEORETICAL (COLTRIMS) VARIABLES //
  //////////////////////////////////////

  // Theoretical (COLTRIMS) reconstructed variables
  double theo_plong()     { return ( ( ( mass*( splat_z-init_z ) )/( tof*c*1e-3 ) ) - Efield/2.*tof*c*1e-6 );      }
  double theo_prad()      { return ( Bfield*splat_radius )/( 2.*33.356*TMath::Abs( TMath::Sin( omega*tof/2. ) ) ); }
  double theo_ptot()      { return TMath::Sqrt( TMath::Power( theo_prad(), 2 ) + TMath::Power( theo_plong(), 2 ) );  }
  double theo_azimuth()   { double output = 0.; output = ( splat_azimuth()*TMath::DegToRad() + ( omega*tof - TMath::Floor( omega*tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() ) /2. ) *TMath::RadToDeg(); if( output > 360. ){ output -= 360.; } return output; }
  double theo_elevation()   {
    double output = -1000000.;
    if( theo_ptot() != 0. ){ output = TMath::ACos( theo_plong() / theo_ptot() ) * TMath::RadToDeg(); }
    return output;
  }
  double theo_energy()    { return 1e3 * ( TMath::Sqrt( TMath::Power( auger_mass, 2 ) + TMath::Power( theo_ptot(), 2 ) ) - auger_mass ); }

  // Theoretical (COLTRIMS) errors  
  double delta_theo_plong(){
    double output = 0.;
    double Lelectron = splat_z-init_z;
    double dplong_dt_deltat = ( - ( mass * Lelectron ) / ( tof*tof*c*1e-3 ) - Efield/2.*c*1e-6 ) * DEtof_MCP ;
    double dplong_dL_deltaL = mass / tof * DELelectron ;
    double dplong_dE_deltaE = -tof / 2. * DEEfield ;
    output = TMath::Sqrt( TMath::Power( dplong_dt_deltat, 2 ) + TMath::Power( dplong_dL_deltaL, 2 ) + TMath::Power( dplong_dE_deltaE, 2 ) );
    return output;
  }
  
  double delta_theo_prad(){
    double output = 0.;
    double wt2 = omega*tof/2.;
    double prad_r = ( Bfield )/( 2.*33.356*TMath::Abs( TMath::Sin( omega*tof/2. ) ) ) ;
    double dprad_dr_deltar = prad_r * DEsplat_radius_MCP ;
    double dprad_dt_deltat = prad_r * ( omega*splat_radius/2. ) * ( DEtof_MCP / TMath::Tan( wt2 ) ) ;
    double dprad_dB_deltaB = prad_r * ( ( splat_radius/Bfield ) - ( ( wt2/Bfield ) / TMath::Tan( wt2 ) ) ) * DEBfield ;
    output = TMath::Sqrt( TMath::Power( dprad_dr_deltar, 2 ) + TMath::Power( dprad_dt_deltat, 2 ) + TMath::Power( dprad_dB_deltaB, 2 ) );
    return output;
  }
  
  double delta_theo_azimuth(){
    double output = 0.;
    double dphi_dtheta_deltatheta = DEsplat_radius_MCP / splat_radius ;
    double dphi_dt_deltat         = -omega / 2. * DEtof_MCP ; 
    double dphi_dB_deltaB         = -( omega * tof ) / Bfield * DEBfield ; 
    output = TMath::Sqrt( TMath::Power( dphi_dtheta_deltatheta, 2 ) + TMath::Power( dphi_dt_deltat, 2 ) + TMath::Power( dphi_dB_deltaB, 2 ) ) * TMath::RadToDeg();
    return output;
  }
  
  double delta_theo_ptot(){
    // δptot = √[ ( ∂ptot/∂px · δpx )² + ( ∂ptot/∂py · δpy )² + ( ∂ptot/∂pz · δpz )² ] = √[ ( px · δpx )² + ( py · δpy )² + ( pz · δpz )² ] / ptot
    // px = pr · cosφ
    // py = pr · sinφ
    // pz = pl
    // δpx = √[ ( ∂px/∂pr · δpr )² + ( ∂px/∂φ · δf )² ] = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ] 
    // δpy = √[ ( ∂py/∂pr · δpr )² + ( ∂py/∂φ · δf )² ] = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    // δpz = δpl
    double output = 0.;
    double theo_phi = theo_azimuth() * TMath::DegToRad();
    double delta_theo_phi = delta_theo_azimuth() * TMath::DegToRad();
    
    double theo_px = theo_prad() * TMath::Cos( theo_phi );    // px = pr · cosφ
    double theo_py = theo_prad() * TMath::Sin( theo_phi );    // py = pr · sinφ
    double theo_pz = theo_plong();                            // pz = pl
    double theo_pt = theo_ptot();
    double deltapx = TMath::Sqrt( TMath::Power( TMath::Cos( theo_phi ) * delta_theo_prad() , 2 ) + TMath::Power( - theo_prad() * TMath::Sin( theo_phi ) * delta_theo_phi , 2 ) );   // δpx = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ]
    double deltapy = TMath::Sqrt( TMath::Power( TMath::Sin( theo_phi ) * delta_theo_prad() , 2 ) + TMath::Power(   theo_prad() * TMath::Cos( theo_phi ) * delta_theo_phi , 2 ) );   // δpy = √[ ( sinφ · δpr )² + (  pr · cosφ · δφ )² ]
    double deltapz = TMath::Sqrt( TMath::Power( delta_theo_plong() , 2 ) );                                                                                                         // δpz = δpl
            output = TMath::Sqrt( TMath::Power( theo_px * deltapx , 2 ) + TMath::Power( theo_py * deltapy , 2 ) + TMath::Power( theo_pz * deltapz , 2 ) ) / theo_pt;
    return output;
  }
  
    
  double delta_theo_elevation(){
    double output = -1000000.;
    
    double  pl   = theo_plong();
    double  pt   = theo_ptot();
    double dpl   = delta_theo_plong();
    double dpt   = delta_theo_ptot();
    double dacos = -1. / TMath::Sqrt( 1. - TMath::Power( pl/pt, 2 ) );
    
    double dth_dpl =  dacos / pt;
    double dth_dpt = -dacos / ( pt * pt );
    
    if( pt != 0. ){ output = TMath::Sqrt( TMath::Power( dth_dpl * dpl , 2 ) + TMath::Power( dth_dpt * dpt , 2 ) ); }
    return output;
  }


  /////////////////////////////////
  // CALCULATED "CALC" VARIABLES //
  /////////////////////////////////

  // Reconstructed prad [keV/c]
  double calc_prad()      { return iso_init_prad->Interpolate( tof, splat_radius ); }
  
  
// ISO_CURVE_THISTO - COMMENTA DA QUI...
  // Reconstructed prad [keV/c] - ntuple method (Lagrange polynomials) with TH3 histogram and vector

// EXPLAINATION OF THE METHOD.
// Given a 3-tuple P{tP, rP, aP} the following happens.
// The method is divided in two phases: I) search of the "cuboid" that contains P; 
// II) identification of the 8 closest points around P through which compute the 
// Lagrange polynomials; III) the actual interpolation via Lagrange pols.

// There are three structures that are accessed for the method:
// * htree: this is a TH3D that spans the 3D space described by tof, splat_radius 
// and splat_azimuth. It is used to quickly associate a "bin number" to any 3-tuple 
// via the method FindBin.
// * rpb: this is a TTree that contains a single branch, which holds a vector of 
// ints. This TTree is built in a way that it has a number of entries equal to the 
// number of bins in htree: the vector of ints in each entry contains the calibration 
// entries that would correspond to that specific htree bin, i.e. a specific range 
// of tofs, splat_radii and splat_azimuths.
// * tthist: this is a TTree that contains the calibration entries for the 
// application of the method. It's the TTree equivalent of the iso_init_* 
// TGraph2D's used for the iso_curves method. tthist contains five branches: the 
// first three are the calibration tof, splat_radius and splat_azimuth; the fourth 
// is the calibration init_prad; and the fifth is the htree bin that would 
// correspond to those calibration {t, r, a}. [Note that the fifth branch is not 
// used in the method, it's an old piece of data that could be eliminated in the 
// next iterations of the method]

  
    double vect_prad(){
      
// I. Search of the cuboid containing P
//   The search mode does the following:
//   1. it identifies the main htree bin that contains P;
//   2. it locates the additional 26 bins that surround this central bin (in other 
//   words, it identifies 27 bins centered on the one that contains P);
//   3. it creates a set, called runs, that will contain the tthist entries contained 
//   in all these bins (read from rpb->runids).
      
    // Central bin containing P
    int binstart = htree->FindBin( tof, splat_radius, splat_azimuth() );
    
    // This takes the bin value and converts it to three ints that describe bins on X, Y and Z
    int bbx, bby, bbz;
    htree->GetBinXYZ(binstart, bbx, bby, bbz);
//     cout << "binstart " << binstart << " bbx = " << bbx << " bby = " << bby << " bbz = " << bbz << endl;
    
    // Identification of the 27 bins around the central bin
    // The bin numbers are stored in the vector binids
    vector<int> *binids = new vector<int>;
    TString selection = "bi < 0 ";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
        for(int k=0; k<3; ++k){
          binids->push_back( htree->GetBin(bbx-1+i, bby-1+j, bbz-1+k) );
        }
      }
    }
//     for(int i=0;i<binids->size();++i){cout << binids->at(i) << endl;}
    
    // Identification of the tthist are already sorted
    set<int> *runs = new set<int>;
    for(int i = 0; i < binids->size(); ++i){
      rpb->GetEntry( binids->at(i) );
      for(int j=0;j<runids->size();++j){
        runs->insert( runids->at(j) );
      }
    }
    
    
//     cout << "runsize " << runs->size() << endl;
//     cout << "runs = ";
//     for(std::set<int>::iterator it = runs->begin(); it != runs->end(); ++it){ cout << *it << " "; } cout << endl;

// II. Identification of cuboid edges
//   Once the various tthist entries that surround P have been identified, it's 
//   time to locate the 8 closest points surrounding P.
//   To do this, what happens is the following:
//   1. for each selected tthist run, calculate the distance between the run and P 
//   ( via √[(trun - tP)² + (rrun - rP)² + (arun - aP)²] );
//   2. for each one of the 8 edges of the cuboid around P, check if the distance 
//   is the lowest measured: this is achieved storing the distances and the runid 
//   into a series of pair<double, int>, and with comparisons between the values of 
//   trun, rrun and arun vs tP, rP and aP;
//   3. once the required 8 tthist runids have been identified, reload them and 
//   store their 3-tuples into service variables that will be used in the last phase.

    // Definition of service variables for the completion of next phase
    double tfn(tof), srn(splat_radius), san(splat_azimuth());   // P quantities tP, rP, aP
    double tf0, tf1, sr0, sr1, sa0, sa1;
    double f000, f100, f010, f001, f110, f101, f011, f111;
    
    // pairs that will contain distance between tthist 3-tuples and P, and tthist runids
    pair<double, int> *d000 = new pair<double, int>; d000->first = 999999; d000->second = -999;
    pair<double, int> *d100 = new pair<double, int>; d100->first = 999999; d100->second = -999;
    pair<double, int> *d010 = new pair<double, int>; d010->first = 999999; d010->second = -999;
    pair<double, int> *d001 = new pair<double, int>; d001->first = 999999; d001->second = -999;
    pair<double, int> *d110 = new pair<double, int>; d110->first = 999999; d110->second = -999;
    pair<double, int> *d101 = new pair<double, int>; d101->first = 999999; d101->second = -999;
    pair<double, int> *d011 = new pair<double, int>; d011->first = 999999; d011->second = -999;
    pair<double, int> *d111 = new pair<double, int>; d111->first = 999999; d111->second = -999;
    
    // Loop on identified tthist runids
    for(std::set<int>::iterator it = runs->begin(); it != runs->end(); ++it){
      tthist->GetEntry( *it );
//       cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
      
      // Distance calculation
      double dist = TMath::Sqrt( TMath::Power(ltf-tfn,2) + TMath::Power(lsr-srn,2) + TMath::Power(lsa-san,2) );
      
      // Distance check for the 8 cuboid edges
      if( dist < d000->first && ltf <= tfn && lsr <= srn && lsa <= san ){ d000->first = dist; d000->second = *it; continue; }
      if( dist < d100->first && ltf >= tfn && lsr <= srn && lsa <= san ){ d100->first = dist; d100->second = *it; continue; }
      if( dist < d010->first && ltf <= tfn && lsr >= srn && lsa <= san ){ d010->first = dist; d010->second = *it; continue; }
      if( dist < d001->first && ltf <= tfn && lsr <= srn && lsa >= san ){ d001->first = dist; d001->second = *it; continue; }
      if( dist < d110->first && ltf >= tfn && lsr >= srn && lsa <= san ){ d110->first = dist; d110->second = *it; continue; }
      if( dist < d101->first && ltf >= tfn && lsr <= srn && lsa >= san ){ d101->first = dist; d101->second = *it; continue; }
      if( dist < d011->first && ltf <= tfn && lsr >= srn && lsa >= san ){ d011->first = dist; d011->second = *it; continue; }
      if( dist < d111->first && ltf >= tfn && lsr >= srn && lsa >= san ){ d111->first = dist; d111->second = *it; continue; }
    }
//     cout << d000->second << " ";
//     cout << d100->second << " ";
//     cout << d010->second << " ";
//     cout << d001->second << " ";
//     cout << d110->second << " ";
//     cout << d101->second << " ";
//     cout << d011->second << " ";
//     cout << d111->second << endl;
    
    // Reloading of the "correct" cuboid edges at the end of the loop, assignment of service variables
    tthist->GetEntry( d000->second ); tf0 = ltf; sr0 = lsr; sa0 = lsa; f000 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d100->second ); tf1 = ltf; sr0 = lsr; sa0 = lsa; f100 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d010->second ); tf0 = ltf; sr1 = lsr; sa0 = lsa; f010 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d001->second ); tf0 = ltf; sr0 = lsr; sa1 = lsa; f001 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d110->second ); tf1 = ltf; sr1 = lsr; sa0 = lsa; f110 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d101->second ); tf1 = ltf; sr0 = lsr; sa1 = lsa; f101 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d011->second ); tf0 = ltf; sr1 = lsr; sa1 = lsa; f011 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    tthist->GetEntry( d111->second ); tf1 = ltf; sr1 = lsr; sa1 = lsa; f111 = lpr; // cout << ltf << " " << lsr << " " << lsa << " " << lpr << endl;
    
// III. Calculation of Lagrange polynomials
//   The (1st order) basis polynomials are constructed as products of the (1st 
//   order) basis elements for the three dimensions (see 
//   https://math.stackexchange.com/questions/2099379/lagrange-polynomial-in-3-d-variable-of-interest-is-a-vector
//   for more info), and these are represented as 
//   TF3 objects. The TF3 parameters are the ones coming from phase II, so that the 
//   1st order Lagrange polynomials are correctly built. As per Lagrange method, they 
//   are then multiplied by the tthist calibration init_prad (the "value of the 
//   function at the point") and summed all together to give the final result fijk.
//   There is a built in safeguard if the result is non-physical, like if it's a 
//   NaN or infinity, and in that case the fijk is artificially set to 999999.
    
    // Definition of the Lagrange polynomials as TF3
    TF3 *l000 = new TF3("l000", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l000->SetNpx(1000); l000->SetNpy(1000); l000->SetNpz(1000);
    TF3 *l100 = new TF3("l100", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l100->SetNpx(1000); l100->SetNpy(1000); l100->SetNpz(1000);
    TF3 *l010 = new TF3("l010", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l010->SetNpx(1000); l010->SetNpy(1000); l010->SetNpz(1000);
    TF3 *l001 = new TF3("l001", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l001->SetNpx(1000); l001->SetNpy(1000); l001->SetNpz(1000);
    TF3 *l110 = new TF3("l110", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l110->SetNpx(1000); l110->SetNpy(1000); l110->SetNpz(1000);
    TF3 *l101 = new TF3("l101", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l101->SetNpx(1000); l101->SetNpy(1000); l101->SetNpz(1000);
    TF3 *l011 = new TF3("l011", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l011->SetNpx(1000); l011->SetNpy(1000); l011->SetNpz(1000);
    TF3 *l111 = new TF3("l111", " (x-[0])/([1]-[0]) * (y-[2])/([3]-[2]) * (z-[4])/([5]-[4]) ", 0, 10, 0, 100, 0, 400  ); l111->SetNpx(1000); l111->SetNpy(1000); l111->SetNpz(1000);
    
    // Assignment of the correct parameters to the TF3's, coming from phase II.
    l000->SetParameters( tf1,tf0,sr1,sr0,sa1,sa0 ); 
    l100->SetParameters( tf0,tf1,sr1,sr0,sa1,sa0 ); 
    l010->SetParameters( tf1,tf0,sr0,sr1,sa1,sa0 ); 
    l001->SetParameters( tf1,tf0,sr1,sr0,sa0,sa1 ); 
    l110->SetParameters( tf0,tf1,sr0,sr1,sa1,sa0 ); 
    l101->SetParameters( tf0,tf1,sr1,sr0,sa0,sa1 ); 
    l011->SetParameters( tf1,tf0,sr0,sr1,sa0,sa1 ); 
    l111->SetParameters( tf0,tf1,sr0,sr1,sa0,sa1 );
    
    // Calculation of the interpolated quantity
    double fijk = l000->Eval(tfn,srn,san)*f000 + l100->Eval(tfn,srn,san)*f100 + l010->Eval(tfn,srn,san)*f010 + l001->Eval(tfn,srn,san)*f001 + l110->Eval(tfn,srn,san)*f110 + l101->Eval(tfn,srn,san)*f101 + l011->Eval(tfn,srn,san)*f011 + l111->Eval(tfn,srn,san)*f111;
    
    // Cleanup
    delete binids; delete runs;
    delete d000; delete l000;
    delete d100; delete l100;
    delete d010; delete l010;
    delete d001; delete l001;
    delete d110; delete l110;
    delete d101; delete l101;
    delete d011; delete l011;
    delete d111; delete l111;
    
    // Output
    if( TMath::IsNaN(fijk) || TMath::Abs(fijk) == TMath::Infinity() ){ fijk = 999999.; }
    return fijk;
  }
// ...A QUI, VEDI DOPO


  // Reconstructed plong [keV/c]
  double calc_plong()     { return iso_init_plong->Interpolate( tof, splat_radius ); }

  // Reconstructed azimuth [deg]
  double calc_azimuth()   { double output = 0.; output = iso_init_azimuth->Interpolate( tof, splat_azimuth() ); if( output > 360. ){ output -= 360.; } return output; }

  // Reconstructed ptot [keV/c]
  double calc_ptot()      { return TMath::Sqrt( TMath::Power( calc_prad(), 2 ) + TMath::Power( calc_plong(), 2 ) ); }
  double calc_ptot_cart()      {
    double pl = calc_plong();
    double pr = calc_prad();
    double az = calc_azimuth() * TMath::DegToRad();
    
    double px = pl;
    double py = pr * TMath::Cos( az );
    double pz = pr * TMath::Sin( az );
    return TMath::Sqrt( TMath::Power( px, 2 ) + TMath::Power( py, 2 ) + TMath::Power( pz, 2 ) );
  }

  // Reconstructed elevation [deg]
  double calc_elevation()   {
    double output = -1000000.;
    double pr = calc_prad();
    double pl = calc_plong();
    double pt = TMath::Sqrt( pr*pr + pl*pl );
    if( pt != 0. ){ output = TMath::ACos( pl / pt ) * TMath::RadToDeg(); }
    return output;
  }

  // Reconstructed energy [keV]
  double calc_energy()    { return 1e3 * ( TMath::Sqrt( TMath::Power( auger_mass, 2 ) + TMath::Power( calc_ptot(), 2 ) ) - auger_mass ); }

  // Reconstructed azimuth ( phi ) [deg] [Ullrich formula]
  double calc_phi()       {
    double val = ( splat_azimuth()*TMath::DegToRad() + ( omega*tof - TMath::Floor( omega*tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2. ) * TMath::RadToDeg();
    if( val < 0. ){ val +=360.; }
    return val;
  }
    
  // Calculated initial value of py and pz, from calculated azimuth
  double calc_pz()  { return  calc_plong(); }
  double calc_px()  { return ( calc_prad() * TMath::Cos( calc_azimuth() * TMath::DegToRad() ) ); }
  double calc_py()  { return ( calc_prad() * TMath::Sin( calc_azimuth() * TMath::DegToRad() ) ); }

  ///////////////////////////
  // UNCERTAINTIES "DELTA" //
  ///////////////////////////

  // Uncertainty on prad [keV/c]
  double delta_calc_prad()        { // ERROR PROPAGATION METHOD: δprad = √[ ( dprad/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double tofepsilon = 1e-6;
    double srepsilon  = 1e-3;
    double dprad_dtof  = ( iso_init_prad->Interpolate( tof+( tofepsilon/2. ), splat_radius ) - iso_init_prad->Interpolate( tof-( tofepsilon/2. ), splat_radius ) ) / tofepsilon;
    double dprad_dsr   = ( iso_init_prad->Interpolate( tof, splat_radius+( srepsilon/2. ) ) - iso_init_prad->Interpolate( tof, splat_radius-( srepsilon/2. ) ) ) / srepsilon;

    output = TMath::Sqrt( TMath::Power( dprad_dtof * DEtof_MCP, 2 ) + TMath::Power( dprad_dsr * DEsplat_radius_MCP, 2 ) );

    return output;
  }
  
  // Uncertainty on plong [keV/c]
  double delta_calc_plong()        { // ERROR PROPAGATION METHOD: δplong = √[ ( dplong/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double tofepsilon = 1e-6;
    double srepsilon  = 1e-3;
    double dplong_dtof  = ( iso_init_plong->Interpolate( tof+( tofepsilon/2. ), splat_radius ) - iso_init_plong->Interpolate( tof-( tofepsilon/2. ), splat_radius ) ) / tofepsilon;
    double dplong_dsr   = ( iso_init_plong->Interpolate( tof, splat_radius+( srepsilon/2. ) ) - iso_init_plong->Interpolate( tof, splat_radius-( srepsilon/2. ) ) ) / srepsilon;

    output = TMath::Sqrt( TMath::Power( dplong_dtof * DEtof_MCP, 2 ) + TMath::Power( dplong_dsr * DEsplat_radius_MCP, 2 ) );

    return output;
  }

  // Uncertainty on phi [deg]
  double delta_calc_phi()          { return ( TMath::Sqrt( TMath::Power( DEsplat_radius_MCP/splat_radius, 2 ) + TMath::Power( 0.5*omega*DEtof_MCP, 2 ) ) ) * TMath::RadToDeg(); }

  // Uncertainty on splat_azimuth [deg]
  double delta_splat_azimuth(){ return ( TMath::Sqrt( TMath::Power( splat_x*DEsplat_radius_MCP, 2 ) + TMath::Power( splat_y*DEsplat_radius_MCP, 2 ) ) / ( splat_radius*splat_radius ) ) * TMath::RadToDeg(); }

  // Uncertainty on calc_azimuth [keV/c]
  double delta_calc_azimuth()        { // ERROR PROPAGATION METHOD: δazimuth = √[ ( dazimuth/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double epsilon = 1e-6;
    double dazimuth_dtof  = ( iso_init_azimuth->Interpolate( tof+( epsilon/2. ), splat_azimuth() ) - iso_init_azimuth->Interpolate( tof-( epsilon/2. ), splat_azimuth() ) ) / epsilon;
    double dazimuth_dsr   = ( iso_init_azimuth->Interpolate( tof, splat_azimuth()+( epsilon/2. ) ) - iso_init_azimuth->Interpolate( tof, splat_azimuth()-( epsilon/2. ) ) ) / epsilon;

    double DEsplat_azimuth_MCP = delta_splat_azimuth();
    output = TMath::Sqrt( TMath::Power( dazimuth_dtof * DEtof_MCP, 2 ) + TMath::Power( dazimuth_dsr * DEsplat_azimuth_MCP, 2 ) );

    return output;
  }
  
  double delta_calc_ptot(){
    // δptot = √[ ( ∂ptot/∂px · δpx )² + ( ∂ptot/∂py · δpy )² + ( ∂ptot/∂pz · δpz )² ] = √[ ( px · δpx )² + ( py · δpy )² + ( pz · δpz )² ] / ptot
    // px = pr · cosφ
    // py = pr · sinφ
    // pz = pl
    // δpx = √[ ( ∂px/∂pr · δpr )² + ( ∂px/∂φ · δf )² ] = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ] 
    // δpy = √[ ( ∂py/∂pr · δpr )² + ( ∂py/∂φ · δf )² ] = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    // δpz = δpl
    double output = 0.;
    double  phi =       calc_azimuth() * TMath::DegToRad();
    double dphi = delta_calc_azimuth() * TMath::DegToRad();
    
    double  px = calc_prad() * TMath::Cos( phi );    // px = pr · cosφ
    double  py = calc_prad() * TMath::Sin( phi );    // py = pr · sinφ
    double  pz = calc_plong();                            // pz = pl
    double  pt = calc_ptot();
    double dpx = TMath::Sqrt( TMath::Power( TMath::Cos( phi ) * delta_calc_prad() , 2 ) + TMath::Power( - calc_prad() * TMath::Sin( phi ) * dphi , 2 ) );   // δpx = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ]
    double dpy = TMath::Sqrt( TMath::Power( TMath::Sin( phi ) * delta_calc_prad() , 2 ) + TMath::Power(   calc_prad() * TMath::Cos( phi ) * dphi , 2 ) );   // δpy = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    double dpz = TMath::Sqrt( TMath::Power( delta_calc_plong() , 2 ) );                                                                                     // δpz = δpl
    output = TMath::Sqrt( TMath::Power( px * dpx , 2 ) + TMath::Power( py * dpy , 2 ) + TMath::Power( pz * dpz , 2 ) ) / pt;
    return output;
  }
  
  double delta_calc_elevation(){
    double output = -1000000.;
    
    double  pl   = calc_plong();
    double  pt   = calc_ptot();
    double dpl   = delta_calc_plong();
    double dpt   = delta_calc_ptot();
    double dacos = -1. / TMath::Sqrt( 1. - TMath::Power( pl/pt, 2 ) );
    
    double dth_dpl =  dacos / pt;
    double dth_dpt = -dacos / ( pt * pt );
    
    if( pt != 0. ){ output = TMath::Sqrt( TMath::Power( dth_dpl * dpl , 2 ) + TMath::Power( dth_dpt * dpt , 2 ) ); }
    return output;
  }

  
  // Uncertainty on the reconstructed energy [keV]
  double delta_calc_energy()  { return 1e3 * ( ( calc_ptot() * delta_calc_ptot() ) / ( TMath::Sqrt( TMath::Power( auger_mass, 2 ) + TMath::Power( calc_ptot(), 2 ) ) ) ); }
  
  
  //////////////////////////////////////////
  // QUANTITIES WITH ERRORS "CALC_*_WERR" //
  //////////////////////////////////////////

  double generate_tof_werr()           { return ra->Gaus( tof,           DEtof_MCP ); }
  double generate_splat_x_werr()       { return ra->Gaus( splat_x,       DEsplat_radius_MCP ); }
  double generate_splat_y_werr()       { return ra->Gaus( splat_y,       DEsplat_radius_MCP ); }
//   double generate_splat_radius_werr()  { return ra->Gaus( splat_radius,  DEsplat_radius_MCP ); }
  double generate_splat_radius_werr()  { return TMath::Sqrt( splat_x_werr*splat_x_werr + splat_y_werr*splat_y_werr ); }

  double tof_werr;
  double splat_x_werr;
  double splat_y_werr;
  double splat_radius_werr;

  // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
  double splat_azimuth_werr()      {
    double val = TMath::ATan2( splat_y_werr, splat_x_werr ) * TMath::RadToDeg();
    if( val < 0. ){ val += 360.; }  // comment for -pi to +pi
    return val;

  }

  // Reconstructed plong [keV/c]
  double calc_plong_werr() {
    return iso_init_plong->Interpolate( tof_werr, splat_radius_werr );
  }

  // Reconstructed plong [keV/c]
  double calc_prad_werr() {
    return iso_init_prad->Interpolate( tof_werr, splat_radius_werr );
  }

  // Reconstructed ptot [keV/c]
  double calc_ptot_werr()      { return TMath::Sqrt( TMath::Power( calc_prad_werr(), 2 ) + TMath::Power( calc_plong_werr(), 2 ) ); }

  // Reconstructed azimuth [deg]
  double calc_azimuth_werr()   {
    return iso_init_azimuth->Interpolate( ( omega*( tof_werr ) - TMath::Floor( omega*( tof_werr )/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2., splat_azimuth_werr() );
  }
  
  // Calculated initial value of px, py and pz, from calculated azimuth
  double calc_pz_werr()  { return  calc_plong_werr(); }
  double calc_px_werr()  { return ( calc_prad_werr() * TMath::Cos( calc_azimuth_werr() * TMath::DegToRad() ) ); }
  double calc_py_werr()  { return ( calc_prad_werr() * TMath::Sin( calc_azimuth_werr() * TMath::DegToRad() ) ); }

  // Uncertainty on plong with errors [keV/c]
  double delta_plong_werr()        {
    vector<double> *pr = new vector<double>;
    
    for( double i=-1.0;i<1.1;i+=0.1 ){
      for( double j=-1.0;j<1.1;j+=0.1 ){
        pr->push_back( iso_init_plong->Interpolate( tof_werr+i*DEtof_MCP, splat_radius_werr+j*DEsplat_radius_MCP ) );
      }
    }
    
    return ( TMath::RMS( pr->begin(), pr->end() ) );
  }
  
  // Uncertainty on prad with errors [keV/c]
  double delta_prad_werr()        {
    vector<double> *pr = new vector<double>;
    double sinomegat = TMath::Sin( omega*tof_werr/2. );
    double deltat = omega/2. * TMath::Cos( omega*tof_werr/2. );
    
    for( double i=-1.0;i<1.1;i+=0.1 ){
      for( double j=-1.0;j<1.1;j+=0.1 ){
        pr->push_back( iso_init_prad->Interpolate( tof_werr+i*DEtof_MCP, splat_radius_werr+j*DEsplat_radius_MCP ) );
      }
    }
    
    return ( TMath::RMS( pr->begin(), pr->end() ) );
  }
  
  // Uncertainty on ptot with errors [keV/c]
  double delta_ptot_werr()         { return ( TMath::Sqrt( TMath::Power( calc_prad_werr()*delta_prad_werr(), 2 ) + TMath::Power( calc_plong_werr()*delta_plong_werr(), 2 ) ) / calc_ptot_werr() ); }
  
  // Uncertainty on phi with errors [deg]
  double delta_phi_werr()          { return ( TMath::Sqrt( TMath::Power( DEsplat_radius_MCP/splat_radius_werr, 2 ) + TMath::Power( 0.5*omega*DEtof_MCP, 2 ) ) ) * TMath::RadToDeg(); }
  
  // Uncertainty on splat_azimuth with errors [deg]
  double delta_splat_azimuth_werr(){ return ( TMath::Sqrt( TMath::Power( splat_x_werr*DEsplat_radius_MCP, 2 ) + TMath::Power( splat_y_werr*DEsplat_radius_MCP, 2 ) ) / ( splat_radius_werr*splat_radius_werr ) ) * TMath::RadToDeg(); }
  
  // Uncertainty on calc_azimuth with errors [keV/c]
  double delta_azimuth_werr()        {
    vector<double> *pr = new vector<double>;
    double omegat = ( omega*tof_werr - TMath::Floor( omega*tof_werr/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2.;
    
    for( double i=-1.0;i<1.1;i+=0.1 ){
      for( double j=-1.0;j<1.1;j+=0.1 ){
        pr->push_back( iso_init_azimuth->Interpolate( omegat+omega*i*DEtof_MCP/2, splat_azimuth_werr()+j*delta_splat_azimuth_werr() ) );
      }
    }
    
    return ( TMath::RMS( pr->begin(), pr->end() ) );
  }
  

  //////////////////
  // CUT CHECKERS //
  //////////////////

  bool cut_checker( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
  {
    bool VIOLATES_ONMCP = false;  // Auger must hit MCP
    bool VIOLATES_SPRAD = false;  // Auger splat_radius must be CX_SR_LO( 3 ) < splat_radius < CX_SR_HI( 60 )
    bool VIOLATES_TOF   = false;  // Auger tof must be < CX_TOF
//     bool VIOLATES_TOFLO = false;  // Low-energy Auger tof must be < CX_TOF_LO
//     bool VIOLATES_TOFPP = false;  // Low-energy Auger with 0.18 < tof< CX_TOF_HI must have splat_radius > pp3( tof ) [this cut preserves backward Augers removing hi-energy e's which could not be reconstructed]
//     bool VIOLATES_TOFHI = false;  // Hi-energy Auger tof must be < CX_TOF_HI
    bool VIOLATES_NODEF = false;  // Auger must not arrive in a node ( this is a tof check )
    bool VIOLATES_CALCP = false;  // Auger's prad and plong must be, ultimately, reconstructable, i.e. != 0.
    double nodelo( 0. ), nodehi( 0. );
    if( onmcp ){ if( !( TMath::Abs( splat_z - ele_MCP_position ) < 0.01 ) )                                                                                                          { VIOLATES_ONMCP = true; } }
    if( sprad ){ if( !( splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH ) )                                                                     { VIOLATES_SPRAD = true; } }
    if( timef ){
      if( !( tof < CX_TOF ) )                                                                                                                     { VIOLATES_TOF = true; }   // There is no more distinction between TOF_LO e TOF_HI
//       if( init_ke < 90. && !( tof < CX_TOF_LO ) )                                                                                              { VIOLATES_TOFLO = true; }   // Remember: we capture ALL the xenons meaning that
//       if( init_ke < 90. && tof > 0.18 && tof < CX_TOF_LO && !( splat_radius > pp3->Eval( tof ) ) )                                               { VIOLATES_TOFPP = true; }   // we are able to determine if we have 2 Augers
//       if( init_ke > 90. && !( tof < CX_TOF_HI ) )                                                                                              { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 ) TEMPORANEAMENTE DISATTIVATO
//       if( init_ke > 90. && !( ( tof < CX_TOF_HI ) || ( tof > CX_TOF_HI_LOW && tof < CX_TOF_HI_HIGH ) ) )                                           { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 )
   }
    if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE ) { VIOLATES_NODEF = true; nodelo = nodes[nnodes]-CX_NODE; nodehi = nodes[nnodes]+CX_NODE; } } }
    if( calcp ){ if( calc_prad() == 0. || calc_plong() == 0. )                                                                                        { VIOLATES_CALCP = true; } }
//     if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOFLO || VIOLATES_TOFPP || VIOLATES_TOFHI || VIOLATES_NODEF || VIOLATES_CALCP ){
    if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOF || VIOLATES_NODEF || VIOLATES_CALCP ){
      if( VIOLATES_ONMCP )  { cout << "Auger not on MCP plane" << endl;                                                                                     }
      if( VIOLATES_SPRAD )  { cout << "Auger splat_radius = " << splat_radius << endl;                                                                      }
      if( VIOLATES_TOF )  { cout << "Auger TOF = " << tof << " [> CX_TOF = " << CX_TOF << "]" << endl;                                                    }
//       if( VIOLATES_TOFLO )  { cout << "Auger TOF = " << tof << " [> CX_TOF_LO = " << CX_TOF_LO << "]" << endl;                                              }
//       if( VIOLATES_TOFPP )  { cout << "Backward Auger, TOF = " << tof << " ; SR = " << splat_radius << " < pp3[" << tof << "] = " << pp3->Eval( tof ) << endl; }
//       if( VIOLATES_TOFHI )  { cout << "Auger TOF = " << tof << " [> CX_TOF_HI = " << CX_TOF_HI << "]" << endl;                                              }
      if( VIOLATES_NODEF )  { cout << "Auger TOF = " << tof << " cut by CX_NODE: " << nodelo << " < " << tof << " < " << nodehi << endl;                    }
      if( VIOLATES_CALCP )  { cout << "Auger momentum not reconstructed: px = " << calc_px() << " ; py = " << calc_py() << " ; pz = " << calc_pz() << endl; }
      return false;
   }
    return true;
  }

  bool silent_cut_checker( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
  {
    if( onmcp ){ if( !( TMath::Abs( splat_z - ele_MCP_position ) < 0.01 ) ){ return false; } }
    if( sprad ){ if( !( splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH ) ){ return false; } }
    if( timef ){
      if( !( tof < CX_TOF ) )                                                         { return false; }
//       if( init_ke < 90. && !( tof < CX_TOF_LO ) )                                                         { return false; }
//       else if( init_ke < 90. && tof > 0.18 && tof < CX_TOF_LO && !( splat_radius > pp3->Eval( tof ) ) )     { return false; }
//       else if( init_ke > 90. && !( tof < CX_TOF_HI ) )                                                { return false; }
//       else if( init_ke > 90. && !( ( tof < CX_TOF_HI ) || ( tof > CX_TOF_HI_LOW && tof < CX_TOF_HI_HIGH ) ) ) { return false; }
   }
    if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE ){ return false; } } }
    if( calcp ){ if( calc_prad() == 0. || calc_plong() == 0. ){ return false; } }
    return true;
  }


















  //////////////////////////////////////////////////
  //    __  __  _____ _____   __   ____     __    //
  //   |  \/  |/ ____|  __ \  \ \ / /\ \   / /    //
  //   | \  / | |    | |__) |  \ V /  \ \_/ /     //
  //   | |\/| | |    |  ___/    > <    \   /      //
  //   | |  | | |____| |       / . \    | |       //
  //   |_|  |_|\_____|_|      /_/ \_\   |_|       //
  //                                              //
  //////////////////////////////////////////////////

  /////////////////////////////////
  // MCP COMPENSATED COORDINATES //
  /////////////////////////////////

  // MCP splat_xy
  double mcp_splat_x()   { return ( round( splat_x * ( 1./DEsplat_radius_MCP ) ) ) / ( 1./DEsplat_radius_MCP ); }
  double mcp_splat_y()   { return ( round( splat_y * ( 1./DEsplat_radius_MCP ) ) ) / ( 1./DEsplat_radius_MCP ); }
  
  /////////////////////////////////
  // ADDITIONAL "WORK" VARIABLES //
  /////////////////////////////////

  // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
  double mcp_splat_azimuth()      {
    double val = TMath::ATan2( mcp_splat_y(), mcp_splat_x() ) * TMath::RadToDeg();
    if( val < 0. ){ val += 360.; }  // comment for -pi to +pi
    return val;
  }
    
  // MCP splat radius: the MCP has holes, this rounds to the MCP resolution
  double mcp_splat_radius()   { double mcpx = mcp_splat_x(); double mcpy = mcp_splat_y(); return TMath::Sqrt( mcpx*mcpx + mcpy*mcpy ); }
  
  
  //////////////////////////////////////
  // THEORETICAL (COLTRIMS) VARIABLES //
  //////////////////////////////////////

  // Theoretical (COLTRIMS) reconstructed variables
  double mcp_theo_plong()     { return ( ( ( mass*( splat_z-init_z ) )/( tof*c*1e-3 ) ) - Efield/2.*tof*c*1e-6 );      }
  double mcp_theo_prad()      { return ( Bfield*mcp_splat_radius() )/( 2.*33.356*TMath::Abs( TMath::Sin( omega*tof/2. ) ) ); }
  double mcp_theo_ptot()      { return TMath::Sqrt( TMath::Power( mcp_theo_prad(), 2 ) + TMath::Power( mcp_theo_plong(), 2 ) );  }
  double mcp_theo_azimuth()   { double output = 0.; output = ( mcp_splat_azimuth()*TMath::DegToRad() + ( omega*tof - TMath::Floor( omega*tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() ) /2. ) *TMath::RadToDeg(); if( output > 360. ){ output -= 360.; } return output; }
  double mcp_theo_elevation() {
    double output = -1000000.;
    if( mcp_theo_ptot() != 0. ){ output = TMath::ACos( mcp_theo_plong() / mcp_theo_ptot() ) * TMath::RadToDeg(); }
    return output;
  }
  
  // Theoretical (COLTRIMS) errors  
  double delta_mcp_theo_plong(){
    double output = 0.;
    double Lelectron = splat_z-init_z;
    double dplong_dt_deltat = ( - ( mass * Lelectron ) / ( tof*tof*c*1e-3 ) - Efield/2.*c*1e-6 ) * DEtof_MCP ;
    double dplong_dL_deltaL = mass / tof * DELelectron ;
    double dplong_dE_deltaE = -tof / 2. * DEEfield ;
    output = TMath::Sqrt( TMath::Power( dplong_dt_deltat, 2 ) + TMath::Power( dplong_dL_deltaL, 2 ) + TMath::Power( dplong_dE_deltaE, 2 ) );
    return output;
  }
  
  double delta_mcp_theo_prad(){
    double output = 0.;
    double wt2 = omega*tof/2.;
    double prad_r = ( Bfield )/( 2.*33.356*TMath::Abs( TMath::Sin( omega*tof/2. ) ) ) ;
    double dprad_dr_deltar = prad_r * DEsplat_radius_MCP ;
    double dprad_dt_deltat = prad_r * ( omega*mcp_splat_radius()/2. ) * ( DEtof_MCP / TMath::Tan( wt2 ) ) ;
    double dprad_dB_deltaB = prad_r * ( ( mcp_splat_radius()/Bfield ) - ( ( wt2/Bfield ) / TMath::Tan( wt2 ) ) ) * DEBfield ;
    output = TMath::Sqrt( TMath::Power( dprad_dr_deltar, 2 ) + TMath::Power( dprad_dt_deltat, 2 ) + TMath::Power( dprad_dB_deltaB, 2 ) );
    return output;
  }
  
  double delta_mcp_theo_azimuth(){
    double output = 0.;
    double dphi_dtheta_deltatheta = DEsplat_radius_MCP / mcp_splat_radius() ;
    double dphi_dt_deltat         = -omega / 2. * DEtof_MCP ; 
    double dphi_dB_deltaB         = -( omega * tof ) / Bfield * DEBfield ; 
    output = TMath::Sqrt( TMath::Power( dphi_dtheta_deltatheta, 2 ) + TMath::Power( dphi_dt_deltat, 2 ) + TMath::Power( dphi_dB_deltaB, 2 ) ) * TMath::RadToDeg();
    return output;
  }
  
  double delta_mcp_theo_ptot(){
    // δptot = √[ ( ∂ptot/∂px · δpx )² + ( ∂ptot/∂py · δpy )² + ( ∂ptot/∂pz · δpz )² ] = √[ ( px · δpx )² + ( py · δpy )² + ( pz · δpz )² ] / ptot
    // px = pr · cosφ
    // py = pr · sinφ
    // pz = pl
    // δpx = √[ ( ∂px/∂pr · δpr )² + ( ∂px/∂φ · δf )² ] = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ] 
    // δpy = √[ ( ∂py/∂pr · δpr )² + ( ∂py/∂φ · δf )² ] = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    // δpz = δpl
    double output = 0.;
    double theo_phi = mcp_theo_azimuth() * TMath::DegToRad();
    double delta_theo_phi = delta_theo_azimuth() * TMath::DegToRad();
    
    double theo_px = mcp_theo_prad() * TMath::Cos( theo_phi );    // px = pr · cosφ
    double theo_py = mcp_theo_prad() * TMath::Sin( theo_phi );    // py = pr · sinφ
    double theo_pz = mcp_theo_plong();                            // pz = pl
    double theo_pt = mcp_theo_ptot();
    double deltapx = TMath::Sqrt( TMath::Power( TMath::Cos( theo_phi ) * delta_mcp_theo_prad() , 2 ) + TMath::Power( - mcp_theo_prad() * TMath::Sin( theo_phi ) * delta_theo_phi , 2 ) );   // δpx = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ]
    double deltapy = TMath::Sqrt( TMath::Power( TMath::Sin( theo_phi ) * delta_mcp_theo_prad() , 2 ) + TMath::Power(   mcp_theo_prad() * TMath::Cos( theo_phi ) * delta_theo_phi , 2 ) );   // δpy = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    double deltapz = TMath::Sqrt( TMath::Power( delta_mcp_theo_plong() , 2 ) );                                                                                                             // δpz = δpl
            output = TMath::Sqrt( TMath::Power( theo_px * deltapx , 2 ) + TMath::Power( theo_py * deltapy , 2 ) + TMath::Power( theo_pz * deltapz , 2 ) ) / theo_pt;
    return output;
  }
  
  double mcp_delta_theo_elevation(){
    double output = -1000000.;
    
    double  pl   = mcp_theo_plong();
    double  pt   = mcp_theo_ptot();
    double dpl   = delta_mcp_theo_plong();
    double dpt   = delta_mcp_theo_ptot();
    double dacos = -1. / TMath::Sqrt( 1. - TMath::Power( pl/pt, 2 ) );
    
    double dth_dpl =  dacos / pt;
    double dth_dpt = -dacos / ( pt * pt );
    
    if( pt != 0. ){ output = TMath::Sqrt( TMath::Power( dth_dpl * dpl , 2 ) + TMath::Power( dth_dpt * dpt , 2 ) ); }
    return output;
  }



  /////////////////////////////////
  // CALCULATED "CALC" VARIABLES //
  /////////////////////////////////

  // Reconstructed plong [keV/c]
  double mcp_calc_prad()      { return iso_init_prad->Interpolate( tof, mcp_splat_radius() ); }

  // Reconstructed plong [keV/c]
  double mcp_calc_plong()     { return iso_init_plong->Interpolate( tof, mcp_splat_radius() ); }

  // Reconstructed azimuth [deg]
  double mcp_calc_azimuth()   { double output = 0.; output = iso_init_azimuth->Interpolate( tof, mcp_splat_azimuth() ); if( output > 360. ){ output -= 360.; } return output; }

  // Reconstructed ptot [keV/c]
  double mcp_calc_ptot()      { return TMath::Sqrt( TMath::Power( mcp_calc_prad(), 2 ) + TMath::Power( mcp_calc_plong(), 2 ) ); }

  // Reconstructed elevation [deg]
  double mcp_calc_elevation()   {
    double output = -1000000.;
    if( mcp_calc_ptot() != 0. ){ output = TMath::ACos( mcp_calc_plong() / mcp_calc_ptot() ) * TMath::RadToDeg(); }
    return output;
  }

  // Reconstructed energy [keV]
  double mcp_calc_energy()    { return 1e3 * ( TMath::Sqrt( TMath::Power( auger_mass, 2 ) + TMath::Power( mcp_calc_ptot(), 2 ) ) - auger_mass ); }

  // Reconstructed azimuth ( phi ) [deg] [Ullrich formula]
  double mcp_calc_phi()       {
    double val = ( mcp_splat_azimuth()*TMath::DegToRad() + ( omega*tof - TMath::Floor( omega*tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2. ) * TMath::RadToDeg();
    if( val < 0. ){ val +=360.; }
    return val;
  }
    
  // Calculated initial value of py and pz, from calculated azimuth
  double mcp_calc_pz()  { return   mcp_calc_plong(); }
  double mcp_calc_px()  { return ( mcp_calc_prad() * TMath::Cos( mcp_calc_azimuth() * TMath::DegToRad() ) ); }
  double mcp_calc_py()  { return ( mcp_calc_prad() * TMath::Sin( mcp_calc_azimuth() * TMath::DegToRad() ) ); }

  ///////////////////////////
  // UNCERTAINTIES "DELTA" //
  ///////////////////////////

  // Uncertainty on prad [keV/c]
  double delta_mcp_calc_prad()        { // ERROR PROPAGATION METHOD: δprad = √[ ( dprad/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double tofepsilon = 1e-6;
    double srepsilon  = 1e-3;
    double dprad_dtof  = ( iso_init_prad->Interpolate( tof+( tofepsilon/2. ), mcp_splat_radius() ) - iso_init_prad->Interpolate( tof-( tofepsilon/2. ), mcp_splat_radius() ) ) / tofepsilon;
    double dprad_dsr   = ( iso_init_prad->Interpolate( tof, mcp_splat_radius()+( srepsilon/2.  ) ) - iso_init_prad->Interpolate( tof, mcp_splat_radius()-( srepsilon/2.  ) ) ) / srepsilon;

    output = TMath::Sqrt( TMath::Power( dprad_dtof * DEtof_MCP, 2 ) + TMath::Power( dprad_dsr * DEsplat_radius_MCP, 2 ) );

    return output;
  }
  
  // Uncertainty on plong [keV/c]
  double delta_mcp_calc_plong()        { // ERROR PROPAGATION METHOD: δplong = √[ ( dplong/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double tofepsilon = 1e-6;
    double srepsilon  = 1e-3;
    double dplong_dtof  = ( iso_init_plong->Interpolate( tof+( tofepsilon/2. ), mcp_splat_radius() ) - iso_init_plong->Interpolate( tof-( tofepsilon/2. ), mcp_splat_radius() ) ) / tofepsilon;
    double dplong_dsr   = ( iso_init_plong->Interpolate( tof, mcp_splat_radius()+( srepsilon/2.  ) ) - iso_init_plong->Interpolate( tof, mcp_splat_radius()-( srepsilon/2.  ) ) ) / srepsilon;

    output = TMath::Sqrt( TMath::Power( dplong_dtof * DEtof_MCP, 2 ) + TMath::Power( dplong_dsr * DEsplat_radius_MCP, 2 ) );

    return output;
  }

  // Uncertainty on phi [deg]
  double delta_mcp_calc_phi()          { return ( TMath::Sqrt( TMath::Power( DEsplat_radius_MCP/mcp_splat_radius(), 2 ) + TMath::Power( 0.5*omega*DEtof_MCP, 2 ) ) ) * TMath::RadToDeg(); }

  // Uncertainty on splat_azimuth [deg]
  double delta_mcp_splat_azimuth()         { return ( TMath::Sqrt( TMath::Power( mcp_splat_x()*DEsplat_radius_MCP, 2 ) + TMath::Power( mcp_splat_y()*DEsplat_radius_MCP, 2 ) ) / ( mcp_splat_radius()*mcp_splat_radius() ) ) * TMath::RadToDeg(); }

  // Uncertainty on calc_azimuth [keV/c]
  double delta_mcp_calc_azimuth()        { // ERROR PROPAGATION METHOD: δazimuth = √[ ( dazimuth/dtof · δtof )² + ( dplong/dtof · δtof )² ] - USE THIS, HAS NO COVARIANCE ELEMENT
    double output = 0.;

    double epsilon = 1e-6;
    double dazimuth_dtof  = ( iso_init_azimuth->Interpolate( tof+( epsilon/2. ), mcp_splat_azimuth() ) - iso_init_azimuth->Interpolate( tof-( epsilon/2. ), mcp_splat_azimuth() ) ) / epsilon;
    double dazimuth_dsr   = ( iso_init_azimuth->Interpolate( tof, mcp_splat_azimuth()+( epsilon/2. ) ) - iso_init_azimuth->Interpolate( tof, mcp_splat_azimuth()-( epsilon/2. ) ) ) / epsilon;

    double DEsplat_azimuth_MCP = delta_mcp_splat_azimuth();
    output = TMath::Sqrt( TMath::Power( dazimuth_dtof * DEtof_MCP, 2 ) + TMath::Power( dazimuth_dsr * DEsplat_azimuth_MCP, 2 ) );

    return output;
  }

  
  double delta_mcp_calc_ptot(){
    // δptot = √[ ( ∂ptot/∂px · δpx )² + ( ∂ptot/∂py · δpy )² + ( ∂ptot/∂pz · δpz )² ] = √[ ( px · δpx )² + ( py · δpy )² + ( pz · δpz )² ] / ptot
    // px = pr · cosφ
    // py = pr · sinφ
    // pz = pl
    // δpx = √[ ( ∂px/∂pr · δpr )² + ( ∂px/∂φ · δf )² ] = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ] 
    // δpy = √[ ( ∂py/∂pr · δpr )² + ( ∂py/∂φ · δf )² ] = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    // δpz = δpl
    double output = 0.;
    double local_phi      =       mcp_calc_azimuth() * TMath::DegToRad();
    double local_deltaphi = delta_mcp_calc_azimuth() * TMath::DegToRad();
    
    double local_px = mcp_calc_prad() * TMath::Cos( local_phi );    // px = pr · cosφ
    double local_py = mcp_calc_prad() * TMath::Sin( local_phi );    // py = pr · sinφ
    double local_pz = mcp_calc_plong();                            // pz = pl
    double local_pt = mcp_calc_ptot();
    double deltapx = TMath::Sqrt( TMath::Power( TMath::Cos( local_phi ) * delta_mcp_calc_prad() , 2 ) + TMath::Power( - mcp_calc_prad() * TMath::Sin( local_phi ) * local_deltaphi , 2 ) );   // δpx = √[ ( cosφ · δpr )² + ( -pr · sinφ · δφ )² ]
    double deltapy = TMath::Sqrt( TMath::Power( TMath::Sin( local_phi ) * delta_mcp_calc_prad() , 2 ) + TMath::Power(   mcp_calc_prad() * TMath::Cos( local_phi ) * local_deltaphi , 2 ) );   // δpy = √[ ( sinφ · δpr )² + ( pr · cosφ · δφ )² ]
    double deltapz = TMath::Sqrt( TMath::Power( delta_mcp_calc_plong() , 2 ) );                                                                                                         // δpz = δpl
    output = TMath::Sqrt( TMath::Power( local_px * deltapx , 2 ) + TMath::Power( local_py * deltapy , 2 ) + TMath::Power( local_pz * deltapz , 2 ) ) / local_pt;
    return output;
  }
  
  double delta_mcp_calc_elevation(){
    double output = -1000000.;
    
    double  pl   = mcp_calc_plong();
    double  pt   = mcp_calc_ptot();
    double dpl   = delta_mcp_calc_plong();
    double dpt   = delta_mcp_calc_ptot();
    double dacos = -1. / TMath::Sqrt( 1. - TMath::Power( pl/pt, 2 ) );
    
    double dth_dpl =  dacos / pt;
    double dth_dpt = -dacos / ( pt * pt );
    
    if( pt != 0. ){ output = TMath::Sqrt( TMath::Power( dth_dpl * dpl , 2 ) + TMath::Power( dth_dpt * dpt , 2 ) ); }
    return output;
  }

  
  // Uncertainty on the reconstructed energy [keV]
  double delta_mcp_calc_energy()  { return 1e3 * ( ( mcp_calc_ptot() * delta_mcp_calc_ptot() ) / ( TMath::Sqrt( TMath::Power( auger_mass, 2 ) + TMath::Power( mcp_calc_ptot(), 2 ) ) ) ); }
  
  
//   //////////////////////////////////////////
//   // QUANTITIES WITH ERRORS "CALC_*_WERR" //
//   //////////////////////////////////////////
// 
//   double generate_tof_werr()           { return ra->Gaus( tof,           DEtof_MCP ); }
//   double generate_mcp_splat_x()_werr()       { return ra->Gaus( mcp_splat_x(),       DEsplat_radius_MCP ); }
//   double generate_mcp_splat_y()_werr()       { return ra->Gaus( mcp_splat_y(),       DEsplat_radius_MCP ); }
// //   double generate_mcp_splat_radius()_werr()  { return ra->Gaus( mcp_splat_radius(),  DEsplat_radius_MCP ); }
//   double generate_mcp_splat_radius()_werr()  { return TMath::Sqrt( mcp_splat_x()_werr*mcp_splat_x()_werr + mcp_splat_y()_werr*mcp_splat_y()_werr ); }
// 
//   double tof_werr;
//   double mcp_splat_x()_werr;
//   double mcp_splat_y()_werr;
//   double mcp_splat_radius()_werr;
// 
//   // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
//   double splat_azimuth_werr()      {
//     double val = TMath::ATan2( mcp_splat_y()_werr, mcp_splat_x()_werr ) * TMath::RadToDeg();
//     if( val < 0. ){ val += 360.; }  // comment for -pi to +pi
//     return val;
// 
//   }
// 
//   // Reconstructed plong [keV/c]
//   double calc_plong_werr() {
//     return iso_init_plong->Interpolate( tof_werr, mcp_splat_radius()_werr );
//   }
// 
//   // Reconstructed plong [keV/c]
//   double calc_prad_werr() {
//     return iso_init_prad->Interpolate( tof_werr, mcp_splat_radius()_werr );
//   }
// 
//   // Reconstructed ptot [keV/c]
//   double calc_ptot_werr()      { return TMath::Sqrt( TMath::Power( calc_prad_werr(), 2 ) + TMath::Power( calc_plong_werr(), 2 ) ); }
// 
//   // Reconstructed azimuth [deg]
//   double calc_azimuth_werr()   {
//     return iso_init_azimuth->Interpolate( ( omega*( tof_werr ) - TMath::Floor( omega*( tof_werr )/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2., splat_azimuth_werr() );
//   }
//   
//   // Calculated initial value of px, py and pz, from calculated azimuth
//   double calc_pz_werr()  { return  calc_plong_werr(); }
//   double calc_px_werr()  { return ( calc_prad_werr() * TMath::Cos( calc_azimuth_werr() * TMath::DegToRad() ) ); }
//   double calc_py_werr()  { return ( calc_prad_werr() * TMath::Sin( calc_azimuth_werr() * TMath::DegToRad() ) ); }
// 
//   // Uncertainty on plong with errors [keV/c]
//   double delta_plong_werr()        {
//     vector<double> *pr = new vector<double>;
//     
//     for( double i=-1.0;i<1.1;i+=0.1 ){
//       for( double j=-1.0;j<1.1;j+=0.1 ){
//         pr->push_back( iso_init_plong->Interpolate( tof_werr+i*DEtof_MCP, mcp_splat_radius()_werr+j*DEsplat_radius_MCP ) );
//       }
//     }
//     
//     return ( TMath::RMS( pr->begin(), pr->end() ) );
//   }
// 
//   // Uncertainty on prad with errors [keV/c]
//   double delta_prad_werr()        {
//     vector<double> *pr = new vector<double>;
//     double sinomegat = TMath::Sin( omega*tof_werr/2. );
//     double deltat = omega/2. * TMath::Cos( omega*tof_werr/2. );
//     
//     for( double i=-1.0;i<1.1;i+=0.1 ){
//       for( double j=-1.0;j<1.1;j+=0.1 ){
//         pr->push_back( iso_init_prad->Interpolate( tof_werr+i*DEtof_MCP, mcp_splat_radius()_werr+j*DEsplat_radius_MCP ) );
//       }
//     }
//     
//     return ( TMath::RMS( pr->begin(), pr->end() ) );
//   }
// 
//   // Uncertainty on ptot with errors [keV/c]
//   double delta_ptot_werr()         { return ( TMath::Sqrt( TMath::Power( calc_prad_werr()*delta_prad_werr(), 2 ) + TMath::Power( calc_plong_werr()*delta_plong_werr(), 2 ) ) / calc_ptot_werr() ); }
// 
//   // Uncertainty on phi with errors [deg]
//   double delta_phi_werr()          { return ( TMath::Sqrt( TMath::Power( DEsplat_radius_MCP/mcp_splat_radius()_werr, 2 ) + TMath::Power( 0.5*omega*DEtof_MCP, 2 ) ) ) * TMath::RadToDeg(); }
// 
//   // Uncertainty on splat_azimuth with errors [deg]
//   double delta_splat_azimuth_werr(){ return ( TMath::Sqrt( TMath::Power( mcp_splat_x()_werr*DEsplat_radius_MCP, 2 ) + TMath::Power( mcp_splat_y()_werr*DEsplat_radius_MCP, 2 ) ) / ( mcp_splat_radius()_werr*mcp_splat_radius()_werr ) ) * TMath::RadToDeg(); }
// 
//   // Uncertainty on calc_azimuth with errors [keV/c]
//   double delta_azimuth_werr()        {
//     vector<double> *pr = new vector<double>;
//     double omegat = ( omega*tof_werr - TMath::Floor( omega*tof_werr/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2.;
//     
//     for( double i=-1.0;i<1.1;i+=0.1 ){
//       for( double j=-1.0;j<1.1;j+=0.1 ){
//         pr->push_back( iso_init_azimuth->Interpolate( omegat+omega*i*DEtof_MCP/2, splat_azimuth_werr()+j*delta_splat_azimuth_werr() ) );
//       }
//     }
//     
//     return ( TMath::RMS( pr->begin(), pr->end() ) );
//   }
  

  //////////////////
  // CUT CHECKERS //
  //////////////////

  bool mcp_cut_checker( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
  {
    bool VIOLATES_ONMCP = false;  // Auger must hit MCP
    bool VIOLATES_SPRAD = false;  // Auger mcp_splat_radius() must be CX_SR_LO( 3 ) < mcp_splat_radius() < CX_SR_HI( 60 )
    bool VIOLATES_TOF   = false;  // Auger tof must be < CX_TOF
//     bool VIOLATES_TOFLO = false;  // Low-energy Auger tof must be < CX_TOF_LO
//     bool VIOLATES_TOFPP = false;  // Low-energy Auger with 0.18 < tof< CX_TOF_HI must have mcp_splat_radius() > pp3( tof ) [this cut preserves backward Augers removing hi-energy e's which could not be reconstructed]
//     bool VIOLATES_TOFHI = false;  // Hi-energy Auger tof must be < CX_TOF_HI
    bool VIOLATES_NODEF = false;  // Auger must not arrive in a node ( this is a tof check )
    bool VIOLATES_CALCP = false;  // Auger's prad and plong must be, ultimately, reconstructable, i.e. != 0.
    double nodelo( 0. ), nodehi( 0. );
    if( onmcp ){ if( !( TMath::Abs( splat_z - ele_MCP_position ) < 0.01 ) )                                                                                 { VIOLATES_ONMCP = true; } }
    if( sprad ){ if( !( mcp_splat_radius() >= CX_SR_LOW && mcp_splat_radius() <= CX_SR_HIGH ) )                                                             { VIOLATES_SPRAD = true; } }
    if( timef ){
      if( !( tof < CX_TOF ) )                                                                                                                               { VIOLATES_TOF = true; }   // There is no more distinction between TOF_LO e TOF_HI
//       if( init_ke < 90. && !( tof < CX_TOF_LO ) )                                                                                              { VIOLATES_TOFLO = true; }   // Remember: we capture ALL the xenons meaning that
//       if( init_ke < 90. && tof > 0.18 && tof < CX_TOF_LO && !( mcp_splat_radius() > pp3->Eval( tof ) ) )                                               { VIOLATES_TOFPP = true; }   // we are able to determine if we have 2 Augers
//       if( init_ke > 90. && !( tof < CX_TOF_HI ) )                                                                                              { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 ) TEMPORANEAMENTE DISATTIVATO
//       if( init_ke > 90. && !( ( tof < CX_TOF_HI ) || ( tof > CX_TOF_HI_LOW && tof < CX_TOF_HI_HIGH ) ) )                                           { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 )
   }
    if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE ) { VIOLATES_NODEF = true; nodelo = nodes[nnodes]-CX_NODE; nodehi = nodes[nnodes]+CX_NODE; } } }
    if( calcp ){ if( mcp_calc_prad() == 0. || mcp_calc_plong() == 0. )                                                                                      { VIOLATES_CALCP = true; } }
//     if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOFLO || VIOLATES_TOFPP || VIOLATES_TOFHI || VIOLATES_NODEF || VIOLATES_CALCP ){
    if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOF || VIOLATES_NODEF || VIOLATES_CALCP ){
      if( VIOLATES_ONMCP )  { cout << "Auger not on MCP plane" << endl;                                                                                                 }
      if( VIOLATES_SPRAD )  { cout << "Auger splat_radius = " << mcp_splat_radius() << endl;                                                                            }
      if( VIOLATES_TOF )    { cout << "Auger TOF = " << tof << " [> CX_TOF = " << CX_TOF << "]" << endl;                                                                }
//       if( VIOLATES_TOFLO )  { cout << "Auger TOF = " << tof << " [> CX_TOF_LO = " << CX_TOF_LO << "]" << endl;                                              }
//       if( VIOLATES_TOFPP )  { cout << "Backward Auger, TOF = " << tof << " ; SR = " << mcp_splat_radius() << " < pp3[" << tof << "] = " << pp3->Eval( tof ) << endl; }
//       if( VIOLATES_TOFHI )  { cout << "Auger TOF = " << tof << " [> CX_TOF_HI = " << CX_TOF_HI << "]" << endl;                                              }
      if( VIOLATES_NODEF )  { cout << "Auger TOF = " << tof << " cut by CX_NODE: " << nodelo << " < " << tof << " < " << nodehi << endl;                                }
      if( VIOLATES_CALCP )  { cout << "Auger momentum not reconstructed: px = " << mcp_calc_px() << " ; py = " << mcp_calc_py() << " ; pz = " << mcp_calc_pz() << endl; }
      return false;
   }
    return true;
  }

  bool mcp_silent_cut_checker( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
  {
    if( onmcp ){ if( !( TMath::Abs( splat_z - ele_MCP_position ) < 0.01 ) )                                                                                 { return false; } }
    if( sprad ){ if( !( mcp_splat_radius() >= CX_SR_LOW && mcp_splat_radius() <= CX_SR_HIGH ) )                                                             { return false; } }
    if( timef ){
      if( !( tof < CX_TOF ) )                                                                                                                               { return false; }
//       if( init_ke < 90. && !( tof < CX_TOF_LO ) )                                                         { return false; }
//       else if( init_ke < 90. && tof > 0.18 && tof < CX_TOF_LO && !( mcp_splat_radius() > pp3->Eval( tof ) ) )     { return false; }
//       else if( init_ke > 90. && !( tof < CX_TOF_HI ) )                                                { return false; }
//       else if( init_ke > 90. && !( ( tof < CX_TOF_HI ) || ( tof > CX_TOF_HI_LOW && tof < CX_TOF_HI_HIGH ) ) ) { return false; }
   }
    if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE ) { return false; } } }
    if( calcp ){ if( mcp_calc_prad() == 0. || mcp_calc_plong() == 0. ){ return false; } }
    return true;
  }

} e;

// TChain that contains all the produced .root data files
TChain *so = new TChain( "tends" );


// Tree initialization routine
void tree_initializer(){
  so->Add( "toutput_roo/*.root" );

  so->SetScanField( 15 );

  so->SetBranchAddress( "ionn",          &e.ionn );
  so->SetBranchAddress( "tof",           &e.tof );
  so->SetBranchAddress( "mass",          &e.mass );
  so->SetBranchAddress( "charge",        &e.charge );
  so->SetBranchAddress( "init_x",        &e.init_x );
  so->SetBranchAddress( "init_y",        &e.init_y );
  so->SetBranchAddress( "init_z",        &e.init_z );
  so->SetBranchAddress( "init_azm",      &e.init_azm );
  so->SetBranchAddress( "init_elev",     &e.init_elev );
  so->SetBranchAddress( "init_dircos",   &e.init_dircos );
  so->SetBranchAddress( "init_vx",       &e.init_vx );
  so->SetBranchAddress( "init_vy",       &e.init_vy );
  so->SetBranchAddress( "init_vz",       &e.init_vz );
  so->SetBranchAddress( "init_vtot",     &e.init_vtot );
  so->SetBranchAddress( "init_vrad",     &e.init_vrad );
  so->SetBranchAddress( "init_vlong",    &e.init_vlong );
  so->SetBranchAddress( "init_px",       &e.init_px );
  so->SetBranchAddress( "init_py",       &e.init_py );
  so->SetBranchAddress( "init_pz",       &e.init_pz );
  so->SetBranchAddress( "init_ptot",     &e.init_ptot );
  so->SetBranchAddress( "init_prad",     &e.init_prad );
  so->SetBranchAddress( "init_plong",    &e.init_plong );
  so->SetBranchAddress( "init_ke",       &e.init_ke );
  so->SetBranchAddress( "splat_x",       &e.splat_x );
  so->SetBranchAddress( "splat_y",       &e.splat_y );
  so->SetBranchAddress( "splat_z",       &e.splat_z );
  so->SetBranchAddress( "splat_azm",     &e.splat_azm );
  so->SetBranchAddress( "splat_elev",    &e.splat_elev );
  so->SetBranchAddress( "splat_dircos",  &e.splat_dircos );
  so->SetBranchAddress( "splat_vx",      &e.splat_vx );
  so->SetBranchAddress( "splat_vy",      &e.splat_vy );
  so->SetBranchAddress( "splat_vz",      &e.splat_vz );
  so->SetBranchAddress( "splat_vtot",    &e.splat_vtot );
  so->SetBranchAddress( "splat_vrad",    &e.splat_vrad );
  so->SetBranchAddress( "splat_vlong",   &e.splat_vlong );
  so->SetBranchAddress( "splat_px",      &e.splat_px );
  so->SetBranchAddress( "splat_py",      &e.splat_py );
  so->SetBranchAddress( "splat_pz",      &e.splat_pz );
  so->SetBranchAddress( "splat_ptot",    &e.splat_ptot );
  so->SetBranchAddress( "splat_prad",    &e.splat_prad );
  so->SetBranchAddress( "splat_plong",   &e.splat_plong );
  so->SetBranchAddress( "splat_ke",      &e.splat_ke );
  so->SetBranchAddress( "splat_radius",  &e.splat_radius );
  so->SetBranchAddress( "run_id",        &e.run_id );
  so->SetBranchAddress( "ontarget",      &e.ontarget );
  so->SetBranchAddress( "onmcpplane",    &e.onmcpplane );

  so->SetBranchAddress( "vec_tof",       &e.vec_tof );
  so->SetBranchAddress( "vec_x",         &e.vec_x );
  so->SetBranchAddress( "vec_y",         &e.vec_y );
  so->SetBranchAddress( "vec_z",         &e.vec_z );
  so->SetBranchAddress( "vec_azm",       &e.vec_azm );
  so->SetBranchAddress( "vec_elev",      &e.vec_elev );
  so->SetBranchAddress( "vec_dircos",    &e.vec_dircos );
  so->SetBranchAddress( "vec_vx",        &e.vec_vx );
  so->SetBranchAddress( "vec_vy",        &e.vec_vy );
  so->SetBranchAddress( "vec_vz",        &e.vec_vz );
  so->SetBranchAddress( "vec_vtot",      &e.vec_vtot );
  so->SetBranchAddress( "vec_vrad",      &e.vec_vrad );
  so->SetBranchAddress( "vec_vlong",     &e.vec_vlong );
  so->SetBranchAddress( "vec_px",        &e.vec_px );
  so->SetBranchAddress( "vec_py",        &e.vec_py );
  so->SetBranchAddress( "vec_pz",        &e.vec_pz );
  so->SetBranchAddress( "vec_ptot",      &e.vec_ptot );
  so->SetBranchAddress( "vec_prad",      &e.vec_prad );
  so->SetBranchAddress( "vec_plong",     &e.vec_plong );
  so->SetBranchAddress( "vec_volt",      &e.vec_volt );
  so->SetBranchAddress( "vec_e",         &e.vec_e );
  so->SetBranchAddress( "vec_ex",        &e.vec_ex );
  so->SetBranchAddress( "vec_ey",        &e.vec_ey );
  so->SetBranchAddress( "vec_ez",        &e.vec_ez );
  so->SetBranchAddress( "vec_b",         &e.vec_b );
  so->SetBranchAddress( "vec_bx",        &e.vec_bx );
  so->SetBranchAddress( "vec_by",        &e.vec_by );
  so->SetBranchAddress( "vec_bz",        &e.vec_bz );
  so->SetBranchAddress( "vec_ke",        &e.vec_ke );
  so->SetBranchAddress( "vec_radius",    &e.vec_radius );

  pp3->SetParameters( CX_TOF_p0, CX_TOF_p1, CX_TOF_p2, CX_TOF_p3 );
  

// ISO_CURVE_THISTO - COMMENTA DA QUI...
// //   // ntuple for Lagrange polynomials interpolation
  tthist->SetBranchAddress("tf", &ltf);
  tthist->SetBranchAddress("sr", &lsr);
  tthist->SetBranchAddress("sa", &lsa);
  tthist->SetBranchAddress("pr", &lpr);
//   tthist->SetBranchAddress("bi", &tbi);
  
  rpb->SetBranchAddress("runs",&runids);
// ...A QUI
  
}

// Returns a pointer to vector<int> containing the track numbers corresponding to the selection parameters
vector<int> *FindTrack( const char* selection ){
  int nentries = so->GetEntries();
  TH1I *hist = new TH1I( "hist","hist",nentries,0,nentries );
  so->Draw( "Entry$>>hist",selection,"goff" );
  vector<int> *runs = new vector<int>;
  for( int i=0;i<hist->GetNbinsX()+1;++i ){
    if( hist->GetBinContent( i ) > 0 ){ runs->push_back( hist->GetBinLowEdge( i ) ); }
  }
  delete hist;
  return runs;
}

  ////////////////////////////////////////////////////////
  //    _____  _    _ _____  ______  __   ____     __   //
  //   |  __ \| |  | |  __ \|  ____| \ \ / /\ \   / /   //
  //   | |__) | |  | | |__) | |__     \ V /  \ \_/ /    //
  //   |  ___/| |  | |  _  /|  __|     > <    \   /     //
  //   | |    | |__| | | \ \| |____   / . \    | |      //
  //   |_|     \____/|_|  \_\______| /_/ \_\   |_|      //
  //                                                    //
  ////////////////////////////////////////////////////////

// Checks if an event passes selected cuts
bool cut_checker( int event, bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
{
  so->GetEntry( event );
  bool VIOLATES_ONMCP = false;  // Auger must hit MCP
  bool VIOLATES_SPRAD = false;  // Auger splat_radius must be CX_SR_LO( 3 ) < splat_radius < CX_SR_HI( 60 )
  bool VIOLATES_TOF   = false;  // Auger tof must be < CX_TOF
//   bool VIOLATES_TOFLO = false;  // Low-energy Auger tof must be < CX_TOF_LO
//   bool VIOLATES_TOFPP = false;  // Low-energy Auger with 0.18 < tof< CX_TOF_HI must have splat_radius > pp3( tof ) [this cut preserves backward Augers removing hi-energy e's which could not be reconstructed]
//   bool VIOLATES_TOFHI = false;  // Hi-energy Auger tof must be < CX_TOF_HI
  bool VIOLATES_NODEF = false;  // Auger must not arrive in a node ( this is a tof check )
  bool VIOLATES_CALCP = false;  // Auger's prad and plong must be, ultimately, reconstructable, i.e. != 0.
  double nodelo( 0. ), nodehi( 0. );
  if( onmcp ){ if( !( e.onmcpplane == true ) )                                                                                                                { VIOLATES_ONMCP = true; } }
  if( sprad ){ if( !( e.splat_radius >= CX_SR_LOW && e.splat_radius <= CX_SR_HIGH ) )                                                                         { VIOLATES_SPRAD = true; } }
  if( timef ){
    if( !( e.tof < CX_TOF ) )                                                                                                                                 { VIOLATES_TOF = true; }   // There is no more distinction between TOF_LO and TOF_HI
//     if( e.init_ke < 90. && !( e.tof < CX_TOF_LO ) )                                                                                              { VIOLATES_TOFLO = true; }   // Remember: we capture ALL the xenons meaning that
//     if( e.init_ke < 90. && e.tof > 0.18 && !( e.tof < CX_TOF_LO && e.splat_radius > pp3->Eval( e.tof ) ) )                                               { VIOLATES_TOFPP = true; }   // we are able to determine if we have 2 Augers
//     if( e.init_ke > 90. && !( e.tof < CX_TOF_HI ) )                                                                                              { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 )
  }
  if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( e.tof >= nodes[nnodes]-CX_NODE && e.tof <= nodes[nnodes]+CX_NODE ) { VIOLATES_NODEF = true; nodelo = nodes[nnodes]-CX_NODE; nodehi = nodes[nnodes]+CX_NODE; } } }
  if( calcp ){ if( e.calc_prad() == 0. || e.calc_plong() == 0. )                                                                                              { VIOLATES_CALCP = true; } }
//   if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOFLO || VIOLATES_TOFPP || VIOLATES_TOFHI || VIOLATES_NODEF || VIOLATES_CALCP ){
  if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOF || VIOLATES_NODEF || VIOLATES_CALCP ){
    if( VIOLATES_ONMCP )  { cout << "Auger not on MCP plane" << endl;                                                                                     }
    if( VIOLATES_SPRAD )  { cout << "Auger splat_radius = " << e.splat_radius << endl;                                                                      }
    if( VIOLATES_TOF )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF = " << CX_TOF << "]" << endl;                                              }
//     if( VIOLATES_TOFLO )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF_LO = " << CX_TOF_LO << "]" << endl;                                              }
//     if( VIOLATES_TOFPP )  { cout << "Backward Auger, TOF = " << e.tof << " ; SR = " << e.splat_radius << " < pp3[" << e.tof << "] = " << pp3->Eval( e.tof ) << endl; }
//     if( VIOLATES_TOFHI )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF_HI = " << CX_TOF_HI << "]" << endl;                                              }
    if( VIOLATES_NODEF )  { cout << "Auger TOF = " << e.tof << " cut by CX_NODE: " << nodelo << " < " << e.tof << " < " << nodehi << endl;                    }
    if( VIOLATES_CALCP )  { cout << "Auger momentum not reconstructed: px = " << e.calc_px() << " ; py = " << e.calc_py() << " ; pz = " << e.calc_pz() << endl; }
    return false;
  }
  return true;
}

// Checks if an event passes selected cuts
bool silent_cut_checker( int event, bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
{
  so->GetEntry( event );
  if( onmcp ){ if( !( e.onmcpplane == true ) ){ return false; } }
  if( sprad ){ if( !( e.splat_radius >= CX_SR_LOW && e.splat_radius <= CX_SR_HIGH ) )                                                                         { return false; } }
  //       if( timef ){ if( !( tof < CX_TOF ) ){ return false; } }    // static TOF cut
  //       if( timef ){ if( !( tof < CX_TOF_LO ) ){ return false; } }    // static TOF cut
  if( timef ){
    if     ( !( e.tof < CX_TOF ) )                                                                                                                            { return false; }
//     if     ( e.init_ke < 90. && !( e.tof < CX_TOF_LO ) )                                                { return false; }
//     else if( e.init_ke < 90. && !( e.tof > 0.18 && e.tof < CX_TOF_LO && e.splat_radius > pp3->Eval( e.tof ) ) ) { return false; }
//     else if( e.init_ke > 90. && !( e.tof < CX_TOF_HI ) )                                                { return false; }
  }
  if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( e.tof >= nodes[nnodes]-CX_NODE && e.tof <= nodes[nnodes]+CX_NODE ) { return false; } } }
  if( calcp ){ if( e.calc_prad() == 0. || e.calc_plong() == 0. )                                                                                              { return false; } }
  return true;
}

// Locates electrons on the three iso_curves
void iso_curve_locator(){
  TCanvas *cele = new TCanvas( "cele","cele",1200,400 ); cele->Divide( 3,1 );

  // iso_init_prad
  cele->cd( 1 ); iso_init_plong->Draw( "CONT1Z" );
  double plong_lox = iso_init_plong->GetXmin();
  double plong_loy = iso_init_plong->GetYmin();
  double plong_hix = iso_init_plong->GetXmax();
  double plong_hiy = iso_init_plong->GetYmax();
  if( plong_lox > CX_TOF ) { plong_lox = CX_TOF   ; }
  if( plong_hix < CX_TOF ) { plong_hix = CX_TOF   ; }
  if( plong_loy > CX_SR_LOW ) { plong_loy = CX_SR_LOW; }
  if( plong_hiy < CX_SR_LOW ) { plong_hiy = CX_SR_LOW; }
  TBox *blong_tof = new TBox( CX_TOF   ,plong_loy,plong_hix,plong_hiy ) ; blong_tof ->SetFillColor( kRed ); blong_tof ->SetFillStyle( 3002 ); blong_tof ->Draw( "SAME" );
  TBox *blong_sr  = new TBox( plong_lox,plong_loy,plong_hix,CX_SR_LOW ) ; blong_sr  ->SetFillColor( kRed ); blong_sr  ->SetFillStyle( 3002 ); blong_sr  ->Draw( "SAME" );
  TBox *blong_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < plong_lox && nodes[nnodes]+CX_NODE < plong_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > plong_hix && nodes[nnodes]+CX_NODE > plong_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < plong_lox && nodes[nnodes]+CX_NODE > plong_lox ){ blong_nodes = new TBox( plong_lox            , plong_loy, nodes[nnodes]+CX_NODE, plong_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < plong_hix && nodes[nnodes]+CX_NODE > plong_hix ){ blong_nodes = new TBox( nodes[nnodes]-CX_NODE, plong_loy, plong_hix            , plong_hiy ); }
    else                                                                              { blong_nodes = new TBox( nodes[nnodes]-CX_NODE, plong_loy, nodes[nnodes]+CX_NODE, plong_hiy ); }
    blong_nodes ->SetFillColor( kRed ); blong_nodes ->SetFillStyle( 3002 ); blong_nodes ->Draw( "SAME" );
  }
  TMarker *me0_l = new TMarker( e.tof, e.splat_radius, 22 );
  me0_l->SetMarkerColor( kBlue );
  me0_l->Draw( "SAME" );

  // iso_init_plong
  cele->cd( 2 ); iso_init_prad->Draw( "CONT1Z" );
  double prad_lox = iso_init_prad->GetXmin();
  double prad_loy = iso_init_prad->GetYmin();
  double prad_hix = iso_init_prad->GetXmax();
  double prad_hiy = iso_init_prad->GetYmax();
  if( prad_lox > CX_TOF ) { prad_lox = CX_TOF   ; }
  if( prad_hix < CX_TOF ) { prad_hix = CX_TOF   ; }
  if( prad_loy > CX_SR_LOW ) { prad_loy = CX_SR_LOW; }
  if( prad_hiy < CX_SR_LOW ) { prad_hiy = CX_SR_LOW; }
  TBox *brad_tof = new TBox( CX_TOF   ,prad_loy,prad_hix,prad_hiy ) ; brad_tof ->SetFillColor( kRed ); brad_tof ->SetFillStyle( 3002 ); brad_tof ->Draw( "SAME" );
  TBox *brad_sr  = new TBox( prad_lox,prad_loy,prad_hix,CX_SR_LOW ) ; brad_sr  ->SetFillColor( kRed ); brad_sr  ->SetFillStyle( 3002 ); brad_sr  ->Draw( "SAME" );
  TBox *brad_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < prad_lox && nodes[nnodes]+CX_NODE < prad_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > prad_hix && nodes[nnodes]+CX_NODE > prad_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < prad_lox && nodes[nnodes]+CX_NODE > prad_lox ){ brad_nodes = new TBox( prad_lox              , prad_loy, nodes[nnodes]+CX_NODE , prad_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < prad_hix && nodes[nnodes]+CX_NODE > prad_hix ){ brad_nodes = new TBox( nodes[nnodes]-CX_NODE , prad_loy, prad_hix              , prad_hiy ); }
    else                                                                            { brad_nodes = new TBox( nodes[nnodes]-CX_NODE , prad_loy, nodes[nnodes]+CX_NODE , prad_hiy ); }
    brad_nodes ->SetFillColor( kRed ); brad_nodes ->SetFillStyle( 3002 ); brad_nodes ->Draw( "SAME" );
  }
  TMarker *me0_r = new TMarker( e.tof, e.splat_radius, 22 );
  me0_r->SetMarkerColor( kBlue );
  me0_r->Draw( "SAME" );

  // iso_init_azimuth
  cele->cd( 3 ); iso_init_azimuth->Draw( "CONT1Z" );
  double azimuth_lox = iso_init_azimuth->GetXmin();
  double azimuth_loy = iso_init_azimuth->GetYmin();
  double azimuth_hix = iso_init_azimuth->GetXmax();
  double azimuth_hiy = iso_init_azimuth->GetYmax();
  if( azimuth_lox > CX_TOF ) { azimuth_lox = CX_TOF   ; }
  if( azimuth_hix < CX_TOF ) { azimuth_hix = CX_TOF   ; }
  if( azimuth_loy > CX_SR_LOW ) { azimuth_loy = CX_SR_LOW; }
  if( azimuth_hiy < CX_SR_LOW ) { azimuth_hiy = CX_SR_LOW; }
  TBox *bazi_tof = new TBox( CX_TOF   ,azimuth_loy,azimuth_hix,azimuth_hiy ) ; bazi_tof ->SetFillColor( kRed ); bazi_tof ->SetFillStyle( 3002 ); bazi_tof ->Draw( "SAME" );
  TBox *bazi_sr  = new TBox( azimuth_lox,azimuth_loy,azimuth_hix,CX_SR_LOW ) ; bazi_sr  ->SetFillColor( kRed ); bazi_sr  ->SetFillStyle( 3002 ); bazi_sr  ->Draw( "SAME" );
  TBox *bazi_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < azimuth_lox && nodes[nnodes]+CX_NODE < azimuth_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > azimuth_hix && nodes[nnodes]+CX_NODE > azimuth_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < azimuth_lox && nodes[nnodes]+CX_NODE > azimuth_lox ){ bazi_nodes = new TBox( azimuth_lox           , azimuth_loy, nodes[nnodes]+CX_NODE, azimuth_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < azimuth_hix && nodes[nnodes]+CX_NODE > azimuth_hix ){ bazi_nodes = new TBox( nodes[nnodes]-CX_NODE , azimuth_loy, azimuth_hix          , azimuth_hiy ); }
    else                                                                                  { bazi_nodes = new TBox( nodes[nnodes]-CX_NODE , azimuth_loy, nodes[nnodes]+CX_NODE, azimuth_hiy ); }
    bazi_nodes ->SetFillColor( kRed ); bazi_nodes ->SetFillStyle( 3002 ); bazi_nodes ->Draw( "SAME" );
  }
  TMarker *me0_a = new TMarker( e.tof, e.splat_azimuth(), 22 );
  me0_a->SetMarkerColor( kBlue );
  me0_a->Draw( "SAME" );
  cout << "e0 = \e[34m▲\e[0m ; e1 = \e[31m▼\e[0m" << endl;
}

// Presents a table containing the number of fired/collected/reconstructed electrons for each energy
void fired_collected_reconstructed( bool b023 = true, bool b034 = true, bool b045 = true, bool b054 = true, bool b065 = true, bool b100 = true, bool b111 = true, bool b122 = true ){
  int n023( 0 ), c023( 0 ), r023( 0 );
  int n034( 0 ), c034( 0 ), r034( 0 );
  int n045( 0 ), c045( 0 ), r045( 0 );
  int n054( 0 ), c054( 0 ), r054( 0 );
  int n065( 0 ), c065( 0 ), r065( 0 );
  int n100( 0 ), c100( 0 ), r100( 0 );
  int n111( 0 ), c111( 0 ), r111( 0 );
  int n122( 0 ), c122( 0 ), r122( 0 );
  
  for( int i=0;i<so->GetEntries();++i ){ so->GetEntry( i );
    if( b023 ){ if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++n023; if( e.ontarget ){ ++c023; } if( e.silent_cut_checker() ){ ++r023; } }  }
    if( b034 ){ if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++n034; if( e.ontarget ){ ++c034; } if( e.silent_cut_checker() ){ ++r034; } }  }
    if( b045 ){ if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++n045; if( e.ontarget ){ ++c045; } if( e.silent_cut_checker() ){ ++r045; } }  }
    if( b054 ){ if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++n054; if( e.ontarget ){ ++c054; } if( e.silent_cut_checker() ){ ++r054; } }  }
    if( b065 ){ if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++n065; if( e.ontarget ){ ++c065; } if( e.silent_cut_checker() ){ ++r065; } }  }
    if( b100 ){ if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++n100; if( e.ontarget ){ ++c100; } if( e.silent_cut_checker() ){ ++r100; } }  }
    if( b111 ){ if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++n111; if( e.ontarget ){ ++c111; } if( e.silent_cut_checker() ){ ++r111; } }  }
    if( b122 ){ if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++n122; if( e.ontarget ){ ++c122; } if( e.silent_cut_checker() ){ ++r122; } }  }
  }
  
            printf( "E ( eV ),Fired,Collected,Frac.Collected,Reconstructed,Rec.Efficiency,Overall Eff.\n"                 );
  if( b023 ){ printf(  "23,%d,%d,%.03f,%d,%.03f,%.03f\n", n023, c023, ( double )c023/n023, r023, ( double )r023/c023, ( double )r023/n023 ); }
  if( b034 ){ printf(  "34,%d,%d,%.03f,%d,%.03f,%.03f\n", n034, c034, ( double )c034/n034, r034, ( double )r034/c034, ( double )r034/n034 ); }
  if( b045 ){ printf(  "45,%d,%d,%.03f,%d,%.03f,%.03f\n", n045, c045, ( double )c045/n045, r045, ( double )r045/c045, ( double )r045/n045 ); }
  if( b054 ){ printf(  "54,%d,%d,%.03f,%d,%.03f,%.03f\n", n054, c054, ( double )c054/n054, r054, ( double )r054/c054, ( double )r054/n054 ); }
  if( b065 ){ printf(  "65,%d,%d,%.03f,%d,%.03f,%.03f\n", n065, c065, ( double )c065/n065, r065, ( double )r065/c065, ( double )r065/n065 ); }
  if( b100 ){ printf( "100,%d,%d,%.03f,%d,%.03f,%.03f\n", n100, c100, ( double )c100/n100, r100, ( double )r100/c100, ( double )r100/n100 ); }
  if( b111 ){ printf( "111,%d,%d,%.03f,%d,%.03f,%.03f\n", n111, c111, ( double )c111/n111, r111, ( double )r111/c111, ( double )r111/n111 ); }
  if( b122 ){ printf( "122,%d,%d,%.03f,%d,%.03f,%.03f\n", n122, c122, ( double )c122/n122, r122, ( double )r122/c122, ( double )r122/n122 ); }
            printf( "%s\n"                            , "https://ozh.github.io/ascii-tables/"                           );
}

// Presents a table containing the number of Augers surviving each cut
void cuts_efficiency( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true ){
  std::array<bool,   5> init_vals { onmcp   , sprad  , timef   , nodec    , calcp     };
  std::array<bool,   5> test_vals { 0       , 0      , 0       , 0        , 0         };
  std::array<string, 5> strings   { "CX_MCP", "CX_SR", "CX_TOF", "CX_NODE", "CX_PXYZ" };
  
  std::array<int, 6> t023 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t034 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t045 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t054 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t065 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t100 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t111 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t122 { 0, 0, 0, 0, 0, 0 };
  
  for( int iter=0;iter<so->GetEntries();++iter ){
    so->GetEntry( iter );
    if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++t023[5]; }
    if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++t034[5]; }
    if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++t045[5]; }
    if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++t054[5]; }
    if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++t065[5]; }
    if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++t100[5]; }
    if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++t111[5]; }
    if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++t122[5]; }
    
    for( int i=0;i<5;++i ){
      
      if ( test_vals[i] == init_vals[i] ){ continue; }
      else{
        test_vals[i] = init_vals[i];
        if( e.silent_cut_checker( test_vals[0], test_vals[1], test_vals[2], test_vals[3], test_vals[4] ) ){
          if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++t023[i]; }
          if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++t034[i]; }
          if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++t045[i]; }
          if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++t054[i]; }
          if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++t065[i]; }
          if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++t100[i]; }
          if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++t111[i]; }
          if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++t122[i]; }
       }
     }
   }
    for( int j=0;j<5;++j ){ test_vals[j] = 0; }
  }
  cout << "E ( eV ),Fired,"       ; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout << strings[i].c_str() << "," ; } } cout << endl;
  cout << " 23," << t023[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t023[i]/t023[5])         << "," ; } } cout << endl;
  cout << " 34," << t034[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t034[i]/t034[5])         << "," ; } } cout << endl;
  cout << " 45," << t045[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t045[i]/t045[5])         << "," ; } } cout << endl;
  cout << " 54," << t054[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t054[i]/t054[5])         << "," ; } } cout << endl;
  cout << " 65," << t065[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t065[i]/t065[5])         << "," ; } } cout << endl;
  cout << "100," << t100[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t100[i]/t100[5])         << "," ; } } cout << endl;
  cout << "111," << t111[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t111[i]/t111[5])         << "," ; } } cout << endl;
  cout << "122," << t122[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    TString::Format("%.03f", (double)t122[i]/t122[5])         << "," ; } } cout << endl;
  cout << "https://ozh.github.io/ascii-tables/"                                                                       << endl;
}









  //////////////////////////////////////////////////
  //    __  __  _____ _____   __   ____     __    //
  //   |  \/  |/ ____|  __ \  \ \ / /\ \   / /    //
  //   | \  / | |    | |__) |  \ V /  \ \_/ /     //
  //   | |\/| | |    |  ___/    > <    \   /      //
  //   | |  | | |____| |       / . \    | |       //
  //   |_|  |_|\_____|_|      /_/ \_\   |_|       //
  //                                              //
  //////////////////////////////////////////////////

// Checks if an event passes selected cuts
bool mcp_cut_checker( int event, bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
{
  so->GetEntry( event );
  bool VIOLATES_ONMCP = false;  // Auger must hit MCP
  bool VIOLATES_SPRAD = false;  // Auger mcp_splat_radius() must be CX_SR_LO( 3 ) < mcp_splat_radius() < CX_SR_HI( 60 )
  bool VIOLATES_TOF   = false;  // Auger tof must be < CX_TOF
//   bool VIOLATES_TOFLO = false;  // Low-energy Auger tof must be < CX_TOF_LO
//   bool VIOLATES_TOFPP = false;  // Low-energy Auger with 0.18 < tof< CX_TOF_HI must have mcp_splat_radius() > pp3( tof ) [this cut preserves backward Augers removing hi-energy e's which could not be reconstructed]
//   bool VIOLATES_TOFHI = false;  // Hi-energy Auger tof must be < CX_TOF_HI
  bool VIOLATES_NODEF = false;  // Auger must not arrive in a node ( this is a tof check )
  bool VIOLATES_CALCP = false;  // Auger's prad and plong must be, ultimately, reconstructable, i.e. != 0.
  double nodelo( 0. ), nodehi( 0. );
  if( onmcp ){ if( !( e.onmcpplane == true ) )                                                                                                                { VIOLATES_ONMCP = true; } }
  if( sprad ){ if( !( e.mcp_splat_radius() >= CX_SR_LOW && e.mcp_splat_radius() <= CX_SR_HIGH ) )                                                             { VIOLATES_SPRAD = true; } }
  if( timef ){
    if( !( e.tof < CX_TOF ) )                                                                                                                                 { VIOLATES_TOF = true; }   // There is no more distinction between TOF_LO and TOF_HI
//     if( e.init_ke < 90. && !( e.tof < CX_TOF_LO ) )                                                                                              { VIOLATES_TOFLO = true; }   // Remember: we capture ALL the xenons meaning that
//     if( e.init_ke < 90. && e.tof > 0.18 && !( e.tof < CX_TOF_LO && e.mcp_splat_radius() > pp3->Eval( e.tof ) ) )                                               { VIOLATES_TOFPP = true; }   // we are able to determine if we have 2 Augers
//     if( e.init_ke > 90. && !( e.tof < CX_TOF_HI ) )                                                                                              { VIOLATES_TOFHI = true; }   // ( init_ke < 90 ) or 1 Auger ( init_ke > 90 )
  }
  if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( e.tof >= nodes[nnodes]-CX_NODE && e.tof <= nodes[nnodes]+CX_NODE ) { VIOLATES_NODEF = true; nodelo = nodes[nnodes]-CX_NODE; nodehi = nodes[nnodes]+CX_NODE; } } }
  if( calcp ){ if( e.mcp_calc_prad() == 0. || e.mcp_calc_plong() == 0. )                                                                                      { VIOLATES_CALCP = true; } }
//   if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOFLO || VIOLATES_TOFPP || VIOLATES_TOFHI || VIOLATES_NODEF || VIOLATES_CALCP ){
  if( VIOLATES_ONMCP || VIOLATES_SPRAD || VIOLATES_TOF || VIOLATES_NODEF || VIOLATES_CALCP ){
    if( VIOLATES_ONMCP )  { cout << "Auger not on MCP plane" << endl;                                                                                                       }
    if( VIOLATES_SPRAD )  { cout << "Auger splat_radius = " << e.mcp_splat_radius() << endl;                                                                                }
    if( VIOLATES_TOF )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF = " << CX_TOF << "]" << endl;                                                                      }
//     if( VIOLATES_TOFLO )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF_LO = " << CX_TOF_LO << "]" << endl;                                              }
//     if( VIOLATES_TOFPP )  { cout << "Backward Auger, TOF = " << e.tof << " ; SR = " << e.mcp_splat_radius() << " < pp3[" << e.tof << "] = " << pp3->Eval( e.tof ) << endl; }
//     if( VIOLATES_TOFHI )  { cout << "Auger TOF = " << e.tof << " [> CX_TOF_HI = " << CX_TOF_HI << "]" << endl;                                              }
    if( VIOLATES_NODEF )  { cout << "Auger TOF = " << e.tof << " cut by CX_NODE: " << nodelo << " < " << e.tof << " < " << nodehi << endl;                                  }
    if( VIOLATES_CALCP )  { cout << "Auger momentum not reconstructed: px = " << e.mcp_calc_px() << " ; py = " << e.mcp_calc_py() << " ; pz = " << e.mcp_calc_pz() << endl; }
    return false;
  }
  return true;
}

// Checks if an event passes selected cuts
bool mcp_silent_cut_checker( int event, bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true )
{
  so->GetEntry( event );
  if( onmcp ){ if( !( e.onmcpplane == true ) ){ return false; } }
  if( sprad ){ if( !( e.mcp_splat_radius() >= CX_SR_LOW && e.mcp_splat_radius() <= CX_SR_HIGH ) )                                                             { return false; } }
  //       if( timef ){ if( !( tof < CX_TOF ) ){ return false; } }    // static TOF cut
  //       if( timef ){ if( !( tof < CX_TOF_LO ) ){ return false; } }    // static TOF cut
  if( timef ){
    if     ( !( e.tof < CX_TOF ) )                                                                                                                            { return false; }
//     if     ( e.init_ke < 90. && !( e.tof < CX_TOF_LO ) )                                                { return false; }
//     else if( e.init_ke < 90. && !( e.tof > 0.18 && e.tof < CX_TOF_LO && e.mcp_splat_radius() > pp3->Eval( e.tof ) ) ) { return false; }
//     else if( e.init_ke > 90. && !( e.tof < CX_TOF_HI ) )                                                { return false; }
  }
  if( nodec ){ for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ if( e.tof >= nodes[nnodes]-CX_NODE && e.tof <= nodes[nnodes]+CX_NODE ) { return false; } } }
  if( calcp ){ if( e.calc_prad() == 0. || e.calc_plong() == 0. )                                                                                              { return false; } }
  return true;
}

// Locates electrons on the three iso_curves
void mcp_iso_curve_locator(){
  TCanvas *cele = new TCanvas( "cele","cele",1200,400 ); cele->Divide( 3,1 );

  // iso_init_prad
  cele->cd( 1 ); iso_init_plong->Draw( "CONT1Z" );
  double plong_lox = iso_init_plong->GetXmin();
  double plong_loy = iso_init_plong->GetYmin();
  double plong_hix = iso_init_plong->GetXmax();
  double plong_hiy = iso_init_plong->GetYmax();
  if( plong_lox > CX_TOF ) { plong_lox = CX_TOF   ; }
  if( plong_hix < CX_TOF ) { plong_hix = CX_TOF   ; }
  if( plong_loy > CX_SR_LOW ) { plong_loy = CX_SR_LOW; }
  if( plong_hiy < CX_SR_LOW ) { plong_hiy = CX_SR_LOW; }
  TBox *blong_tof = new TBox( CX_TOF   ,plong_loy,plong_hix,plong_hiy ) ; blong_tof ->SetFillColor( kRed ); blong_tof ->SetFillStyle( 3002 ); blong_tof ->Draw( "SAME" );
  TBox *blong_sr  = new TBox( plong_lox,plong_loy,plong_hix,CX_SR_LOW ) ; blong_sr  ->SetFillColor( kRed ); blong_sr  ->SetFillStyle( 3002 ); blong_sr  ->Draw( "SAME" );
  TBox *blong_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < plong_lox && nodes[nnodes]+CX_NODE < plong_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > plong_hix && nodes[nnodes]+CX_NODE > plong_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < plong_lox && nodes[nnodes]+CX_NODE > plong_lox ){ blong_nodes = new TBox( plong_lox            , plong_loy, nodes[nnodes]+CX_NODE, plong_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < plong_hix && nodes[nnodes]+CX_NODE > plong_hix ){ blong_nodes = new TBox( nodes[nnodes]-CX_NODE, plong_loy, plong_hix            , plong_hiy ); }
    else                                                                              { blong_nodes = new TBox( nodes[nnodes]-CX_NODE, plong_loy, nodes[nnodes]+CX_NODE, plong_hiy ); }
    blong_nodes ->SetFillColor( kRed ); blong_nodes ->SetFillStyle( 3002 ); blong_nodes ->Draw( "SAME" );
  }
  TMarker *me0_l = new TMarker( e.tof, e.mcp_splat_radius(), 22 );
  me0_l->SetMarkerColor( kBlue );
  me0_l->Draw( "SAME" );

  // iso_init_plong
  cele->cd( 2 ); iso_init_prad->Draw( "CONT1Z" );
  double prad_lox = iso_init_prad->GetXmin();
  double prad_loy = iso_init_prad->GetYmin();
  double prad_hix = iso_init_prad->GetXmax();
  double prad_hiy = iso_init_prad->GetYmax();
  if( prad_lox > CX_TOF ) { prad_lox = CX_TOF   ; }
  if( prad_hix < CX_TOF ) { prad_hix = CX_TOF   ; }
  if( prad_loy > CX_SR_LOW ) { prad_loy = CX_SR_LOW; }
  if( prad_hiy < CX_SR_LOW ) { prad_hiy = CX_SR_LOW; }
  TBox *brad_tof = new TBox( CX_TOF   ,prad_loy,prad_hix,prad_hiy ) ; brad_tof ->SetFillColor( kRed ); brad_tof ->SetFillStyle( 3002 ); brad_tof ->Draw( "SAME" );
  TBox *brad_sr  = new TBox( prad_lox,prad_loy,prad_hix,CX_SR_LOW ) ; brad_sr  ->SetFillColor( kRed ); brad_sr  ->SetFillStyle( 3002 ); brad_sr  ->Draw( "SAME" );
  TBox *brad_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < prad_lox && nodes[nnodes]+CX_NODE < prad_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > prad_hix && nodes[nnodes]+CX_NODE > prad_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < prad_lox && nodes[nnodes]+CX_NODE > prad_lox ){ brad_nodes = new TBox( prad_lox              , prad_loy, nodes[nnodes]+CX_NODE , prad_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < prad_hix && nodes[nnodes]+CX_NODE > prad_hix ){ brad_nodes = new TBox( nodes[nnodes]-CX_NODE , prad_loy, prad_hix              , prad_hiy ); }
    else                                                                            { brad_nodes = new TBox( nodes[nnodes]-CX_NODE , prad_loy, nodes[nnodes]+CX_NODE , prad_hiy ); }
    brad_nodes ->SetFillColor( kRed ); brad_nodes ->SetFillStyle( 3002 ); brad_nodes ->Draw( "SAME" );
  }
  TMarker *me0_r = new TMarker( e.tof, e.mcp_splat_radius(), 22 );
  me0_r->SetMarkerColor( kBlue );
  me0_r->Draw( "SAME" );

  // iso_init_azimuth
  cele->cd( 3 ); iso_init_azimuth->Draw( "CONT1Z" );
  double azimuth_lox = iso_init_azimuth->GetXmin();
  double azimuth_loy = iso_init_azimuth->GetYmin();
  double azimuth_hix = iso_init_azimuth->GetXmax();
  double azimuth_hiy = iso_init_azimuth->GetYmax();
  if( azimuth_lox > CX_TOF ) { azimuth_lox = CX_TOF   ; }
  if( azimuth_hix < CX_TOF ) { azimuth_hix = CX_TOF   ; }
  if( azimuth_loy > CX_SR_LOW ) { azimuth_loy = CX_SR_LOW; }
  if( azimuth_hiy < CX_SR_LOW ) { azimuth_hiy = CX_SR_LOW; }
  TBox *bazi_tof = new TBox( CX_TOF   ,azimuth_loy,azimuth_hix,azimuth_hiy ) ; bazi_tof ->SetFillColor( kRed ); bazi_tof ->SetFillStyle( 3002 ); bazi_tof ->Draw( "SAME" );
  TBox *bazi_sr  = new TBox( azimuth_lox,azimuth_loy,azimuth_hix,CX_SR_LOW ) ; bazi_sr  ->SetFillColor( kRed ); bazi_sr  ->SetFillStyle( 3002 ); bazi_sr  ->Draw( "SAME" );
  TBox *bazi_nodes;
  for( int nnodes=0;nnodes<sizeof( nodes )/sizeof( double );++nnodes ){
    if( nodes[nnodes]-CX_NODE < azimuth_lox && nodes[nnodes]+CX_NODE < azimuth_lox ){ continue; }
    if( nodes[nnodes]-CX_NODE > azimuth_hix && nodes[nnodes]+CX_NODE > azimuth_hix ){ continue; }
    if      ( nodes[nnodes]-CX_NODE < azimuth_lox && nodes[nnodes]+CX_NODE > azimuth_lox ){ bazi_nodes = new TBox( azimuth_lox           , azimuth_loy, nodes[nnodes]+CX_NODE, azimuth_hiy ); }
    else if ( nodes[nnodes]-CX_NODE < azimuth_hix && nodes[nnodes]+CX_NODE > azimuth_hix ){ bazi_nodes = new TBox( nodes[nnodes]-CX_NODE , azimuth_loy, azimuth_hix          , azimuth_hiy ); }
    else                                                                                  { bazi_nodes = new TBox( nodes[nnodes]-CX_NODE , azimuth_loy, nodes[nnodes]+CX_NODE, azimuth_hiy ); }
    bazi_nodes ->SetFillColor( kRed ); bazi_nodes ->SetFillStyle( 3002 ); bazi_nodes ->Draw( "SAME" );
  }
  TMarker *me0_a = new TMarker( e.tof, e.splat_azimuth(), 22 );
  me0_a->SetMarkerColor( kBlue );
  me0_a->Draw( "SAME" );
  cout << "e0 = \e[34m▲\e[0m ; e1 = \e[31m▼\e[0m" << endl;
}

// Presents a table containing the number of fired/collected/reconstructed electrons for each energy
void mcp_fired_collected_reconstructed( bool b023 = true, bool b034 = true, bool b045 = true, bool b054 = true, bool b065 = true, bool b100 = true, bool b111 = true, bool b122 = true ){
  int n023( 0 ), c023( 0 ), r023( 0 );
  int n034( 0 ), c034( 0 ), r034( 0 );
  int n045( 0 ), c045( 0 ), r045( 0 );
  int n054( 0 ), c054( 0 ), r054( 0 );
  int n065( 0 ), c065( 0 ), r065( 0 );
  int n100( 0 ), c100( 0 ), r100( 0 );
  int n111( 0 ), c111( 0 ), r111( 0 );
  int n122( 0 ), c122( 0 ), r122( 0 );
  
  for( int i=0; i<so->GetEntries(); ++i ){ so->GetEntry( i );
    if( b023 ){ if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++n023; if( e.ontarget ){ ++c023; } if( e.mcp_silent_cut_checker() ){ ++r023; } }  }
    if( b034 ){ if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++n034; if( e.ontarget ){ ++c034; } if( e.mcp_silent_cut_checker() ){ ++r034; } }  }
    if( b045 ){ if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++n045; if( e.ontarget ){ ++c045; } if( e.mcp_silent_cut_checker() ){ ++r045; } }  }
    if( b054 ){ if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++n054; if( e.ontarget ){ ++c054; } if( e.mcp_silent_cut_checker() ){ ++r054; } }  }
    if( b065 ){ if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++n065; if( e.ontarget ){ ++c065; } if( e.mcp_silent_cut_checker() ){ ++r065; } }  }
    if( b100 ){ if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++n100; if( e.ontarget ){ ++c100; } if( e.mcp_silent_cut_checker() ){ ++r100; } }  }
    if( b111 ){ if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++n111; if( e.ontarget ){ ++c111; } if( e.mcp_silent_cut_checker() ){ ++r111; } }  }
    if( b122 ){ if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++n122; if( e.ontarget ){ ++c122; } if( e.mcp_silent_cut_checker() ){ ++r122; } }  }
  }
  
              printf( "E ( eV ),Fired,Collected,Frac.Collected,Reconstructed,Rec.Efficiency,Overall Eff.\n"                               );
  if( b023 ){ printf(  "23,%d,%d,%.03f,%d,%.03f,%.03f\n", n023, c023, ( double )c023/n023, r023, ( double )r023/c023, ( double )r023/n023 ); }
  if( b034 ){ printf(  "34,%d,%d,%.03f,%d,%.03f,%.03f\n", n034, c034, ( double )c034/n034, r034, ( double )r034/c034, ( double )r034/n034 ); }
  if( b045 ){ printf(  "45,%d,%d,%.03f,%d,%.03f,%.03f\n", n045, c045, ( double )c045/n045, r045, ( double )r045/c045, ( double )r045/n045 ); }
  if( b054 ){ printf(  "54,%d,%d,%.03f,%d,%.03f,%.03f\n", n054, c054, ( double )c054/n054, r054, ( double )r054/c054, ( double )r054/n054 ); }
  if( b065 ){ printf(  "65,%d,%d,%.03f,%d,%.03f,%.03f\n", n065, c065, ( double )c065/n065, r065, ( double )r065/c065, ( double )r065/n065 ); }
  if( b100 ){ printf( "100,%d,%d,%.03f,%d,%.03f,%.03f\n", n100, c100, ( double )c100/n100, r100, ( double )r100/c100, ( double )r100/n100 ); }
  if( b111 ){ printf( "111,%d,%d,%.03f,%d,%.03f,%.03f\n", n111, c111, ( double )c111/n111, r111, ( double )r111/c111, ( double )r111/n111 ); }
  if( b122 ){ printf( "122,%d,%d,%.03f,%d,%.03f,%.03f\n", n122, c122, ( double )c122/n122, r122, ( double )r122/c122, ( double )r122/n122 ); }
            printf( "%s\n"                            , "https://ozh.github.io/ascii-tables/"                           );
}

// Presents a table containing the number of Augers surviving each cut
void mcp_cuts_efficiency( bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true, bool calcp = true ){
  std::array<bool,   5> init_vals { onmcp   , sprad  , timef   , nodec    , calcp     };
  std::array<bool,   5> test_vals { 0       , 0      , 0       , 0        , 0         };
  std::array<string, 5> strings   { "CX_MCP", "CX_SR", "CX_TOF", "CX_NODE", "CX_PXYZ" };
  
  std::array<int, 6> t023 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t034 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t045 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t054 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t065 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t100 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t111 { 0, 0, 0, 0, 0, 0 };
  std::array<int, 6> t122 { 0, 0, 0, 0, 0, 0 };
  
  for( int iter=0;iter<so->GetEntries();++iter ){
    so->GetEntry( iter );
    if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++t023[5]; }
    if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++t034[5]; }
    if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++t045[5]; }
    if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++t054[5]; }
    if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++t065[5]; }
    if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++t100[5]; }
    if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++t111[5]; }
    if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++t122[5]; }
    
    for( int i=0;i<5;++i ){
      
      if ( test_vals[i] == init_vals[i] ){ continue; }
      else{
        test_vals[i] = init_vals[i];
        if( e.mcp_silent_cut_checker( test_vals[0], test_vals[1], test_vals[2], test_vals[3], test_vals[4] ) ){
          if( e.init_ke >  22.9 && e.init_ke <  23.1 ){ ++t023[i]; }
          if( e.init_ke >  33.9 && e.init_ke <  34.1 ){ ++t034[i]; }
          if( e.init_ke >  44.9 && e.init_ke <  45.1 ){ ++t045[i]; }
          if( e.init_ke >  53.9 && e.init_ke <  54.1 ){ ++t054[i]; }
          if( e.init_ke >  64.9 && e.init_ke <  65.1 ){ ++t065[i]; }
          if( e.init_ke >  99.9 && e.init_ke < 100.1 ){ ++t100[i]; }
          if( e.init_ke > 110.9 && e.init_ke < 111.1 ){ ++t111[i]; }
          if( e.init_ke > 121.9 && e.init_ke < 122.1 ){ ++t122[i]; }
       }
     }
   }
    for( int j=0;j<5;++j ){ test_vals[j] = 0; }
  }
  cout << "E ( eV ),Fired,"       ; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout << strings[i].c_str() << "," ; } } cout << endl;
  cout << " 23," << t023[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t023[i]         << "," ; } } cout << endl;
  cout << " 34," << t034[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t034[i]         << "," ; } } cout << endl;
  cout << " 45," << t045[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t045[i]         << "," ; } } cout << endl;
  cout << " 54," << t054[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t054[i]         << "," ; } } cout << endl;
  cout << " 65," << t065[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t065[i]         << "," ; } } cout << endl;
  cout << "100," << t100[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t100[i]         << "," ; } } cout << endl;
  cout << "111," << t111[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t111[i]         << "," ; } } cout << endl;
  cout << "122," << t122[5] << ","; for( int i=0;i<5;++i ){ if( init_vals[i] ){ cout <<    t122[i]         << "," ; } } cout << endl;
  cout << "https://ozh.github.io/ascii-tables/"                                                                       << endl;
}
