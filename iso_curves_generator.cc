///////////////////////////////////////////////////////////////////////////////////
// iso_curves_generator.c (vX.0) - Author: Francesco Granato
// 2019-06-18
//
// NOTE: the initial version of this program has been written...a long ago, but I
// didn't comment it. Today 20190618 I'm writing the comments and the description.
//
// This program opens a .root file which has been generated with the increasing
// angle and energy .fly2 file. It plots the associated iso-curves applying the
// cuts described in the program itself.
//
// X.0  - Initial draft
///////////////////////////////////////////////////////////////////////////////////

{
  bool REGULAR = false;

  REGULAR = true;

  // Selection of the runid file
  long runid;
  if(REGULAR) { runid =  64501111; }  // isos-LoRes_HiAccept RealB8 ForwBack
  else        { runid = 1461602; }  // INVERSE isos-LoRes_HiAccept RealB8 ForwBack

//   int runid = 6248462;  // isos-LoRes_HiAccept IdealB8 ForwBack
//   int runid = 4952819;  // isos-HiRes_LoAccept IdealB8 ForwBack
//   vector<int> *tr2 = FindTrack( TString::Format("run_id == %d && init_dircos > -0.5",runid) );
  vector<int> *tr2 = FindTrack( TString::Format("init_ke > %d",0) );
//   vector<int> *tr2 = FindTrack( TString::Format("run_id == %d",runid) );

  // Definition of the TGraph2D's for the iso-curves
  TGraph2D *gprad = new TGraph2D();   // init_prad
  TGraph2D *gplng = new TGraph2D();   // init_plong
  TGraph2D *gazma = new TGraph2D();   // init_azimuth
  
  TGraph2D *gelev = new TGraph2D();   // TEST init_elevation
  
  TGraph2D *gsprd = new TGraph2D();   // TEST splat_radius
  TGraph2D *gprwt = new TGraph2D();   // TEST init_prwt
  TGraph2D *gpxxx = new TGraph2D();   // TEST init_px
  TGraph2D *gpyyy = new TGraph2D();   // TEST init_py
  TGraph2D *gpzzz = new TGraph2D();   // TEST init_pz
  int nprad(0), nplng(0), nazma(0);
  int nelev(0); // TEST pwrt
  
  int nprwt(0); // TEST pwrt
  int nsprd(0); // TEST SPLAT_RADIUS
  int npxxx(0), npyyy(0), npzzz(0); // TEST PX, PY, PZ
  
  
  // TTree FOR LINEAR INTERPOLATION
  double ttf, tsr, tsa, tpr;
  int tbi;
  TTree *tthist = new TTree("tthist","tthist");
  tthist->Branch("tf", &ttf);
  tthist->Branch("sr", &tsr);
  tthist->Branch("sa", &tsa);
  tthist->Branch("pr", &tpr);
  tthist->Branch("bi", &tbi);
  
  vector< vector<int>* > *vcont = new vector< vector<int>* >;
  vector<int> *runs;

  for(int i=0; i<(htree->GetNbinsX()+2)*(htree->GetNbinsY()+2)*(htree->GetNbinsZ()+2); ++i){
    runs = new vector<int>;
    vcont->push_back(runs);
  }
  
// // // //   TTree *ttrunids = new TTree("");  WORKING ON THIS!!!
  
//   double tfmin(0.086), tfmax(0.206);
//   double srmin(0.000), srmax(60.00);
//   double samin(0.000), samax(360.0);
//   TH3D *htree = new TH3D( "htree","htree", 40, tfmin, tfmax, 30, srmin, srmax, 180, samin, samax );

  // Naming and formatting of TGraph2D's
  // iso_init_prad
  gprad->SetName("iso_init_prad");
  gprad->SetTitle( TString::Format("Contour lines for iso-init_prad - %ld; TOF (#mus); splat_radius (mm); init_prad (keV/c)", runid) );
//   gprad->SetTitle( TString::Format("Contour lines for iso-init_prad - %d; sin(#frac{#omegat}{2}); splat_radius (mm); init_prad (keV/c)", runid) );
  gprad->GetZaxis()->SetLabelOffset(0.0025);
  gprad->GetZaxis()->SetTitleSize(0.02);
  gprad->GetZaxis()->SetLabelSize(0.02);

  // iso_init_plong
  gplng->SetName("iso_init_plong");
  gplng->SetTitle( TString::Format("Contour lines for iso-init_plong - %ld; TOF (#mus); splat_radius (mm); init_plong (keV/c)", runid) );
  gplng->GetZaxis()->SetLabelOffset(0.0025);
  gplng->GetZaxis()->SetTitleSize(0.02);
  gplng->GetZaxis()->SetLabelSize(0.02);

  // iso_init_azimuth
  gazma->SetName("iso_init_azimuth");
  gazma->SetTitle( TString::Format("Contour lines for iso-init_azimuth - %ld; TOF (#mus); splat_azimuth (deg); init_azimuth (deg)", runid) );
//   gazma->SetTitle( TString::Format("Contour lines for iso-init_azimuth - %d; (#omegat - floor(#frac{#omegat}{2#pi})#upoint2#pi)/2; splat_azimuth (deg); init_azimuth (deg)", runid) );
  gazma->GetZaxis()->SetLabelOffset(0.0025);
  gazma->GetZaxis()->SetTitleSize(0.02);
  gazma->GetZaxis()->SetLabelSize(0.02);

  // TEST iso_init_elev
  gelev->SetName("iso_init_elevation");
  gelev->SetTitle( TString::Format("Contour lines for iso-init_elevation - %ld; TOF (#mus); splat_radius (mm); init_elevation (deg)", runid) );
  gelev->GetZaxis()->SetLabelOffset(0.0025);
  gelev->GetZaxis()->SetTitleSize(0.02);
  gelev->GetZaxis()->SetLabelSize(0.02);
  
//   // TEST iso_init_prad_wt
//   gprwt->SetName("iso_init_prwt");
// //   gprad->SetTitle( TString::Format("Contour lines for iso-init_prad - %d; TOF (#mus); splat_radius (mm); init_prad (keV/c)", runid) );
//   gprwt->SetTitle( TString::Format("Contour lines for iso-init_prwt - %d; sin(#frac{#omega t}{2}); splat_radius (mm); init_prad (keV/c)", runid) );
//   gprwt->GetZaxis()->SetLabelOffset(0.0025);
//   gprwt->GetZaxis()->SetTitleSize(0.02);
//   gprwt->GetZaxis()->SetLabelSize(0.02);
// 
//   // TEST iso_splat_radius
//   gsprd->SetName("iso_splat_radius");
//   gsprd->SetTitle( TString::Format("Contour lines for iso-splat_radius - %d; TOF (#mus); splat_azimuth (deg); splat_radius (mm)", runid) );
//   gsprd->GetZaxis()->SetLabelOffset(0.0025);
//   gsprd->GetZaxis()->SetTitleSize(0.02);
//   gsprd->GetZaxis()->SetLabelSize(0.02);
// 
//   // TEST iso_init_px/py/pz
//   gpxxx->SetName("iso_init_px");
//   gpxxx->SetTitle( TString::Format("Contour lines for iso-init_px - %d; TOF (#mus); splat_azimuth (deg); init_px (keV/c)", runid) );
//   gpxxx->GetZaxis()->SetLabelOffset(0.0025);
//   gpxxx->GetZaxis()->SetTitleSize(0.02);
//   gpxxx->GetZaxis()->SetLabelSize(0.02);
// 
//   gpyyy->SetName("iso_init_py");
//   gpyyy->SetTitle( TString::Format("Contour lines for iso-init_py - %d; TOF (#mus); splat_azimuth (deg); init_py (keV/c)", runid) );
//   gpyyy->GetZaxis()->SetLabelOffset(0.0025);
//   gpyyy->GetZaxis()->SetTitleSize(0.02);
//   gpyyy->GetZaxis()->SetLabelSize(0.02);
// 
//   gpzzz->SetName("iso_init_pz");
//   gpzzz->SetTitle( TString::Format("Contour lines for iso-init_pz - %d; TOF (#mus); splat_azimuth (deg); init_pz (keV/c)", runid) );
//   gpzzz->GetZaxis()->SetLabelOffset(0.0025);
//   gpzzz->GetZaxis()->SetTitleSize(0.02);
//   gpzzz->GetZaxis()->SetLabelSize(0.02);



  // The iso-curve for init_azimuth requires rescaling the init_azimuth values by
  // +360 depending if their value is <200°.
  // Graphically, there are diagonal zones in the unmodified iso-curve for azimuth
  // that correspond to the gap 360°-0° in the plotting process.
  // The funciton that describes the diagonal line of separation between the various
  // zones depends on the value of the node, and from a value, delta_azimuth, that
  // has been empirically determined to be 179 (I suspect 180° would work too). The
  // final formula is:
  //
  //    -4000*x + p0 + (nn-init_node)*delta_azimuth)
  //
  // where p0 is a starting value empirically determined. nn goes from init_node up
  // to (init_node+active_nodes): these values have to be determined in case there's
  // a change in the magnetic field and the nodes situation changes. In general, the
  // value of init_node is the first one which begins to be relevant and present in
  // the init_azimuth plot; active_nodes is the number of nodes that appear in the
  // aforementioned plot. A good way to obtain these values is the following:
  // * obtain the values for the nodes with the Mathematica code;
  // * plot iso_init_azimuth WITHOUT any azimuth cut;
  // * take note of what are the initial node present on the plot, the number of
  //   active nodes, and the value of p0, and correct the values accordingly.

  // for cut azimuth
  double delta_azimuth = 179.;

  vector<double> *vnodes = new vector<double>;
//   if(REGULAR) { vnodes->insert( vnodes->end() , {0.0446548, 0.0893097, 0.133965, 0.178619, 0.223274, 0.267929, 0.312584, 0.357239, 0.401894, 0.446548, 0.491203, 0.535858, 0.580513, 0.625168, 0.669823, 0.714477, 0.759132, 0.803787, 0.848442, 0.893097, 0.937752} ); }   // REGULAR potentials
  if(REGULAR) { vnodes->insert( vnodes->end() , {0., 0.1255, 0.1663, 0.2032, 0.2475, 10000} ); }   // Manually inserted for COMSOL generated Bfieldmap
  else        { vnodes->insert( vnodes->end() , {0, 0.1345, 0.1755, 0.2203, 0.2675, 0.3145, 0.3600} ); }   // INVERSE potentials

  vector<TF1*> *f_azimuth = new vector<TF1*>;
  int nnodes = vnodes->size();
  int init_node(0), active_nodes(0);
  double p0_base(0), p1_base(0), p0(0), p1(0);
//   if(REGULAR) { init_node = 2; active_nodes = 7; p0_base =  712; p1 = -4000; }   // REGULAR potentials
  if(REGULAR) { init_node = 1; active_nodes = 4; p0_base =  720; p1_base = -4335; }   // Manually inserted for COMSOL generated Bfieldmap - Lowering p0_base shifts the "lines" down without changing the angle (and that should be enough)
  else        { init_node = 1; active_nodes = 5; p0_base = -356; p1 =  4000; }   // INVERSE potentials
//   double p0s[] = {0.,   719.,   898,  1074.,  1251.5};  // NoPMT
//   double p1s[] = {0., -4335., -4337, -4336., -4336.5};  // NoPMT
  double p0s[] = {0.,   723.,   903,  1083.,  1251.5};  // WithPMT
  double p1s[] = {0., -4330., -4337, -4336., -4336.5};  // WithPMT

  for(int nn=init_node;nn<(init_node+active_nodes);++nn){
//     p0 = p0_base-(nn-1)*0.5 - (nn-init_node) * TMath::Power(-1,REGULAR) * delta_azimuth - nn; // (-1)^REGULAR == [-1 if REGULAR == true, +1 if REGULAR == false]
//     p1 = p1_base - (nn-1)*0.5;
    
    p0 = p0s[nn];
    p1 = p1s[nn];
    TF1 *dummy = new TF1("dummy",   TString::Format("%f*x + %f", p1, p0),   vnodes->at(nn-1),vnodes->at(nn));
    f_azimuth->push_back(dummy);
  }

  // CUTS APPLICATION
  // If one of these cuts is not satisfied, the event is rejected
  bool bex = false;   // control bool
  
  
  
  
  // PROGRESS COUNTER
  
  // OUTSIDE THE LOOP
  // PRE-COUNTER START
  int niter = tr2->size();   // This has to be set to the number of items you are iterating over
  int numb = 0;                   // Variable used to count steps
  int step_delta, step;           // step_delta is used to select the counting speed (see formula); step determines when a new printing of progress is displayed
  if(niter > 10.){
    step_delta = 0; // smaller values (even < 0) = faster counter; larger values = slower counter
    step = TMath::Power(10,step_delta + TMath::Floor(TMath::Log10(niter)));    // Event counter
  }
  // THREE WAYS TO COUNT TIME DIFFERENCES
  // ROOT and std::chrono give exactly the same results
  // clock() gives the time based on CPU clock ticks
  
  // ROOT time counting ( . $ROOTSYS/bin/thisroot.sh && # include "TTimeStamp.h" )
  TTimeStamp *rttime_counter = new TTimeStamp();      // ROOT structure to count time
  rttime_counter->Set();                              // Sets the time counter to now
  double rttime_start = rttime_counter->AsDouble();   // Converts the value of the structure into a double
  double rttime_previous = rttime_start;              // Sets the value of "time_previous" to now, for zeroth iteration
  // PRE-COUNTER STOP
  
  
  
  
  
  
  
  
  int tthist_entries = 0;

  for(int i=0;i<tr2->size();++i){
//   for(int i=0;i<5000;++i){
    
    // COUNTER START
    if(niter > 10.){
      if((numb)%step == 0 || numb == niter-1){
        // ROOT
        rttime_counter->Set();
        double rttime_now = rttime_counter->AsDouble();   // Sets a new time counter to now
        double rttelapsed = rttime_now-rttime_start;      // Computes elapsed time
        double rtdeltat = rttime_now-rttime_previous;     // Computes Δt between current and previous step
        rttime_previous = rttime_now;                     // Sets value of time_previous to now for next iteration
        printf( "%d/%d\t%.2fs elapsed\t Δt = %.2fs\n", numb, niter, rttelapsed, rtdeltat );
      }
    }
    ++numb;
    // COUNTER STOP
    
    
    so->GetEntry(tr2->at(i));
//     if(REGULAR) { if(e.init_elevation() > 120.){continue;} }  // events with init_elevation > 120 are rejected tout court
    if(REGULAR) { if(e.init_elevation() > 113.){continue;} }  // TEMPORARY!!!
    else        { if(e.init_elevation() <  90.){continue;} }  // events with init_elevation > 120 are rejected tout court


    if(!(TMath::Abs(e.splat_z - ele_MCP_position) < 0.01)){bex = true;}
//     if(!(e.onmcpplane == true)){bex = true;}
    if(!(e.splat_radius <= CX_SR_HIGH)){bex = true;}  // NOTE: we don't cut on CX_SR_LOW because those events are rejected with the software cuts: the iso-curve MUST contain them
    if(REGULAR) { if(!(e.tof < CX_TOF_LO)){bex = true;} }  // CX_TOF needed only for REGULAR potentials

    if(bex == true){bex = false; continue;}

    // Filling of TGraph2D's
//     AAAAAAAAAAAAAAAAAAAAAAA
    gprad->SetPoint(nprad,e.tof,e.splat_radius,e.init_prad) ; ++nprad;

//     tt->Fill(e.tof,e.splat_radius,e.splat_azimuth(),e.init_prad);
    ttf = e.tof;
    tsr = e.splat_radius;
    tsa = e.splat_azimuth();
    tpr = e.init_prad;
    tbi = htree->FindBin(ttf, tsr, tsa);
    tthist->Fill();
//     cout << ttf << " " << tsr << " " << tsa << " " << tbi << endl;
    vcont->at( tbi )->push_back(tthist_entries);
    ++tthist_entries;
    

    
    
//     AAAAAAAAAAAAAAAAAAAAAAA
    gplng->SetPoint(nplng,e.tof,e.splat_radius,e.init_plong); ++nplng;
    
    // TEST 
//     gelev->SetPoint(nelev,e.tof,e.splat_radius,e.init_elevation()); ++nelev;

    // Uncomment next line to fill gazma without any azimuth cut (when, for example, you have changed the B field)
//     gazma->SetPoint(nazma,e.tof,e.splat_azimuth(),e.init_azimuth())  ; ++nazma;


    // Filling of gazma: the position in the TOF-splat_azimuth plot determines if it's
    // below or over the appropriate diagonal line, and eventually the value of
    // init_azimuth is corrected by adding +360°.
    bool nodevar = false;
    for(int nn=init_node;nn<(init_node+active_nodes);++nn){
      if(e.tof >= vnodes->at(nn-1) && e.tof < vnodes->at(nn) && e.splat_azimuth() > f_azimuth->at(nn-init_node)->Eval(e.tof) && e.init_azimuth() < 356.5 ){
        nodevar = true; continue;
      }
    }
//     nodevar = false;
    if(nodevar) { gazma->SetPoint(nazma,e.tof,e.splat_azimuth(),e.init_azimuth()+360) ; ++nazma; }
    else        { gazma->SetPoint(nazma,e.tof,e.splat_azimuth(),e.init_azimuth())     ; ++nazma; }


    // TEST SPLAT_RADIUS
//     if(e.init_elevation() > 90.){ gsprd->SetPoint(nazma,e.tof,e.splat_azimuth(),e.init_azimuth())  ; ++nsprd; }

    // TEST PX, PY
//     gpzzz->SetPoint(npzzz,e.tof,e.splat_radius,e.init_pz)  ; ++npzzz;
//     gpxxx->SetPoint(npxxx,e.init_azimuth(),e.init_elevation(),e.init_py)  ; ++npxxx;
//     gpyyy->SetPoint(npyyy,e.init_azimuth(),e.init_elevation(),e.init_pz)  ; ++npyyy;

  }

  // Plotting of iso-curves
  int precision = 500;

//   TCanvas *celev = new TCanvas("celev","celev",800,800);
//   gelev->SetNpx(precision);  gelev->SetNpy(precision);
//   gelev->Draw("CONT1Z");
//   celev->Update();

  TCanvas *cprad = new TCanvas("cprad","cprad",800,800);
  gprad->SetNpx(precision);  gprad->SetNpy(precision);
  gprad->Draw("CONT1Z");
  cprad->Update();

  TCanvas *cplng = new TCanvas("cplng","cplng",800,800);
  gplng->SetNpx(precision);  gplng->SetNpy(precision);
  gplng->Draw("CONT1Z");
  cplng->Update();

  TCanvas *cazma = new TCanvas("cazma","cazma",800,800);
  gazma->SetNpx(precision);  gazma->SetNpy(precision);
  gazma->Draw("CONT1Z");
  cazma->Update();

//   TCanvas *cprwt = new TCanvas("cprwt","cprwt",800,800);
//   gprwt->SetNpx(precision);  gprwt->SetNpy(precision);
//   gprwt->Draw("CONT1Z");
//   cprwt->Update();

//   TCanvas *csprd = new TCanvas("csprd","csprd",800,800);
//   gsprd->SetNpx(precision);  gsprd->SetNpy(precision);
//   gsprd->Draw("CONT1Z");
//   csprd->Update();

//   TCanvas *cpxxx = new TCanvas("cpxxx","cpxxx",800,800);
//   gpxxx->SetNpx(precision);  gpxxx->SetNpy(precision);
//   gpxxx->Draw("CONT1Z");
//   cpxxx->Update();
//
//   TCanvas *cpyyy = new TCanvas("cpyyy","cpyyy",800,800);
//   gpyyy->SetNpx(precision);  gpyyy->SetNpy(precision);
//   gpyyy->Draw("CONT1Z");
//   cpyyy->Update();
//
//   TCanvas *cpzzz = new TCanvas("cpzzz","cpzzz",800,800);
//   gpzzz->SetNpx(precision);  gpzzz->SetNpy(precision);
//   gpzzz->Draw("CONT1Z");
//   cpzzz->Update();


  // OUTPUT: after creating the iso-curves it's time to save them. Copy-paste in ROOT
  // the appropriate lines to save the TGraph2D's in a suitable file.
  // Remember that the introduction message in rootlogon.C reads and prints the
  // description TNamed iso_comment. The text is colored depending if using high resolution
  // or low resolution.
  // The color is given with " \033[1;32m " [colors text in GREEN] and " \033[1;31m "
  // [colors text in RED]. It wants " \033[0m " to be closed.

//     TNamed *aa = new TNamed("iso_comment","real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TNamed *aa = new TNamed("iso_comment","INVERSE potentials, real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TFile *f = new TFile("INVERSE_isos_realB8_Forward_LoRes_HiAcc_4100V_tof_360deg_lessthan120degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","real B field (8 gauss), rings in 316LN - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TFile *f = new TFile("isos_realB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan120degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, phi step 1 deg, theta step 1 deg, KE 1-150 step 1 eV - real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, phi step 1 deg, theta step 1 deg, ptot 3-15 step 0.1 keV - real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TFile *f = new TFile("isos_GEANT4phi1_theta1_ptot3-13_0.1_realB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan120degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, 181 uniform theta, 361 uniform phi, KE 1-150 step 1 eV - real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TFile *f = new TFile("isos_GEANT4_181uniformtheta1_361uniformphi_KE1-150_1eV_realB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan120degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, 1801 uniform theta, phi = 0, KE 1-150 step 0.1 eV - real B field (8 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 120 deg");
//     TFile *f = new TFile("isos_GEANT4_1801uniformtheta1_phi0_KE1-150_0.1eV_realB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan120degcut.root","recreate");

//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, phi step 1 deg, theta step 1 deg, KE 1-150 step 1 eV - COMSOL real B field (~8.5 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 113 deg");
//     TFile *f = new TFile("isos_GEANT4phi1_theta1_KE1-150_1eV_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, phi step 1 deg, theta step 0.5 deg, Range 5eV around Standard KE, dE 0.25eV  - COMSOL real B field (~8.5 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 113 deg");
//     TFile *f = new TFile("isos_GEANT4phi1_theta0.5_StandardErange5eV_dE025eV_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root","recreate");
//     TNamed *aa = new TNamed("iso_comment","GEANT4 simulation, phi step 1 deg, theta step 1 deg, KE 1-150 step 1 eV, WithPMT 450mm - COMSOL real B field (~8.5 gauss) - \033[1;31mHiAcceptance field (~0.4100 V/mm)\033[0m, ForwBack hemispheres, 360 deg with cuts, only less than 113 deg");
//     TFile *f = new TFile("isos_GEANT4phi1_theta1_KE1-150_1eV_WithPMT450mm_COMSOLrealB8_ForBack_LoRes_HiAcc_4100V_tof_316LN_360deg_lessthan113degcut.root","recreate");
//     gprad->Write();
//     gplng->Write();
//     gazma->Write();
//     aa->Write();
//     TTree *tthistogram = (TTree*)tthist->Clone();
//     TTree *run_per_bin = new TTree("run_per_bin","run_per_bin");
//     run_per_bin->Branch("runs",&runs);
//     for(int i=0;i<vcont->size();++i){runs = vcont->at(i);run_per_bin->Fill();}
//     tthistogram->Write();
//     run_per_bin->Write();

// .q

}
