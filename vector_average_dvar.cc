{
//   so->SetBranchStatus("*",0);
//   so->SetBranchStatus("init_px"     ,1);
//   so->SetBranchStatus("init_py"     ,1);
//   so->SetBranchStatus("init_pz"     ,1);
//   so->SetBranchStatus("init_ptot"   ,1);
//   so->SetBranchStatus("init_ke"     ,1);
//   so->SetBranchStatus("splat_x"     ,1);
//   so->SetBranchStatus("splat_y"     ,1);
//   so->SetBranchStatus("splat_z"     ,1);
//   so->SetBranchStatus("splat_radius",1);
//   so->SetBranchStatus("tof"         ,1);
  
  int n023 = 12609; n023 = 0;
  int n034 = 12727; n034 = 0;
  int n045 = 12517; n045 = 0;
  int n054 = 12441; n054 = 0;
  int n065 = 12630; n065 = 0;
  int n100 = 12385; n100 = 0;
  int n111 = 12415; n111 = 0;
  int n122 = 12276; n122 = 0;
  
  so->SetBranchStatus("*",0);
  so->SetBranchStatus("init_ke"     ,1);
  for(int i=0;i<so->GetEntries();++i){so->GetEntry(i);
               if(TMath::Abs(e.init_ke- 23)<0.1){ ++n023; }
          else if(TMath::Abs(e.init_ke- 34)<0.1){ ++n034; }
          else if(TMath::Abs(e.init_ke- 45)<0.1){ ++n045; }
          else if(TMath::Abs(e.init_ke- 54)<0.1){ ++n054; }
          else if(TMath::Abs(e.init_ke- 65)<0.1){ ++n065; }
          else if(TMath::Abs(e.init_ke-100)<0.1){ ++n100; }
          else if(TMath::Abs(e.init_ke-111)<0.1){ ++n111; }
          else if(TMath::Abs(e.init_ke-122)<0.1){ ++n122; }
  }
  so->SetBranchStatus("*",1);
  
  gStyle->SetOptStat(1); // Draw histogram with statbox
  int dn = tout->GetEntries();
  
  
  vector<double> *vvar023 = new vector<double>;
  vector<double> *vvar034 = new vector<double>;
  vector<double> *vvar045 = new vector<double>;
  vector<double> *vvar054 = new vector<double>;
  vector<double> *vvar065 = new vector<double>;
  vector<double> *vvar100 = new vector<double>;
  vector<double> *vvar111 = new vector<double>;
  vector<double> *vvar122 = new vector<double>;
  
  std::cout.width(11);
  
  // PRE-COUNTER START
  int niter = tout->GetEntries();   // This has to be set to the number of items you are iterating over
  int numb = 0;                   // Variable used to count steps
  int step_delta, step;           // step_delta is used to select the counting speed (see formula); step determines when a new printing of progress is displayed
  if(niter > 10.){
    step_delta = 0; // smaller values (even < 0) = faster counter; larger values = slower counter
    step = TMath::Power(10,step_delta + TMath::Floor(TMath::Log10(niter)));    // Event counter
  }
  
  // ROOT time counting ( . $ROOTSYS/bin/thisroot.sh && # include "TTimeStamp.h" )
  TTimeStamp *rttime_counter = new TTimeStamp();      // ROOT structure to count time
  rttime_counter->Set();                              // Sets the time counter to now
  double rttime_start = rttime_counter->AsDouble();   // Converts the value of the structure into a double
  double rttime_previous = rttime_start;              // Sets the value of "time_previous" to now, for zeroth iteration
  // PRE-COUNTER STOP
  
  

  for(int i=0;i<tout->GetEntries();++i){
    tout->GetEntry(i);
    
    // COUNTER START
    if(niter > 10.){
      if((numb)%step == 0){
        // ROOT
        rttime_counter->Set();
        double rttime_now = rttime_counter->AsDouble();   // Sets a new time counter to now
        double rttelapsed = rttime_now-rttime_start;      // Computes elapsed time
        double rtdeltat = rttime_now-rttime_previous;     // Computes Δt between current and previous step
        rttime_previous = rttime_now;                     // Sets value of time_previous to now for next iteration
        printf( "%d/%d (%.02f%%)\t%.2fs elapsed\t Δt = %.2fs\n", numb, niter, (double)numb/niter*100., rttelapsed, rtdeltat );
      }
    }
    ++numb;
    // COUNTER STOP
    
//     if(!(e.silent_cut_checker())){continue;}
    double initke = var[2];
    double dvar   = var[1];
    
               if(TMath::Abs(initke- 23)<0.1){ if(TMath::Abs(dvar) < 10){ vvar023->push_back(dvar); } }
          else if(TMath::Abs(initke- 34)<0.1){ if(TMath::Abs(dvar) < 10){ vvar034->push_back(dvar); } }
          else if(TMath::Abs(initke- 45)<0.1){ if(TMath::Abs(dvar) < 10){ vvar045->push_back(dvar); } }
          else if(TMath::Abs(initke- 54)<0.1){ if(TMath::Abs(dvar) < 10){ vvar054->push_back(dvar); } }
          else if(TMath::Abs(initke- 65)<0.1){ if(TMath::Abs(dvar) < 10){ vvar065->push_back(dvar); } }
          else if(TMath::Abs(initke-100)<0.1){ if(TMath::Abs(dvar) < 10){ vvar100->push_back(dvar); } }
          else if(TMath::Abs(initke-111)<0.1){ if(TMath::Abs(dvar) < 10){ vvar111->push_back(dvar); } }
          else if(TMath::Abs(initke-122)<0.1){ if(TMath::Abs(dvar) < 10){ vvar122->push_back(dvar); } }
  }
  
  vector<double> *avgs = new vector<double>;
  vector<double> *RMSs = new vector<double>;
  
  double avg023 = TMath::Mean( vvar023->begin(), vvar023->end() );  double rms023 = TMath::RMS( vvar023->begin(), vvar023->end() ); 
  double avg034 = TMath::Mean( vvar034->begin(), vvar034->end() );  double rms034 = TMath::RMS( vvar034->begin(), vvar034->end() );
  double avg045 = TMath::Mean( vvar045->begin(), vvar045->end() );  double rms045 = TMath::RMS( vvar045->begin(), vvar045->end() );
  double avg054 = TMath::Mean( vvar054->begin(), vvar054->end() );  double rms054 = TMath::RMS( vvar054->begin(), vvar054->end() );
  double avg065 = TMath::Mean( vvar065->begin(), vvar065->end() );  double rms065 = TMath::RMS( vvar065->begin(), vvar065->end() );
  double avg100 = TMath::Mean( vvar100->begin(), vvar100->end() );  double rms100 = TMath::RMS( vvar100->begin(), vvar100->end() );
  double avg111 = TMath::Mean( vvar111->begin(), vvar111->end() );  double rms111 = TMath::RMS( vvar111->begin(), vvar111->end() );
  double avg122 = TMath::Mean( vvar122->begin(), vvar122->end() );  double rms122 = TMath::RMS( vvar122->begin(), vvar122->end() );

  
  if(vvar023->size() > 0) { avgs->push_back( avg023 ); RMSs->push_back( rms023 ); }
  if(vvar034->size() > 0) { avgs->push_back( avg034 ); RMSs->push_back( rms034 ); }
  if(vvar045->size() > 0) { avgs->push_back( avg045 ); RMSs->push_back( rms045 ); }
  if(vvar054->size() > 0) { avgs->push_back( avg054 ); RMSs->push_back( rms054 ); }
  if(vvar065->size() > 0) { avgs->push_back( avg065 ); RMSs->push_back( rms065 ); }
  if(vvar100->size() > 0) { avgs->push_back( avg100 ); RMSs->push_back( rms100 ); }
  if(vvar111->size() > 0) { avgs->push_back( avg111 ); RMSs->push_back( rms111 ); }
  if(vvar122->size() > 0) { avgs->push_back( avg122 ); RMSs->push_back( rms122 ); }
  
  double tot_avg = TMath::Mean( avgs->begin(), avgs->end() );
  std::sort( RMSs->begin(), RMSs->end() ); double tot_RMS = RMSs->at( RMSs->size()-1 );
//   double tot_RMS = TMath::Mean( RMSs->begin(), RMSs->end() );
  tot_avg = 0;
  tot_RMS = 0.12;
  
  TH1D *hvar023 = new TH1D("hvar023", TString::Format( "Distribution of %s -  23 eV", vary->GetTitle() ), 100, avg023 - 5.*rms023, avg023 + 5.*rms023 ); for(int i=0;i<vvar023->size(); ++i ){ hvar023->Fill( vvar023->at(i) ); }
  TH1D *hvar034 = new TH1D("hvar034", TString::Format( "Distribution of %s -  34 eV", vary->GetTitle() ), 100, avg034 - 5.*rms034, avg034 + 5.*rms034 ); for(int i=0;i<vvar034->size(); ++i ){ hvar034->Fill( vvar034->at(i) ); }
  TH1D *hvar045 = new TH1D("hvar045", TString::Format( "Distribution of %s -  45 eV", vary->GetTitle() ), 100, avg045 - 5.*rms045, avg045 + 5.*rms045 ); for(int i=0;i<vvar045->size(); ++i ){ hvar045->Fill( vvar045->at(i) ); }
  TH1D *hvar054 = new TH1D("hvar054", TString::Format( "Distribution of %s -  54 eV", vary->GetTitle() ), 100, avg054 - 5.*rms054, avg054 + 5.*rms054 ); for(int i=0;i<vvar054->size(); ++i ){ hvar054->Fill( vvar054->at(i) ); }
  TH1D *hvar065 = new TH1D("hvar065", TString::Format( "Distribution of %s -  65 eV", vary->GetTitle() ), 100, avg065 - 5.*rms065, avg065 + 5.*rms065 ); for(int i=0;i<vvar065->size(); ++i ){ hvar065->Fill( vvar065->at(i) ); }
  TH1D *hvar100 = new TH1D("hvar100", TString::Format( "Distribution of %s - 100 eV", vary->GetTitle() ), 100, avg100 - 5.*rms100, avg100 + 5.*rms100 ); for(int i=0;i<vvar100->size(); ++i ){ hvar100->Fill( vvar100->at(i) ); }
  TH1D *hvar111 = new TH1D("hvar111", TString::Format( "Distribution of %s - 111 eV", vary->GetTitle() ), 100, avg111 - 5.*rms111, avg111 + 5.*rms111 ); for(int i=0;i<vvar111->size(); ++i ){ hvar111->Fill( vvar111->at(i) ); }
  TH1D *hvar122 = new TH1D("hvar122", TString::Format( "Distribution of %s - 122 eV", vary->GetTitle() ), 100, avg122 - 5.*rms122, avg122 + 5.*rms122 ); for(int i=0;i<vvar122->size(); ++i ){ hvar122->Fill( vvar122->at(i) ); }
//   TH1D *hvar023 = new TH1D("hvar023", TString::Format( "Distribution of %s -  23 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar023->size(); ++i ){ hvar023->Fill( vvar023->at(i) ); }
//   TH1D *hvar034 = new TH1D("hvar034", TString::Format( "Distribution of %s -  34 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar034->size(); ++i ){ hvar034->Fill( vvar034->at(i) ); }
//   TH1D *hvar045 = new TH1D("hvar045", TString::Format( "Distribution of %s -  45 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar045->size(); ++i ){ hvar045->Fill( vvar045->at(i) ); }
//   TH1D *hvar054 = new TH1D("hvar054", TString::Format( "Distribution of %s -  54 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar054->size(); ++i ){ hvar054->Fill( vvar054->at(i) ); }
//   TH1D *hvar065 = new TH1D("hvar065", TString::Format( "Distribution of %s -  65 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar065->size(); ++i ){ hvar065->Fill( vvar065->at(i) ); }
//   TH1D *hvar100 = new TH1D("hvar100", TString::Format( "Distribution of %s - 100 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar100->size(); ++i ){ hvar100->Fill( vvar100->at(i) ); }
//   TH1D *hvar111 = new TH1D("hvar111", TString::Format( "Distribution of %s - 111 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar111->size(); ++i ){ hvar111->Fill( vvar111->at(i) ); }
//   TH1D *hvar122 = new TH1D("hvar122", TString::Format( "Distribution of %s - 122 eV", vary->GetTitle() ), 100, tot_avg - 5.*tot_RMS, tot_avg + 5.*tot_RMS ); for(int i=0;i<vvar122->size(); ++i ){ hvar122->Fill( vvar122->at(i) ); }
  
  hvar023->SetLineColor(kBlue); hvar023->SetLineWidth(2); hvar023->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar023->GetYaxis()->SetTitle( "Counts" );
  hvar034->SetLineColor(kBlue); hvar034->SetLineWidth(2); hvar034->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar034->GetYaxis()->SetTitle( "Counts" );
  hvar045->SetLineColor(kBlue); hvar045->SetLineWidth(2); hvar045->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar045->GetYaxis()->SetTitle( "Counts" );
  hvar054->SetLineColor(kBlue); hvar054->SetLineWidth(2); hvar054->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar054->GetYaxis()->SetTitle( "Counts" );
  hvar065->SetLineColor(kBlue); hvar065->SetLineWidth(2); hvar065->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar065->GetYaxis()->SetTitle( "Counts" );
  hvar100->SetLineColor(kRed ); hvar100->SetLineWidth(2); hvar100->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar100->GetYaxis()->SetTitle( "Counts" );
  hvar111->SetLineColor(kRed ); hvar111->SetLineWidth(2); hvar111->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar111->GetYaxis()->SetTitle( "Counts" );
  hvar122->SetLineColor(kRed ); hvar122->SetLineWidth(2); hvar122->GetXaxis()->SetTitle( TString::Format( "%s %s", vary->GetTitle(), unity->GetTitle() ) ); hvar122->GetYaxis()->SetTitle( "Counts" );
  
  if(vvar023->size() > 0) { TCanvas *c023 = new TCanvas( "c023", "c023", 600, 600 ); c023->SetLogy(1); c023->SetGridy(1); hvar023->Draw(); }
  if(vvar034->size() > 0) { TCanvas *c034 = new TCanvas( "c034", "c034", 600, 600 ); c034->SetLogy(1); c034->SetGridy(1); hvar034->Draw(); }
  if(vvar045->size() > 0) { TCanvas *c045 = new TCanvas( "c045", "c045", 600, 600 ); c045->SetLogy(1); c045->SetGridy(1); hvar045->Draw(); }
  if(vvar054->size() > 0) { TCanvas *c054 = new TCanvas( "c054", "c054", 600, 600 ); c054->SetLogy(1); c054->SetGridy(1); hvar054->Draw(); }
  if(vvar065->size() > 0) { TCanvas *c065 = new TCanvas( "c065", "c065", 600, 600 ); c065->SetLogy(1); c065->SetGridy(1); hvar065->Draw(); }
  if(vvar100->size() > 0) { TCanvas *c100 = new TCanvas( "c100", "c100", 600, 600 ); c100->SetLogy(1); c100->SetGridy(1); hvar100->Draw(); }
  if(vvar111->size() > 0) { TCanvas *c111 = new TCanvas( "c111", "c111", 600, 600 ); c111->SetLogy(1); c111->SetGridy(1); hvar111->Draw(); }
  if(vvar122->size() > 0) { TCanvas *c122 = new TCanvas( "c122", "c122", 600, 600 ); c122->SetLogy(1); c122->SetGridy(1); hvar122->Draw(); }

//   TCanvas *ctot = new TCanvas( "ctot", "ctot", 1900, 900);
//   ctot->Divide(4,2);
//   ctot->cd(1); ctot->cd(1)->SetLogy(1); ctot->cd(1)->SetGridy(1); hvar023->Draw();
//   ctot->cd(2); ctot->cd(2)->SetLogy(1); ctot->cd(2)->SetGridy(1); hvar034->Draw();
//   ctot->cd(3); ctot->cd(3)->SetLogy(1); ctot->cd(3)->SetGridy(1); hvar045->Draw();
//   ctot->cd(4); ctot->cd(4)->SetLogy(1); ctot->cd(4)->SetGridy(1); hvar054->Draw();
//   ctot->cd(5); ctot->cd(5)->SetLogy(1); ctot->cd(5)->SetGridy(1); hvar065->Draw();
//   ctot->cd(6); ctot->cd(6)->SetLogy(1); ctot->cd(6)->SetGridy(1); hvar100->Draw();
//   ctot->cd(7); ctot->cd(7)->SetLogy(1); ctot->cd(7)->SetGridy(1); hvar111->Draw();
//   ctot->cd(8); ctot->cd(8)->SetLogy(1); ctot->cd(8)->SetGridy(1); hvar122->Draw();
  
  cout << endl;
  cout << "nofld (tot = " << vvar023->size() + vvar034->size() + vvar045->size() + vvar054->size() + vvar065->size() + vvar100->size() + vvar111->size() + vvar122->size() << ")"                                     << endl 
                          << " 23 eV\tn.e- = " << vvar023->size() << " / " << n023 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar023->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar023->GetRMS() ) << " " << unity->GetTitle() << endl
                          << " 34 eV\tn.e- = " << vvar034->size() << " / " << n034 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar034->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar034->GetRMS() ) << " " << unity->GetTitle() << endl
                          << " 45 eV\tn.e- = " << vvar045->size() << " / " << n045 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar045->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar045->GetRMS() ) << " " << unity->GetTitle() << endl
                          << " 54 eV\tn.e- = " << vvar054->size() << " / " << n054 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar054->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar054->GetRMS() ) << " " << unity->GetTitle() << endl
                          << " 65 eV\tn.e- = " << vvar065->size() << " / " << n065 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar065->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar065->GetRMS() ) << " " << unity->GetTitle() << endl
                          << "100 eV\tn.e- = " << vvar100->size() << " / " << n100 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar100->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar100->GetRMS() ) << " " << unity->GetTitle() << endl
                          << "111 eV\tn.e- = " << vvar111->size() << " / " << n111 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar111->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar111->GetRMS() ) << " " << unity->GetTitle() << endl
                          << "122 eV\tn.e- = " << vvar122->size() << " / " << n122 << "\tavg " << vary->GetTitle() << " = " << TString::Format( "%.08f", hvar122->GetMean() ) << " " << unity->GetTitle() << " ; rms = " << TString::Format( "%.08f", hvar122->GetRMS() ) << " " << unity->GetTitle() << endl << endl;

                          
  printf( "E (eV),Fired,Collected,<%s> %s,RMS <%s> %s\n", vary->GetTitle(), unity->GetTitle(), vary->GetTitle(), unity->GetTitle() );
  printf(  "23,%d,%zu,%.03f,%.03f\n", n023, vvar023->size(), avg023, rms023 );
  printf(  "34,%d,%zu,%.03f,%.03f\n", n034, vvar034->size(), avg034, rms034 );
  printf(  "45,%d,%zu,%.03f,%.03f\n", n045, vvar045->size(), avg045, rms045 );
  printf(  "54,%d,%zu,%.03f,%.03f\n", n054, vvar054->size(), avg054, rms054 );
  printf(  "65,%d,%zu,%.03f,%.03f\n", n065, vvar065->size(), avg065, rms065 );
  printf( "100,%d,%zu,%.03f,%.03f\n", n100, vvar100->size(), avg100, rms100 );
  printf( "111,%d,%zu,%.03f,%.03f\n", n111, vvar111->size(), avg111, rms111 );
  printf( "122,%d,%zu,%.03f,%.03f\n", n122, vvar122->size(), avg122, rms122 );
  printf( "%s\n"                   , "https://ozh.github.io/ascii-tables/" );
  
  cout << endl;
  
  cout << TString::Format( "%zu", vvar023->size() ) << endl;
  cout << TString::Format( "%zu", vvar034->size() ) << endl;
  cout << TString::Format( "%zu", vvar045->size() ) << endl;
  cout << TString::Format( "%zu", vvar054->size() ) << endl;
  cout << TString::Format( "%zu", vvar065->size() ) << endl;
  cout << TString::Format( "%zu", vvar100->size() ) << endl;
  cout << TString::Format( "%zu", vvar111->size() ) << endl;
  cout << TString::Format( "%zu", vvar122->size() ) << endl << endl;
                                                    
  cout << TString::Format( "%d", n023 ) << endl;
  cout << TString::Format( "%d", n034 ) << endl;
  cout << TString::Format( "%d", n045 ) << endl;
  cout << TString::Format( "%d", n054 ) << endl;
  cout << TString::Format( "%d", n065 ) << endl;
  cout << TString::Format( "%d", n100 ) << endl;
  cout << TString::Format( "%d", n111 ) << endl;
  cout << TString::Format( "%d", n122 ) << endl << endl;
                                                    
  cout << TString::Format( "%.08f", hvar023->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar034->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar045->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar054->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar065->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar100->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar111->GetMean() ) << endl;
  cout << TString::Format( "%.08f", hvar122->GetMean() ) << endl << endl;
  
  cout << TString::Format( "%.08f", hvar023->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar034->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar045->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar054->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar065->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar100->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar111->GetRMS() ) << endl;
  cout << TString::Format( "%.08f", hvar122->GetRMS() ) << endl << endl;
}
