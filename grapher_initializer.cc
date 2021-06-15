{
  cout << endl;
  cout << "\e[35mLatest files:\e[0m" << endl;
  gROOT->ProcessLine(".! ls -A1tr *.root");
  cout << endl;
  
  string infile;
  cout << "Select infile\e[93m" << endl;
  std::getline( std::cin, infile);
  cout << "\e[0m" << endl;

  TFile *fout = new TFile( infile.c_str() );
  TNamed *varx         = (TNamed*)fout->Get("varx");
  TNamed *vary         = (TNamed*)fout->Get("vary");
  TNamed *varz         = (TNamed*)fout->Get("varz");
  TNamed *unitx        = (TNamed*)fout->Get("unitx");
  TNamed *unity        = (TNamed*)fout->Get("unity");
  TNamed *unitz        = (TNamed*)fout->Get("unitz");
  TNamed *cutscode     = (TNamed*)fout->Get("cutscode");
//   TString treename = TString::Format( "tout_0-%s-%s-%s-%s", varx->GetTitle(), vary->GetTitle(), varz->GetTitle(), cutscode->GetTitle() );    // For older files, use this
  TString treename = TString::Format( "tout_%s-%s-%s-%s", varx->GetTitle(), vary->GetTitle(), varz->GetTitle(), cutscode->GetTitle() );   // NOTE THAT THIS IS RUN_ID AGNOSTIC
  TTree *tout = (TTree*)fout->Get( treename.Data() );
  TNamed *axis_x_title = (TNamed*)fout->Get("axis_x_title");
  TNamed *axis_y_title = (TNamed*)fout->Get("axis_y_title");
  TNamed *axis_z_title = (TNamed*)fout->Get("axis_z_title");
  TNamed *thisto_title = (TNamed*)fout->Get("thisto_title");

  double var[3];
  long long outn;
  tout->SetBranchAddress(varx->GetTitle(), &var[0]);
  tout->SetBranchAddress(vary->GetTitle(), &var[1]);
  tout->SetBranchAddress(varz->GetTitle(), &var[2]);
  tout->SetBranchAddress("outn", &outn);

  cout << "TFile pointer: fout" << endl;
  cout << "TTree pointer: tout" << endl;

  cout << endl;
  cout << "Found following branches -> associated variable" << endl;
  cout << " - " << varx->GetTitle() << "\t\t-> var[0]" << endl;
  cout << " - " << vary->GetTitle() << "\t\t-> var[1]" << endl;
  cout << " - " << varz->GetTitle() << "\t\t-> var[2]" << endl;
  cout << endl;
  cout << "Entry number accessible via \'outn\' variable" << endl;
  cout << endl;
  cout << "To quickly draw the original plot:" << endl << endl;
  cout << "TCanvas *ccout = new TCanvas(\"ccout\",\"ccout\",800,800);TPad *pout = new TPad(\"pout\",\"\",0,0,1,1);pout->Draw();pout->SetRightMargin(0.15);pout->SetLeftMargin(0.12);pout->cd();" << endl;
  cout << TString::Format("tout->Draw(TString::Format(\"%s\", \"%s\",\"%s\",\"%s\"),\"\",\"COLZ\");","%s:%s:%s >>thisto", vary->GetTitle(), varx->GetTitle(), varz->GetTitle()) << endl;
  cout << "thisto->GetXaxis()->SetTitle(axis_x_title->GetTitle()); thisto->GetXaxis()->SetLabelSize(0.02);" << endl;
  cout << "thisto->GetYaxis()->SetTitle(axis_y_title->GetTitle()); thisto->GetYaxis()->SetLabelSize(0.02);" << endl;
  cout << "thisto->GetZaxis()->SetTitle(axis_z_title->GetTitle()); thisto->GetZaxis()->SetLabelSize(0.02); thisto->GetZaxis()->SetTitleOffset(1.3);" << endl;
  cout << "thisto            ->SetTitle(thisto_title->GetTitle());" << endl;
  cout << endl;
//   gStyle->SetOptStat(0); // Draw histogram without statbox
  gStyle->SetPalette(kRainBow);
//   gStyle->SetPalette(kGreenRedViolet);
//   gStyle->SetPalette(kInvertedDarkBodyRadiator);
}
