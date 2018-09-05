//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TLorentzVector.h>

//#include "rp_aperture_config.h"

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include "TRandom3.h"
#include <fstream>

#include "statusbar.h"

//#include "pid_tools.h"
//#include "shared_fill_info.h"

#include <vector> 

#define MASS_MU 0.1057 // GeV
#define MASS_E  0.000511 // GeV
#define MASS_P  0.938272029 // GeV
#define pi 3.14159265359
#define ECM 13000.0 // GeV

const float ns_to_s_ = 1e-9;
const float m_to_cm_ = 1e2;

using namespace std;

bool showBar = false;

TFile* fOut;
TTree* trout;

//Particle flow meaning
//0 unknown
//1 charged (trk)
//2 electron (ch+ECAL)
//3 muon (trk+muon chamber)
//4 gamma (ecal and no trk)
//5 neutral hadron (HCal and no tracks)
//6 HF hadronics
//7 HF electromagnetic
///eos/totem/data/cmstotem/2015/90m/Merged/4510/9985

//#define preTS2 0

int preTS2;
int getXangle(int run,int lumi, const char* filename);

struct FillInfo
{
  unsigned int fillNumber;
  unsigned int runMin, runMax;
  string alignmentTag;

  FillInfo(unsigned int _fn=0, unsigned int _rmin=0, unsigned int _rmax=0, const string &at=""):fillNumber(_fn), runMin(_rmin), runMax(_rmax), alignmentTag(at){}
};

//----------------------------------------------------------------------------------------------------

struct FillInfoCollection : public vector<FillInfo> 
{
  FillInfo FindByrun(unsigned int run) const
  {
    for (const auto it : *this)
    {
      if (it.runMin <= run && it.runMax >= run)
	return it;
    }
    for (const auto it : *this)
    {
      if (it.runMin <= 300000  && it.runMax >= 300000)
      {
	return it;
      }
    }
    printf("THIS INSTEAD SHOULD NOT HAPPEN. EXIT\n");
    exit(1);
    return FillInfo();
  }
};

FillInfoCollection fillInfoCollection;

void InitFillInfoCollection()
{
  fillInfoCollection.push_back(FillInfo(4947, 273725, 273730, "period1_physics_margin/fill_4947/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4953, 274094, 274094, "period1_physics_margin/fill_4953/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4961, 274198, 274200, "period1_physics_margin/fill_4961/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274240, 274241, "period1_physics_margin/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274244, 274244, "period1_physics/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4976, 274282, 274286, "period1_physics_margin/fill_4976/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4985, 274387, 274388, "period1_physics/fill_4985/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4988, 274420, 274422, "period1_physics/fill_4988/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4990, 274440, 274443, "period1_physics/fill_4990/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5005, 274954, 274959, "period1_physics/fill_5005/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5013, 274966, 274971, "period1_physics/fill_5013/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5017, 274998, 275001, "period1_physics/fill_5017/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5020, 275059, 275074, "period1_physics/fill_5020/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5021, 275124, 275125, "period1_physics/fill_5021/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5024, 275282, 275293, "period1_physics/fill_5024/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5026, 275309, 275311, "period1_physics/fill_5026/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5027, 275319, 275338, "period1_physics/fill_5027/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5028, 275344, 275345, "period1_physics/fill_5028/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5029, 275370, 275371, "period1_physics/fill_5029/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5030, 275375, 275376, "period1_physics/fill_5030/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5038, 275656, 275659, "period1_physics/fill_5038/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5043, 275757, 275783, "period1_physics/fill_5043/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5045, 275828, 275847, "period1_physics/fill_5045/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5048, 275886, 275890, "period1_physics/fill_5048/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5052, 275911, 275931, "period1_physics/fill_5052/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5261, 279760, 279767, "period1_physics/fill_5261/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5264, 279794, 279794, "period1_physics/fill_5264/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5265, 279823, 279823, "period1_physics/fill_5265/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5266, 279841, 279841, "period1_physics/fill_5266/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5267, 279844, 279865, "period1_physics/fill_5267/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5274, 279931, 279931, "period1_physics/fill_5274/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5275, 279966, 279966, "period1_physics/fill_5275/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5276, 279975, 279975, "period1_physics/fill_5276/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5277, 279993, 280024, "period1_physics/fill_5277/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5279, 280187, 280194, "period1_physics/fill_5279/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5287, 280327, 280364, "period1_physics/fill_5287/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5288, 280383, 280385, "period1_physics/fill_5288/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(8000, 300000, 300001, "period1_physics/fill_8000/2017_01_17"));

  fillInfoCollection.push_back(FillInfo(5393, 282730, 282735, "TODO"));
  fillInfoCollection.push_back(FillInfo(5401, 282920, 282924, "TODO"));
  fillInfoCollection.push_back(FillInfo(5405, 283040, 283043, "TODO"));
  fillInfoCollection.push_back(FillInfo(5406, 283049, 283067, "TODO"));
  fillInfoCollection.push_back(FillInfo(5418, 283305, 283308, "TODO"));
  fillInfoCollection.push_back(FillInfo(5421, 283353, 283359, "TODO"));
  fillInfoCollection.push_back(FillInfo(5423, 283407, 283416, "TODO"));
  fillInfoCollection.push_back(FillInfo(5424, 283453, 283453, "TODO"));
  fillInfoCollection.push_back(FillInfo(5427, 283478, 283481, "TODO"));
  fillInfoCollection.push_back(FillInfo(5433, 283548, 283560, "TODO"));
  fillInfoCollection.push_back(FillInfo(5437, 283672, 283685, "TODO"));
  fillInfoCollection.push_back(FillInfo(5439, 283820, 283835, "TODO"));
  fillInfoCollection.push_back(FillInfo(5441, 283863, 283865, "TODO"));
  fillInfoCollection.push_back(FillInfo(5442, 283876, 283878, "TODO"));
  fillInfoCollection.push_back(FillInfo(5443, 283884, 283885, "TODO"));
  fillInfoCollection.push_back(FillInfo(5446, 283933, 283934, "TODO"));
  fillInfoCollection.push_back(FillInfo(5448, 283946, 283964, "TODO"));
  fillInfoCollection.push_back(FillInfo(5450, 284006, 284014, "TODO"));
  fillInfoCollection.push_back(FillInfo(5451, 284025, 284044, "TODO"));
}

struct AlignmentResult
{
  double sh_x, sh_x_unc;		// mm
  double sh_y, sh_y_unc;		// mm

  AlignmentResult(double _sh_x=0., double _sh_x_unc=0., double _sh_y=0., double _sh_y_unc=0.) :
    sh_x(_sh_x), sh_x_unc(_sh_x_unc), sh_y(_sh_y), sh_y_unc(_sh_y_unc)
  {
  }

  void Write(FILE *f) const
  {
    fprintf(f, "sh_x=%.3f,sh_x_unc=%.3f,sh_y=%.3f,sh_y_unc=%.3f\n", sh_x, sh_x_unc, sh_y, sh_y_unc);
  }
};

struct AlignmentResults : public map<unsigned int, AlignmentResult>
{
  void Write(FILE *f) const
  {
    for (auto &p : *this)
    {
      fprintf(f, "id=%u,", p.first);
      p.second.Write(f);
    }
  }

  int Add(char *line)
  {
    bool idSet = false;
    unsigned int id = 0;
    AlignmentResult result;

    // loop over entries separated by ","
    char *p = strtok(line, ",");
    while (p != NULL)
    {
      // isolate key and value strings
      char *pe = strstr(p, "=");
      if (pe == NULL)
      {
	printf("ERROR in AlignmentResults::Add > entry missing = sign: %s.\n", p);
	return 2;
      }

      char *s_key = p;

      p = strtok(NULL, ",");
      *pe = 0;
      char *s_val = pe+1;

      // interprete keys
      if (strcmp(s_key, "id") == 0)
      {
	idSet = true;
	id = atoi(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_x") == 0)
      {
	result.sh_x = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_x_unc") == 0)
      {
	result.sh_x_unc = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_y") == 0)
      {
	result.sh_y = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_y_unc") == 0)
      {
	result.sh_y_unc = atof(s_val);
	continue;
      }
      printf("ERROR in AlignmentResults::Add > unknown key: %s.\n", s_key);
      return 3;
    }

    if (!idSet)
    {
      printf("ERROR in AlignmentResults::Add > id not set on the following line:\n%s.\n", line);
      return 4;
    }
    insert({id, result});
    return 0;
  }

  double GETX(int rpnumber) const
  {
    auto ait = find(rpnumber);
    double toret=0;
    if (ait == end())
    {
      printf("ERROR in AlignmentResults::Apply > no alignment data for RP %u.\n", rpnumber);
      exit(1);
    }
    else
    {
      toret=ait->second.sh_x;
    }
    return toret;
  }

};

struct AlignmentResultsCollection : public map<string, AlignmentResults>
{
  int Write(const string &fn) const
  {
    FILE *f = fopen(fn.c_str(), "w");
    if (!f)
      return -1;
    Write(f);
    fclose(f);
    return 0;
  }

  void Write(FILE *f) const
  {
    for (auto &p : *this)
    {
      fprintf(f, "\n[%s]\n", p.first.c_str());
      p.second.Write(f);
    }
  }

  int Load(const string &fn)
  {
    FILE *f = fopen(fn.c_str(), "r");
    if (!f)
      return -1;
    return Load(f);
    fclose(f);
  }

  int Load(FILE *f)
  {
    string label = "unknown";
    AlignmentResults block;

    while (!feof(f))
    {
      char line[300];
      char *success = fgets(line, 300, f);
      if (success == NULL)
	break;
      if (line[0] == '\n')
	continue;
      if (line[0] == '[')
      {
	if (block.size() > 0)
	  insert({label, block});
	block.clear();
	char *end = strstr(line, "]");
	if (end == NULL)
	{
	  printf("ERROR in AlignmentResultsCollection::Load > missing closing ].\n");
	  return 1;
	}
	*end = 0;
	label = line+1;
	continue;
      }
      if (block.Add(line) != 0)
	return 2;
    }
    if (block.size() > 0)
      insert({label, block});
    return 0;
  }
};

TRandom3 *Rgenerator;

double Rndm(){
  return ( double(rand())/RAND_MAX );
}

class Particle{

  public:
    Particle();
    Particle(double);
    double   px, py, pz, E, m, pAbs, pt, p[4]; 
    void     p4(double, double, double, double);
    void     boost(Particle);
    void     boost2(double px_, double py_, double pz_, double E_);
    void     twoBodyDecay(Particle&, Particle&);
    void     print();
};

Particle::Particle(){
  px = py = pz = E = m = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}
Particle::Particle(double mass){
  m = mass;
  px = py = pz = E = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}

//
//*** Sets components of 4-momentum vector -------------------------------------
//
void Particle::p4(double momx, double momy, double momz, double energy){
  // components of 4-momenta
  px = p[0] = momx;
  py = p[1] = momy;
  pz = p[2] = momz;
  E  = p[3] = energy;
  // transverse momentum and the magnitude of the space momentum
  pt      = sqrt(momx*momx + momy*momy);
  pAbs    = sqrt(momx*momx + momy*momy + momz*momz);
}

//
//*** Prints 4-vector ----------------------------------------------------------
//
void Particle::print(){
  cout << "(" << p[0] <<",\t" << p[1] <<",\t"<< p[2] <<",\t"<< p[3] << ")" << endl;
}

//
//*** Evaluates 4-vector of decay product in the rest frame of parent. ---------
//
void Particle::twoBodyDecay(Particle& prod1, Particle& prod2) {

  double m1 = prod1.m;
  double m2 = prod2.m;

  // check if particle m can decay
  if( m < m1+m2 ){
    cout << "Particle::twoBodyDecay  parent mass is less than sum of products." 
      << endl;
    prod1.p4 ( 0., 0., 0., m1);
    prod2.p4 ( 0., 0., 0., m2);
    return;
  }

  // CM energies and momentum
  double e1 = (m*m + m1*m1 - m2*m2) / (2.0*m);
  double e2 = (m*m - m1*m1 + m2*m2) / (2.0*m);
  double P  = sqrt(e1*e1 - m1*m1);

  // Isotropic random angles
  double theta = acos( 2.0*Rndm() - 1.0 );
  double phi   = 2.0*M_PI *Rndm();

  double pX = P*sin(theta)*cos(phi);
  double pY = P*sin(theta)*sin(phi);
  double pZ = P*cos(theta);

  // set the 4-momenta
  prod1.p4(  pX,  pY,  pZ, e1 );
  prod2.p4( -pX, -pY, -pZ, e2 );
}

//
//*** Lorentz Boost  -----------------------------------------------------------
//
void Particle::boost(Particle parent){

  // beta and gamma values
  double betax = parent.px / parent.E;
  double betay = parent.py / parent.E;
  double betaz = parent.pz / parent.E;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}

void Particle::boost2(double px_, double py_, double pz_, double E_){

  // beta and gamma values
  double betax = px_ / E_;
  double betay = py_ / E_;
  double betaz = pz_ / E_;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}

void DarkMatterSearchLeadingMuonsFinn(string year, vector<string> const& fileNames, string const& outputFileNameTree = "TTreeOutput.root", const Int_t nevt_max = 100)
{ 

  bool isMC  = false;
  bool verbose = false;

  int nMuonCand=0;
  int nLocalProtCand=0;

  string treeName = "SlimmedNtuple";
  map<string,TH1F*> histosTH1F;

  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;

  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;
  vector<TString>* vfiles = new vector<TString>; 

  for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file)
  {
    vfiles->push_back( fileNames[idx_file] );
    std::cout<<"Loaded Input file: "<< fileNames[idx_file]<<std::endl;
  }

  // Declaration of tree and its branches variables
  TTree* tree = NULL;

  std::vector<float> *muon_pt = 0;
  std::vector<float> *muon_eta = 0;
  std::vector<float> *muon_phi = 0;
  std::vector<float> *muon_px = 0;
  std::vector<float> *muon_py = 0;
  std::vector<float> *muon_pz = 0;
  std::vector<float> *muon_e = 0;
  std::vector<float> *muon_charge = 0;
  std::vector<float> *muon_iso = 0;

  std::vector<float> *rp_tracks_x = 0;
  std::vector<float> *rp_tracks_y = 0;
  std::vector<uint32_t> *rp_tracks_detId = 0;

  TBranch *bmuon_pt = 0;
  TBranch *bmuon_eta = 0;
  TBranch *bmuon_phi = 0;
  TBranch *bmuon_px = 0;
  TBranch *bmuon_py = 0;
  TBranch *bmuon_pz = 0;
  TBranch *bmuon_e = 0;
  TBranch *bmuon_charge = 0;
  TBranch *bmuon_iso = 0;

  TBranch *brp_tracks_x = 0;
  TBranch *brp_tracks_y = 0;
  TBranch *brp_tracks_detId = 0;

  TBranch *bispps = 0;
  TBranch *brun = 0;
  TBranch *blumiblock = 0;
  TBranch *bevent = 0;

  bool ispps=false;
  int run=0;
  int lumiblock=0;
  Long64_t event=0;

  double pixelcoordX[2][150];
  double pixelcoordY[2][150];

  int NpixelLL=0;
  int NpixelRR=0;
  int pixelLL=0;
  int pixelRR=0;

  int i_tot = 0 , nevt_tot = 0; int indicatorefile=-1;

  for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max_corr ; ++itfiles){

    std::cout<<"File loop:"<<indicatorefile<<std::endl;
    TFile* file = TFile::Open(*itfiles,"READ");
    tree = (TTree*) file->Get("demo/SlimmedNtuple");
    indicatorefile++;

    int nev = int(tree->GetEntries());
    nevt_tot += nev;
    cout <<"The current file ("<< indicatorefile <<") has " << nev << " entries : "<< endl << *itfiles << endl;

    //  to TTree ----------------------------------------------------------------------
    tree->SetBranchAddress("run",&run, &brun);
    tree->SetBranchAddress("event", &event, &bevent);
    tree->SetBranchAddress("lumiblock", &lumiblock, &blumiblock);
    tree->SetBranchAddress("muon_pt",&muon_pt, &bmuon_pt);
    tree->SetBranchAddress("muon_px",&muon_px); 
    tree->SetBranchAddress("muon_py",&muon_py, &bmuon_py);    
    tree->SetBranchAddress("muon_pz",&muon_pz, &bmuon_pz);    
    tree->SetBranchAddress("muon_eta",&muon_eta, &bmuon_eta);
    tree->SetBranchAddress("muon_phi",&muon_phi, &bmuon_phi);
    tree->SetBranchAddress("muon_e",&muon_e, &bmuon_e);
    tree->SetBranchAddress("muon_charge",&muon_charge, &bmuon_charge);    
    tree->SetBranchAddress("muon_iso",&muon_iso, &bmuon_iso);
    tree->SetBranchAddress("rp_tracks_x", &rp_tracks_x, &brp_tracks_x);
    tree->SetBranchAddress("rp_tracks_y", &rp_tracks_y, &brp_tracks_y);
    tree->SetBranchAddress("rp_tracks_detId", &rp_tracks_detId, &brp_tracks_detId);
    tree->SetBranchAddress("ispps",&ispps, &bispps);

    //  to TTree ----------------------------------------------------------------------
    //  load alignment collection
    AlignmentResultsCollection alignmentCollection;
    alignmentCollection.Load("./analysis_ctpps_alignment/shared_alignment/collect_alignments.out");
    AlignmentResults *alignments = NULL;
    InitFillInfoCollection();

    double px_1= 0.;
    double py_1 = 0.;
    double pz_1= 0.;
    double pE_1 = 0.;
    double px_2 = 0.;
    double py_2 = 0.;
    double pz_2  = 0.;
    double pE_2= 0.;
    double px_3 = 0.;
    double py_3 = 0.;
    double pz_3  = 0.;
    double pE_3= 0.;
    double px_4 = 0.;
    double py_4 = 0.;
    double pz_4  = 0.;
    double pE_4= 0.;
    int chargeVal=0;

    //calculate total PT of the system
    double LeadingMuonPt=0.;
    double SecondMuonPt=0.;
    int LeadingMuonIso=0;
    int SecondMuonIso=0;
    double Muons_pxtot=0.;
    double Muons_pytot=0.;
    double Muons_pztot=0.;
    double Muons_pttot=0.;
    double pairMass=0.;
    double pairPhi=0.;
    double pairEta=0.;
    double pairpT=0.;
    double acoplan=0.;
    double Muons_Etot=0.;
    double proton1_pz=0.;
    double proton2_pz=0.;
    double protonmass=0.;
    double missingmass=0.;
    double missingmassLorentz=0.;
    double xi1=0.;
    double xi2=0.;
    int angle = -1;
    double pixelXarm0=0.;
    double pixelXarm1=0.;
    double pixelYarm0=0.;
    double pixelYarm1=0.;

    int therunnumber=0;
    int thelumisection=0;
    int theevent=0;
    bool theispps=false;

    TString ttreefilename = outputFileNameTree;

    fOut = new TFile(ttreefilename, "RECREATE");
    fOut->cd();

    trout = new TTree("Events", "Events");
    trout->Branch("ispps",&theispps,"ispps/B");
    trout->Branch("run",&therunnumber,"run/I");
    trout->Branch("lumiblock",&thelumisection,"lumiblock/I");
    trout->Branch("event",&theevent,"event/I");
    trout->Branch("missingmass",&missingmass,"missingmass/D");
    trout->Branch("missingmassLorentz",&missingmassLorentz,"missingmassLorentz/D");
    trout->Branch("nMuonCand",&nMuonCand,"nMuonCand/I");
    trout->Branch("LeadingMuonPt",&LeadingMuonPt,"LeadingMuonPt/D");
    trout->Branch("LeadingMuonIso",&LeadingMuonIso,"LeadingMuonIso/D");
    trout->Branch("SecondMuonPt",&SecondMuonPt,"SecondMuonPt/D");
    trout->Branch("SecondMuonIso",&SecondMuonIso,"SecondMuonIso/D");
    trout->Branch("acoplan",&acoplan,"acoplan/D");
    trout->Branch("chargeVal",&chargeVal,"chargeVal/I");
    trout->Branch("nLocalProtCand",&nLocalProtCand,"nLocalProtCand/I");
    trout->Branch("proton1_pz",&proton1_pz,"proton1_pz/D");
    trout->Branch("proton2_pz",&proton2_pz,"proton2_pz/D");
    trout->Branch("xi1",&xi1,"xi1/D");
    trout->Branch("xi2",&xi2,"xi2/D");
    trout->Branch("protonmass",&protonmass,"protonmass/D");
    trout->Branch("pairMass",&pairMass,"pairMass/D");
    trout->Branch("pairEta",&pairEta,"pairEta/D");
    trout->Branch("pairPhi",&pairPhi,"pairPhi/D");
    trout->Branch("pairpT",&pairpT,"pairpT/D");
    trout->Branch("pixelXarm0",&pixelXarm0,"pixelXarm0/D");
    trout->Branch("pixelXarm1",&pixelXarm1,"pixelXarm1/D");
    trout->Branch("pixelYarm0",&pixelYarm0,"pixelYarm0/D");
    trout->Branch("pixelYarm1",&pixelYarm1,"pixelYarm1/D");
    trout->Branch("angle",&angle,"angle/I");

    bool debug = false;

    for(int unsigned i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot)
    {

      if(showBar){
	if (i_evt==0) {
	  std::cout << "" << std::endl;
	  std::cout<< "Status Bar" << std::endl;
	}
	loadBar(i_evt, nev);
      }

      tree->GetEntry(i_evt);

      for(int i=0;i<150;i++){
	pixelcoordX[0][i]=0;
	pixelcoordX[1][i]=0;
	pixelcoordY[0][i]=0;
	pixelcoordY[1][i]=0;
      }

      chargeVal=0;
      NpixelLL=0;
      NpixelRR=0;
      pixelLL=0;
      pixelRR=0;
      Muons_pxtot=0.;
      Muons_pytot=0.;
      Muons_pztot=0.;
      Muons_pttot=0.;
      LeadingMuonPt=0;
      SecondMuonPt=0;
      LeadingMuonIso=0;
      SecondMuonIso=0;
      pairMass=0.;
      pairPhi=0.;
      pairEta=0.;
      pairpT=0.;
      acoplan=0.;
      Muons_Etot=0.;
      proton1_pz=0.;
      proton2_pz=0.;
      protonmass=0.;
      missingmass=0.;
      missingmassLorentz=0.;
      xi1=0;
      xi2=0;
      angle = -1;
      pixelXarm0=0;
      pixelXarm1=0;
      pixelYarm0=0;
      pixelYarm1=0;

      nLocalProtCand = 0;
      NpixelLL=0;
      NpixelRR=0;

      for(std::vector<int>::size_type s = 0; s != rp_tracks_detId->size(); s++) {

	static const uint32_t startbit = 24;
	static const uint32_t mask = 0x1;
	uint32_t pixelID = (uint32_t(rp_tracks_detId->at(s))>>startbit) & mask;

	static const uint32_t startStationBit = 22;
	static const uint32_t maskStation = 0x3;
	uint32_t stationID = (uint32_t(rp_tracks_detId->at(s))>>startStationBit) & maskStation;

	static const uint32_t startRPBit = 19;
	static const uint32_t maskRP = 0x7;
	uint32_t rpID = (uint32_t(rp_tracks_detId->at(s))>>startRPBit) & maskRP;

	if(pixelID==0 && (stationID==2) && rpID==3){
	  pixelcoordX[0][NpixelLL]=rp_tracks_x->at(s); 
	  pixelcoordY[0][NpixelLL]=rp_tracks_y->at(s); 
	  if(debug){
	    std::cout << "Map(X): " <<  pixelcoordX[0][NpixelLL] << ", Track X [mm]: " << rp_tracks_x->at(s) << ", Arm: " << pixelID << ", stationID: "<< stationID << ", rpID: " << rpID << std::endl;
	    std::cout << "Map(Y): " <<  pixelcoordY[0][NpixelLL] << ", Track Y [mm]: " << rp_tracks_y->at(s) << ", Arm: " << pixelID << ", stationID: "<< stationID << ", rpID: " << rpID << std::endl;
	  }
	  NpixelLL++;
	  pixelLL=s;
	}
	if(pixelID==1 && (stationID==2) && rpID==3){
	  pixelcoordX[1][NpixelRR]=rp_tracks_x->at(s); 
	  pixelcoordY[1][NpixelRR]=rp_tracks_y->at(s); 
	  if(debug){
	    std::cout << "Map(X): " <<  pixelcoordX[1][NpixelRR] << ", Track X [mm]: " << rp_tracks_x->at(s) << ", Arm: " << pixelID << ", stationID: "<< stationID << ", rpID: " << rpID << std::endl;
	    std::cout << "Map(Y): " <<  pixelcoordY[1][NpixelRR] << ", Track Y [mm]: " << rp_tracks_y->at(s) << ", Arm: " << pixelID << ", stationID: "<< stationID << ", rpID: " << rpID << std::endl;
	  }
	  NpixelRR++;
	  pixelRR=s;
	}
      }
      nLocalProtCand=NpixelLL+NpixelRR;

      if(nLocalProtCand>0)
      {

	if(nLocalProtCand==2)
	{

	  //calculate total PT of the system
	  LeadingMuonPt=0.;
	  SecondMuonPt=0.;
	  LeadingMuonIso=0.;
	  SecondMuonIso=0.;
	  Muons_pxtot=0.;
	  Muons_pytot=0.;
	  Muons_pztot=0.;
	  Muons_pttot=0.;
	  Muons_Etot=0.;

	  theispps=false;
	  therunnumber=0;
	  thelumisection=0;
	  theevent=0;

	  TLorentzVector p1;
	  TLorentzVector p2;

	  double dispersionL=-9999;
	  double dispersionR=-9999;

	  double Dispersion[2][200];
	  Dispersion[0][120]=9.145*10;
	  Dispersion[0][130]=8.591*10;
	  Dispersion[0][140]=8.226*10;
	  Dispersion[1][120]=7.291*10;
	  Dispersion[1][130]=6.621*10;
	  Dispersion[1][140]=6.191*10;

	  double relcorr=0;
	  preTS2 ? relcorr=(42.05):relcorr=(42.2);

	  if(preTS2==1){
	    angle=getXangle(run,  lumiblock , "./inp/new/xangle_tillTS2_STABLEBEAMS_CLEANUP.csv");
	  }

	  if(preTS2==0) angle=getXangle(run,  lumiblock , "./inp/xangle_afterTS2_cleanup.csv");
	  if(angle==120 || angle ==130 || angle==140){dispersionL=Dispersion[0][angle];dispersionR=Dispersion[1][angle];}
	  if(angle !=120 && angle != 130 && angle != 140){
	    dispersionL = Dispersion[0][140] - 4.6*((angle-140.)/10.) ;
	    dispersionR = Dispersion[1][140] - 4.6*((angle-140.)/10.) ;
	  }

	  if(angle>0){
	    if(NpixelLL==1){
	      for(int k=0;k<NpixelLL;k++){
		xi1=(pixelcoordX[0][k]-relcorr)/dispersionL;
		proton1_pz=(1-(xi1))*(ECM/2.);
		pixelXarm0=pixelcoordX[0][k];
		pixelYarm0=pixelcoordY[0][k];
		p1.SetPxPyPzE(0.,0., ECM*xi1/2., ECM*xi1/2.);
	      }
	    }else{
	      pixelXarm0=0.;
	      pixelYarm0=0.;
	    }
	    if(NpixelRR==1){
	      for(int k=0;k<NpixelRR;k++){
		xi2=(pixelcoordX[1][k]-relcorr)/dispersionR;
		proton2_pz=-(1-(xi2))*(ECM/2.);
		pixelXarm1=pixelcoordX[1][k];
		pixelYarm1=pixelcoordY[1][k];
		p2.SetPxPyPzE(0.,0., -ECM*xi2/2., ECM*xi2/2.);
	      }
	    }else{
	      pixelXarm1=0.;
	      pixelYarm1=0.;
	    }
	    protonmass=ECM*sqrt((  (xi1)*(xi2) ));
	  }

	  TLorentzVector muon0;
	  TLorentzVector muon1;
	  TLorentzVector muonpair;
	  TLorentzVector ppZsystem;

	  if(NpixelRR==1 && NpixelLL==1){

	    nMuonCand=muon_pt->size();
	    if(nMuonCand>1){

	      theispps=ispps;
	      therunnumber=run;
	      thelumisection=lumiblock;
	      theevent=event;

	      LeadingMuonPt=muon_pt->at(0);
	      SecondMuonPt=muon_pt->at(1);
	      LeadingMuonIso=muon_iso->at(0);
	      SecondMuonIso=muon_iso->at(1);
	      muon0.SetPtEtaPhiE(muon_pt->at(0), muon_eta->at(0), muon_phi->at(0), muon_e->at(0));
	      muon1.SetPtEtaPhiE(muon_pt->at(1), muon_eta->at(1), muon_phi->at(1), muon_e->at(1));
	      muonpair = muon0+muon1;

	      Muons_Etot=sqrt(muon_e->at(0)*muon_e->at(0)+MASS_MU*MASS_MU) + sqrt(muon_e->at(1)*muon_e->at(1)+MASS_MU*MASS_MU);
	      Muons_pxtot = muon0.Px() + muon1.Px();
	      Muons_pytot = muon0.Py() + muon1.Py();
	      Muons_pztot = muon0.Pz() + muon1.Pz();
	      chargeVal=(int)muon_charge->at(0)+(int)muon_charge->at(1);

	      pairMass=muonpair.M();
	      pairEta=muonpair.Eta();
	      pairPhi=muonpair.Phi();
	      pairpT=muonpair.Pt();

	      double dphi = fabs(muon0.Phi()-muon1.Phi());
	      dphi = (dphi<pi) ? dphi : 2.*pi-dphi;
	      acoplan = 1 - fabs(dphi)/pi;

	      missingmass=sqrt((ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))*(ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))-(Muons_pxtot)*(Muons_pxtot)-(Muons_pytot)*(Muons_pytot)- (Muons_pztot+proton1_pz+proton2_pz)*(Muons_pztot+proton1_pz+proton2_pz));

	      ppZsystem = p1 + p2 - muonpair;
	      missingmassLorentz = ppZsystem.M(); 

	      fOut->cd();
	      trout->Fill();

	    }
	  }
	}
      }
    }

    file->Close();
    trout->Write();

  } // End of loop over files

  fOut->cd();
  fOut->Close();

}

//----------------------------------------------------------------------------------------------------


int main(int argc, char** argv)
{

  std::cout << "Main is called with " << argc << " arguments:" <<" Expected: int maxevent, string extrastr, era B,C,G"<< std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout << argv[i] << std::endl;
  }
  int maxevent=atoi(argv[1]);
  string extrastr=argv[2];
  string era_str=argv[3];
  string year_str=argv[4];

  std::cout<<"CALLING MAIN.."<<std::endl;
  std::vector<string> allnamesWithPrefixFinal;

  if(era_str=="B"||era_str=="C"||era_str=="D") preTS2 = 1;
  if(era_str=="E"||era_str=="F") preTS2 = 0;
  if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/dmf/Finn/DoubleMuC_SlimmedNtupleMerged.root");
  if(era_str=="D")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/dmf/Finn/DoubleMuD_SlimmedNtupleMerged.root");
  if(era_str=="E")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/dmf/Finn/DoubleMuE_SlimmedNtupleMerged.root");
  if(era_str=="F")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/dmf/Finn/DoubleMuF_SlimmedNtupleMerged.root");

  string filenameTree="/afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/DarkMatterSearch/Analyser/TTree"+era_str+"_"+extrastr+".root"; 
  DarkMatterSearchLeadingMuonsFinn(year_str, allnamesWithPrefixFinal, filenameTree, maxevent);

  return 0;
}

int getXangle(int run, int lumi, const char* filename)
{
  TString drun;
  TString dfill;
  TString dlumi;
  TString temp;
  int Xangle=-1;
  TString runs;runs.Form("%d", run);
  TString lumis;lumis.Form("%d", lumi);
  ifstream F;
  F.open((const char*)filename);
  //F.clear();
  //F.seekg(0, ios::beg);
  int counter=0;
  if(F){
    while (!F.eof())
    {
      F>>drun;
      F>>dfill;
      F>>dlumi;
      F>>temp;
      F>>Xangle;
      //if(runfill==runlumi && )
      if( runs == drun &&  lumis==dlumi )
      {
	//  cout << "Read: "<< run<<" " << lumi<<"    -----> "<<Xangle << endl;   
	break;
	//return Xangle;
      }
      //    prevrunfill=runfill;
      //    prevXangle=Xangle;
    }
  }//endif
  else cout << "[!] gerXangle():Error reading from file" << endl;
  if(F.eof())
  {
    cout << "[getXangle() warning:] No Xangle data for this run found!" <<endl;
    F.close();
    return -1;
  }
  else  {F.close(); //cout << "~~~~~~~ xangle: "<< Xangle << endl;
    return Xangle;
  }
}
