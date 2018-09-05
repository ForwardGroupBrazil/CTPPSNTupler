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
  FillInfo FindByRun(unsigned int run) const
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

void DarkMatterSearchLeadingMuonsJets(string year, vector<string> const& fileNames, string const& outputFileNameTree = "TTreeOutput.root",const Int_t nevt_max = 100)
{ 

  bool isMC  = false;
  bool verbose = false;

  UInt_t nMuonCand;
  int nJetCand;
  int nLocalProtCand;

  string treeName;
  treeName = (!isMC) ? "ggll_aod/ntp1" : "ggll_aod/ntp1";
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

  unsigned short LocalProtCand_arm[100];
  unsigned short LocalProtCand_side[100];
  int HLT_Accept[10];
  double Pair_mass[100];
  double Pair_eta[100];
  double Pair_phi[100];
  double Pair_dphi[100];
  double Pair_pt[100];
  double ClosestExtraTrack_vtxdxyz[100];

  double ExtraTrack_z[6500];
  double ExtraTrack_px[6500];
  double ExtraTrack_py[6500];
  double ExtraTrack_pz[6500];

  double MuonCand_vtxz[100];
  double JetCand_vz[100];
  double KalmanVertexCand_z[100];

  std::vector<double> VExtraTrack_z;
  std::vector<double> VExtraTrack_px;
  std::vector<double> VExtraTrack_py;
  std::vector<double> VExtraTrack_pz;

  double JetCand_eta[100];
  double JetCand_phi[100];
  double JetCand_pt[100];
  double JetCand_e[100];

  double MuonCand_pt[100];
  int MuonCand_istight[100];
  int MuonCand_isglobal[100];
  double MuonCand_eta[100];
  double MuonCand_phi[100];
  double MuonCand_e[100];
  int MuonCand_charge[100];
  double Etmiss;

  double LocalProtCand_x[100];
  double LocalProtCand_z[100];

  UInt_t nRPPixelTrack=0;
  double RPPixelTrack_x[100];
  double RPPixelTrack_y[100];
  int RPPixelTrack_arm[100];

  UInt_t Run=0;
  UInt_t LumiSection=0;
  UInt_t EventNum=0;

  double pixelcoordX[2][150];
  double pixelcoordY[2][150];

  double RPTimingTrack_Time[800];
  int RPTimingTrack_arm[800];
  UInt_t nRPTimingTrack = 0;

  int i_tot = 0 , nevt_tot = 0; int indicatorefile=-1;

  for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max_corr ; ++itfiles){
    std::cout<<"File loop:"<<indicatorefile<<std::endl;
    TFile* file = TFile::Open(*itfiles,"READ");
    indicatorefile++;

    // Access TTree from current file
    tree = (TTree*) file->Get( treeName.c_str() );

    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file ("<< indicatorefile <<") has " << nev << " entries : "<< endl << *itfiles << endl;
    std::vector<unsigned int> rp_list;

    //  to TTree ----------------------------------------------------------------------
    tree->SetBranchAddress("nMuonCand",&nMuonCand);
    tree->SetBranchAddress("nJetCand",&nJetCand);
    tree->SetBranchAddress("Run",&Run);  

    int NpixelLL;
    int NpixelRR;
    int pixelLL;
    int pixelRR;

    if(year=="2016"){
      tree->SetBranchAddress("nLocalProtCand",&nLocalProtCand);
      tree->SetBranchAddress("LocalProtCand_x",&LocalProtCand_x);
      tree->SetBranchAddress("LocalProtCand_z",&LocalProtCand_z);
    }else if (year=="2017"){
      tree->SetBranchAddress("HLT_Accept", &HLT_Accept);
      tree->SetBranchAddress("nRPPixelTrack", &nRPPixelTrack);
      tree->SetBranchAddress("RPPixelTrack_x", &RPPixelTrack_x);
      tree->SetBranchAddress("RPPixelTrack_y", &RPPixelTrack_y);
      tree->SetBranchAddress("RPPixelTrack_arm", &RPPixelTrack_arm);
      tree->SetBranchAddress("Pair_mass", &Pair_mass);
      tree->SetBranchAddress("Pair_pt", &Pair_pt);
      tree->SetBranchAddress("Pair_eta", &Pair_eta);
      tree->SetBranchAddress("Pair_phi", &Pair_phi);
      tree->SetBranchAddress("Pair_dphi", &Pair_dphi);
      tree->SetBranchAddress("ClosestExtraTrack_vtxdxyz", &ClosestExtraTrack_vtxdxyz);
      tree->SetBranchAddress("nRPTimingTrack", &nRPTimingTrack);
      tree->SetBranchAddress("RPTimingTrack_arm",&RPTimingTrack_arm);
      tree->SetBranchAddress("RPTimingTrack_Time",&RPTimingTrack_Time);
      tree->SetBranchAddress("ExtraTrack_z",&ExtraTrack_z);
      tree->SetBranchAddress("ExtraTrack_px",&ExtraTrack_px);
      tree->SetBranchAddress("ExtraTrack_py",&ExtraTrack_py);
      tree->SetBranchAddress("ExtraTrack_pz",&ExtraTrack_pz);
      tree->SetBranchAddress("Etmiss",&Etmiss);
    }else{
      exit (EXIT_FAILURE);
    }
    tree->SetBranchAddress("EventNum", &EventNum);
    tree->SetBranchAddress("LumiSection", &LumiSection);
    tree->SetBranchAddress("KalmanVertexCand_z",&KalmanVertexCand_z);
    tree->SetBranchAddress("MuonCand_vtxz",&MuonCand_vtxz);
    tree->SetBranchAddress("MuonCand_pt",&MuonCand_pt);
    tree->SetBranchAddress("MuonCand_eta",&MuonCand_eta);
    tree->SetBranchAddress("MuonCand_phi",&MuonCand_phi);
    tree->SetBranchAddress("MuonCand_e",&MuonCand_e);
    tree->SetBranchAddress("MuonCand_charge",&MuonCand_charge);    
    tree->SetBranchAddress("MuonCand_isglobal",&MuonCand_isglobal);    
    tree->SetBranchAddress("MuonCand_istight",&MuonCand_istight);
    tree->SetBranchAddress("JetCand_vz",&JetCand_vz);
    tree->SetBranchAddress("JetCand_pt",&JetCand_pt);
    tree->SetBranchAddress("JetCand_eta",&JetCand_eta);
    tree->SetBranchAddress("JetCand_phi",&JetCand_phi);
    tree->SetBranchAddress("JetCand_e",&JetCand_e); 

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

    double alignX=0.;
    double alignX_2=0.;
    double alignX_3=0.;
    double alignX_102=0.;
    double alignX_103=0.;
    int chargeVal=0;
    int compactRpID=0;

    //starting loop over events, stops when reached end of file or nevt_max
    double fakedispersion_mm=80.;

    //calculate total PT of the system
    double LeadingMuonPt=0.;
    double LeadingJetPt=0.;
    double SecondMuonPt=0.;
    double SecondJetPt=0.;
    double LeadingMuonVtxZ=0.;
    double LeadingJetVtxZ=0.;
    double SecondMuonVtxZ=0.;
    double SecondJetVtxZ=0.;
    double VtxZ=0.;
    double VtxZpps=0.;

    int LeadingMuonTightID=0;
    int SecondMuonTightID=0;
    double MeT = 0.;
    double Muons_pxtot=0.;
    double Muons_pytot=0.;
    double Muons_pztot=0.;
    double Muons_pttot=0.;
    double dileptonMass=0.;
    double dileptonPhi=0.;
    double dileptonEta=0.;
    double dileptonpT=0.;
    double dijetMass=0.;
    double dijetPhi=0.;
    double dijetEta=0.;
    double dijetpT=0.;
    double acoplan=0.;
    double closestExtraTrack=0.;

    int nExtraTrack0pt5mm = 0;
    int nExtraTrack1mm = 0;
    double Muons_Etot=0.;
    double align_trk0=0.;double align_trk1=0.;
    double corrected_X_0=0.;
    double corrected_X_1=0.;
    double proton1_pz=0.;
    double proton2_pz=0.;
    double protonmass=0.;
    double missingmass=0.;
    double missingmassLorentzdilepton=0.;
    double missingmassLorentzdijet=0.;
    double missingmassLorentz=0.;
    double xi1=0.;
    double xi2=0.;
    int angle = -1;
    int triggerA = 0;
    int triggerB = 0;
    double pixelXarm0=0.;
    double pixelXarm1=0.;
    double pixelYarm0=0.;
    double pixelYarm1=0.;

    int therunnumber=0;
    int thelumisection=0;
    int theevent=0;

    double AverageTimeSec45 = 0.;
    double AverageTimeSec56 = 0.;

    TString ttreefilename = outputFileNameTree;
    fOut = new TFile(ttreefilename, "RECREATE");
    fOut->cd();

    trout = new TTree("Events", "Events");
    trout->Branch("Run",&therunnumber,"Run/i");
    trout->Branch("LumiSection",&thelumisection,"LumiSection/I");
    trout->Branch("EventNum",&theevent,"EventNum/I");
    trout->Branch("missingmass",&missingmass,"missingmass/D");
    trout->Branch("missingmassLorentzdilepton",&missingmassLorentzdilepton,"missingmassLorentzdilepton/D");
    trout->Branch("missingmassLorentzdijet",&missingmassLorentzdijet,"missingmassLorentzdijet/D");
    trout->Branch("missingmassLorentz",&missingmassLorentz,"missingmassLorentz/D");
    trout->Branch("nJetCand",&nJetCand,"nJetCand/I");
    trout->Branch("nMuonCand",&nMuonCand,"nMuonCand/i");
    trout->Branch("LeadingJetPt",&LeadingJetPt,"LeadingJetPt/D");
    trout->Branch("SecondJetPt",&SecondJetPt,"SecondJetPt/D");
    trout->Branch("LeadingJetVtxZ",&LeadingJetVtxZ,"LeadingJetVtxZ/D");
    trout->Branch("SecondJetVtxZ",&SecondJetVtxZ,"SecondJetVtxZ/D");
    trout->Branch("LeadingMuonPt",&LeadingMuonPt,"LeadingMuonPt/D");
    trout->Branch("LeadingMuonVtxZ",&LeadingMuonVtxZ,"LeadingMuonVtxZ/D");
    trout->Branch("LeadingMuonTightID",&LeadingMuonTightID,"LeadingMuonTightID/I");
    trout->Branch("SecondMuonPt",&SecondMuonPt,"SecondMuonPt/D");
    trout->Branch("SecondMuonVtxZ",&SecondMuonVtxZ,"SecondMuonVtxZ/D");
    trout->Branch("SecondMuonTightID",&SecondMuonTightID,"SecondMuonTightID/I");
    trout->Branch("acoplan",&acoplan,"acoplan/D");
    trout->Branch("VtxZ",&VtxZ,"VtxZ/D");
    trout->Branch("chargeVal",&chargeVal,"chargeVal/I");
    trout->Branch("nLocalProtCand",&nLocalProtCand,"nLocalProtCand/I");
    trout->Branch("proton1_pz",&proton1_pz,"proton1_pz/D");
    trout->Branch("proton2_pz",&proton2_pz,"proton2_pz/D");
    trout->Branch("xi1",&xi1,"xi1/D");
    trout->Branch("xi2",&xi2,"xi2/D");
    trout->Branch("protonmass",&protonmass,"protonmass/D");
    trout->Branch("triggerA",&triggerA,"triggerA/I");
    trout->Branch("triggerB",&triggerB,"triggerB/I");
    trout->Branch("dileptonMass",&dileptonMass,"dileptonMass/D");
    trout->Branch("dileptonEta",&dileptonEta,"dileptonEta/D");
    trout->Branch("dileptonPhi",&dileptonPhi,"dileptonPhi/D");
    trout->Branch("dileptonpT",&dileptonpT,"dileptonpT/D");
    trout->Branch("dijetMass",&dijetMass,"dijetMass/D");
    trout->Branch("dijetEta",&dijetEta,"dijetEta/D");
    trout->Branch("dijetPhi",&dijetPhi,"dijetPhi/D");
    trout->Branch("dijetpT",&dijetpT,"dijetpT/D");
    if(year=="2017"){
      trout->Branch("MeT",&MeT,"MeT/D");
      trout->Branch("pixelXarm0",&pixelXarm0,"pixelXarm0/D");
      trout->Branch("pixelXarm1",&pixelXarm1,"pixelXarm1/D");
      trout->Branch("pixelYarm0",&pixelYarm0,"pixelYarm0/D");
      trout->Branch("pixelYarm1",&pixelYarm1,"pixelYarm1/D");
      trout->Branch("VtxZpps",&VtxZpps,"VtxZpps/D");
      trout->Branch("closestExtraTrack",&closestExtraTrack,"closestExtraTrack/D");
      trout->Branch("nRPTimingTrack",&nRPTimingTrack, "nRPTimingTrack/I");
      trout->Branch("ExtraTrack_z",&VExtraTrack_z);
      trout->Branch("ExtraTrack_px",&VExtraTrack_px);
      trout->Branch("ExtraTrack_py",&VExtraTrack_py);
      trout->Branch("ExtraTrack_pz",&VExtraTrack_pz);
      trout->Branch("AverageTimeSec45",&AverageTimeSec45, "AverageTimeSec45/D");
      trout->Branch("AverageTimeSec56",&AverageTimeSec56, "AverageTimeSec56/D");
      trout->Branch("angle",&angle,"angle/I");
    }

    bool debug = false;

    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot)
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
      triggerA=0;
      triggerB=0;
      Muons_pxtot=0.;
      Muons_pytot=0.;
      Muons_pztot=0.;
      Muons_pttot=0.;
      LeadingMuonPt=0;
      LeadingJetPt=0;
      SecondMuonPt=0;
      SecondJetPt=0;
      LeadingMuonVtxZ=0;
      SecondMuonVtxZ=0;
      LeadingJetVtxZ=0;
      SecondJetVtxZ=0;
      MeT=0;
      VtxZ=0;
      VtxZpps=0;
      LeadingMuonTightID=0;
      SecondMuonTightID=0;
      dileptonMass=0.;
      dileptonPhi=0.;
      dileptonEta=0.;
      dileptonpT=0.;
      dijetMass=0.;
      dijetPhi=0.;
      dijetEta=0.;
      dijetpT=0.;
      acoplan=0.;
      closestExtraTrack=0.;
      Muons_Etot=0.;
      align_trk0=0.;
      align_trk1=0.;
      corrected_X_0=0.;
      corrected_X_1=0.;
      proton1_pz=0.;
      proton2_pz=0.;
      protonmass=0.;
      missingmass=0.;
      missingmassLorentzdilepton=0.;
      missingmassLorentzdijet=0.;
      missingmassLorentz=0.;
      xi1=0;
      xi2=0;
      angle = -1;
      pixelXarm0=0;
      pixelXarm1=0;
      pixelYarm0=0;
      pixelYarm1=0;

      // HLT Menu: 'HLT_DoubleEle33_CaloIdL_MW_v*','HLT_Ele27_WPTight_Gsf_v*'
      if(HLT_Accept[0] == 1) triggerA = 1;
      if(HLT_Accept[1] == 1 || HLT_Accept[2] == 1) triggerB = 1;

      if(year=="2017"){
	nLocalProtCand = 0;
	NpixelLL=0;
	NpixelRR=0;
	for(int s=0;s<nRPPixelTrack;s++)
	{
	  if(RPPixelTrack_arm[s]==0){
	    pixelcoordX[0][NpixelLL]=RPPixelTrack_x[s]; 
	    pixelcoordY[0][NpixelLL]=RPPixelTrack_y[s]; 
	    NpixelLL++; 
	    pixelLL=s; 
	    //nLocalProtCand++;
	  }
	  if(RPPixelTrack_arm[s]==1){
	    pixelcoordX[1][NpixelRR]=RPPixelTrack_x[s]; 
	    pixelcoordY[1][NpixelRR]=RPPixelTrack_y[s]; 
	    NpixelRR++;
	    pixelRR=s; 
	    //nLocalProtCand++;
	  }
	  if(debug){ 
	    std::cout << "pixel X [mm]: " << RPPixelTrack_x[s] << ", ARM: " << RPPixelTrack_arm[s] << std::endl;
	    std::cout << "pixel Y [mm]: " << RPPixelTrack_y[s] << ", ARM: " << RPPixelTrack_arm[s] << std::endl;
	  }
	}
	nLocalProtCand=NpixelLL+NpixelRR;
      }

      if(debug) std::cout << "nLocalProtCand: " << nLocalProtCand << std::endl;

      if(year=="2016"){

	const auto &fillInfo = fillInfoCollection.FindByRun(Run);
	if(fillInfo.fillNumber==8000)
	  continue;

	alignX=0.;
	alignX_2=0.;
	alignX_3=0.;
	alignX_102=0.;
	alignX_103=0.;

	const auto alignment_it = alignmentCollection.find(fillInfo.alignmentTag);
	if (alignment_it == alignmentCollection.end())
	{
	  printf("ERROR: no alignment for tag '%s'.\n", fillInfo.alignmentTag.c_str());
	  continue;
	}
	else
	{
	  alignments = &alignment_it->second;	      
	  alignX_2=alignments->GETX(2);
	  alignX_3=alignments->GETX(3);
	  alignX_102=alignments->GETX(102);
	  alignX_103=alignments->GETX(103);
	}
      }

      int anyproton=0;

      if(nLocalProtCand>0)
      {

	if(debug) std::cout << "Analyzing event with at least one proton..." << std::endl;
	anyproton++;

	if(nLocalProtCand==2)
	{

	  fakedispersion_mm=80.;

	  //calculate total PT of the system
	  LeadingMuonPt=0.;
	  LeadingJetPt=0.;
	  SecondMuonPt=0.;
	  SecondJetPt=0.;
	  LeadingMuonVtxZ=0.;
	  SecondMuonVtxZ=0.;
	  LeadingJetVtxZ=0.;
	  SecondJetVtxZ=0.;
	  MeT=0.;
	  VtxZ=0.;
	  VtxZpps=0.;
	  LeadingMuonTightID=0.;
	  SecondMuonTightID=0.;
	  Muons_pxtot=0.;
	  Muons_pytot=0.;
	  Muons_pztot=0.;
	  Muons_pttot=0.;
	  Muons_Etot=0.;
	  align_trk0=0.; align_trk1=0.;

	  therunnumber=0;
	  thelumisection=0;
	  theevent=0;

	  TLorentzVector p1;
	  TLorentzVector p2;

	  if(year=="2016"){
	    if(LocalProtCand_z[0]>0)
	      align_trk0=alignX_2;
	    else
	      align_trk0=alignX_102;
	    if(LocalProtCand_z[1]>0)
	      align_trk1=alignX_2;
	    else
	      align_trk1=alignX_102;
	    corrected_X_0=align_trk0+1000.*LocalProtCand_x[0];
	    corrected_X_1=align_trk1+1000.*LocalProtCand_x[1];
	  }

	  if(year=="2017"){

	    angle = -1;
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
	      angle=getXangle(Run,  LumiSection , "./inp/new/xangle_tillTS2_STABLEBEAMS_CLEANUP.csv");
	    }

	    if(preTS2==0) angle=getXangle(Run,  LumiSection , "./inp/xangle_afterTS2_cleanup.csv");
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
	  }else{
	    proton1_pz=(1-(corrected_X_0/fakedispersion_mm))*(ECM/2.);
	    proton2_pz=-(1-(corrected_X_1/fakedispersion_mm))*(ECM/2.);
	    xi1 = corrected_X_0/fakedispersion_mm;
	    xi2 = corrected_X_1/fakedispersion_mm;
	    p1.SetPxPyPzE(0.,0., ECM*xi1/2., ECM*xi1/2.);
	    p2.SetPxPyPzE(0.,0., -ECM*xi2/2., ECM*xi2/2.);
	    protonmass=ECM*sqrt((  (corrected_X_0/fakedispersion_mm)*(corrected_X_1/fakedispersion_mm) ));
	  }

	  TLorentzVector muon0;
	  TLorentzVector muon1;
	  TLorentzVector muonpair;

	  TLorentzVector jet0;
	  TLorentzVector jet1;
	  TLorentzVector jetpair;

	  TLorentzVector ppZsystem;

	  bool protonsAccept = false;

	  if(year=="2017"){
	    if(NpixelRR==1 && NpixelLL==1) protonsAccept = true;
	  }

	  if(year=="2016"){
	    if((LocalProtCand_z[0]>0 && LocalProtCand_z[1]<0) || (LocalProtCand_z[0]<0 && LocalProtCand_z[1]>0)) protonsAccept = true;
	  }

	  if(protonsAccept){

	    if(nMuonCand>1 && nJetCand>1){

	      therunnumber=Run;
	      thelumisection=LumiSection;
	      theevent=EventNum;

	      LeadingJetPt = JetCand_pt[0];
	      SecondJetPt = JetCand_pt[1];
	      LeadingJetVtxZ = JetCand_vz[0];
	      SecondJetVtxZ = JetCand_vz[1];
	      LeadingMuonPt=MuonCand_pt[0];
	      SecondMuonPt=MuonCand_pt[1];
	      LeadingMuonTightID=MuonCand_istight[0];
	      SecondMuonTightID=MuonCand_istight[1];
	      LeadingMuonVtxZ=MuonCand_vtxz[0];
	      SecondMuonVtxZ=MuonCand_vtxz[1];

	      VtxZ=KalmanVertexCand_z[0];

	      muon0.SetPtEtaPhiE(MuonCand_pt[0], MuonCand_eta[0], MuonCand_phi[0], MuonCand_e[0]);
	      muon1.SetPtEtaPhiE(MuonCand_pt[1], MuonCand_eta[1], MuonCand_phi[1], MuonCand_e[1]);
	      muonpair = muon0+muon1;

	      jet0.SetPtEtaPhiE(JetCand_pt[0], JetCand_eta[0], JetCand_phi[0], JetCand_e[0]);
	      jet1.SetPtEtaPhiE(JetCand_pt[1], JetCand_eta[1], JetCand_phi[1], JetCand_e[1]);
	      jetpair = jet0+jet1;

	      Muons_Etot=sqrt(MuonCand_e[0]*MuonCand_e[0]+MASS_MU*MASS_MU) + sqrt(MuonCand_e[1]*MuonCand_e[1]+MASS_MU*MASS_MU);
	      Muons_pxtot = muon0.Px() + muon1.Px();
	      Muons_pytot = muon0.Py() + muon1.Py();
	      Muons_pztot = muon0.Pz() + muon1.Pz();
	      chargeVal=(int)MuonCand_charge[0]+(int)MuonCand_charge[1];

	      dileptonMass=muonpair.M();
	      dileptonEta=muonpair.Eta();
	      dileptonPhi=muonpair.Phi();
	      dileptonpT=muonpair.Pt();

	      dijetMass=jetpair.M();
	      dijetEta=jetpair.Eta();
	      dijetPhi=jetpair.Phi();
	      dijetpT=jetpair.Pt();

	      double dphi = fabs(muon0.Phi()-muon1.Phi());
	      dphi = (dphi<pi) ? dphi : 2.*pi-dphi;
	      acoplan = 1 - fabs(Pair_dphi[0])/pi;

	      missingmass=sqrt((ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))*(ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))-(Muons_pxtot)*(Muons_pxtot)-(Muons_pytot)*(Muons_pytot)- (Muons_pztot+proton1_pz+proton2_pz)*(Muons_pztot+proton1_pz+proton2_pz));

	      ppZsystem = p1 + p2 - muonpair;
	      missingmassLorentzdilepton = ppZsystem.M(); 

	      ppZsystem = p1 + p2 - jetpair;
	      missingmassLorentzdijet = ppZsystem.M();

	      ppZsystem = p1 + p2 - (jetpair + muonpair);
	      missingmassLorentz = ppZsystem.M();

	      if(year=="2017"){

		closestExtraTrack=ClosestExtraTrack_vtxdxyz[0];
		MeT = Etmiss;

		double TotalTimeSec45 = 0.;
		double TotalTimeSec56 = 0.;

		int NSec45 = 0;
		int NSec56 = 0;

		for(int k=0; k<nRPTimingTrack; k++){
		  if(RPTimingTrack_arm[k]==0){
		    TotalTimeSec45 += RPTimingTrack_Time[k];
		    NSec45++;
		  }
		  if(RPTimingTrack_arm[k]==1){
		    TotalTimeSec56 += RPTimingTrack_Time[k];
		    NSec56++;
		  }
		}

		VExtraTrack_z.clear();
		VExtraTrack_px.clear();
		VExtraTrack_py.clear();
		VExtraTrack_pz.clear();
		for(int k=0; k<21; k++){
		  VExtraTrack_z.push_back(ExtraTrack_z[k]);
		  VExtraTrack_px.push_back(ExtraTrack_px[k]);
		  VExtraTrack_py.push_back(ExtraTrack_py[k]);
		  VExtraTrack_pz.push_back(ExtraTrack_pz[k]);
		}

		double deltaT = 0;

		if(NSec45>0) AverageTimeSec45 = TotalTimeSec45/NSec45;
		else AverageTimeSec45 = 0;
		if(NSec56>0) AverageTimeSec56 = TotalTimeSec56/NSec56;
		else AverageTimeSec56 = 0;
		if(NSec45 > 0 && NSec56 > 0) deltaT = AverageTimeSec45-AverageTimeSec56;

		VtxZpps=deltaT * ns_to_s_ * TMath::C() * 0.5 * m_to_cm_;

	      }

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

  if(year_str=="2016"){
    preTS2 = 0;
    std::cout << "Analyzing 2016 data..." << std::endl;
    if(era_str=="B")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016B_23Sep2016v3_FullFinal_pT40.root");
    if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016C_23Sep2016v1_FullFinal_pT40.root");
    if(era_str=="G")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016G_23Sep2016v1_FullFinal_pT40.root");
  } else if(year_str=="2017"){
    std::cout << "Analyzing 2017 data..." << std::endl;
    if(era_str=="B"||era_str=="C"||era_str=="D") preTS2 = 1;
    if(era_str=="E"||era_str=="F") preTS2 = 0;
    if(era_str=="B")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonB.root");
    if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonC.root");
    if(era_str=="D")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonD.root");
    if(era_str=="E")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonE.root");
    if(era_str=="F")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonF.root");
  }else{
    std::cout << "Please, put the right year." << std::endl;
    return 0;
  }

  string filenameTree="/afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/DarkMatterSearch/Analyser/TTree"+era_str+"_"+extrastr+".root";
  DarkMatterSearchLeadingMuonsJets(year_str, allnamesWithPrefixFinal, filenameTree, maxevent);

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
