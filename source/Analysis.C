#include <map>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TRandom.h>
#include <TMath.h>
#include <string.h>
#include <TError.h>
#include <sstream>
#include <fstream>

using namespace std;

class WaveForms {
public:
  WaveForms();
  ~WaveForms();
  vector<vector<Bool_t>> hasSignal; // To let us know if there was a signal that fired
  vector<Bool_t> hasHit; // To let us know if there was any signal other than FPNs
  vector<Bool_t> hasFPN; // To let us know if there was any FPN signal
  vector<Bool_t> doneFPN; // To let us know if there was any FPN signal averaged
  vector<vector<Bool_t>> isOverflow; // To let us know if the signal was overflow.
  Bool_t hasOverflow; // To let us know if the event has overflow channels.,
  vector<vector<Int_t>> isDecay;
  Bool_t hasImplant;
  Bool_t hasDecay;
  UInt_t frameIdx;
  UInt_t coboIdx;
  UInt_t asadIdx;
  vector<vector<vector<UInt_t>>> waveform; // To save the waveform for each Aget, Channel and Bucket
  vector<vector<vector<Int_t>>> corrwaveform; // To save the corrected waveform by the averaged FPN waveform
  vector<vector<UInt_t>> fpnwaveform; // To save the waveform for each Aget, FPN Channel and Bucket
  vector<vector<UInt_t>> energy; // Digitized energy value
  vector<vector<UInt_t>> time; // time value
  vector<vector<UInt_t>> background; // background value
  vector<vector<Double_t>> PSDIntegral; // To save ratio of peak to integral.
  vector<vector<Double_t>> PSDRatio; // To save ratio of peak to integral.
};

class MM_Track {
public:
  MM_Track();
  ~MM_Track();
  UInt_t hasTrack; // total number of hits in track
  UInt_t hasTrackChain; // total number of chain hits in track
  UInt_t hasTrackStrip; // total number of strip hits in track
  Bool_t hasProjectile; // incoming projectile flag 
  Bool_t hasEjectile; // outgoing ejectile flag 
  Bool_t hasRecoil; // outgoing recoil flag 
  Bool_t hasOverflow; // overflow flag
  Bool_t hasImplant; // implant flag
  Bool_t hasDecay; // decay flag
  UInt_t decay[150][170]; // decay flag by pixel X and Y, 0=implant, 1=decay
  UInt_t pixel[150][170]; // hit pattern by pixel X and Y
  UInt_t energy[150][170]; // Digitized energy value
  UInt_t time[150][170]; // time value
  Double_t avgposx[170]; // Average x position value
  Double_t avgposy[170]; // Average y position value
  Double_t sumenergy[170]; // Average energy value
  Double_t sum2penergy[4]; // Sum of all 2p energies, 0=bc,1=bl,2=br,4=all
  UInt_t coloridx[150][170]; // color value
  Double_t posx[150][170]; // Position X by pixel X
  Double_t posy[150][170]; // Position Y by pixel Y
  Double_t posz[150][170]; // Position Z by time and drift velocity
  Double_t posxerr[150][170]; // Position X Error by pixel X
  Double_t posyerr[150][170]; // Position Y Error by pixel Y
  Double_t poszerr[150][170]; // Position Z Error by time and drift velocity
};

class MapChanToMM {
public:
  MapChanToMM();
  ~MapChanToMM();
  UInt_t pxidx[4][4][64];
  UInt_t pxidy[4][4][64];
  Double_t pxposx[4][4][64];
  Double_t pxposy[4][4][64];
};

WaveForms::WaveForms() {
}

WaveForms::~WaveForms() {
}

MM_Track::MM_Track() {
}

MM_Track::~MM_Track() {
}

MapChanToMM::MapChanToMM() {
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<64;k++){
        pxidx[i][j][k]=-1;
        pxidy[i][j][k]=-1;
      }
    }
  }
}

MapChanToMM::~MapChanToMM() {
}

void InitWaveforms();
void ResetWaveforms();
void Init();
void RootROpenFile(string & inputFileName, string & inputTreeName);
void RootRInit();
void RootRInitWaveforms();
void RootReadEvent();
void RootRResetWaveforms();
void RootRReset();
void RootRCloseFile();
void RootHOpenFile(string & outputFileName);
void RootHWrite();
void RootHCloseFile();
void SetHWrite(int flag);
int GetHWrite();
void SetRootConverter(int flag);
void SetDrawWaveform(int flag);
void SetScaler(int flag);
void Set2pMode(int flag);
void SetBucketSize(int BucketSize);
void ReadMapChanToMM(string filename);
void ReadResponseWaveform(string filename);
void SetResponseWaveform();
Bool_t IsResponseSample(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
void SetResponseSample(Int_t type, Int_t evtno, Int_t buckcut, Int_t buckwidth, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t rep, Int_t iter, Int_t boost);
void GetMaxResponseWaveform();
void GetSigmaResponseWaveform();
Double_t ShaperF_GET(Double_t *x, Double_t *p);
Double_t ShaperF_MSCF(Double_t *x, Double_t *p);
void GetAverageFPN(Int_t cobo, Int_t asad, Int_t aget);
void GetCorrWaveform(Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
void GetEnergyTime(Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
Int_t GetEnergybyFitWaveform(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t maxAt, Int_t peakAt);
//Int_t GetEnergybyFitWaveform(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t peakAt);
void DrawWaveForm(Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
void FillTrack();
void ResetHitPattern();
void DrawHitPattern(UInt_t cobo);
void DrawdEvsE();
void CleanTrack();
void FilldEvsE();
void DrawTrack();
void HoughTransform();
void GetTrackPosYLimit();
void GetXYZTrack();
void InitTrack();
void ResetTrack();
void ResetTrackHist();
void ResetdEvsE();

static unsigned int triggerrate[2];
static int maxcobo;
static int maxasad;
static WaveForms* waveforms;
static WaveForms* rwaveforms[2];
static MM_Track* mm_tracks;
static MapChanToMM* mapchantomm;
static TH1D* hWaveForm[32];
static TH1D* hCorrWaveForm[32];
static TH2D* hWaveFormbyEvent[16];
static TH2D* hCorrWaveFormbyEvent[16];
static TSpectrum* sCorrWaveForm[32];
static TH2D* hMM_TrackAll;
static TH2D* hMM_Track[16];
static TH2D* hMM_TrackDecay[16];
static TGraph2D* gMM_TrackDecay[16];
static TH2D* hMM_TrackvsE[16];
static TH2D* hMM_TrackPosAll;
static TH2D* hMM_TrackPos[16];
static TGraph2D* gMM_TrackDecayPos[16];
static TH2D* hMM_TrackdEvsEALL[3];
static TH2D* hMM_TrackdEvsE[16][3];
static TH2D* hMM_TrackPosHough[16][3];
static TH2D* hMM_TrackXZ[16];
static TH2D* hMM_TrackYZ[16];
static TH2D* hMM_TrackPosXYAll;
static TH2D* hMM_TrackPosXYTAll;
static TH2D* hMM_TrackPosXZAll;
static TH2D* hMM_TrackPosYZAll;
static TH2D* hMM_TrackPosXY[16];
static TH2D* hMM_TrackPosXZ[16];
static TH2D* hMM_TrackPosYZ[16];
static TH2D* hMM_TimevsPxIDX[16];
static TH2D* hMM_TimevsPxIDXPosAll;
static TH2D* hMM_TimevsPxIDXPos[16];
static TH2D* hMM_TimevsPxIDXPosHough[16][3];
static TH2D* hMM_TimevsPxIDY[16];
static TH2D* hMM_TimevsPxIDYPosAll;
static TH2D* hMM_TimevsPxIDYPos[16];
static TH2D* hMM_TimevsPxIDYPosHough[16][3];
static TH2D* hMM_EnergyvsPxIDYALL;
static TH2D* hMM_EnergyvsPxIDY[16];
static TH1I* hMM_SumEnergyvsPxIDY[16];
static TH2D* hMM_PSDIntegralvsMax[16];
static TH2D* hMM_PSDRatiovsMax[16];
static TH2D* hMM_PSDRatio[16];
static TH1D* hMM_D2PTime;
static TH1D* hMM_Sum2pEnergy[4];
static TH1D* hMM_Time[16];
static TH1D* hMM_Energy[16];
static TH1D* hGET_EALL[32];
static TH1D* hGET_E[32][68];
static TH1I* hGET_HitPattern;
static TH2I* hGET_EHitPattern2D;
static TH2I* hGET_THitPattern2D;
static TH2D* hGET_EHitPattern[16];
static TH2D* hGET_THitPattern[16];
static TH2D* hGET_ERHitPattern[16];
static TH1D* hMM_a0;
static TH1D* hMM_a1;
static TH1D* hMM_b0;
static TH1D* hMM_b1;
static TRandom* rpos;
static int bucketmax;
static int readtype;
static int enableroot;
static int enabledraww;
static int enablescaler;
static int enablehwrite;
static int enableupdatefast;
static int hasDrawn[16];
static UInt_t weventIdx;
static UInt_t reventIdx;
static UInt_t rd2ptime;
static UInt_t FirsteventIdx;
static UInt_t LasteventIdx;
static Bool_t IsFirstevent;
static UInt_t IsTrig[16];
static int evtcounter;
static int goodevtcounter;
static int badevtcounter;
static int printed=0;
static int hitcounter;
static int prevmidx;
static int mutantcounter;
static vector<UInt_t> mevtidx; // eventidx value
static vector<UInt_t> md2ptime; // d2ptime value
static double evtno;
static double scaler1;
static double scaler1start;
static double scaler1end;
static double scaler2;
static double scaler2start;
static double scaler2end;
static double scaler3;
static double scaler3start;
static double scaler3end;
static double scaler4;
static double scaler4start;
static double scaler4end;
static double scaler5;
static double scaler5start;
static double scaler5end;
static double driftv;
static double rx;
static double ry;
static double rz;
static double rt;
static UInt_t mm_minenergy;
static UInt_t mm_maxenergy;
static UInt_t mm_mintime;
static UInt_t mm_maxtime;
static UInt_t si_minenergy;
static UInt_t si_maxenergy;
static TFile* fInputFile;
static TTree* fInputTree;
static TFile* fOutputFile;
static UInt_t fInputFileSize;
static TTree* fOutputTree;
static TFile* fhOutputFile;
static Int_t fNumberEvents;
static Float_t fTimePerBin;
static Int_t wGETMul;
static Int_t wGETHit;
static Int_t wGETEventIdx;
static Int_t wGETD2PTime;
static Int_t wGETFrameNo[4352];
static Int_t wGETDecayNo[4352];
static Float_t wGETTime[4352];
static Float_t wGETEnergy[4352];
static Int_t wGETCobo[4352];
static Int_t wGETAsad[4352];
static Int_t wGETAget[4352];
static Int_t wGETChan[4352];
static Int_t wGETWaveformX[4352][512];
static Int_t wGETWaveformY[4352][512];
static Int_t rGETMul;
static Int_t rGETHit;
static Int_t rGETEventIdx;
static Int_t rGETD2PTime;
static Int_t rGETFrameNo[4352];
static Int_t rGETDecayNo[4352];
static Float_t rGETTime[4352];
static Float_t rGETEnergy[4352];
static Int_t rGETCobo[4352];
static Int_t rGETAsad[4352];
static Int_t rGETAget[4352];
static Int_t rGETChan[4352];
static Int_t rGETWaveformX[4352][512];
static Int_t rGETWaveformY[4352][512];
static Double_t posenergysum[170];
static Double_t energysum[170];
static Int_t sumcounter[170];
static Int_t hMM_TrackCounter[3];
static Int_t hMM_TrackdECounter[3];
static Double_t hMM_dE[3];
static Double_t hMM_E[3];
static int csi_channo[9]={2,7,10,16,19,25,28,33,36};
static Double_t trackposymin[3];
static Double_t trackposymax[3];
static ostringstream scalerout;
static string mapChanToMM = "mapchantomm.txt";
static string rffilename = "responsewaveform.txt";
static string inrfname = "example.mfm.root";
static string inrtname = "TEvent";
static string outrfname = "histogram.root";
static int BucketSize = 512;
static int ScalerMode = 0;
static int DrawWaveform = 0;
static int RootOutput = 1;
static int maxrespsample = 0;
static TF1* fResponse[10];
static TH1D* hResponse[10];
static Double_t response[10][512];
static Int_t responsesample[10][7];
static UInt_t maxresponse[10];
static UInt_t responsesigma[10];
static UInt_t deconvrep[10];
static UInt_t deconviter[10];
static UInt_t deconvboost[10];
static Double_t corrwaveformdec[512];
static Double_t corrwaveformfit[512];

void Analysis(){
  SetBucketSize(BucketSize);
  Init();
  SetScaler(ScalerMode);
  ReadMapChanToMM(mapChanToMM);
  //ReadResponseWaveform(rffilename);
  SetResponseSample(0,7847,150,150,0,2,1,38,5,5,1); // MM (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
  SetResponseSample(1,7848,100,200,1,1,0,2,30,30,1); // IC (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
  //SetResponseSample(2,7847,100,200,1,1,0,2,20,20,1); // IC (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
  SetDrawWaveform(DrawWaveform);
  SetHWrite(RootOutput);
  RootROpenFile(inrfname, inrtname);
  RootReadEvent();
  RootRCloseFile();
  if(GetHWrite()==1){
    RootHOpenFile(outrfname);
    RootHWrite();
    RootHCloseFile();
    cout << "Data Analysis completed." << endl;
  }
}

void Init() {
  enableroot = 0;
  enablehwrite = 0;
  enableupdatefast = 0;
  maxasad = 4;
  FirsteventIdx = 0;
  IsFirstevent = true;
  LasteventIdx = 0;
  //mm_minenergy = 100;
  //mm_maxenergy = 4096;
  mm_minenergy = 200; //for AllChannelReading
  //mm_minenergy = 0;
  mm_maxenergy = 4096;
  //mm_mintime = 130;
  //mm_maxtime = 380;
  mm_mintime = 0;
  mm_maxtime = 512;
  si_minenergy = 0;
  si_maxenergy = 4000;
  driftv = 0.0025; // mm/ns, for P5 at 100 Torr and -800V
  //driftv = 0.017; // mm/ns, for CO2 at 20 Torr and -800V
  fTimePerBin = 40.; // 40 ns per bin
  scalerout.clear();
  scalerout.str("");

  mapchantomm = new MapChanToMM();
  for(int i=0; i<maxasad ; i++) {
    for(int j=0; j<4 ; j++) {
      for(int k=0; k<64 ; k++) {
	mapchantomm->pxidx[i][j][k]=0;
	mapchantomm->pxidy[i][j][k]=0;
      }
    }
  }

  for(int i=0;i<2;i++){
    rwaveforms[i] = new WaveForms();
  }
  mm_tracks = new MM_Track();
  
  RootRInitWaveforms();

  hMM_D2PTime = new TH1D(Form("hMM_D2PTime"),Form("hMM_D2PTime;dTime (10ns/bin);Counts"),3000000,0,30000000);
  for(int i=0;i<3;i++){
    hMM_Sum2pEnergy[i] = new TH1D(Form("hMM_Sum2pEnergy_%d",i),Form("hMM_Sum2pEnergy_%d;Energy (ADC Channel);Counts",i),1000000,0,1000000);
  }
  hMM_Sum2pEnergy[3] = new TH1D(Form("hMM_Sum2pEnergy_All"),Form("hMM_Sum2pEnergy_All;Energy (ADC Channel);Counts"),1000000,0,1000000);

  for(int cobo=0; cobo<2 ; cobo++) {
    for(int i=0; i<maxasad ; i++) {
      for(int j=0; j<4 ; j++) {
        hWaveForm[cobo*maxasad*4+i*4+j] = new TH1D(Form("hWaveForm_%d_%d_%d",cobo,i,j),Form("hWaveForm_%d_%d_%d;Bucket No;ADC Channel",cobo,i,j),512,0,512);
        hCorrWaveForm[cobo*maxasad*4+i*4+j] = new TH1D(Form("hCorrWaveForm_%d_%d_%d",cobo,i,j),Form("hCorrWaveForm_%d_%d_%d;Bucket No;ADC Channel",cobo,i,j),512,0,512);
        sCorrWaveForm[cobo*maxasad*4+i*4+j] = new TSpectrum();
      }
    }
  }
  for(int i=0; i<16 ; i++) {
    hWaveFormbyEvent[i] = new TH2D(Form("hWaveFormbyEvent_%d",i),Form("hWaveFormbyEvent_%d;Bucket No;ADC Channel",i),512,0,512,512,0,4096);
    hCorrWaveFormbyEvent[i] = new TH2D(Form("hCorrWaveFormbyEvent_%d",i),Form("hCorrWaveFormbyEvent_%d;Bucket No;ADC Channel",i),512,0,512,512,0,4096);
  }

  // Histograms for particle tracks in Micromega
  hMM_a0 = new TH1D(Form("hMM_a0"),Form("hMM_a0;a0 (1);Counts"),20000,-100,100);
  hMM_a1 = new TH1D(Form("hMM_a1"),Form("hMM_a1;a1 (mm);Counts"),20000,-100,100);
  hMM_b0 = new TH1D(Form("hMM_b0"),Form("hMM_b0;b0 (ns/mm);Counts"),20000,-100,100);
  hMM_b1 = new TH1D(Form("hMM_b1"),Form("hMM_b1;b1 (ns);Counts"),20000,-100,100);
  hMM_TrackPosAll = new TH2D(Form("hMM_TrackPosAll"),Form("All Event Tracks Position in the Micromega;X (mm);Y (mm)"),300,-150,150,700,-300,400);
  hMM_TrackAll = new TH2D(Form("hMM_TrackAll"),Form("All Event Tracks in the Micromega;Cell ID X;Cell ID Y"),150,0,150,170,0,170);
  hMM_TrackPosXYAll = new TH2D(Form("hMM_TrackPosXYAll"),Form("All Event Tracks Pos Y vs Pos X;X (mm);Y (mm)"),300,-150,150,700,-300,400);
  hMM_TrackPosXYTAll = new TH2D(Form("hMM_TrackPosXYTAll"),Form("All Event Tracks Pos Y vs Pos X by Time;X (mm);Y (mm)"),300,-150,150,700,-300,400);
  hMM_TrackPosXZAll = new TH2D(Form("hMM_TrackPosXZAll"),Form("All Event Tracks Pos Z vs Pos X;X (mm);Z (mm)"),300,-150,150,130,0,130);
  hMM_TrackPosYZAll = new TH2D(Form("hMM_TrackPosYZAll"),Form("All Event Tracks Pos Z vs Pos Y;Y (mm);Z (mm)"),700,-300,400,130,0,130);
  hMM_TimevsPxIDXPosAll = new TH2D(Form("hMM_TimevsPxIDXPosAll"),Form("All Events Time vs Position X in the Micromega;X (mm);Time (%fns/bin)",fTimePerBin),300,-150,150,512,0,512);
  hMM_TimevsPxIDYPosAll = new TH2D(Form("hMM_TimevsPxIDYPosAll"),Form("All Events Time vs Position Y in the Micromega;Y (mm);Time (%fns/bin)",fTimePerBin),700,-300,400,512,0,512);
  hMM_EnergyvsPxIDYALL = new TH2D(Form("hMM_EnergyvsPxIDYALL"),Form("All Event Energy vs Pixel Y in the Central Pad Micromega;Y (pixel);Energy (ADC Channel)"),150,0,150,4000,0,4000);
  for(int i=0;i<16;i++){
    hMM_Track[i] = new TH2D(Form("hMM_Track_%d",i),Form("Single Event Track in the Micromega;X (pixel);Y (pixel)"),150,0,150,170,0,170);
    hMM_TrackDecay[i] = new TH2D(Form("hMM_TrackDecay_%d",i),Form("Single Event Decay Track in the Micromega;X (pixel);Y (pixel)"),150,0,150,170,0,170);
    hMM_TrackPos[i] = new TH2D(Form("hMM_TrackPos_%d",i),Form("Single Event Track Position in the Micromega;X (mm);Y (mm)"),300,-150,150,700,-300,400);
    gMM_TrackDecay[i] = new TGraph2D();
    gMM_TrackDecay[i]->SetName(Form("gMM_TrackDecay_%d",i));
    gMM_TrackDecay[i]->SetTitle(Form("Single Event Track in the Micromega;X (pixel);Y (pixel)"));
    gMM_TrackDecayPos[i] = new TGraph2D();
    gMM_TrackDecayPos[i]->SetName(Form("gMM_TrackDecayPos_%d",i));
    gMM_TrackDecayPos[i]->SetTitle(Form("Single Event Track Position in the Micromega;X (mm);Y (mm)"));
    gMM_TrackDecayPos[i]->GetXaxis()->SetRangeUser(-150,150);
    gMM_TrackDecayPos[i]->GetYaxis()->SetRangeUser(-300,400);
    //gMM_TrackDecayPos[i]->GetZaxis()->SetRangeUser(-300,400);
    hMM_TrackvsE[i] = new TH2D(Form("hMM_TrackvsE_%d",i),Form("Single Event Track vs E in the Micromega;X (pixel);Y (pixel)"),150,0,150,170,0,170);
    for(int j=0;j<3;j++){
      hMM_TrackPosHough[i][j] = new TH2D(Form("hMM_TrackPosHough_%d_%d",i,j),Form("Single Event Hough Space from TrackPos_%d_%d;Theta (deg);Radius (mm)",i,j),1800,0,180,600,-300,300);
      //hMM_TrackdEvsE[i][j] = new TH2D(Form("hMM_TrackdEvsE_%d_%d",i,j),Form("Single Event dE vs E %d_%d;E (ch);dE (ch)",i,j),512,0,4096,512,0,4096);
    }
    hMM_TrackXZ[i] = new TH2D(Form("hMM_TrackXZ_%d",i),Form("Single Event Track Pixel Z vs Pixel X;X (pixel);Z (pixel)"),150,0,150,512,0,512);
    hMM_TrackYZ[i] = new TH2D(Form("hMM_TrackYZ_%d",i),Form("Single Event Track Pixel Z vs Pixel Y;Y (pixel);Z (pixel)"),170,0,170,512,0,512);
    hMM_TrackPosXY[i] = new TH2D(Form("hMM_TrackPosXY_%d",i),Form("Single Event Track Pos Y vs Pos X;X (mm);Y (mm)"),300,-150,150,700,-300,400);
    hMM_TrackPosXZ[i] = new TH2D(Form("hMM_TrackPosXZ_%d",i),Form("Single Event Track Pos Z vs Pos X;X (mm);Z (mm)"),300,-150,150,130,0,130);
    hMM_TrackPosYZ[i] = new TH2D(Form("hMM_TrackPosYZ_%d",i),Form("Single Event Track Pos Z vs Pos Y;Y (mm);Z (mm)"),700,-300,400,130,0,130);
    hMM_TimevsPxIDX[i] = new TH2D(Form("hMM_TimevsPxIDX_%d",i),Form("Single Event Time vs Pixel X in the Micromega;X (pixel);Time (%fns/bin)",fTimePerBin),150,0,150,512,0,512);
    hMM_TimevsPxIDXPos[i] = new TH2D(Form("hMM_TimevsPxIDXPos_%d",i),Form("Single Event Time vs Position X in the Micromega;X (mm);Time (%fns/bin)",fTimePerBin),300,-150,150,512,0,512);
    for(int j=0;j<3;j++){
      hMM_TimevsPxIDXPosHough[i][j] = new TH2D(Form("hMM_TimevsPxIDXPosHough_%d_%d",i,j),Form("Single Event Hough Space from Time vs Pos X %d_%d;Theta (deg);Radius (mm)",i,j),1800,0,180,1000,-500,500);
      hMM_TimevsPxIDYPosHough[i][j] = new TH2D(Form("hMM_TimevsPxIDYPosHough_%d_%d",i,j),Form("Single Event Hough Space from Time vs Pos Y %d_%d;Theta (deg);Radius (mm)",i,j),1800,0,180,1000,-400,600);
    }
    hMM_TimevsPxIDY[i] = new TH2D(Form("hMM_TimevsPxIDY_%d",i),Form("Single Event Time vs Pixel Y in the Micromega;Y (pixel);Time (%fns/bin)",fTimePerBin),170,0,170,512,0,512);
    hMM_TimevsPxIDYPos[i] = new TH2D(Form("hMM_TimevsPxIDYPos_%d",i),Form("Single Event Time vs Position Y in the Micromega;Y (mm);Time (%fns/bin)",fTimePerBin),700,-300,400,512,0,512);
    hMM_EnergyvsPxIDY[i] = new TH2D(Form("hMM_EnergyvsPxIDY_%d",i),Form("Single Event Energy vs Pixel Y in the Micromega;Y (pixel);Energy (ADC Channel)"),150,0,150,4000,0,4000);
    hMM_SumEnergyvsPxIDY[i] = new TH1I(Form("hMM_SumEnergyvsPxIDY_%d",i),Form("Single Event Sum Energy vs Pixel Y in the Micromega;Y (pixel);Sum Energy (ADC Channel)"),150,0,150);
    hMM_Time[i] = new TH1D(Form("hMM_Time_%d",i),Form("Single Event Time_%d;Time (%fns/bin);Counts",i,fTimePerBin),512,0,512);
    hMM_Energy[i] = new TH1D(Form("hMM_Energy_%d",i),Form("Single Event Energy_%d;Energy (ADC Channel);Counts",i),4000,0,4000);
    for(int j=0;j<3;j++){
    }
    for(int j=0;j<3;j++){
    }
    hMM_PSDIntegralvsMax[i] = new TH2D(Form("hMM_PSDIntegralvsMax_%d",i),Form("PSD Integral vs Amplitude for Single Event in the Micromega;Amplitude (ch);Integral (ch)"),512,0,4096,1000,0,200000);
    hMM_PSDRatiovsMax[i] = new TH2D(Form("hMM_PSDRatiovsMax_%d",i),Form("PSD Ratio vs Amplitude for Single Event in the Micromega;Amplitude (ch);Ratio"),512,0,4096,1000,0,1);
    hMM_PSDRatio[i] = new TH2D(Form("hMM_PSDRatio_%d",i),Form("PSD Integral vs Cases for Single Event in the Micromega;Case No.;PSDRatio"),10,0,10,1000,0,1);
  }
  for(int j=0;j<3;j++){
      hMM_TrackdEvsEALL[j] = new TH2D(Form("hMM_TrackdEvsE_ALL_%d",j),Form("All Event dE vs E %d;E (ch);dE (ch)",j),512,0,4096,512,0,4096);
  }
  
  // Histograms for the GET electronics raw channels
  hGET_HitPattern = new TH1I(Form("hGET_HitPattern"),Form("hGET_HitPattern"),4000,0,4000);
  hGET_EHitPattern2D = new TH2I(Form("hGET_EHitPattern2D"),Form("hGET_EHitPattern2D"),2500,0,2500,512,0,4096);
  hGET_THitPattern2D = new TH2I(Form("hGET_THitPattern2D"),Form("hGET_THitPattern2D"),2500,0,2500,512,0,512);
  for(int i=0; i<16 ; i++) {
    hGET_EHitPattern[i] = new TH2D(Form("hGET_EHitPattern_%d",i),Form("hGET_EHitPattern_%d",i),2000,0,2000,512,0,4096);
    hGET_THitPattern[i] = new TH2D(Form("hGET_THitPattern_%d",i),Form("hGET_THitPattern_%d",i),2000,0,2000,512,0,512);
    hGET_ERHitPattern[i] = new TH2D(Form("hGET_ERHitPattern_%d",i),Form("hGET_ERHitPattern_%d",i),2000,0,2000,100,0,1);
  }
  
  for(int cobo=0; cobo<2 ; cobo++) {
    for(int i=0; i<maxasad ; i++) {
      for(int j=0; j<4 ; j++) {
        hGET_EALL[cobo*16+i*4+j] = new TH1D(Form("hGET_EALL_%d_%d_%d",cobo,i,j),Form("hGET_EALL_%d_%d_%d",cobo,i,j),4000,0,4000);
        hGET_EALL[cobo*16+i*4+j]->SetTitle(Form("hGET_EALL_%d_%d_%d;ADC Channel;Counts",cobo,i,j));
        for(int k=0; k<68; k++) {
          hGET_E[cobo*16+i*4+j][k] = new TH1D(Form("hGET_E_%d_%d_%d_%d",cobo,i,j,k),Form("hGET_E_%d_%d_%d_%d",cobo,i,j,k),4000,0,4000);
          hGET_E[cobo*16+i*4+j][k]->SetTitle(Form("hGET_E_%d_%d_%d_%d;ADC Channel;Counts",cobo,i,j,k));
          if(k==11 || k==22 || k==45 || k==56){
            hGET_E[cobo*16+i*4+j][k]->SetTitle(Form("hGET_E_%d_%d_%d_%d (FPN);ADC Channel;Counts",cobo,i,j,k));
          }else{
          }
        }
      }
    }
  }
  InitTrack();
}

void InitWaveforms() {
  waveforms->waveform.resize(maxasad*4);
  waveforms->hasSignal.resize(maxasad*4);
  waveforms->hasHit.resize(maxasad*4); //default set to false
  waveforms->hasFPN.resize(maxasad*4); //default set to false
  waveforms->doneFPN.resize(maxasad*4); //default set to false
  for(int i=0;i<maxasad;i++){
    for(int j=0;j<4;j++){
      waveforms->waveform[i*4+j].resize(68);
      waveforms->hasSignal[i*4+j].resize(68); //default set to false
      for(int k=0;k<68;k++){
	waveforms->waveform[i*4+j][k].resize(bucketmax);
      }
    }
  }
}

void ResetWaveforms() {
  for(int i=0;i<maxasad;i++){
    for(int j=0;j<4;j++){
      if(waveforms->hasHit[i*4+j] || waveforms->hasFPN[i*4+j]){
        for(int k=0;k<68;k++){
          if(waveforms->hasSignal[i*4+j][k]){
            waveforms->hasSignal[i*4+j][k] = false;
	    fill(waveforms->waveform[i*4+j][k].begin(),
		      waveforms->waveform[i*4+j][k].end(),0);
          }
        }
      }
      waveforms->hasHit[i*4+j] = false;
      waveforms->hasFPN[i*4+j] = false;
      waveforms->doneFPN[i*4+j] = false;
    }
  }
}

void RootReadEvent() {
  Int_t frameIdx = 0;
  Int_t coboIdx = 0;
  Int_t asadIdx = 0;
  Int_t agetIdx = 0;
  Int_t chanIdx = 0;
  UInt_t MMCoboId = 0;
  UInt_t percent = 0;
  UInt_t Entries = 1;
  UInt_t oldbucketmax = bucketmax;
  ostringstream s4out;
  ostringstream s5out;
  ofstream f4out;
  ofstream f5out;
  RootRInit();
  fNumberEvents = fInputTree->GetEntries();
  if(fNumberEvents>0){
    Entries = fNumberEvents;
    percent = Entries/10;
  }
  SetResponseWaveform();
  GetMaxResponseWaveform();
  GetSigmaResponseWaveform();
  fNumberEvents=30;
  for(int i=0;i<fNumberEvents;i++){
    cout << i << " " << fNumberEvents << endl;
    fInputFile->cd();
    RootRInit();
    fInputTree->GetEntry(i);
    reventIdx = rGETEventIdx;
    LasteventIdx = reventIdx;
    if(IsFirstevent){
      FirsteventIdx = reventIdx;
      evtcounter=0;
      goodevtcounter=0;
      badevtcounter=0;
      IsFirstevent=false;
      RootRResetWaveforms();
    }else{
      evtcounter++;
      RootRResetWaveforms();
    }
    if((reventIdx-FirsteventIdx)%100==0){
      if(reventIdx==FirsteventIdx){
        cout << Form("Starting from Event No. %d.",reventIdx) << endl;
      }else{
        cout << Form("Upto Event No. %d Processed (from %d)..",reventIdx, FirsteventIdx) << endl;
      }
    }

    //cout << rGETEventIdx << " " << rGETD2PTime << " " << rGETMul << endl;
    for(Int_t i=0; i<rGETMul; i++){
      coboIdx = rGETCobo[i];
      asadIdx = rGETAsad[i];
      agetIdx = rGETAget[i];
      chanIdx = rGETChan[i];
      if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)){
       	rwaveforms[coboIdx]->hasFPN[asadIdx*4+agetIdx] = true;
      }
      //cout << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << rwaveforms[coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx] << endl;
      for(int j=0;j<bucketmax;j++){
        rwaveforms[coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
      }
    }
    //cout << "Done with reading events from root tree" << endl;
    
    for(Int_t i=0; i<rGETMul; i++){
      frameIdx = rGETFrameNo[i];
      coboIdx = rGETCobo[i];
      asadIdx = rGETAsad[i];
      agetIdx = rGETAget[i];
      chanIdx = rGETChan[i];
      rwaveforms[coboIdx]->frameIdx = frameIdx;
      if(chanIdx!=11 && chanIdx!=22 && chanIdx!=45 && chanIdx!=56){ // We want to skip the FPN channels.
        GetAverageFPN(coboIdx,asadIdx,agetIdx);
        GetCorrWaveform(coboIdx,asadIdx,agetIdx,chanIdx);
        GetEnergyTime(coboIdx,asadIdx,agetIdx,chanIdx);
      }
      //cout << i << " " << rGETMul << ", Done with energy/time" << endl;
      if(enabledraww==1){
        DrawWaveForm(coboIdx,asadIdx,agetIdx,chanIdx);
        //cout << "Done with drawing waveform" << endl;
      }
    }
    if(enabledraww==1){
      hWaveFormbyEvent[evtcounter%16]->SetTitle(Form("hWaveFormbyEvent(EvtNo=%d);ADC Channel;Counts [D2PTime=%d usec]",reventIdx,(rd2ptime/1000)));
      hCorrWaveFormbyEvent[evtcounter%16]->SetTitle(Form("hCorrWaveFormbyEvent(EvtNo=%d);ADC Channel;Counts [ D2PTime=%d usec]",reventIdx,(rd2ptime/1000)));
    }
    ResetHitPattern();
    DrawHitPattern(MMCoboId);
    //WaveformShapeFilter(MMCoboId);
    //DrawPSDFilter(MMCoboId);
    //cout << "Done with filters" << endl;

    FillTrack();
    //cout << "Done with filling tracks" << endl;
    
    if(mm_tracks->hasTrack>0){
      //if((!mm_tracks->hasOverflow)||
      if((!mm_tracks->hasOverflow) && (mm_tracks->hasTrack>0)){
        goodevtcounter++;
        CleanTrack();
        cout << "Done with cleaning tracks" << endl;
        FilldEvsE();
        cout << "Done with filling dE vs E" << endl;
        ResetTrackHist();
        cout << "Done with resetting track histograms" << endl;
        DrawTrack();
        cout << "Done with drawing tracks" << endl;
        DrawdEvsE();
        cout << "Done with drawing dEvsE" << endl;
        HoughTransform();
        cout << "Done with Hough Transformation" << endl;
        GetTrackPosYLimit();
        cout << "Done with getting track posy min/max" << endl;
        GetXYZTrack();
        cout << "Done with getting track" << endl;
        ResetTrack();
        ResetdEvsE();
        cout << "Done with resetting tracks" << endl;
      }else{
        badevtcounter++;
        ResetTrack();
        ResetdEvsE();
      }
    }else{
      badevtcounter++;
      ResetTrack();
    }
    cout << "Done with finding tracks" << endl;
  }
  bucketmax = oldbucketmax;
}

//We will get the averaged FPN waveform for each Aget.
void GetAverageFPN(Int_t cobo, Int_t asad, Int_t aget){
  if(rwaveforms[cobo]->hasFPN[asad*4+aget] && !rwaveforms[cobo]->doneFPN[asad*4+aget]){
    for(UInt_t buck=0;buck<bucketmax;buck++){
      rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck] = rwaveforms[cobo]->waveform[asad*4+aget][11][buck];
      rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[cobo]->waveform[asad*4+aget][22][buck];
      rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[cobo]->waveform[asad*4+aget][45][buck];
      rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[cobo]->waveform[asad*4+aget][56][buck];
      rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck] = rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck]/4;
    }
    rwaveforms[cobo]->doneFPN[asad*4+aget] = true;
  }
}

void DrawWaveForm(Int_t cobo, Int_t asad, Int_t aget, Int_t chan){
  Int_t background = 0;
  if(chan==11 || chan==22 || chan==45 || chan==56) return; // We want to skip the FPN channels.
  if(!rwaveforms[cobo]->hasSignal[asad*4+aget][chan]) return; // skip signals that did not fire

  //if(cobo>=0 && ((rwaveforms[cobo]->energy[asad*4+aget][chan]>mm_minenergy && rwaveforms[cobo]->energy[asad*4+aget][chan]<mm_maxenergy && rwaveforms[cobo]->time[asad*4+aget][chan]>mm_mintime && rwaveforms[cobo]->time[asad*4+aget][chan]<mm_maxtime) || ())){
  if(cobo>=0){
    for(int i=0;i<=512;i++) hWaveForm[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
    for(int i=0;i<=512;i++) hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);

    background = rwaveforms[cobo]->background[asad*4+aget][chan];
    for(UInt_t buck=0;buck<bucketmax;buck++){
      if(cobo>=0){
        hWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[cobo]->waveform[asad*4+aget][chan][buck]);
        hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]-background);
        hWaveFormbyEvent[evtcounter%16]->Fill(buck,rwaveforms[cobo]->waveform[asad*4+aget][chan][buck]);
        hCorrWaveFormbyEvent[evtcounter%16]->Fill(buck,rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]-background);
      }
    }
    hWaveForm[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hWaveForm_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
    hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hCorrWaveForm_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
    UInt_t nfound = 0;
    //UInt_t nfound = sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Search(hCorrWaveForm[cobo*maxasad*4+asad*4+aget],2,"",0.4);
    if(nfound>100){
      rwaveforms[cobo]->hasSignal[asad*4+aget][chan] = 0;
    }
  }
}

void GetCorrWaveform(Int_t cobo, Int_t asad, Int_t aget, Int_t chan){
  Int_t minValue = 100000;
  for(UInt_t buck=0;buck<bucketmax;buck++){
    // Find the maximum
    if(cobo==0) {
      //rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
      rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
    }else if(cobo==1) {
      rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
      //rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
    }
    if(minValue>rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]){
      minValue = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
    }
  }

  if(minValue<0){
    for(UInt_t buck=0;buck<bucketmax;buck++){
      rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] - minValue;
    }
  }
}

void GetEnergyTime(Int_t cobo, Int_t asad, Int_t aget, Int_t chan){
  Int_t maxValue = -100000;
  Int_t maxValueBucket = 0;
  Int_t maxValuedec = -100000;
  Int_t maxValueBucketdec = 0;
  Int_t minValue = 100000;
  Int_t corrwaveformintegral = 0;
  Double_t psdratio = 0;
  Int_t backgroundSum = 0;
  Int_t backgroundSamples = 0;
  Int_t background = 0;
  UInt_t type;

  for(int i=0;i<512;i++){
   corrwaveformdec[i]=0;
   corrwaveformfit[i]=0;
  }

  for(UInt_t buck=0;buck<bucketmax;buck++){
    corrwaveformdec[buck] = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
  }

  for(UInt_t buck=0;buck<bucketmax;buck++){
    // Find the maximum
    if(maxValue<rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]){
      maxValue = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
      maxValueBucket = buck;
    }
  }

/*
  if(printed==0 && cobo==0 && asad==3 && aget==1 && chan==20){
    for(UInt_t buck=0;buck<bucketmax;buck++){
      cout << (buck-150) << " " << rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] << endl;
    }
    printed++;
  }
*/

  if(cobo==0){
    type=0;
  }else if(cobo==1){
    if(asad==0){
      type=1;
    }else if(asad==1){
      type=1;
/*
      for(int i=0;i<512;i++) hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(i,corrwaveformdec[i]);
      UInt_t nfound = sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Search(hCorrWaveForm[cobo*maxasad*4+asad*4+aget],1,"",0.4);
      if(nfound<2){
        type=1;
      }else{
        type=2;
      }
      for(int i=0;i<=512;i++) hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
*/
    }
  }
  sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Deconvolution(corrwaveformdec,response[type],bucketmax,deconvrep[type],deconviter[type],deconvboost[type]);
  for(UInt_t buck=maxresponse[type]*1.2;buck<bucketmax;buck++){
    // Find the maximum
    if(maxValuedec<corrwaveformdec[buck]) {
      maxValuedec = corrwaveformdec[buck];
      maxValueBucketdec = buck;
    }
  }
  maxValuedec = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][maxValueBucketdec];
  Int_t dBucket = TMath::Abs(maxValueBucket-maxValueBucketdec);
  if(maxValuedec>300 && dBucket>2 && dBucket<100){
    //maxValuedec = GetEnergybyFitWaveform(type,cobo,asad,aget,chan,maxValueBucketdec);
  //if(reventIdx==7861 && cobo==0 && asad==2 && aget==1 && chan==42){
    cout << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " " << type << " - ";
    cout << "(" << maxValue << "," << maxValueBucket << ") ";
    cout << ": (" << maxValuedec << "," << maxValueBucketdec << ") ";
    cout << "= (" << (maxValue-maxValuedec) << "," << (maxValueBucket-maxValueBucketdec) << ") "<<endl;
  //}
    maxValuedec = GetEnergybyFitWaveform(type,cobo,asad,aget,chan,maxValueBucket,maxValueBucketdec);
  }
  maxValue = maxValuedec;
  maxValueBucket = maxValueBucketdec;


 // 7861 0 2 1 42
  if(reventIdx==7861 && cobo==0 && asad==2 && aget==1 && chan==42){
    cout << reventIdx << " " << type << " " << cobo << " " << asad << " " << aget << " " << chan << " - ";
    cout << "(" << maxValue << "," << maxValueBucket << ") ";
    cout << ": (" << maxValuedec << "," << maxValueBucketdec << ") ";
    cout << "= (" << (maxValue-maxValuedec) << "," << (maxValueBucket-maxValueBucketdec) << ") "<<endl;
    for(UInt_t buck=0;buck<bucketmax;buck++){
      hWaveForm[0]->Fill(buck,corrwaveformdec[buck]);
      hWaveForm[1]->Fill(buck,rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]);
      hWaveForm[3]->Fill(buck,corrwaveformfit[buck]);
      hWaveForm[2]->Fill(buck,response[type][buck]);
      //cout << buck << " " << rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] << " " << corrwaveformdec[buck] << " ";
      //cout << (rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]-corrwaveformdec[buck]) << " ";
      //cout << response[cobo][buck] << endl;
    }
    hWaveForm[0]->SetLineColor(kRed);
    hWaveForm[1]->SetLineColor(kBlue);
    hWaveForm[2]->SetLineColor(kBlack);
    hWaveForm[3]->SetLineColor(kGreen);
    hWaveForm[4]->SetLineColor(kMagenta);
    hWaveForm[0]->SetOption("hist");
  }
/*
*/

  UInt_t bgmin=50;
  UInt_t bgmax=bucketmax-50;
  for(UInt_t buck=0;buck<bucketmax;buck++){
    if(maxValueBucket>bgmin && maxValueBucket<bgmax){
      if(buck>(bgmin-30) && buck<(bgmin-20)){
        backgroundSum += rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
        backgroundSamples++;
      }else if((buck>(bgmax+20) && buck<(bgmax+30))){
        backgroundSum += rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
        backgroundSamples++;
      }
    }else if(maxValueBucket<bgmin){
      if((buck>bgmax && buck<(bgmax+20))){
        backgroundSum += rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
        backgroundSamples++;
      }
    }else if(maxValueBucket>bgmax){
      if(buck>(bgmin-20) && buck<bgmin){
        backgroundSum += rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck];
        backgroundSamples++;
      }
    }
  }

  background = (backgroundSamples>0) ? round(backgroundSum/backgroundSamples) : 0;
  //background = 50;

  corrwaveformintegral=0;
  maxValue = maxValuedec-background;
  maxValueBucket = maxValueBucketdec;
  for(UInt_t buck=0;buck<bucketmax;buck++){
    corrwaveformintegral+=(rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]-background);
  }
  psdratio = ((Double_t)maxValue/(Double_t)corrwaveformintegral);
  //if(cobo==0 && maxValue>mm_minenergy && maxValue<mm_maxenergy && maxValueBucket>mm_mintime && maxValueBucket<mm_maxtime){
  if(cobo>=0){
    rwaveforms[cobo]->hasHit[asad*4+aget] = true;
    rwaveforms[cobo]->hasSignal[asad*4+aget][chan] = true;
    rwaveforms[cobo]->energy[asad*4+aget][chan] = maxValue;
    rwaveforms[cobo]->time[asad*4+aget][chan] = maxValueBucket;
    rwaveforms[cobo]->background[asad*4+aget][chan] = background;
    rwaveforms[cobo]->PSDIntegral[asad*4+aget][chan] = corrwaveformintegral;
    rwaveforms[cobo]->PSDRatio[asad*4+aget][chan] = psdratio;
    hGET_E[cobo*maxasad*4+asad*4+aget][chan]->Fill(rwaveforms[cobo]->energy[asad*4+aget][chan]);
    hGET_EALL[cobo*maxasad*4+asad*4+aget]->Fill(rwaveforms[cobo]->energy[asad*4+aget][chan]);
    hGET_EHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,rwaveforms[cobo]->energy[asad*4+aget][chan]);
    hGET_THitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,rwaveforms[cobo]->time[asad*4+aget][chan]);
  }
}

Int_t GetEnergybyFitWaveform(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t maxAt, Int_t peakAt){
  //TH1D* hWaveform;
  //TH1D* hScaledResponse;
  Int_t maxValuedec = 0;
  Int_t bucketmin = peakAt - maxresponse[type];
  Int_t bucketgap = peakAt - maxAt;
  Int_t bucketmax = peakAt + responsesigma[type]*1.5;

  Double_t ratio = rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][maxAt-1] / response[type][maxresponse[type]-(peakAt-maxAt)-1];
  Double_t chisq=0;
  Double_t ndof=0;
  Double_t minchisq=1E9;
  Double_t tryratio=0;
  Double_t bestratio=ratio;
/*
  hWaveform = new TH1D("hWaveform","hWaveform",512,0,512);
  for(int buck=bucketmin;buck<bucketmax;buck++){
    hWaveform->Fill(buck-bucketmin,rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]);
  }
*/
  Int_t maxstep=40;
  Double_t stepsize=0.05;
  for(int i=0;i<maxstep;i++){
    chisq = 0;
    for(int buck=0;buck<bucketmax;buck++){
      if((buck+bucketmin)<bucketmax){
        if(buck<maxresponse[type]-bucketgap || buck>maxresponse[type]+bucketgap){
          tryratio = ratio*(0.5+i*stepsize);
          chisq += TMath::Power(rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck+bucketmin] - response[type][buck]*tryratio,2);
          ndof += rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck+bucketmin];
        }
      }
    }
/*
    hScaledResponse = (TH1D*) hResponse[type]->Clone();
    hScaledResponse->Scale(ratio*(0.8+i*stepsize));
    chisq = hWaveform->Chi2Test(hScaledResponse,"CHI2");
*/
    chisq = chisq/ndof;
    if(minchisq>chisq){
      minchisq = chisq;
      bestratio = tryratio;
    }
    if(reventIdx==7861 && cobo==0 && asad==2 && aget==1 && chan==42){
      cout << i << " " << tryratio << " " << chisq << " " << minchisq << " " << bestratio << endl;
    }
    //delete hScaledResponse;
  }
  for(int buck=0;buck<bucketmax;buck++){
    if((buck+bucketmin)<bucketmax){
      corrwaveformfit[buck+bucketmin] = response[type][buck]*bestratio;
    }
  }
/*
  if(reventIdx==7847 && cobo==0 && asad==2 && aget==1 && chan==34){
    for(int buck=0;buck<bucketmax;buck++){
      if((buck+bucketmin)<bucketmax){
        hWaveForm[3]->Fill(buck+bucketmin,response[type][buck]*bestratio);
      }
    }
  }
*/
  if(response[type][maxresponse[type]]*bestratio>0 && response[type][maxresponse[type]]*bestratio>100000){
    maxValuedec = (Int_t)response[type][maxresponse[type]]*bestratio;
  }else{
    maxValuedec = 0;
  }
  //hWaveform->SetTitle(Form("minchisq=%f, bestratio=%f, maxValuedec=%d",minchisq,bestratio,maxValuedec));
  //delete hWaveform;
  return maxValuedec;
}
/*
Int_t GetEnergybyFitWaveform(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t peakAt){
  TF1* fWaveform;
  TH1D* hWaveform;
  Int_t maxValuedec = 0;
  Int_t bucketmin = peakAt - responsesigma[type]*1.5;
  Int_t bucketmax = peakAt + responsesigma[type]*1.5;

  if(type==0){
    //fWaveform = new TF1("fWaveform", ShaperF_GET, bucketmin, bucketmax, 6);
    fWaveform = new TF1("fWaveform", ShaperF_GET, 0, 512, 6);
    fWaveform->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
    fWaveform->SetParameters(10, 500, peakAt, responsesigma[type], 3, 0.2);
    fWaveform->FixParameter(3, responsesigma[type]);
    fWaveform->SetParLimits(1, 0, 10000);
  }else if(type>0) {
    fWaveform = new TF1("fWaveform", ShaperF_MSCF, bucketmin, bucketmax, 6);
    fWaveform->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
    fWaveform->SetParameters(10, 2000, peakAt, responsesigma[type], 3, 0.2);
    fWaveform->FixParameter(3, responsesigma[type]);
    fWaveform->SetParLimits(1, 0, 10000);
  }
  hWaveform = new TH1D("hWaveform","hWaveform",512,0,512);
  for(int buck=0;buck<512;buck++){
    hWaveform->Fill(buck,rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck]);
  }
  //hWaveform->Fit("fWaveform","RN"); //Q=Quiet and N=Do not draw
  hWaveform->Fit("fWaveform","N"); //Q=Quiet and N=Do not draw
  //maxValuedec = fWaveform->GetParameter(1);
  maxValuedec = fWaveform->Eval(peakAt);
  if(reventIdx==7847 && cobo==0 && asad==2 && aget==1 && chan==34){
    for(int buck=0;buck<512;buck++){
      cout << buck << " " << fWaveform->Eval(buck) << " " << rwaveforms[cobo]->corrwaveform[asad*4+aget][chan][buck] << endl;
      hWaveForm[3]->Fill(buck,fWaveform->Eval(buck));
    }
  }
  delete hWaveform;
  delete fWaveform;
  return maxValuedec;
}
*/

void ResetHitPattern() {
  for(int i=0;i<2000;i++){
    for(int j=0;j<512;j++){
      hGET_THitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
      hGET_EHitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
    }
    for(int j=0;j<100;j++){
      hGET_ERHitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
}

void DrawHitPattern(UInt_t cobo) {
  Double_t psdratio = 0;
  for(UInt_t asad=0; asad<maxasad; asad++) {
    for(UInt_t aget=0; aget<4; aget++) {
      if(rwaveforms[cobo]->hasHit[asad*4+aget]){
        for(UInt_t chan=0; chan<68; chan++) {
          if(!rwaveforms[cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire
          psdratio = rwaveforms[cobo]->PSDRatio[asad*4+aget][chan];
          hGET_EHitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,rwaveforms[cobo]->energy[asad*4+aget][chan]);
          hGET_THitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,rwaveforms[cobo]->time[asad*4+aget][chan]);
          hGET_ERHitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,psdratio);
        }
        hGET_EHitPattern[goodevtcounter%16]->SetTitle(Form("hGET_EHitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;ADC Channel (ch)",goodevtcounter%16,reventIdx));
        hGET_THitPattern[goodevtcounter%16]->SetTitle(Form("hGET_THitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;ADC Channel (ch)",goodevtcounter%16,reventIdx));
        hGET_ERHitPattern[goodevtcounter%16]->SetTitle(Form("hGET_ERHitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;Ratio",goodevtcounter%16,reventIdx));
      }
    }
  }
}

void WaveformShapeFilter(UInt_t cobo) {
  Int_t maxValue = -100000;
  UInt_t maxValueBucket = 0;
  Int_t projybegin,projyend;
  Int_t eentriesfh,emeanfh,elimitfh,eintegralfh;
  Int_t tentriesfh,tmeanfh,tlimitfh,tintegralfh;
  Int_t eentriessh,emeansh,elimitsh,eintegralsh;
  Int_t tentriessh,tmeansh,tlimitsh,tintegralsh;
  Int_t bin1,bin2;
  ostringstream s1out;
  ostringstream s2out;
  ostringstream s3out;
  ofstream f1out;
  ofstream f2out;
  ofstream f3out;
  for(UInt_t asad=0; asad<maxasad; asad++) {
    for(UInt_t aget=0; aget<4; aget++) {
      if(!rwaveforms[cobo]->hasHit[asad*4+aget]) continue; // skip aget that did not fire
      if((asad==2||asad==3)&&(aget==0||aget==1)){
        projybegin=asad*4*100+aget*100;
        projyend=asad*4*100+aget*100+32;
        eentriesfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
        emeanfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
        elimitfh = 50;
        bin1 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeanfh-elimitfh);
        bin2 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeanfh+elimitfh);
        eintegralfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
        tentriesfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
        tmeanfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
        tlimitfh = 2;
        bin1 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeanfh-tlimitfh);
        bin2 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeanfh+tlimitfh);
        tintegralfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
        projybegin=asad*4*100+aget*100+33;
        projyend=asad*4*100+aget*100+67;
        eentriessh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
        emeansh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
        elimitsh = 50;
        bin1 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeansh-elimitsh);
        bin2 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeansh+elimitsh);
        eintegralsh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
        tentriessh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
        tmeansh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
        tlimitsh = 2;
        bin1 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeansh-tlimitsh);
        bin2 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeansh+tlimitsh);
        tintegralsh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
      }
      for(UInt_t chan=0; chan<68; chan++) {
	if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
        if(!rwaveforms[cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

	maxValue = rwaveforms[cobo]->energy[asad*4+aget][chan];
	maxValueBucket = rwaveforms[cobo]->time[asad*4+aget][chan];
        if((asad==2||asad==3)&&(aget==0||aget==1) &&
	(!rwaveforms[cobo]->isOverflow[asad*4+aget][chan])){
	  if(chan<33 && eintegralfh*3>eentriesfh && tintegralfh*3>tentriesfh){
	    if( maxValue>=(emeanfh-elimitfh) &&
	       	maxValue<=(emeanfh+elimitfh) &&
	       	maxValueBucket>=(tmeanfh-tlimitfh) &&
	       	maxValueBucket<=(tmeanfh+tlimitfh)){
	          rwaveforms[cobo]->hasSignal[asad*4+aget][chan] = 0;
	          //continue; // skip fake signals by saturated signal
		  //s2out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 2" << endl;
            }
	  }else if(chan>32 && eintegralsh*3>eentriessh && tintegralsh*3>tentriessh){
	    if( maxValue>=(emeansh-elimitsh) &&
	       	maxValue<=(emeansh+elimitsh) &&
	       	maxValueBucket>=(tmeansh-tlimitsh) &&
	       	maxValueBucket<=(tmeansh+tlimitsh)){
	          rwaveforms[cobo]->hasSignal[asad*4+aget][chan] = 0;
	          //continue; // skip fake signals by saturated signal
		  //s2out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 2" << endl;
	    }
	  }
	}else if((cobo==0)&&(asad==0||asad==1)&&(!rwaveforms[cobo]->isOverflow[asad*4+aget][chan])){
	  //s1out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 1" << endl;
	}
      }
    }
  }
  /*
  if(s1out.str().size()!=0){
    f1out.open("goodlabels.txt", std::ofstream::out|std::ofstream::app);
    f1out << s1out.str() << endl;
    f1out.close();
    s1out.clear();
    s1out.str("");
  }
  if(s2out.str().size()!=0){
    f2out.open("badlabels.txt", std::ofstream::out|std::ofstream::app);
    f2out << s2out.str() << endl;
    f2out.close();
    s2out.clear();
    s2out.str("");
  }
  */
}

void DrawPSDFilter(UInt_t cobo) {
  Double_t ratiocut=0.0397+0.00412; // centroid + 1sigma for good data
  Int_t maxValue = -100000;
  Int_t corrwaveformintegral = 0;
  Double_t psdratio = 0;
  Int_t background = 0;
  for(UInt_t asad=0; asad<maxasad; asad++) {
    for(UInt_t aget=0; aget<4; aget++) {
      if(!rwaveforms[cobo]->hasHit[asad*4+aget]) continue; // skip aget that did not fire
      for(UInt_t chan=0; chan<68; chan++) {
	if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
        if(!rwaveforms[cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

        if((asad==2||asad==3)&&(aget==0||aget==1)&&(rwaveforms[cobo]->hasOverflow)){
	  maxValue = rwaveforms[cobo]->energy[asad*4+aget][chan];
	  corrwaveformintegral = rwaveforms[cobo]->PSDIntegral[asad*4+aget][chan];
	  psdratio = rwaveforms[cobo]->PSDRatio[asad*4+aget][chan];
	  //((maxValue<300 && psdratio>0.05) || (maxValue>300 && psdratio>ratiocut))){
          if(rwaveforms[cobo]->isOverflow[asad*4+aget][chan]){
	    hMM_PSDIntegralvsMax[0]->Fill(maxValue,corrwaveformintegral);
	    hMM_PSDRatiovsMax[0]->Fill(maxValue,psdratio);
	    hMM_PSDRatio[0]->Fill(0.0,psdratio);
	  }else{
	    hMM_PSDIntegralvsMax[1]->Fill(maxValue,corrwaveformintegral);
	    hMM_PSDRatiovsMax[1]->Fill(maxValue,psdratio);
	    hMM_PSDRatio[0]->Fill(1.0,psdratio);
	  }
        }else if((asad==2||asad==3)&&(aget==0||aget==1)){
	  maxValue = rwaveforms[cobo]->energy[asad*4+aget][chan];
	  background = rwaveforms[cobo]->background[asad*4+aget][chan];
	  corrwaveformintegral = rwaveforms[cobo]->PSDIntegral[asad*4+aget][chan];
	  psdratio = rwaveforms[cobo]->PSDRatio[asad*4+aget][chan];
	  hMM_PSDIntegralvsMax[2]->Fill(maxValue,corrwaveformintegral);
	  hMM_PSDRatiovsMax[2]->Fill(maxValue,psdratio);
	  hMM_PSDRatio[0]->Fill(2.0,psdratio);
        }
      }
    }
  }
}

void FillTrack() {
  UInt_t dchan = 0;
  UInt_t spxidx = 0;
  UInt_t spxidy = 0;
  int rgidx=0;

  for(UInt_t cobo=0; cobo<2; cobo++) {
    for(UInt_t asad=0; asad<maxasad; asad++) {
      for(UInt_t aget=0; aget<4; aget++) {
        if(!rwaveforms[cobo]->hasHit[asad*4+aget]) continue; // skip agets that did not have any fire
        for(UInt_t chan=0; chan<68; chan++) {
	  if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
          if(!rwaveforms[cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

          if(cobo==0 && rwaveforms[cobo]->energy[asad*4+aget][chan]>mm_minenergy && rwaveforms[cobo]->energy[asad*4+aget][chan]<mm_maxenergy && rwaveforms[cobo]->time[asad*4+aget][chan]>mm_mintime && rwaveforms[cobo]->time[asad*4+aget][chan]<mm_maxtime){
	    if(chan<11) {
	        dchan = chan;
	    }else if(chan>11 && chan<22) {
	        dchan = chan - 1;
	    }else if(chan>22 && chan<45) {
	        dchan = chan - 2;
	    }else if(chan>45 && chan<56) {
	        dchan = chan - 3;
	    }else if(chan>56) {
	        dchan = chan - 4;
	    }

            //if((asad==2||asad==3)&&(aget==0||aget==1)&&(rwaveforms[cobo]->hasOverflow[asad*4+aget])){
            if(rwaveforms[cobo]->hasOverflow){
	      mm_tracks->hasOverflow=true;
	    }
	    mm_tracks->hasTrack+=1;
	    spxidx = mapchantomm->pxidx[asad][aget][dchan];
	    spxidy = mapchantomm->pxidy[asad][aget][dchan];
	    if(spxidx!=-1 && spxidy!=-1){
	      mm_tracks->pixel[spxidx][spxidy]+=1;
	      mm_tracks->time[spxidx][spxidy]=rwaveforms[cobo]->time[asad*4+aget][chan];
	      mm_tracks->energy[spxidx][spxidy]=rwaveforms[cobo]->energy[asad*4+aget][chan];
	      mm_tracks->posx[spxidx][spxidy]=spxidx*3.4-67*3.4 + 3.4/2;
	      mm_tracks->posy[spxidx][spxidy]=spxidy*1.7 + 1.7/2;
	      mm_tracks->posxerr[spxidx][spxidy]=3.4/2;
	      mm_tracks->posyerr[spxidx][spxidy]=1.7/2;
	      mm_tracks->poszerr[spxidx][spxidy]=fTimePerBin*driftv/2;
	      rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[spxidx][spxidy];
	      ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[spxidx][spxidy];
	      rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[spxidx][spxidy];
	      //if(spxidx>63 && spxidx<68) cout << spxidx << " " << spxidy << ", posx=" <<  mm_tracks->posx[spxidx][spxidy] << ", posy=" << mm_tracks->posy[spxidx][spxidy] << endl;
	      if(asad<2) mm_tracks->coloridx[spxidx][spxidy]=1;
	      else mm_tracks->coloridx[spxidx][spxidy]=2;

	      if(mm_tracks->energy[spxidx][spxidy]>0){
	        posenergysum[spxidy]+=(mm_tracks->posx[spxidx][spxidy]+rx)*(mm_tracks->energy[spxidx][spxidy]);
	        energysum[spxidy]+=mm_tracks->energy[spxidx][spxidy];
	        sumcounter[spxidy]++;
	      }
	    }
	    if(spxidx==-1){ //strip
	      mm_tracks->hasTrackStrip+=1;
	      if(asad==2){
  	        for(int i=0;i<64;i++){
	          mm_tracks->pixel[i][spxidy]+=1;
	          mm_tracks->time[i][spxidy]=rwaveforms[cobo]->time[asad*4+aget][chan];
	          mm_tracks->energy[i][spxidy]=rwaveforms[cobo]->energy[asad*4+aget][chan];
	          mm_tracks->posx[i][spxidy]=i*1.7-64*1.7-3*3.4 + 1.7/2;
	          mm_tracks->posy[i][spxidy]=spxidy*1.7 + 1.7/2;
	          mm_tracks->posxerr[i][spxidy]=1.7/2;
	          mm_tracks->posyerr[i][spxidy]=1.7/2;
	          mm_tracks->poszerr[i][spxidy]=fTimePerBin*driftv/2;
	          mm_tracks->coloridx[i][spxidy]=3;
	          rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[i][spxidy];
	          ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[i][spxidy];
	          rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[i][spxidy];
	        }
	        rgidx=1;
	        hMM_TrackCounter[rgidx]++;
	        hMM_dE[rgidx]+=rwaveforms[cobo]->energy[asad*4+aget][chan];
	      }else if(asad==3){
  	        for(int i=0;i<64;i++){
	          mm_tracks->pixel[i+70][spxidy]+=1;
	          mm_tracks->time[i+70][spxidy]=rwaveforms[cobo]->time[asad*4+aget][chan];
	          mm_tracks->energy[i+70][spxidy]=rwaveforms[cobo]->energy[asad*4+aget][chan];
	          mm_tracks->posx[i+70][spxidy]=i*1.7 + 3*3.4 + 1.7/2;
	          mm_tracks->posy[i+70][spxidy]=spxidy*1.7 + 1.7/2;
	          mm_tracks->posxerr[i+70][spxidy]=1.7/2;
	          mm_tracks->posyerr[i+70][spxidy]=1.7/2;
	          mm_tracks->poszerr[i+70][spxidy]=fTimePerBin*driftv/2;
	          mm_tracks->coloridx[i+70][spxidy]=3;
	          rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[i+70][spxidy];
	          ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[i+70][spxidy];
	          rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[i+70][spxidy];
	        }
	        rgidx=2;
	        hMM_TrackCounter[rgidx]++;
	        hMM_dE[rgidx]+=rwaveforms[cobo]->energy[asad*4+aget][chan];
	      }
	    }
	    if(spxidy==-1){ //chain
	      mm_tracks->hasTrackChain+=1;
	      for(int i=0;i<64;i++){
	        mm_tracks->pixel[spxidx][i*2+1]+=1;
	        mm_tracks->time[spxidx][i*2+1]=rwaveforms[cobo]->time[asad*4+aget][chan];
	        mm_tracks->energy[spxidx][i*2+1]=rwaveforms[cobo]->energy[asad*4+aget][chan];
	        if(spxidx<64) mm_tracks->posx[spxidx][i*2+1]=(spxidx-64)*1.7 - 3*3.4 + 1.7/2;
	        else mm_tracks->posx[spxidx][i*2+1]=(spxidx-70)*1.7 + 3*3.4 + 1.7/2;
	        mm_tracks->posy[spxidx][i*2+1]=(i*2+1)*1.7;
	        mm_tracks->posxerr[spxidx][i*2+1]=1.7/2;
	        mm_tracks->posyerr[spxidx][i*2+1]=1.7/2;
	        mm_tracks->poszerr[spxidx][i*2+1]=fTimePerBin*driftv/2;
	        mm_tracks->coloridx[spxidx][i*2+1]=4;
	        rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[spxidx][i*2+1];
	        ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[spxidx][i*2+1];
	        rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[spxidx][i*2+1];
	      }
	    }
	  }
        }
      }
    }
  }
}

void CleanTrack(){
  if(mm_tracks->hasTrack>0){
    for(int i=0;i<150;i++){ // pixel_x
      for(int j=0;j<140;j++){ // pixel_y
	if((i<64 || i>=70) && mm_tracks->pixel[i][j]>0){
	  if(mm_tracks->pixel[i][j+1]>0){
	    //if(mm_tracks->time[i][j]<=((mm_tracks->time[i][j+1])+1) &&
	       //mm_tracks->time[i][j]>=((mm_tracks->time[i][j+1])-1)){ //for slow drift velocity
	    if(mm_tracks->time[i][j]==mm_tracks->time[i][j+1]){ //for fast drift velocity
	      mm_tracks->pixel[i][j] = mm_tracks->pixel[i][j];
	    }else{
	      mm_tracks->pixel[i][j] = 0;
	    }
	  }else{
	    mm_tracks->pixel[i][j] = 0;
	  }
	}
      }
    }
  }
}

void FillDecayFlag() {
  UInt_t dchan = 0;
  UInt_t spxidx = 0;
  UInt_t spxidy = 0;
  for(UInt_t cobo=0; cobo<2; cobo++) {
    for(UInt_t asad=0; asad<maxasad; asad++) {
      for(UInt_t aget=0; aget<4; aget++) {
        for(UInt_t chan=0; chan<68; chan++) {
          if(rwaveforms[cobo]->isDecay[asad*4+aget][chan]==1){
	    //cout << "FDF:" << cobo << " " << asad << " " << aget << " " << chan << " " << rwaveforms[cobo]->isDecay[asad*4+aget][chan] << endl;
	    mm_tracks->hasDecay=true;
	    if(chan<11) {
	        dchan = chan;
	    }else if(chan>11 && chan<22) {
	        dchan = chan - 1;
	    }else if(chan>22 && chan<45) {
	        dchan = chan - 2;
	    }else if(chan>45 && chan<56) {
	        dchan = chan - 3;
	    }else if(chan>56) {
	        dchan = chan - 4;
	    }

	    spxidx = mapchantomm->pxidx[asad][aget][dchan];
	    spxidy = mapchantomm->pxidy[asad][aget][dchan];
	    if(spxidx!=-1 && spxidy!=-1 && mm_tracks->pixel[spxidx][spxidy]>0){
	      mm_tracks->decay[spxidx][spxidy]=1;
	    }
	    if(spxidx==-1){ //strip
	      if(asad==2){
  	        for(int i=0;i<64;i++){
	          if(mm_tracks->pixel[i][spxidy]>0){
		    mm_tracks->decay[i][spxidy]=1;
	          }
	        }
	      }else if(asad==3){
  	        for(int i=0;i<64;i++){
	          if(mm_tracks->pixel[i+70][spxidy]>0){
		    mm_tracks->decay[i+70][spxidy]=1;
	          }
	        }
	      }
	    }
	    if(spxidy==-1){ //chain
	      for(int i=0;i<64;i++){
	        if(mm_tracks->pixel[spxidx][i*2+1]>0){
		  mm_tracks->decay[spxidx][i*2+1]=1;
	        }
	      }
	    }
	  }
        }
      }
    }
  }
}

void GetTrackPosYLimit(){
  int rgidx=0;
  if(mm_tracks->hasTrack>0){
    for(rgidx=0;rgidx<3;rgidx++){
      trackposymin[rgidx] = 10000;
      trackposymax[rgidx] = -10000;
    }
    for(int i=0;i<150;i++){
      for(int j=0;j<170;j++){
	if(mm_tracks->pixel[i][j]>0){
	  if(i>63 && i<71) rgidx=0;
	  if(i<64) rgidx=1;
	  if(i>70) rgidx=2;
	  if(trackposymin[rgidx]>mm_tracks->posy[i][j]) trackposymin[rgidx] = mm_tracks->posy[i][j];
	  if(trackposymax[rgidx]<mm_tracks->posy[i][j]) trackposymax[rgidx] = mm_tracks->posy[i][j];
        }
      }
    }
  }
}

void HoughTransform(){
  double theta;
  double costheta;
  double sintheta;
  double radiusxy;
  double radiusxt;
  double radiusyt;
  int rgidx=0;
  int avgflag=0;
  if(mm_tracks->hasTrack>0){
  //if(mm_tracks->hasTrack>0 && mm_tracks->hasTrackChain>10 && mm_tracks->hasTrackStrip>10){
    for(int i=0;i<150;i++){
      for(int j=0;j<170;j++){
	if(mm_tracks->pixel[i][j]>0){
	  if(i>63 && i<71) rgidx=0;
	  if(i<64) rgidx=1;
	  if(i>70) rgidx=2;
	  for(theta=0;theta<180;theta+=0.1){
	    costheta = TMath::Cos(theta*TMath::DegToRad());
	    sintheta = TMath::Sin(theta*TMath::DegToRad());
	    if(rgidx==0){
	      radiusxy = (mm_tracks->avgposx[j])*costheta+(mm_tracks->avgposy[j])*sintheta;
	    }else if(rgidx>0){
	      radiusxy = (mm_tracks->posx[i][j])*costheta+(mm_tracks->posy[i][j])*sintheta;
	    }
	    if(radiusxy>-300 && radiusxy<300) hMM_TrackPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusxy);
	    radiusxt = (mm_tracks->posx[i][j])*costheta+(mm_tracks->time[i][j])*sintheta;
	    if(radiusxt>-400 && radiusxt<600) hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusxt);
	    radiusyt = (mm_tracks->posy[i][j])*costheta+(mm_tracks->time[i][j])*sintheta;
	    if(radiusyt>-400 && radiusyt<600) hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusyt);
	  }
	}
      }
    }
  }
}

void GetXYZTrack(){
  int maxbinxy=0;
  int maxbinxt=0;
  int maxbinyt=0;
  int xbin=0;
  int ybin=0;
  int zbin=0;
  double thetaxy=0;
  double thetaxt=0;
  double thetayt=0;
  double costhetaxy=0;
  double sinthetaxy=0;
  double costhetaxt=0;
  double sinthetaxt=0;
  double costhetayt=0;
  double sinthetayt=0;
  double radiusxy=0;
  double radiusxt=0;
  double radiusyt=0;
  double stddevxy=0;
  double minstddevxy=100000;
  double minthetayt=0;
  double stddevyt=0;
  double minstddevyt=100000;
  double minthetaxt=0;
  double stddevxt=0;
  double minstddevxt=100000;
  int entries=0;
  double maxr=0;
  double maxtheta=0;
  int binmax=0;
  int maximum=0;
  //int minentries[3]={20,10,10};
  int minentries[3]={0,0,0};
  if(mm_tracks->hasTrack>0){
  //if(mm_tracks->hasTrack>0 && mm_tracks->hasTrackChain>10 && mm_tracks->hasTrackStrip>10){
    for(int rgidx=0;rgidx<3;rgidx++){
      minstddevxy=100000;
      minstddevyt=100000;
      minstddevxt=100000;
      entries = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",450,451)->GetEntries();
      if(entries>minentries[rgidx]){
        for(int i=0;i<1800;i++){
          if(hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetEntries()==entries){
            stddevxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
	      if(minstddevxy>=stddevxy){
	      minstddevxy = stddevxy;
              thetaxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
              ybin = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
              radiusxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(ybin);
            }
  
            stddevxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
            if(minstddevxt>stddevxt){
	      minstddevxt = stddevxt;
              thetaxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
              xbin = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
              radiusxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(xbin);
            }
            stddevyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
            if(minstddevyt>stddevyt){
	      minstddevyt = stddevyt;
              thetayt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
              ybin = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
              radiusyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(ybin);
            }
          }
        }
  
        costhetaxy = TMath::Cos(thetaxy*TMath::DegToRad());
        sinthetaxy = TMath::Sin(thetaxy*TMath::DegToRad());
        double ifzeroxy = sinthetaxy*costhetaxy;
	if(ifzeroxy==0){
	  thetaxy+=0.1;
          costhetaxy = TMath::Cos(thetaxy*TMath::DegToRad());
          sinthetaxy = TMath::Sin(thetaxy*TMath::DegToRad());
	}
        maxbinxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetMaximumBin();
        costhetaxt = TMath::Cos(thetaxt*TMath::DegToRad());
        sinthetaxt = TMath::Sin(thetaxt*TMath::DegToRad());
        double ifzeroxt = sinthetaxt*costhetaxt;
	if(ifzeroxt==0){
	  thetaxt+=0.1;
          costhetaxt = TMath::Cos(thetaxt*TMath::DegToRad());
          sinthetaxt = TMath::Sin(thetaxt*TMath::DegToRad());
	}
        maxbinyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetMaximumBin();
        costhetayt = TMath::Cos(thetayt*TMath::DegToRad());
        sinthetayt = TMath::Sin(thetayt*TMath::DegToRad());
        double ifzeroyt = sinthetayt*costhetayt;
	if(ifzeroyt==0){
	  thetayt+=0.1;
          costhetayt = TMath::Cos(thetayt*TMath::DegToRad());
          sinthetayt = TMath::Sin(thetayt*TMath::DegToRad());
	}
        int i = 0; // pixel_x
        double posxt = 0; //pos_xt 
        double posx = 0; //pos_x 
        double posy = 0; //pos_y
        double posz = 0; //pos_z
        double timemm = 0; // time_mm
        double timez = 0; // time_z
        double timesi = 0; // time_si
	int mmcounts=0;
	int maxmmcounts=0;
        posx = (radiusxy+300*sinthetaxy)/costhetaxy;
	if(rgidx==0 && TMath::Abs(posx)>5){
          for(int i=65;i<69;i++){ // pixel_x
	    mmcounts=0;
            for(int j=0;j<128;j++){ // pixel_y
	      if(mm_tracks->pixel[i][j]>0) mmcounts++;
            }
	    if(mmcounts>maxmmcounts){
	      radiusxy = (i-67)*3.4+3.4/2;
	      maxmmcounts = mmcounts;
            }
          }
	  thetaxy=0;
	  costhetaxy = 1;
	  sinthetaxy = 0;
	}
        for(int j=-300;j<400;j++){ // pos_y (mm)
          //posxt = (radiusxt - sinthetaxt*(radiusyt-j*costhetayt)/sinthetayt)/costhetaxt;
          posx = (radiusxy-j*sinthetaxy)/costhetaxy;
          hMM_a0->Fill(costhetaxy/sinthetaxy*-1);
          hMM_a1->Fill(radiusxy/sinthetaxy);
          timemm = (radiusyt-j*costhetayt)/sinthetayt*fTimePerBin;
          hMM_b0->Fill(costhetayt/sinthetayt*-1)*fTimePerBin;
          hMM_b1->Fill(radiusyt/sinthetayt)*fTimePerBin;
          timez = timemm - timesi;
          posz = 70 - (timez * driftv);
          //if(j==0) cout << a0 << " " << a1 << " " << posx << " " << b0 << " " << b1 << " " << timemm << " " << timesi << " " << timez << " " << posz << endl;
	  if(j>=-164 && j<-100) hMM_TrackPosXZ[(164+j)]->Fill(posx,posz);
          if(j>=-205 && j<-195) hMM_TrackPosXZAll->Fill(posx,posz);
          for(int k=0;k<10;k++){
       	    //hMM_TrackPos[goodevtcounter%16]->Fill(posxt,j);
       	    //hMM_TrackPosXY[goodevtcounter%16]->Fill(posxt,j);
          }
	  if(j<=trackposymax[rgidx]){
	    hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
            hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
            hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);
	    /*
	    if(j<0 && ((rgidx==0) || (posx<0 && rgidx==1) || (posx>0 && rgidx==2))){
	      hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
              hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
              hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);
	      
	      //if(posxt<150){
              //  hMM_TrackPosXYTAll->Fill(posxt,j);
              //  timemm = (radiusxt-posxt*costhetaxt)/sinthetaxt;
	      //  hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(posxt,timemm);
	      //}
	    }else if(j>=0 && ((posx>=-10.2 && posx<=10.2 && rgidx==0) || (posx<-10.2 && rgidx==1) || (posx>10.2 && rgidx==2))){
              hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
              hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
              hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);
	      //if(posxt<150){
              //  hMM_TrackPosXYTAll->Fill(posxt,j);
              //  timemm = (radiusxt-posxt*costhetaxt)/sinthetaxt;
	      //  hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(posxt,timemm);
	      //}
	    }
	    */
	  }
          hMM_TrackPosXYAll->Fill(posx,j);
          hMM_TrackPosYZAll->Fill(j,posz);
          timemm = (radiusyt-j*costhetayt)/sinthetayt;
	  hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(j,timemm);
          //hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez);
          //hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez);
          //hMM_TrackXZ[goodevtcounter%16]->Fill(i,posz);
          //hMM_TrackYZ[goodevtcounter%16]->Fill(j,posz);
        }
      }
    }
  }
}
  
void FilldEvsE(){
  if(mm_tracks->hasTrack>0){
    for(int j=0;j<170;j++){
      if(sumcounter[j]>0){
        mm_tracks->avgposx[j] = posenergysum[j]/energysum[j];
        mm_tracks->sumenergy[j] = energysum[j];
      }
    }
    for(int i=64;i<71;i++){
      for(int j=0;j<170;j++){
	if(mm_tracks->pixel[i][j]>0){
	  mm_tracks->avgposy[j] = mm_tracks->posy[i][j];
        }
      }
    }
  }
}

void Sum2pEnergy(){
  int rgidx=0;
  if(mm_tracks->hasTrack>0 && (mm_tracks->hasDecay)){
    for(int i=0;i<150;i++){
      for(int j=0;j<170;j++){
        //Micromega
	if(mm_tracks->decay[i][j]>0){
          if(i>63 && i<71) rgidx=0;
          if(i<64) rgidx=1;
          if(i>70) rgidx=2;
	  mm_tracks->sum2penergy[rgidx] += mm_tracks->energy[i][j];
	  mm_tracks->sum2penergy[3] += mm_tracks->energy[i][j];
	} 
      } 
    } 
    //cout << "sum = " << mm_tracks->sum2penergy << endl;;
    for(int i=0;i<4;i++){
      if(mm_tracks->sum2penergy[i]>0){
        hMM_Sum2pEnergy[i]->Fill(mm_tracks->sum2penergy[i]);
      } 
    } 
    if(mm_tracks->sum2penergy[3]>0){
        //hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X [SumE(bc,bl,br,all)=%d,%d,%d,%d];Cell ID Y [EventNo=%d, goodevtcounter=%d]",mm_tracks->sum2penergy[0],mm_tracks->sum2penergy[1],mm_tracks->sum2penergy[2],mm_tracks->sum2penergy[3],reventIdx,goodevtcounter));
        hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X;Cell ID Y [EventNo=%d, goodevtcounter=%d]",reventIdx,goodevtcounter));
    } 
  } 
} 

void DrawTrack(){
  int rgidx=0;
  if(mm_tracks->hasTrack>0){
    for(int i=0;i<150;i++){
      for(int j=0;j<170;j++){
        //Micromega
	if(mm_tracks->pixel[i][j]>0){
          if(i>63 && i<71) rgidx=0;
          if(i<64) rgidx=1;
          if(i>70) rgidx=2;
	  hMM_Time[goodevtcounter%16]->Fill(mm_tracks->time[i][j]);
	  hMM_Energy[goodevtcounter%16]->Fill(mm_tracks->energy[i][j]);
	  hMM_TrackAll->Fill(i,j);
	  hMM_TrackPosAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry);
	  hMM_TimevsPxIDXPosAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j]);
	  hMM_TimevsPxIDYPosAll->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j]);
	  hMM_TrackvsE[goodevtcounter%16]->Fill(i,j,mm_tracks->energy[i][j]);
	  hMM_Track[goodevtcounter%16]->Fill(i,j,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	  hMM_TrackPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	  if(i<64 || i>69){
	    if(j%2!=0){ //chain
	      hMM_TimevsPxIDX[goodevtcounter%16]->Fill(i,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	      hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    }else{ //strip
	      hMM_TimevsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	      hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    }
	  }else{
	    hMM_TimevsPxIDX[goodevtcounter%16]->Fill(i,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    hMM_TimevsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
	    hMM_EnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j],mm_tracks->coloridx[i][j]);
	    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j]);
	  }
	  if(i>=64 && i<70){
	    hMM_EnergyvsPxIDYALL->Fill(j,mm_tracks->energy[i][j]);
	  }
	}
      }
    }
/*
    for(int i=64;i<71;i++){
      for(int j=0;j<170;j++){
	if(mm_tracks->pixel[i][j]>0){
	  hMM_TrackPos[goodevtcounter%16]->Fill(mm_tracks->avgposx[j],mm_tracks->posy[i][j],20);
        }
      }
    }
*/
  }
}
void DrawTrack2pMode(){
  int rgidx=0;
  int gtrackidx=0;
  double posz=0;
  double timez=0;
  double poszoffset=88; //tb=32 is the beam position=0
  if(mm_tracks->hasTrack>0){
    gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx++,-150,-300,-150);
    gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx++,150,400,150);
    for(int i=0;i<150;i++){
      for(int j=0;j<170;j++){
        //Micromega
	if(mm_tracks->pixel[i][j]>0){
          if(i>63 && i<71) rgidx=0;
          if(i<64) rgidx=1;
          if(i>70) rgidx=2;
	  if(mm_tracks->posx[i][j]>-150 && mm_tracks->posx[i][j]<150){
            posz = (poszoffset-(mm_tracks->time[i][j]-mm_tracks->decay[i][j]*256))*fTimePerBin*driftv;
	    timez = bucketmax - (mm_tracks->time[i][j]-mm_tracks->decay[i][j]*256);
	    //cout << "time=" << mm_tracks->time[i][j] << ", posz=" << posz << endl;
	    hMM_TrackPosXYAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry);
	    hMM_TrackPosYZAll->Fill(mm_tracks->posy[i][j]+ry,posz+rz);
	    hMM_TrackPosXZAll->Fill(mm_tracks->posx[i][j]+rx,posz+rz);
	    if(i<64 || i>69){
	      if(j%2!=0){ //chain
	        hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	        hMM_TrackPosXZ[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      }else{ //strip
	        hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	        hMM_TrackPosYZ[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      }
	    }else{
	      hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      hMM_TrackPosXZ[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      hMM_TrackPosYZ[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	    }
	    if(mm_tracks->decay[i][j]==1){
	      hMM_TrackDecay[goodevtcounter%16]->Fill(i,j,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
	      gMM_TrackDecay[goodevtcounter%16]->SetPoint(gtrackidx,i,j,(int)(posz/1.7));
	      gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx,mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry,posz+rz);
	      gtrackidx++;
	    }
	  }
	}
      }
    }
  }
}

void InitTrack(){
  rpos = new TRandom();
  ResetTrack();
}

void ResetdEvsE(){
  for(int i=0;i<3;i++){
    hMM_TrackCounter[i]=0;
    hMM_dE[i]=0;
    hMM_E[i]=0;
  }
}

void DrawdEvsE(){
  int rgidx=0;
  for(int j=112;j<128;j++){
    if(mm_tracks->sumenergy[j]>0){
      hMM_TrackCounter[rgidx]++;
      hMM_dE[rgidx]+=mm_tracks->sumenergy[j];
    }
  }
  for(rgidx=0;rgidx<3;rgidx++){
    if(hMM_E[rgidx]>0 && hMM_dE[rgidx]>0){
      hMM_dE[rgidx]=hMM_dE[rgidx]/hMM_TrackCounter[rgidx];
      hMM_TrackdEvsEALL[rgidx]->Fill(hMM_E[rgidx],hMM_dE[rgidx]);
    }
  }
}

void ResetTrack(){
  rx = 0;
  ry = 0;
  rz = 0;
  rt = 0;
  mm_tracks->hasTrack=0;
  mm_tracks->hasTrackChain=0;
  mm_tracks->hasTrackStrip=0;
  mm_tracks->hasProjectile=0;
  mm_tracks->hasEjectile=0;
  mm_tracks->hasRecoil=0;
  mm_tracks->hasDecay=false;
  for(int i=0;i<4;i++) mm_tracks->sum2penergy[i] = 0;
  for(int i=0;i<150;i++){
    for(int j=0;j<170;j++){
      mm_tracks->decay[i][j]=0;
      mm_tracks->pixel[i][j]=0;
      mm_tracks->energy[i][j]=0;
      mm_tracks->time[i][j]=0;
      mm_tracks->coloridx[i][j]=0;
    }
  }
  for(int j=0;j<170;j++){
    mm_tracks->avgposx[j] = 0;
    mm_tracks->sumenergy[j] = 0;
    posenergysum[j]=0;
    energysum[j]=0;
    sumcounter[j]=0;
  }
}

void ResetTrackHist(){
  gMM_TrackDecay[goodevtcounter%16]->Clear();
  gMM_TrackDecayPos[goodevtcounter%16]->Clear();
  for(int i=0;i<150;i++){
    for(int j=0;j<170;j++){
      hMM_Track[goodevtcounter%16]->SetBinContent(i,j,0);
      hMM_TrackvsE[goodevtcounter%16]->SetBinContent(i,j,0);
    }
    for(int j=0;j<512;j++){
      hMM_TimevsPxIDX[goodevtcounter%16]->SetBinContent(i,j,0);
      hMM_TimevsPxIDY[goodevtcounter%16]->SetBinContent(i,j,0);
      hMM_TrackXZ[goodevtcounter%16]->SetBinContent(i,j,0);
    }
    for(int j=0;j<4000;j++){
      hMM_EnergyvsPxIDY[goodevtcounter%16]->SetBinContent(i,j,0);
    }
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetBinContent(i,0);
  }
  for(int i=0;i<700;i++){
    for(int j=0;j<512;j++){
      hMM_TimevsPxIDYPos[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
  for(int i=0;i<300;i++){
    for(int j=0;j<512;j++){
      hMM_TimevsPxIDXPos[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
  for(int i=0;i<170;i++){
    for(int j=0;j<512;j++){
      hMM_TrackYZ[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
  for(int i=0;i<300;i++){
    for(int j=0;j<700;j++){
      hMM_TrackPos[goodevtcounter%16]->SetBinContent(i,j,0);
      hMM_TrackPosXY[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
  for(int i=0;i<700;i++){
    for(int j=0;j<400;j++){
      hMM_TrackPosYZ[goodevtcounter%16]->SetBinContent(i,j,0);
    }
  }
  for(int rgidx=0;rgidx<3;rgidx++){
    for(int i=0;i<512;i++){
      for(int j=0;j<512;j++){
	//hMM_TrackdEvsE[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
      }
    }
    for(int i=0;i<1800;i++){
      for(int j=0;j<600;j++){
        hMM_TrackPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
      }
      for(int j=0;j<1000;j++){
        hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
        hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
      }
    }
  }
  hMM_Track[goodevtcounter%16]->SetTitle(Form("Single Event Track by channels;Cell ID X [D2PTime=%d usec];Cell ID Y [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TrackDecay[goodevtcounter%16]->SetTitle(Form("Single Event Decay Track by channels;Cell ID X [D2PTime=%d usec];Cell ID Y [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TrackPos[goodevtcounter%16]->SetTitle(Form("Single Event Track by Position;Pos X (mm) [D2PTime=%d usec];Pos Y (mm) [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TrackPosXY[goodevtcounter%16]->SetTitle(Form("Single Event Track Pos Y vs Pos X;Pos X (mm) [D2PTime=%d usec];Pos Y (mm) [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TimevsPxIDX[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel X;Cell ID X [D2PTime=%d usec];Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TimevsPxIDXPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos X;Pos X (mm) [D2PTime=%d usec];Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TimevsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel Y;Cell ID Y [D2PTime=%d usec];Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_TimevsPxIDYPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos Y;Pos Y (mm) [D2PTime=%d usec];Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_EnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Energy vs Pixel Y;Cell ID Y [D2PTime=%d usec];Energy [FrameNo~%d, EventNo=%d, goodevtcounter=%d]",(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter));
  hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Sum Energy vs Pixel Y;Cell ID Y;Energy [EventNo=%d]",reventIdx));
  hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetFillStyle(3001);
  hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetFillColor(kMagenta+3);
  hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetLineColor(1);
  hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetOption("hist");
  for(int j=0;j<512;j++){
    hMM_Time[goodevtcounter%16]->SetBinContent(j,0);
  }
  for(int j=0;j<4000;j++){
    hMM_Energy[goodevtcounter%16]->SetBinContent(j,0);
  }
  hMM_Time[goodevtcounter%16]->SetTitle(Form("Single Event Time;Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d, D2PTime=%d usec]",(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter,(rd2ptime/1000)));
  hMM_Energy[goodevtcounter%16]->SetTitle(Form("Single Event Energy;Energy [FrameNo~%d, EventNo=%d, goodevtcounter=%d, D2PTime=%d usec]",(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter,(rd2ptime/1000)));
  /*
  for(int j=0;j<170;j++){
    pMM_Track[goodevtcounter%16]->SetBinContent(j,0);
  }
  */
}

void SetScaler(int flag){
  enablescaler = flag;
}

void RootROpenFile(string & inputFileName, string & inputTreeName){
    size_t f = inputFileName.find("online");
    if(f!=std::string::npos) gErrorIgnoreLevel = kFatal;
    fInputFile = new TFile(inputFileName.c_str(),"read");
    fInputTree = (TTree*)fInputFile->Get(inputTreeName.c_str());
    fInputFileSize = fInputFile->GetSize();
    if(f!=std::string::npos) gErrorIgnoreLevel = kError;
}

void RootRInit(){
  fInputTree->SetBranchAddress("mmMul",&rGETMul);
  fInputTree->SetBranchAddress("mmHit",&rGETHit);
  fInputTree->SetBranchAddress("mmEventIdx",&rGETEventIdx);
  fInputTree->SetBranchAddress("mmD2PTime",&rGETD2PTime);
  fInputTree->SetBranchAddress("mmFrameNo",rGETFrameNo);
  fInputTree->SetBranchAddress("mmDecayNo",rGETDecayNo);
  fInputTree->SetBranchAddress("mmCobo",rGETCobo);
  fInputTree->SetBranchAddress("mmAsad",rGETAsad);
  fInputTree->SetBranchAddress("mmAget",rGETAget);
  fInputTree->SetBranchAddress("mmChan",rGETChan);
  fInputTree->SetBranchAddress("mmTime",rGETTime);
  fInputTree->SetBranchAddress("mmEnergy",rGETEnergy);
  fInputTree->SetBranchAddress("mmWaveformX",rGETWaveformX);
  fInputTree->SetBranchAddress("mmWaveformY",rGETWaveformY);
  /*
  TBranch *branch;
  branch = fInputTree->GetBranch("mmMul");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmHit");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmEventIdx");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmD2PTime");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmFrameNo");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmDecayNo");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmCobo");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmAsad");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmAget");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmChan");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmTime");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmEnergy");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmWaveformX");
  branch->SetAutoDelete(kTRUE);
  branch = fInputTree->GetBranch("mmWaveformY");
  branch->SetAutoDelete(kTRUE);
  */
  rGETMul = 0;
  rGETHit = 0;
  rGETEventIdx = 0;
  rGETD2PTime = 0;
  for(int i=0;i<4352;i++){
    rGETFrameNo[i] = 0;
    rGETDecayNo[i] = 0;
    rGETEnergy[i] = 0;
    rGETCobo[i] = 0;
    rGETAsad[i] = 0;
    rGETAget[i] = 0;
    rGETChan[i] = 0;
    for(int j=0;j<512;j++){
      rGETWaveformX[i][j] = 0;
      rGETWaveformY[i][j] = 0;
    }
  }
  RootRInitWaveforms();
}

void RootRInitWaveforms() {
  for(int cobo=0;cobo<2;cobo++){
    rwaveforms[cobo]->hasSignal.resize(maxasad*4);
    rwaveforms[cobo]->hasHit.resize(maxasad*4);
    rwaveforms[cobo]->hasFPN.resize(maxasad*4);
    rwaveforms[cobo]->doneFPN.resize(maxasad*4);
    rwaveforms[cobo]->isOverflow.resize(maxasad*4);
    rwaveforms[cobo]->hasOverflow = false;
    rwaveforms[cobo]->isDecay.resize(maxasad*4);
    rwaveforms[cobo]->hasImplant = false;
    rwaveforms[cobo]->hasDecay = false;
    rwaveforms[cobo]->coboIdx = 0;
    rwaveforms[cobo]->asadIdx = 0;
    rwaveforms[cobo]->waveform.resize(maxasad*4);
    rwaveforms[cobo]->corrwaveform.resize(maxasad*4);
    rwaveforms[cobo]->fpnwaveform.resize(maxasad*4);
    rwaveforms[cobo]->energy.resize(maxasad*4);
    rwaveforms[cobo]->time.resize(maxasad*4);
    rwaveforms[cobo]->background.resize(maxasad*4);
    rwaveforms[cobo]->PSDIntegral.resize(maxasad*4);
    rwaveforms[cobo]->PSDRatio.resize(maxasad*4);
    for(int i=0;i<maxasad;i++){
      for(int j=0;j<4;j++){
	rwaveforms[cobo]->hasSignal[i*4+j].resize(68);
	rwaveforms[cobo]->isOverflow[i*4+j].resize(68);
	rwaveforms[cobo]->isDecay[i*4+j].resize(68);
        rwaveforms[cobo]->waveform[i*4+j].resize(68);
        rwaveforms[cobo]->corrwaveform[i*4+j].resize(68);
	rwaveforms[cobo]->energy[i*4+j].resize(68);
	rwaveforms[cobo]->time[i*4+j].resize(68);
	rwaveforms[cobo]->background[i*4+j].resize(68);
	rwaveforms[cobo]->PSDIntegral[i*4+j].resize(68);
	rwaveforms[cobo]->PSDRatio[i*4+j].resize(68);
        for(int k=0;k<68;k++){
	  rwaveforms[cobo]->waveform[i*4+j][k].resize(bucketmax);
	  rwaveforms[cobo]->corrwaveform[i*4+j][k].resize(bucketmax);
        }
	rwaveforms[cobo]->fpnwaveform[i*4+j].resize(bucketmax);
      }
    }
  }
}

void RootRResetWaveforms() {
  for(int cobo=0;cobo<2;cobo++){
    rwaveforms[cobo]->hasOverflow = false;
    rwaveforms[cobo]->hasImplant = false;
    rwaveforms[cobo]->hasDecay = false;
    for(int i=0;i<maxasad;i++){
      for(int j=0;j<4;j++){
        if(rwaveforms[cobo]->hasHit[i*4+j] || rwaveforms[cobo]->hasFPN[i*4+j]){
          for(int k=0;k<68;k++){
	    if(!rwaveforms[cobo]->hasSignal[i*4+j][k]) continue;
            rwaveforms[cobo]->hasSignal[i*4+j][k] = false;
            rwaveforms[cobo]->isOverflow[i*4+j][k] = false;
            rwaveforms[cobo]->isDecay[i*4+j][k] = 0;
            rwaveforms[cobo]->energy[i*4+j][k] = 0.0;
            rwaveforms[cobo]->time[i*4+j][k] = 0.0;
            rwaveforms[cobo]->background[i*4+j][k] = 0.0;
            rwaveforms[cobo]->PSDIntegral[i*4+j][k] = 0.0;
            rwaveforms[cobo]->PSDRatio[i*4+j][k] = 0.0;
	    fill(rwaveforms[cobo]->waveform[i*4+j][k].begin(),
	         rwaveforms[cobo]->waveform[i*4+j][k].end(),0);
	    fill(rwaveforms[cobo]->corrwaveform[i*4+j][k].begin(),
	         rwaveforms[cobo]->corrwaveform[i*4+j][k].end(),0);
          }
        }
        if(rwaveforms[cobo]->hasFPN[i*4+j]){
	  fill(rwaveforms[cobo]->fpnwaveform[i*4+j].begin(),
	       rwaveforms[cobo]->fpnwaveform[i*4+j].end(),0);
        }
        rwaveforms[cobo]->hasHit[i*4+j] = false;
        rwaveforms[cobo]->hasFPN[i*4+j] = false;
        rwaveforms[cobo]->doneFPN[i*4+j] = false;
      }
    }
  }

  if(enabledraww==1){
    for(int i=0;i<512;i++){
      for(int j=0;j<512;j++){
        hWaveFormbyEvent[evtcounter%16]->SetBinContent(i,j,0);
        hCorrWaveFormbyEvent[evtcounter%16]->SetBinContent(i,j,0);
      }
    }
  }
}

void RootRReset(){
  rGETHit = 0;
  rGETEventIdx = 0;
  rGETD2PTime = 0;
  for(int i=0;i<=rGETMul;i++){
    rGETFrameNo[i] = 0;
    rGETDecayNo[i] = 0;
    rGETEnergy[i] = 0;
    rGETCobo[i] = 0;
    rGETAsad[i] = 0;
    rGETAget[i] = 0;
    rGETChan[i] = 0;
    for(int j=0;j<512;j++){
      rGETWaveformX[i][j] = 0;
      rGETWaveformY[i][j] = 0;
    }
  }
  rGETMul = 0;
}

void RootRCloseFile(){
  fInputFile->cd();
  fInputFile->Close();
}

void RootHOpenFile(string & outputFileName){
  fhOutputFile = new TFile(outputFileName.c_str(),"recreate");
}

void RootHWrite(){
  fhOutputFile->cd();
 for(int i=0; i<16 ; i++) {
   hWaveFormbyEvent[i]->Write();
   hCorrWaveFormbyEvent[i]->Write();
 }

 // Histograms for particle tracks in Micromega
 hMM_TrackPosAll->Write();
 hMM_TrackAll->Write();
 hMM_TrackPosXYAll->Write();
 hMM_TrackPosXYTAll->Write();
 hMM_TrackPosXZAll->Write();
 hMM_TrackPosYZAll->Write();
 hMM_TimevsPxIDXPosAll->Write();
 hMM_TimevsPxIDYPosAll->Write();
 hMM_EnergyvsPxIDYALL->Write();
 //tdir->mkdir("hMM_EnergyvsPxIDYALL,Form("Tracks_Energy2D"));
 for(int i=0;i<16;i++){
   hMM_Track[i]->Write();
   hMM_TrackPos[i]->Write();
   hMM_TrackvsE[i]->Write();
   for(int j=0;j<3;j++){
     hMM_TrackPosHough[i][j]->Write();
     //hMM_TrackdEvsE[i][j]->Write();
   }
   hMM_TrackXZ[i]->Write();
   hMM_TrackYZ[i]->Write();
   hMM_TrackPosXY[i]->Write();
   hMM_TrackPosXZ[i]->Write();
   hMM_TrackPosYZ[i]->Write();
   hMM_TimevsPxIDX[i]->Write();
   hMM_TimevsPxIDXPos[i]->Write();
   for(int j=0;j<3;j++){
     hMM_TimevsPxIDXPosHough[i][j]->Write();
     hMM_TimevsPxIDYPosHough[i][j]->Write();
   }
   hMM_SumEnergyvsPxIDY[i]->Write();
   hMM_TimevsPxIDY[i]->Write();
 }
 for(int j=0;j<3;j++){
     hMM_TrackdEvsEALL[j]->Write();
 }

 // Histograms for the GET electronics raw channels
 hGET_HitPattern->Write();
 hGET_EHitPattern2D->Write();
 hGET_THitPattern2D->Write();
}

void RootHCloseFile(){
  fhOutputFile->Close();
}

void SetHWrite(int flag){
  enablehwrite = flag;
}

int GetHWrite(){
  return enablehwrite;
}
void SetBucketSize(int BucketSize){
  bucketmax = BucketSize;
}

void SetRootConverter(int flag){
  enableroot = flag;
}

void SetDrawWaveform(int flag){
  enabledraww = flag;
}

void ReadMapChanToMM(string filename){
  ifstream MapEData;
  UInt_t asadid;
  UInt_t agetid;
  UInt_t dchanid;
  UInt_t spxidx;
  UInt_t spxidy;

  MapEData.open(filename.data());
  if(MapEData.fail()==true){
    cerr<<"The ChanToMM_Map file wasn't opened!"<<endl;
  }
  while(MapEData.good()){
    MapEData >> asadid >> agetid >> dchanid >> spxidx >> spxidy;
    mapchantomm->pxidx[asadid][agetid][dchanid]=spxidx;
    mapchantomm->pxidy[asadid][agetid][dchanid]=spxidy;
  }
  MapEData.close();
}

void ReadResponseWaveform(string filename){
  ifstream ResponseData;
  UInt_t type;
  UInt_t timebucket;
  UInt_t amplitude;

  ResponseData.open(filename.data());
  if(ResponseData.fail()==true){
    cerr<<"The ResponseData file wasn't opened!"<<endl;
  }
  while(ResponseData.good()){
    ResponseData >> type >> timebucket >> amplitude;
    response[type][timebucket] = amplitude;
    
  }
  ResponseData.close();
}

void SetResponseSample(Int_t type, Int_t evtno, Int_t buckcut, Int_t buckwidth, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t rep, Int_t iter, Int_t boost){
  responsesample[type][0] = evtno;
  responsesample[type][1] = cobo;
  responsesample[type][2] = asad;
  responsesample[type][3] = aget;
  responsesample[type][4] = chan;
  responsesample[type][5] = buckcut;
  responsesample[type][6] = buckwidth;
  deconvrep[type] = rep;
  deconviter[type] = iter;
  deconvboost[type] = boost;
  if(maxrespsample<=type) maxrespsample = type;
}

Bool_t IsResponseSample(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan){
  if(responsesample[type][1] == cobo && responsesample[type][2] == asad && responsesample[type][3] == aget && responsesample[type][4] == chan) return true;
  return false;
}

void SetResponseWaveform(){
  Int_t coboIdx = 0;
  Int_t asadIdx = 0;
  Int_t agetIdx = 0;
  Int_t chanIdx = 0;
  Int_t buckbegin;
  Int_t buckend;
  Int_t maxreventIdx=0;
  for(int j=0;j<=maxrespsample;j++){
      if(maxreventIdx <= responsesample[j][0]){
        maxreventIdx = responsesample[j][0];
      }
  }
  for(int i=0;i<fNumberEvents;i++){
    fInputFile->cd();
    RootRInit();
    fInputTree->GetEntry(i);
    reventIdx = rGETEventIdx;
    if(reventIdx > maxreventIdx) return;
    for(int j=0;j<=maxrespsample;j++){
      if(responsesample[j][0] == reventIdx){
        for(Int_t i=0; i<rGETMul; i++){
          coboIdx = rGETCobo[i];
          asadIdx = rGETAsad[i];
          agetIdx = rGETAget[i];
          chanIdx = rGETChan[i];
          if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)){
           	rwaveforms[coboIdx]->hasFPN[asadIdx*4+agetIdx] = true;
          }
          for(int j=0;j<bucketmax;j++){
            rwaveforms[coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
          }
        }
    
        for(Int_t i=0; i<rGETMul; i++){
          coboIdx = rGETCobo[i];
          asadIdx = rGETAsad[i];
          agetIdx = rGETAget[i];
          chanIdx = rGETChan[i];
          if(chanIdx!=11 && chanIdx!=22 && chanIdx!=45 && chanIdx!=56){ // We want to skip the FPN channels.
            if(IsResponseSample(j,coboIdx,asadIdx,agetIdx,chanIdx)){
              buckbegin = responsesample[j][5];
              buckend = responsesample[j][5] + responsesample[j][6];
              GetAverageFPN(coboIdx,asadIdx,agetIdx);
              GetCorrWaveform(coboIdx,asadIdx,agetIdx,chanIdx);
              for(int buck=buckbegin;buck<buckend;buck++){
                response[j][buck-buckbegin] = rwaveforms[coboIdx]->corrwaveform[asadIdx*4+agetIdx][chanIdx][buck];
              }
            }
          }
        }
      }
    }
  }
}

void GetMaxResponseWaveform(){
  Int_t maxamplitude;
  for(int i=0;i<=maxrespsample;i++){
    maxamplitude=-10000;
    for(int buck=0;buck<512;buck++){
      if(maxamplitude<response[i][buck]){
        maxamplitude = response[i][buck];
        maxresponse[i] = buck;
      }
    }
  }
}

void GetSigmaResponseWaveform(){
  for(int i=0;i<=maxrespsample;i++){
    hResponse[i] = new TH1D(Form("hResponse_%d",i),Form("hResponse_%d",i),512,0,512);
    if(i==0){
      fResponse[i] = new TF1(Form("fResponse_%d",i), ShaperF_GET, 0, 512, 6);
      fResponse[i]->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
      fResponse[i]->SetParameters(10, 500, 80, 15, 3, 0.2);
      fResponse[i]->SetParLimits(1, 0, 10000);
      fResponse[i]->SetParLimits(3, 0, 50);
    }else if(i>0) {
      fResponse[i] = new TF1(Form("fResponse_%d",i), ShaperF_MSCF, 0, 512, 6);
      fResponse[i]->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
      fResponse[i]->SetParameters(10, 2000, 80, 25, 3, 0.2);
      fResponse[i]->SetParLimits(1, 0, 10000);
      fResponse[i]->SetParLimits(3, 0, 50);
    }
/*
    fResponse[i] = new TF1(Form("fResponse_%d",i),conv, 0,512);
    for(int ii=0; ii<p2[0]; ii++) {
      f_conv->FixParameter(ii, p2[ii]);
    }
*/

    for(int buck=0;buck<512;buck++){
      hResponse[i]->Fill(buck,response[i][buck]);
    }
    hResponse[i]->Fit(Form("fResponse_%d",i),"QN"); //Q=Quiet and N=Do not draw
    responsesigma[i] = fResponse[i]->GetParameter(3);
    cout << "Sigma[" << i << "] = " << responsesigma[i] << endl;
  }
}

Double_t ShaperF_GET(Double_t *x, Double_t *p) {
   Double_t semiGaus = p[1]*22.68113723*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);
   return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
/*
  Double_t val=0.;
  for(int i=0; i<(p[0]-2)/2; i++) {
    if(!(x[0]<p[2*i+2+1] || x[0]>1000000.)) 
      val += p[2*i+2] * 22.68113723 * exp(-3.*(x[0]-p[2*i+2+1])/p[1]) * sin((x[0]-p[2*i+2+1])/p[1]) * pow((x[0]-p[2*i+2+1])/p[1], 3);
  }
  return(val);
*/
}

Double_t ShaperF_MSCF(Double_t *x, Double_t *p) {
   Double_t semiGaus = p[1]*21.928*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);
   return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
}
