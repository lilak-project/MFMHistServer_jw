//#define DB_LOCATION
#include <stdio.h>      /* for printf() and fprintf() */
#include <stdlib.h>     /* for malloc */
#include <ncurses.h>    /* for getch() */
#include <sys/socket.h> /* for socket(), connect(), send(), and recv() */
#include <arpa/inet.h>  /* for sockaddr_in and inet_addr() */
#include <string.h>     /* for memset() */
#include <unistd.h>     /* for close() */
#include <errno.h>
#include "TH1.h"
#include "TH1.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include <stdlib.h>
#include <unistd.h>
#include "mfm/FrameBuilder.h"
#include "HistServer.h"
#include "json/json.h"
#include <sstream>
#include <fstream>
#include <TSystem.h>
#include <thread>

using namespace std;

static HistServer* histServer;
static HistServer* convServer;
static int rsize_buffer=512;
static string inrfname;
static int maxinrfidx;
static const int maxfileno=2000;
static string listinrfname[maxfileno];
static string inrtname;
static string inoutrfname;
static string listinoutrfname[maxfileno];
static string maxinoutrfidx;
static string inoutrtname;
static char keyinput='r';
static thread tgetkey;
void getkey (){
  while(1){
    keyinput = cin.get();
  }
}
static thread tupdatehistogram;
void updatehistogram (){
  int prevcounter=0;
  int counter=0;
  bool IsConverting=true;
  while(keyinput=='r'){
    if(inrfname=="list"){
      for(int i=0;i<maxinrfidx;i++){
	IsConverting=true;
        prevcounter=0;
	inrtname = "TEvent";
        cout << "Reading TTree(" << inrtname << ") and TFile(" << listinrfname[i] << ").." << endl;
        while(IsConverting){
          histServer->RootROpenFile(listinrfname[i], inrtname);
	  if(histServer->GetReadMode()==2){
            counter = (histServer->GetRootFileSize())/(rsize_buffer*1000);
            if(counter>20 && counter>prevcounter){
              prevcounter = counter;
              histServer->RootReadEvent();
              //histServer->RootFindEvent();
	    }
	  }else if(histServer->GetReadMode()==4){
            counter = (histServer->GetRootFileSize())/(rsize_buffer);
            if(counter>20 && counter>prevcounter){
              prevcounter = counter;
              histServer->RootReadHistEvent();
	    }
	  }
	  if(!(histServer->RootRRecovered())){
	    IsConverting=false;
            cout << "Sleep for 1s to start reading a next root file." << endl;
            sleep(1);
            histServer->RootRCloseFile();
	  }else{
            histServer->RootRCloseFile();
	  }
        }
      }
      inrfname="done";
    }else if(inrfname=="done"){
      cout << "Data Analysis completed. Sleeping for a minute." << endl;
      sleep(60);
    }else{
      cout << "Reading TTree(" << inrtname << ") and TFile(" << inrfname << ").." << endl;
      while(IsConverting){
        histServer->RootROpenFile(inrfname, inrtname);
	if(histServer->GetReadMode()==2){
          counter = (histServer->GetRootFileSize())/(rsize_buffer*1000);
	  cout << "counter=" << counter<< endl;
	  cout << "prevcounter=" << counter<< endl;
          if(counter>prevcounter){
            prevcounter = counter;
            histServer->RootReadEvent();
            //histServer->RootFindEvent();
	  }
	}else if(histServer->GetReadMode()==4){
          counter = (histServer->GetRootFileSize())/(rsize_buffer);
          if(counter>prevcounter){
            histServer->RootReadHistEvent();
	  }
        }
	if(!(histServer->RootRRecovered())) IsConverting=false;
        histServer->RootRCloseFile();
      }
      inrfname="done";
    }
  }
}

static thread tupdateenergytime;
void updateenergytime (){
  int prevcounter=0;
  int counter=0;
  bool IsConverting=true;
  while(keyinput=='r'){
    if(inrfname=="list"){
      for(int i=0;i<maxinrfidx;i++){
	IsConverting=true;
        prevcounter=0;
        cout << "Reading TTree(" << inrtname << ") and TFile(" << listinrfname[i] << ").." << endl;
	cout << "Writing TTree(" << inoutrtname << ") and TFile(" << listinoutrfname[i] << ").." << endl;
        convServer->Init(3,0);
        while(IsConverting){
          convServer->RootROpenFile(listinrfname[i], inrtname);
          convServer->RootRWOpenFile(listinoutrfname[i], inoutrtname);
          counter = (convServer->GetRootFileSize())/(rsize_buffer*1000);
          if(counter>20 && counter>prevcounter){
            prevcounter = counter;
            convServer->RootReadWriteEvent();
	  }
	  if(!(convServer->RootRRecovered())){
	    IsConverting=false;
            cout << "Sleep for 3s to start reading a next root file." << endl;
            sleep(3);
            convServer->RootRCloseFile();
            convServer->RootRWCloseFile();
	  }else{
            convServer->RootRCloseFile();
            convServer->RootRWCloseFile();
	  }
        }
      }
      inrfname="done";
    }else if(inrfname=="done"){
      cout << "Data Analysis completed." << endl;
      keyinput='q';
    }else{
      cout << "Reading TTree(" << inrtname << ") and TFile(" << inrfname << ").." << endl;
      while(IsConverting){
        convServer->RootROpenFile(inrfname, inrtname);
        convServer->RootRWOpenFile(inoutrfname, inoutrtname);
        counter = (convServer->GetRootFileSize())/(rsize_buffer*1000);
        if(counter>prevcounter){
          prevcounter = counter;
          convServer->RootReadWriteEvent();
        }
	if(!(convServer->RootRRecovered())) IsConverting=false;
        convServer->RootRCloseFile();
        convServer->RootRWCloseFile();
      }
      inrfname="done";
    }
  }
  delete convServer;
  return;
}


string _executeShellCommand(string command) {
  char buffer[256];
  string result = "";
  const char * cmd = command.c_str();
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw runtime_error("popen() failed!");

  try {
    while (!feof(pipe)) 
      if (fgets(buffer, 128, pipe) != NULL)
        result += buffer;
  } catch (...) {
      pclose(pipe);
      throw;
  }
  pclose(pipe);
  return result;
}

bool checkport(int usPort){
  string shellCommand;
  string pcPort;
  //sprintf(shellCommand, "netstat -lntu | awk '{print $4}' | grep ':' | cut -d \":\" -f 2 | sort | uniq | grep %hu", usPort);
  shellCommand = "netstat -an | grep -i " + to_string(usPort);
  pcPort = ":"+to_string(usPort);

  string output =  _executeShellCommand(shellCommand);
  //cout << shellCommand << endl;
  //cout << output << endl;

  if(output.find(pcPort) == string::npos)
    return false;
  else
    return true;
}

int main (int argc, char **argv)
{
  int begin, end;
  int maxit = 0;
  int currit = 0;
  int size_file = 0;
  int percent = 0;

  string configFileName;
  string cfgFileName;
  string smode;
  if(argc<3) {
    cout << "Usage: MFMHistServer configfile mode (mode=0 for reading mfm and plot, mode=1 for mfm2rootwf, mode=2 for reading rootwf, mode=3 for rootwf2rootet and mode=4 for reading rootet)" << endl;
    return 0;
  } else {
    configFileName = string(argv[1]);
    smode = string(argv[2]);
  }

  //Removing comments
  cfgFileName = configFileName+".nc";
  ifstream configStream(configFileName.c_str());
  ofstream ocfgStream(cfgFileName.c_str());
  string strcomment;
  string configline;
  size_t pos = 0;
  //while(configStream.good()){
  //  configStream >> configline;
  while(getline(configStream,configline)){
    strcomment="//";
    pos = configline.find(strcomment);
    if(pos!=std::string::npos){
      configline.erase(pos, configline.length()-pos);
    }
    strcomment="$";
    pos = configline.find(strcomment);
    if (pos!=std::string::npos){
      configline.erase(pos, configline.length()-pos);
    }
    strcomment="#";
    pos = configline.find(strcomment);
    if (pos!=std::string::npos){
      configline.erase(pos, configline.length()-pos);
    }
    strcomment=", ";
    pos = configline.find(strcomment);
    if (pos!=std::string::npos){
      configline.erase(pos+1, 1);
    }
    if(configline.length()>0){
      ocfgStream << configline << endl;
    }
  }
  configStream.close();
  ocfgStream.close();

  int mode = atoi(smode.c_str());
  Json::Value config;
  ifstream icfgStream(cfgFileName.c_str());
  icfgStream >> config;

  string sreadtype = config["ReadType"].asString();
  int readtype = atoi(sreadtype.c_str());
  string mfmfilename = config["MFMFileName"].asString();
  string outputFileName = config["OutputFileName"].asString();
  string snumfiles = config["NumberofFiles"].asString();
  int numfiles = atoi(snumfiles.c_str());
  int currfiles = 1;
  string watcherIP = config["watcherIP"].asString();
  string swatcherPort = config["watcherPort"].asString();
  string sConverterPort = config["ConverterPort"].asString();
  string sCoBoServerPort = config["CoBoServerPort"].asString();
  string sMutantServerPort = config["MutantServerPort"].asString();
  int watcherPort = atoi(swatcherPort.c_str());
  int converterPort = atoi(sConverterPort.c_str());
  int CoBoServerPort = atoi(sCoBoServerPort.c_str());
  int MutantServerPort = atoi(sMutantServerPort.c_str());
  string sBucketSize = config["BucketSize"].asString();
  int BucketSize = atoi(sBucketSize.c_str());
  string mapChanToMM = config["ChanToMMMapFileName"].asString();
  string mapChanToSi = config["ChanToSiMapFileName"].asString();
  string sRootConvert = config["RootConvertEnable"].asString();
  string sEnergyMethod = config["EnergyFindingMethod"].asString();
  int energymethod = atoi(sEnergyMethod.c_str());
  string sReadRW = config["ReadResponseWaveformFlag"].asString();
  int readrw = atoi(sReadRW.c_str());
  string rwfilename = config["ResponseWaveformFileName"].asString();
  int RootConvert = atoi(sRootConvert.c_str());
  string sScalerMode = config["ScalerMode"].asString();
  int ScalerMode = atoi(sScalerMode.c_str());
  string sd2pMode = config["2pMode"].asString();
  int d2pMode = atoi(sd2pMode.c_str());
  string supdatefast = config["UpdateFast"].asString();
  int updatefast = atoi(supdatefast.c_str());
  string sIgnoreMM = config["IgnoreMicromegas"].asString();
  //int IgnoreMM = atoi(sIgnoreMM.c_str());
  int IgnoreMM = 1;
  string sDrawWaveform = config["DrawWaveformEnable"].asString();
  int DrawWaveform = atoi(sDrawWaveform.c_str());
  string sCleanTrack = config["CleanTrackEnable"].asString();
  int cleantrack = atoi(sCleanTrack.c_str());
  string sDrawTrack = config["DrawTrackEnable"].asString();
  int DrawTrack = atoi(sDrawTrack.c_str());
  string sSkipEvents = config["SkipEvents"].asString();
  int SkipEvents = atoi(sSkipEvents.c_str());
  string sfirstEventNo = config["firstEventNo"].asString();
  int firstEventNo=0;
  string goodEventList;
  if(SkipEvents==1) firstEventNo = atoi(sfirstEventNo.c_str());
  if(SkipEvents==2) goodEventList = config["firstEventNo"].asString();

  ifstream in;
  string infname;
  string minfname;
  string outrfname = outputFileName;
  string outrtname;
  string slistinfname;
  string listinfname[maxfileno];
  int maxinfidx=0;
  int infidx=0;

  if(mode==0){
    while(checkport(watcherPort)){
      cout << "The Histogram Server Port No. " << watcherPort << " is being used. Checking " << (watcherPort+1) << "." << endl;
      watcherPort++;
    }
    cout << "The Histogram Server Port No. " << watcherPort << " is a good one. Please use this port number in vigru." << endl;
    convServer = new HistServer(watcherPort);
    convServer->SetBucketSize(BucketSize);
    convServer->Init(mode,d2pMode);
    convServer->SetReadMode(mode);
    convServer->SetReadType(readtype);
    convServer->SetScaler(ScalerMode);
    convServer->Set2pMode(d2pMode);
    convServer->SetHistMode();
    convServer->SetUpdateSpeed(updatefast);
    convServer->ReadMapChanToMM(mapChanToMM);
    convServer->ReadMapChanToSi(mapChanToSi);
    convServer->ReadMapChanToX6();
    convServer->SetIgnoreMM(IgnoreMM);
    convServer->SetDrawWaveform(DrawWaveform);
    convServer->SetCleanTrack(cleantrack);
    convServer->SetDrawTrack(DrawTrack);
    convServer->SetSkipEvents(SkipEvents);
    if(SkipEvents==1) convServer->SetfirstEventNo(firstEventNo);
    if(SkipEvents==2) convServer->ReadGoodEventList(goodEventList);
    convServer->SetEnergyMethod(energymethod);
    convServer->SetReadRF(readrw);
    if(readrw!=1){
      cout << "ReadResponseWaveformFlag should be set to 1." << endl;
      return 0;
    }
  }else if(mode==1){
    convServer = new HistServer(converterPort);
    convServer->SetBucketSize(BucketSize);
    convServer->Init(mode,d2pMode);
    convServer->SetReadMode(mode);
    convServer->SetReadType(readtype);
    convServer->SetScaler(ScalerMode);
    convServer->Set2pMode(d2pMode);
    convServer->SetUpdateSpeed(updatefast);
  }else if(mode==2){
    while(checkport(watcherPort)){
      cout << "The Histogram Server Port No. " << watcherPort << " is being used. Checking " << (watcherPort+1) << "." << endl;
      watcherPort++;
    }
    cout << "The Histogram Server Port No. " << watcherPort << " is a good one. Please use this port number in vigru." << endl;
    histServer = new HistServer(watcherPort);
    histServer->SetBucketSize(BucketSize);
    histServer->Init(mode,d2pMode);
    histServer->SetReadMode(mode);
    histServer->SetReadType(readtype);
    histServer->SetScaler(ScalerMode);
    histServer->Set2pMode(d2pMode);
    histServer->SetHistMode();
    histServer->SetUpdateSpeed(updatefast);
    histServer->ReadMapChanToMM(mapChanToMM);
    histServer->ReadMapChanToSi(mapChanToSi);
    histServer->ReadMapChanToX6();
    histServer->SetIgnoreMM(IgnoreMM);
    histServer->SetDrawWaveform(DrawWaveform);
    histServer->SetCleanTrack(cleantrack);
    histServer->SetDrawTrack(DrawTrack);
    histServer->SetSkipEvents(SkipEvents);
    if(SkipEvents==1) histServer->SetfirstEventNo(firstEventNo);
    if(SkipEvents==2) histServer->ReadGoodEventList(goodEventList);
    histServer->SetEnergyMethod(energymethod);
    histServer->SetReadRF(readrw);
    if(readrw>0){
      histServer->ReadResponseWaveform(rwfilename);
    }else{
      histServer->SetResponseSample(0,295270,90,60,0,0,0,37,15,15,1); // Si (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
      histServer->SetResponseSample(2,38509,100,220,1,1,0,16,20,20,1); // IC (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
    }
  }else if(mode==3){
    convServer = new HistServer(converterPort-10);
    convServer->SetBucketSize(BucketSize);
    convServer->Init(mode,d2pMode);
    convServer->SetReadMode(mode);
    convServer->SetReadType(readtype);
    convServer->SetScaler(ScalerMode);
    convServer->Set2pMode(d2pMode);
    convServer->SetUpdateSpeed(updatefast);
    convServer->SetEnergyMethod(energymethod);
    convServer->SetReadRF(readrw);
    if(readrw>0){
      convServer->ReadResponseWaveform(rwfilename);
    }else{
      convServer->SetResponseSample(0,295270,90,60,0,0,0,37,15,15,1); // Si (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
      convServer->SetResponseSample(2,38509,100,220,1,1,0,16,20,20,1); // IC (type,evtno,buckcut,buckwidth,cobo,asad,aget,chan,decrep,deciter,decboost)
    }
  }else if(mode==4){
    while(checkport(watcherPort)){
      cout << "The Histogram Server Port No. " << watcherPort << " is being used. Checking " << (watcherPort+1) << "." << endl;
      watcherPort++;
    }
    cout << "The Histogram Server Port No. " << watcherPort << " is a good one. Please use this port number in vigru." << endl;
    histServer = new HistServer(watcherPort);
    histServer->SetBucketSize(BucketSize);
    histServer->Init(mode,d2pMode);
    histServer->SetReadMode(mode);
    histServer->SetHistMode();
    histServer->SetReadType(readtype);
    histServer->SetScaler(ScalerMode);
    histServer->Set2pMode(d2pMode);
    histServer->SetUpdateSpeed(updatefast);
    histServer->ReadMapChanToMM(mapChanToMM);
    histServer->ReadMapChanToSi(mapChanToSi);
    histServer->ReadMapChanToX6();
    histServer->SetIgnoreMM(IgnoreMM);
    histServer->SetDrawWaveform(DrawWaveform);
    histServer->SetCleanTrack(cleantrack);
    histServer->SetDrawTrack(DrawTrack);
    histServer->SetSkipEvents(SkipEvents);
    if(SkipEvents==1) histServer->SetfirstEventNo(firstEventNo);
    if(SkipEvents==2) histServer->ReadGoodEventList(goodEventList);
  }

  if(readtype==11){
    string delimiter = "/";
    size_t pos = mfmfilename.find(delimiter);
    if (pos==std::string::npos){
      cout << "You should input an experiment name and runfile name separated by '/' in mfmfilename." << endl;
      cout << "The runfile name should follow a format of 'run_#### (e.g. run_0102)." << endl;
      return 0;
    }else{
      string expname = mfmfilename.substr(0,pos);
      string runno = mfmfilename.substr(pos+1,mfmfilename.length()-pos);
      pos = runno.find("_");
      if (pos==std::string::npos){
        cout << "The runfile name should follow a format of 'run_#### (e.g. run_0123)." << endl;
	return 0;
      }
      string runpath="/mnt/CRIBdisk2/o14apf17/ganacq_manip/"+expname+"/acquisition/run/";
      //string runpath="/hdfs/data/"+expname+"/acquisition/run/";
      cout << runpath << " " << expname << " " << runno << endl;
      TSystemDirectory dir(runpath.c_str(),runpath.c_str());
      TList *files = dir.GetListOfFiles();
      files->Sort();
      if(files){
	TSystemFile *file;
	TString fname;
	TIter next(files);
        ofstream ofrunlist("temprunlist.txt",std::ofstream::out);
	while((file=(TSystemFile*)next())){
	  fname = file->GetName();
	  if(mode==0 || mode==1 || mode==2 || mode==3){
	    if(!file->IsDirectory() && fname.BeginsWith("run") && fname.Contains(runno.c_str(),TString::kIgnoreCase) && !fname.Contains(".root",TString::kIgnoreCase)){
	       ofrunlist << runpath << fname.Data() << endl;
	    }
	  }else if(mode==4){
	    if(!file->IsDirectory() && fname.BeginsWith("run") && fname.Contains(runno.c_str(),TString::kIgnoreCase) && fname.Contains("_et.root",TString::kIgnoreCase)){
	      fname = fname.ReplaceAll("_et.root","");
	       ofrunlist << runpath << fname.Data() << endl;
	    }
	  }
	}
	ofrunlist.close();
      }	
      readtype=3;
      mfmfilename="temprunlist.txt";
    }
  }
  if(readtype==3 || readtype==13){
    infname=mfmfilename;
    cout << "MFM list file name: " << infname << endl;
    in.open(infname.data());
    if(in.fail()==true){
      cerr<<"The MFM list file wasn't opened!"<<endl;
      return 0;
    }
    while(in.good()){
      in >> slistinfname;
      listinfname[maxinfidx]=slistinfname;
      maxinfidx++;
    }
    in.close();
    maxinfidx -= 2;
  }

#ifdef DB_LOCATION
  cout<< __FILE__ << " +" << __LINE__ << " #After Init" << endl;
#endif

  //if (1) {
  //  outrfname = "/home/cens-alpha-00/ejungwoo/mfm_converter/test1.root";
  //}
  //if(ScalerMode==0){
  //  outrfname = "/home/cens-alpha-00/MFMHistServer/online2.root";
  //  inoutrfname = "/home/cens-alpha-00/MFMHistServer/online2_et.root";
  //  //outrfname = "/data/grgroup/gr01/MFMHistServer/online/online.root";
  //  //inoutrfname = "/data/grgroup/gr01/MFMHistServer/online/online_et.root";
  //}else{
  //  outrfname = "/home/cens-alpha-00/MFMHistServer/scaler.root";
  //  inoutrfname = "/home/cens-alpha-00/MFMHistServer/scaler_et.root";
  //  //outrfname = "/data/grgroup/gr01/MFMHistServer/online/scaler.root";
  //  //inoutrfname = "/data/grgroup/gr01/MFMHistServer/online/scaler_et.root";
  //}
  outrtname = "TEvent";
  inoutrtname = "TEvent";

#ifdef DB_LOCATION
  cout<< __FILE__ << " +" << __LINE__ << " #After file name configuration" << endl;
#endif

  int sock;                    /* Socket descriptor */
  struct sockaddr_in servAddr; /* server address */
  int size_buffer;
  int size_to_recv;
  char *buffer;
  int received_data;
  int total_received_data;
  int coboflag=0;
  int onlcnt=0;
  int oflcnt=0;
  cout<<"MODE: "<<mode<<"\t\tREAD TYPE:\t"<<readtype<<"\tROOTCONVERT:\t"<<RootConvert<<endl;
  if(mode==0 || mode==1){
#ifdef DB_LOCATION
    cout<< __FILE__ << " +" << __LINE__ << " #mode 1" << endl;
#endif
    //cout<<"root file name="<< outrfname << ", root tree name=" << outrtname << endl;
    if(RootConvert==0) convServer->RootWOpenFile(outrfname, outrtname);
    while(keyinput=='r'){
      if(readtype==0){
#ifdef DB_LOCATION
      cout<< __FILE__ << " +" << __LINE__ << " #while(r) readtype==0" << endl;
#endif
        if(RootConvert==0 && onlcnt==0){
          convServer->RootWInit();
        }
        sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
        if (sock == 0) {
          printf ("socket call failed\n");
          return 1;
        }
   
        memset(&servAddr, 0, sizeof(servAddr));
        servAddr.sin_family = AF_INET;             /* Internet address family */
        servAddr.sin_addr.s_addr = inet_addr(watcherIP.c_str());   /* Server IP address */
        if(ScalerMode==0){
          if(coboflag!=0 || (coboflag==0 && size_buffer==0)){
            servAddr.sin_port = htons(MutantServerPort); /* Server port */
            coboflag=0;
          }else{
            servAddr.sin_port = htons(CoBoServerPort); /* Server port */
            coboflag++;
          }
        }else{
            servAddr.sin_port = htons(MutantServerPort); /* Server port */
        }
        if (connect(sock, (struct sockaddr *) &servAddr, sizeof(servAddr)) < 0) {
          printf ("connect call failed\n");
          return 2;
        }
   
	recv (sock, &size_buffer, sizeof (int), 0);
        buffer = (char *) malloc (size_buffer);
        total_received_data = 0;
        size_to_recv = size_buffer;
   
  
        while (total_received_data < size_buffer) {
          received_data = recv (sock, buffer, size_to_recv, 0);
          total_received_data += received_data;
          if (received_data == -1) {
            printf ("errno %d\n",errno);
            perror (NULL);
            return 3;
          }
   
          buffer += received_data;
          size_to_recv -= received_data;
   
        }
        close (sock);
   
        if(size_buffer>4) convServer->addDataChunk(buffer-total_received_data,buffer);
  
        if(onlcnt>50){
          cout << "Reading MFM online..." << endl;
          onlcnt=0;
        }
        onlcnt++;
   
      }else if(readtype==1 || readtype==2 || readtype==10){
#ifdef DB_LOCATION
      cout<< __FILE__ << " +" << __LINE__ << " #while(r) readtype==1 2 10" << endl;
#endif
        if(readtype==1 || readtype==2){
          if(numfiles==1){
            infname=mfmfilename;
            if(readtype==1) readtype=4;
            if(readtype==2) readtype=0;
          }else if(currfiles<numfiles){
            if(currfiles==1)  infname=mfmfilename;
            else  infname=mfmfilename+Form(".%d",(currfiles-1));
            currfiles++;
          }else if(currfiles==numfiles){
            infname=mfmfilename+Form(".%d",(currfiles-1));
            currfiles++;
            if(readtype==1) readtype=4;
            if(readtype==2) readtype=0;
          }
        }else if(readtype==10){
          infname=mfmfilename;
          size_t f = infname.find("/run/");
          infname.replace(f, std::string("/run/").length(), "/run/mutant/");
          f = infname.find("s.");
	  if (f!=std::string::npos){
            infname.replace(f+1, (infname.length()-f), "");
	  }
          cout << infname << endl;
          readtype=14;
        }
        cout << "MFM file name: " << infname << endl;
        in.open(infname.c_str(),std::ios::binary | std::ios::in);
        if(!in) {
          std::cout << "Could not open input file!" << std::endl; 
          return 0;
        }
	if(mode==0){
          convServer->RootWInit();
	}else if(mode==1){
#ifdef DB_LOCATION
          cout<< __FILE__ << " +" << __LINE__ << " #while(r) readtype==1 START  RootConvert=" << RootConvert << endl;
#endif
          if(RootConvert==0){
            convServer->RootWInit();
          }
          if(RootConvert==1){
            //outrfname = infname+".root";
            //outrfname = "test2.root";
	    //outrfname.replace(5,9,"CRIBdisk");
	    //outrfname.replace(5,5,"disk01");
            cout<<"root file name="<< outrfname << ", root tree name=" << outrtname << endl;
            convServer->RootWOpenFile(outrfname, outrtname);
            convServer->RootWInit();
            convServer->SetRootConverter(RootConvert);
          }
#ifdef DB_LOCATION
          cout<< __FILE__ << " +" << __LINE__ << " #while(r) readtype==1 END" << endl;
#endif
        }
	if(readtype==10){
          size_buffer = 32;
	}else{
          size_buffer = 512;
	}
        buffer = (char *) malloc (size_buffer);
        currit = 0;
        while(!in.eof()) {
          //std::cout << currit << endl;
          in.read(buffer,size_buffer);
          currit++;
        }
        maxit = currit;
        currit = 0;
        percent = maxit/10;
        in.close();
  
        if(d2pMode==1){
#ifdef DB_LOCATION
          cout<< __FILE__ << " +" << __LINE__ << " # d2pMode" << endl;
#endif
          minfname=mfmfilename;
          size_t f = minfname.find("/run/");
          minfname.replace(f, std::string("/run/").length(), "/run/mutant/");
          f = minfname.find("s.");
	  if (f!=std::string::npos){
            minfname.replace(f+1, (minfname.length()-f), "");
	  }
          in.open(minfname.c_str(),std::ios::binary | std::ios::in);
          while(!in.eof()) {
            in.read(buffer,size_buffer);
            if(!in.eof()) {
              try {
                convServer->addDataChunk(buffer,buffer+size_buffer);
              }catch (const std::exception& e){
                cout << e.what() << endl;
              }
            }else if(in.gcount()>0) {
              convServer->addDataChunk(buffer,buffer+in.gcount());
            }
          }
          in.close();
          cout << "Reading Mutant " << minfname << " done." << endl;
	}
	size_t const matrixSize = 4*68*512*sizeof(double);
	//size_t const matrixSize = 512;
	char *buffer = (char *) malloc (matrixSize);
	cout<<"READ BUFFERS, matrixSize was "<< 4*68*512*sizeof(double) <<endl;
        in.open(infname.c_str(),std::ios::binary | std::ios::in);
#ifdef DB_LOCATION
          cout<< __FILE__ << " +" << __LINE__ << " # before while(!in.eof())" << endl;
#endif
          while(!in.eof()) {
#ifdef DB_LOCATION
              cout<< __FILE__ << " +" << __LINE__ << " # INWHILE in while(!in.eof()) " << currit << endl;
#endif

              //in.read(buffer,size_buffer);
              int filebuffer=0;
              in.seekg(filebuffer, std::ios_base::cur);
              in.read(buffer,matrixSize);
              filebuffer += matrixSize;
              currit++;
              if(!in.eof()) {
                  try {
                      //convServer->addDataChunk(buffer,buffer+size_buffer);
                      convServer->addDataChunk(buffer,buffer+matrixSize);
                  }catch (const std::exception& e){
                      cout << e.what() << endl;
                      return 0;
                  }
                  if(oflcnt>50){
                      oflcnt=0;
                  }
                  oflcnt++;
                  if(currit%percent==0) cout << Form("%d%% (%d/%d) Processed..",(100*currit/maxit),currit,maxit) << endl;
              }
              else if(in.gcount()>0)
              {
                  try {
                      convServer->addDataChunk(buffer,buffer+in.gcount());
                  }catch (const std::exception& e){
                      cout << e.what() << endl;
                      return 0;
                  }
                  if(oflcnt>50){
                      oflcnt=0;
                  }
                  oflcnt++;
                  currit=maxit;
                  cout << Form("100%s (%d/%d) Processed..","%",currit,maxit) << endl;
                  if(mode==1){
                      //if(RootConvert==1) convServer->RootWCloseFile();
                      convServer->RootWCloseFile();
                  }
              }
          }
#ifdef DB_LOCATION
          cout<< __FILE__ << " +" << __LINE__ << " # end of while(!in.eof())" << endl;
#endif
        in.close();
        if(mode==1){
          if(readtype==14){
            delete convServer;
            return 0;
          }else if(RootConvert==1 && (currit==maxit)){
            if(numfiles==1 || currfiles>numfiles){
              delete convServer;
              return 0;
	    }
          }else if(readtype==4 && (currit==maxit)){
            delete convServer;
            return 0;
          }
        }
      }else if(readtype==3 || readtype==13){
        infname=listinfname[infidx];
        cout << "MFM file name: " << infname << endl;
        in.open(infname.c_str(),std::ios::binary | std::ios::in);
        if(!in) {
          std::cout << "Could not open input file!" << std::endl; 
          return 0;
        }
	if(mode==0){
          convServer->RootWInit();
	}else if(mode==1){
          if(RootConvert==0){
            convServer->RootWInit();
          }
          if(RootConvert==1){
            outrfname = infname+".root";
            if(atoi(&outrfname[13])>=1) outrfname.replace(13,1,"");
            cout << "outrfname: " << outrfname << endl;
	    //outrfname.replace(5,5,"disk01");
            convServer->RootWOpenFile(outrfname, outrtname);
            convServer->RootWInit();
            convServer->SetRootConverter(RootConvert);
          }
        }
        buffer = (char *) malloc (size_buffer);
        size_buffer = 512;
        currit = 0;
        while(!in.eof()) {
          in.read(buffer,size_buffer);
          currit++;
        }
        maxit = currit;
        currit = 0;
        percent = maxit/10;
        in.close();
  
        if(d2pMode==1){
	  minfname=listinfname[infidx];
          size_t f = minfname.find("/run/");
          minfname.replace(f, std::string("/run/").length(), "/run/mutant/");
          f = minfname.find("s.");
	  if (f!=std::string::npos){
            minfname.replace(f+1, (minfname.length()-f), "");
	  }
          in.open(minfname.c_str(),std::ios::binary | std::ios::in);
          while(!in.eof()) {
            in.read(buffer,size_buffer);
            if(!in.eof()) {
              try {
                convServer->addDataChunk(buffer,buffer+size_buffer);
              }catch (const std::exception& e){
                cout << e.what() << endl;
              }
            }else if(in.gcount()>0) {
              convServer->addDataChunk(buffer,buffer+in.gcount());
            }
          }
          in.close();
          cout << "Reading Mutant " << minfname << " done." << endl;
	}

        in.open(infname.c_str(),std::ios::binary | std::ios::in);
        while(!in.eof()) {
          in.read(buffer,size_buffer);
          currit++;
          if(!in.eof()) {
            convServer->addDataChunk(buffer,buffer+size_buffer);
            if(oflcnt>50){
              oflcnt=0;
            }
            oflcnt++;
            //if(currit%percent==0 && maxit!=0) cout << Form("Reading %d/%d file: %d%% (%d/%d) Processed..",(infidx+1),(maxinfidx+1),(100*currit/maxit),currit,maxit) << endl;
          } else if(in.gcount()>0) {
            convServer->addDataChunk(buffer,buffer+in.gcount());
            if(oflcnt>50){
              oflcnt=0;
            }
            oflcnt++;
            currit=maxit;
            cout << Form("Reading %d/%d file: 100%s (%d/%d) Processed..",(infidx+1),(maxinfidx+1),"%",currit,maxit) << endl;
	    if(mode==1){
              //if(RootConvert==1) convServer->RootWCloseFile();
              convServer->RootWCloseFile();
            }
          }
        }
        in.close();
        if(readtype==3 && infidx==maxinfidx){
          readtype=4;
        }else if(readtype==13 && infidx==maxinfidx){
          readtype=14;
        }else{
          infidx++;
        }
      }
      if(mode==1){
        if(readtype==14){
          delete convServer;
          return 0;
        }else if(readtype==4 && RootConvert==1){
          delete convServer;
          return 0;
        }
      }
    }
  
#ifdef DB_LOCATION
    cout<< __FILE__ << " +" << __LINE__ << " #end of mode 1" << endl;
#endif
    delete convServer;
    return 0;
  }else if(mode==2){
    if(readtype==0){
      inrfname = outrfname;
      inrtname = outrtname;
      tupdatehistogram=thread(updatehistogram);
    }else if(readtype==1 || readtype==2){
      if(RootConvert==0){
        inrfname = outrfname;
      }
      if(RootConvert==1){
        if(numfiles==1){
          inrfname = mfmfilename+".root";
	  //inrfname.replace(5,5,"disk01");
        }else{
          for(currfiles=0;currfiles<numfiles;currfiles++){
            if(currfiles==0) inrfname=mfmfilename+".root";
            else inrfname=mfmfilename+Form(".%d.root",currfiles);
	    //inrfname.replace(5,5,"disk01");
            listinrfname[currfiles]=inrfname;
          }
	  maxinrfidx=numfiles;
	  inrfname = "list";
        }
      }
      inrtname = outrtname;
      tupdatehistogram=thread(updatehistogram);
    }else if(readtype==3 || readtype==13){
      if(RootConvert==0){
        inrfname = outrfname;
      }
      if(RootConvert==1){
        inrfname = "list";
        for(infidx=0;infidx<(maxinfidx+1);infidx++){
	  //listinfname[infidx].replace(5,5,"disk01");
          listinrfname[infidx]=listinfname[infidx]+".root";
        }
        maxinrfidx=maxinfidx+1;
      }
      inrtname = outrtname;
      tupdatehistogram=thread(updatehistogram);
    }

    while(keyinput=='r'){
    }
  }else if(mode==3){
    if(readtype==0){
      inrfname = outrfname;
      inrtname = outrtname;
      tupdateenergytime=thread(updateenergytime);
    }else if(readtype==1 || readtype==2){
      if(RootConvert==0){
        inrfname = outrfname;
        inrtname = outrtname;
      }
      if(RootConvert==1){
        if(numfiles==1){
          inrfname = mfmfilename+".root";
          inoutrfname = mfmfilename+"_et.root";
	  //inrfname.replace(5,5,"disk01");
        }else{
          for(currfiles=0;currfiles<numfiles;currfiles++){
            if(currfiles==0) inrfname=mfmfilename+".root";
            else inrfname=mfmfilename+Form(".%d.root",currfiles);
	    //inrfname.replace(5,5,"disk01");
            listinrfname[currfiles]=inrfname;
            if(currfiles==0) inoutrfname=mfmfilename+"_et.root";
            else inoutrfname=mfmfilename+Form(".%d_et.root",currfiles);
	    //inrfname.replace(5,5,"disk01");
            listinoutrfname[currfiles]=inoutrfname;
          }
	  maxinrfidx=numfiles;
	  inrfname = "list";
	  inoutrfname = "list";
        }
      }
      inrtname = outrtname;
      tupdateenergytime=thread(updateenergytime);
    }else if(readtype==3 || readtype==13){
      if(RootConvert==0){
        inrfname = outrfname;
        inrtname = outrtname;
      }
      if(RootConvert==1){
        inrfname = "list";
        inoutrfname = "list";
        for(infidx=0;infidx<(maxinfidx+1);infidx++){
	  //listinfname[infidx].replace(5,5,"disk01");
          listinrfname[infidx]=listinfname[infidx]+".root";
          listinoutrfname[infidx]=listinfname[infidx]+"_et.root";
        }
        maxinrfidx=maxinfidx+1;
      }
      inrtname = outrtname;
      tupdateenergytime=thread(updateenergytime);
    }

    while(keyinput=='r'){
    }
    return 0;
  }else if(mode==4){
    if(readtype==0){
      inrfname = inoutrfname;
      inrtname = inoutrtname;
      tupdatehistogram=thread(updatehistogram);
    }else if(readtype==1 || readtype==2){
      if(RootConvert==0){
        inrfname = inoutrfname;
      }
      if(RootConvert==1){
        if(numfiles==1){
          inrfname = mfmfilename+"_et.root";
	  //inrfname.replace(5,5,"disk01");
        }else{
          for(currfiles=0;currfiles<numfiles;currfiles++){
            if(currfiles==0) inrfname=mfmfilename+".root";
            else inrfname=mfmfilename+Form(".%d_et.root",currfiles);
	    //inrfname.replace(5,5,"disk01");
            listinrfname[currfiles]=inrfname;
          }
	  maxinrfidx=numfiles;
	  inrfname = "list";
        }
      }
      inrtname = outrtname;
      tupdatehistogram=thread(updatehistogram);
    }else if(readtype==3 || readtype==13){
      if(RootConvert==0){
        inrfname = inoutrfname;
        inrtname = inoutrtname;
      }
      if(RootConvert==1){
        inrfname = "list";
        for(infidx=0;infidx<(maxinfidx+1);infidx++){
	  //listinfname[infidx].replace(5,5,"disk01");
          listinrfname[infidx]=listinfname[infidx]+"_et.root";
        }
        maxinrfidx=infidx+1;
      }
      inrtname = outrtname;
      tupdatehistogram=thread(updatehistogram);
    }

    while(keyinput=='r'){
    }
  }
}
