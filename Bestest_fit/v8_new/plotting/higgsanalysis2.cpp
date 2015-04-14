#ifndef scatter_cpp
#define scatter_cpp

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <istream>
#include <stdio.h>
#include <sstream>

#include "higgsanalysis.h"
#include "plotting.h"
//#include "style.h"

using namespace std;

int main(int argc, char* argv[])
{
	bool degen = false;
	
	if ( argc > 1 ) {
		string degtext;
		for ( int i=1; i<argc; i++) {
			degtext = argv[i];
			if ( degtext == "deg" || degtext == "degen" ) {
				degen = true;
			}
		}
	}
	
	
	
	ifstream datafile;
	
	if ( degen ) {
		datafile.open("Bestest_fit_deg_out.txt");
	}
	else {
		datafile.open("Bestest_fit_nb_out.txt");
	}
	
	
	
	string line="x";
	stringstream ss;
	int start=0,end=0;
	
	double ggFyySMATLAS=0, VBFyySMATLAS=0, VHyySMATLAS=0, ttHyySMATLAS=0;
	double ggFbbSMATLAS=0, VBFbbSMATLAS=0, VHbbSMATLAS=0, ttHbbSMATLAS=0;
	double ggFtataSMATLAS=0, VBFtataSMATLAS=0, VHtataSMATLAS=0, ttHtataSMATLAS=0;
	double ggFZZSMATLAS=0, VBFZZSMATLAS=0, VHZZSMATLAS=0, ttHZZSMATLAS=0;
	double ggFWWSMATLAS=0, VBFWWSMATLAS=0, VHWWSMATLAS=0, ttHWWSMATLAS=0;
	double ggFZySMATLAS=0, VBFZySMATLAS=0, VHZySMATLAS=0, ttHZySMATLAS=0;
	double ggFyySMCMS=0, VBFyySMCMS=0, VHyySMCMS=0, ttHyySMCMS=0;
	double ggFbbSMCMS=0, VBFbbSMCMS=0, VHbbSMCMS=0, ttHbbSMCMS=0;
	double ggFtataSMCMS=0, VBFtataSMCMS=0, VHtataSMCMS=0, ttHtataSMCMS=0;
	double ggFZZSMCMS=0, VBFZZSMCMS=0, VHZZSMCMS=0, ttHZZSMCMS=0;
	double ggFWWSMCMS=0, VBFWWSMCMS=0, VHWWSMCMS=0, ttHWWSMCMS=0;
	double ggFZySMCMS=0, VBFZySMCMS=0, VHZySMCMS=0, ttHZySMCMS=0;
	
	double muCMSyy1=0.78;
	double muCMSZZ1=0.91, muCMSZZ2=1.22;
	double muCMSWW1=0.76;
	double muCMSbb1=1.3;
	double muCMStata1=1.1, muCMStata2=1.4;
	double emuCMSyy1=sqrt(0.28*0.28+0.26*0.26)/2.0;
	double emuCMSZZ1=sqrt(0.30*0.30+0.24*0.24)/2.0, emuCMSZZ2=sqrt(0.84*0.84+0.57*0.57)/2.0;
	double emuCMSWW1=0.21;
	double emuCMSbb1=sqrt(0.7*0.7+0.6*0.6)/2.0;
	double emuCMStata1=0.4, emuCMStata2=0.6;
	
	double muATLASyy1=1.65;
	double muATLASZZ1=1.7, muATLASZZ2=1.5;
	double muATLASWW1=1.01, muATLASWW2=1.66, muATLASWW3=0.82;
	double muATLASbb1=-0.4;
	double muATLAStata1=0.7;
	double emuATLASyy1=sqrt(0.34*0.34+0.30*0.30)/2.0;
	double emuATLASZZ1=sqrt(0.5*0.5+0.4*0.4)/2.0, emuATLASZZ2=0.4;
	double emuATLASWW1=0.31, emuATLASWW2=0.79, emuATLASWW3=0.36;
	double emuATLASbb1=1.0;
	double emuATLAStata1=0.7;
	
	double muCMSyy1blh, muCMSZZ1blh, muCMSZZ2blh, muCMSWW1blh, muCMSbb1blh, muCMStata1blh, muCMStata2blh, muCMSZy1blh;
	double muATLASyy1blh, muATLASZZ1blh, muATLASWW1blh, muATLASWW2blh, muATLASWW3blh, muATLASbb1blh, muATLAStata1blh, muATLASZy1blh;
	vector<double> chisqCMS;
	vector<double> chisqATLAS;
	vector<double> chisq;
	
	double minATLAS, minCMS, min;
	int CMSlt1, CMSlt2, CMSlt3, CMSgt3;
	int ATLASlt1, ATLASlt2, ATLASlt3, ATLASgt3;
	int lt1, lt2, lt3, gt3;
	int iCMSlt1, iCMSlt2, iCMSlt3, iCMSgt3;
	int iATLASlt1, iATLASlt2, iATLASlt3, iATLASgt3;
	int ilt1, ilt2, ilt3, igt3;
	
	vector<double> rchisqCMS;
	vector<double> rchisqATLAS;
	vector<double> rchisq;
	
	double rminATLAS, rminCMS, rmin;
	int rCMSlt1, rCMSlt2, rCMSlt3, rCMSgt3;
	int rATLASlt1, rATLASlt2, rATLASlt3, rATLASgt3;
	int rlt1, rlt2, rlt3, rgt3;
	int riCMSlt1, riCMSlt2, riCMSlt3, riCMSgt3;
	int riATLASlt1, riATLASlt2, riATLASlt3, riATLASgt3;
	int rilt1, rilt2, rilt3, rigt3;
	
	vector<double> fl;
	vector<double> fh;
	vector<double> mh0;
	vector<double> ma0;
	vector<double> mhh;
	vector<double> tB;
	vector<double> tG;
	vector<double> t12;
	vector<double> t13;
	vector<double> sigfact;
	vector<double> kG;
	vector<double> kY;
	vector<double> kL;
	vector<double> rDSgg;
	vector<double> rSgg;
	vector<double> rDSyy;
	vector<double> rSyy;
	vector<double> iDSgg;
	vector<double> iSgg;
	vector<double> iDSyy;
	vector<double> iSyy;
	vector<double> DSgg;
	vector<double> Sgg;
	vector<double> DSyy;
	vector<double> Syy;
	vector<double> Cyy;
	vector<double> Cgg;
	vector<double> CZZ;
	vector<double> CWW;
	vector<double> Cbb;
	vector<double> Ctata;
	
	vector<double> deltam;
	
	vector<double> vmuCMSyy1blh;
	vector<double> vmuCMSZZ1blh;
	vector<double> vmuCMSZZ2blh;
	vector<double> vmuCMSWW1blh;
	vector<double> vmuCMSbb1blh;
	vector<double> vmuCMStata1blh;
	vector<double> vmuCMStata2blh;
	vector<double> vmuCMSZy1blh;
	
	vector<double> vmuATLASyy1blh;
	vector<double> vmuATLASZZ1blh;
	vector<double> vmuATLASWW1blh;
	vector<double> vmuATLASWW2blh;
	vector<double> vmuATLASWW3blh;
	vector<double> vmuATLASbb1blh;
	vector<double> vmuATLAStata1blh;
	vector<double> vmuATLASZy1blh;
	
	vector<double> hmuATLASyy1blh;
	vector<double> hmuATLASZZ1blh;
	vector<double> hmuATLASWW1blh;
	vector<double> hmuATLASWW2blh;
	vector<double> hmuATLASWW3blh;
	vector<double> hmuATLASbb1blh;
	vector<double> hmuATLAStata1blh;
	vector<double> hmuATLASZy1blh;
	vector<double> AmuATLASyy1blh;
	vector<double> AmuATLAStata1blh;
	vector<double> AmuATLASZy1blh;
	vector<double> HmuATLASyy1blh;
	vector<double> HmuATLASZZ1blh;
	vector<double> HmuATLASWW1blh;
	vector<double> HmuATLASWW2blh;
	vector<double> HmuATLASWW3blh;
	vector<double> HmuATLASbb1blh;
	vector<double> HmuATLAStata1blh;
	vector<double> HmuATLASZy1blh;
	
	double flmin=9999, flmax=0;
	double fhmin=9999, fhmax=0;
	double mh0min=9999, mh0max=0;
	double ma0min=9999, ma0max=0;
	double deltammin=9999, deltammax=0;
	double tBmin=9999, tBmax=0;
	double tGmin=9999, tGmax=0;
	double DSggmin=9999, DSggmax=0;
	double Sggmin=9999, Sggmax=0;
	double DSyymin=9999, DSyymax=0;
	double Syymin=9999, Syymax=0;
	double Cyymin=9999, Cyymax=0;
	double Cggmin=9999, Cggmax=0;
	double CWWmin=9999, CWWmax=0;
	double CZZmin=9999, CZZmax=0;
	double Cbbmin=9999, Cbbmax=0;
	double Ctatamin=9999, Ctatamax=0;
	double muyymin=9999, muyymax=0;
	double muZZmin=9999, muZZmax=0;
	double muWWmin=9999, muWWmax=0;
	double mubbmin=9999, mubbmax=0;
	double mutatamin=9999, mutatamax=0;
	double muZymin=9999, muZymax=0;
	double mhhmin=9999, mhhmax=0;
	
	
	vector<double> ggFhyyATLAS;
	vector<double> VBFhyyATLAS;
	vector<double> VHhyyATLAS;
	vector<double> ttHhyyATLAS;
	vector<double> ggFhbbATLAS;
	vector<double> VBFhbbATLAS;
	vector<double> VHhbbATLAS;
	vector<double> ttHhbbATLAS;
	vector<double> ggFhtataATLAS;
	vector<double> VBFhtataATLAS;
	vector<double> VHhtataATLAS;
	vector<double> ttHhtataATLAS;
	vector<double> ggFhZZATLAS;
	vector<double> VBFhZZATLAS;
	vector<double> VHhZZATLAS;
	vector<double> ttHhZZATLAS;
	vector<double> ggFhWWATLAS;
	vector<double> VBFhWWATLAS;
	vector<double> VHhWWATLAS;
	vector<double> ttHhWWATLAS;
	vector<double> ggFhZyATLAS;
	vector<double> VBFhZyATLAS;
	vector<double> VHhZyATLAS;
	vector<double> ttHhZyATLAS;
	
	vector<double> hggFhyyATLAS;
	vector<double> hVBFhyyATLAS;
	vector<double> hVHhyyATLAS;
	vector<double> httHhyyATLAS;
	vector<double> hggFhbbATLAS;
	vector<double> hVBFhbbATLAS;
	vector<double> hVHhbbATLAS;
	vector<double> httHhbbATLAS;
	vector<double> hggFhtataATLAS;
	vector<double> hVBFhtataATLAS;
	vector<double> hVHhtataATLAS;
	vector<double> httHhtataATLAS;
	vector<double> hggFhZZATLAS;
	vector<double> hVBFhZZATLAS;
	vector<double> hVHhZZATLAS;
	vector<double> httHhZZATLAS;
	vector<double> hggFhWWATLAS;
	vector<double> hVBFhWWATLAS;
	vector<double> hVHhWWATLAS;
	vector<double> httHhWWATLAS;
	vector<double> hggFhZyATLAS;
	vector<double> hVBFhZyATLAS;
	vector<double> hVHhZyATLAS;
	vector<double> httHhZyATLAS;
	
	vector<double> AggFhyyATLAS;
	vector<double> AttHhyyATLAS;
	vector<double> AggFhbbATLAS;
	vector<double> AttHhbbATLAS;
	vector<double> AggFhtataATLAS;
	vector<double> AttHhtataATLAS;
	vector<double> AggFhZyATLAS;
	vector<double> AttHhZyATLAS;
	
	vector<double> HggFhyyATLAS;
	vector<double> HVBFhyyATLAS;
	vector<double> HVHhyyATLAS;
	vector<double> HttHhyyATLAS;
	vector<double> HggFhbbATLAS;
	vector<double> HVBFhbbATLAS;
	vector<double> HVHhbbATLAS;
	vector<double> HttHhbbATLAS;
	vector<double> HggFhtataATLAS;
	vector<double> HVBFhtataATLAS;
	vector<double> HVHhtataATLAS;
	vector<double> HttHhtataATLAS;
	vector<double> HggFhZZATLAS;
	vector<double> HVBFhZZATLAS;
	vector<double> HVHhZZATLAS;
	vector<double> HttHhZZATLAS;
	vector<double> HggFhWWATLAS;
	vector<double> HVBFhWWATLAS;
	vector<double> HVHhWWATLAS;
	vector<double> HttHhWWATLAS;
	vector<double> HggFhZyATLAS;
	vector<double> HVBFhZyATLAS;
	vector<double> HVHhZyATLAS;
	vector<double> HttHhZyATLAS;
	
	
	vector<double> ggFhyyCMS;
	vector<double> VBFhyyCMS;
	vector<double> VHhyyCMS;
	vector<double> ttHhyyCMS;
	vector<double> ggFhbbCMS;
	vector<double> VBFhbbCMS;
	vector<double> VHhbbCMS;
	vector<double> ttHhbbCMS;
	vector<double> ggFhtataCMS;
	vector<double> VBFhtataCMS;
	vector<double> VHhtataCMS;
	vector<double> ttHhtataCMS;
	vector<double> ggFhZZCMS;
	vector<double> VBFhZZCMS;
	vector<double> VHhZZCMS;
	vector<double> ttHhZZCMS;
	vector<double> ggFhWWCMS;
	vector<double> VBFhWWCMS;
	vector<double> VHhWWCMS;
	vector<double> ttHhWWCMS;
	vector<double> ggFhZyCMS;
	vector<double> VBFhZyCMS;
	vector<double> VHhZyCMS;
	vector<double> ttHhZyCMS;
	
	vector<double> Htest;
	
	double CMSyy7=5.1;
	double CMSZZ7=5.1;
	double CMSWW7=4.9;
	double CMSbb7=5.0;
	double CMStata7=4.9;
	double CMSyy8=19.6;
	double CMSZZ8=19.6;
	double CMSWW8=19.5;
	double CMSbb8=12.1;
	double CMStata8=19.4;
	double CMSZy7=5.0;
	double CMSZy8=19.6;
	double ATLASyy7=4.8;
	double ATLASZZ7=4.6;
	double ATLASWW7=4.6;
	double ATLASbb7=4.7;
	double ATLAStata7=4.6;
	double ATLASyy8=20.7;
	double ATLASZZ8=20.7;
	double ATLASWW8=20.7;
	double ATLASbb8=13.0;
	double ATLAStata8=13.0;
	double ATLASZy7=4.6;
	double ATLASZy8=20.7;
	
	int iline = 0;
	if ( datafile.is_open() ) {
		while ( !datafile.eof() ) {
			getline(datafile,line);
			if ( line.find("CMS",0) != string::npos ) {
				//cout << "Entering CMS" << endl;
				start = 4;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFyySMCMS += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFyySMCMS += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMCMS += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMCMS += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHyySMCMS += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFyySMCMS += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFyySMCMS += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMCMS += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMCMS += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHyySMCMS += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFbbSMCMS += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFbbSMCMS += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMCMS += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMCMS += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHbbSMCMS += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFbbSMCMS += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFbbSMCMS += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMCMS += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMCMS += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHbbSMCMS += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFtataSMCMS += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFtataSMCMS += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMCMS += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMCMS += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHtataSMCMS += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFtataSMCMS += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFtataSMCMS += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMCMS += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMCMS += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHtataSMCMS += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZZSMCMS += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZZSMCMS += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMCMS += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMCMS += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZZSMCMS += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZZSMCMS += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZZSMCMS += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMCMS += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMCMS += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZZSMCMS += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFWWSMCMS += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFWWSMCMS += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMCMS += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMCMS += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHWWSMCMS += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFWWSMCMS += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFWWSMCMS += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMCMS += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMCMS += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHWWSMCMS += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZySMCMS += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZySMCMS += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMCMS += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMCMS += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZySMCMS += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZySMCMS += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZySMCMS += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMCMS += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMCMS += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZySMCMS += CMSZy8*atof(line.substr(start,end-start).c_str());
			}
			if ( line.find("ATLAS",0) != string::npos ) {
				//cout << "Entering ATLAS" << endl;
				start = 6;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFyySMATLAS += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFyySMATLAS += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMATLAS += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMATLAS += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHyySMATLAS += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFyySMATLAS += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFyySMATLAS += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMATLAS += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHyySMATLAS += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHyySMATLAS += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFbbSMATLAS += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFbbSMATLAS += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMATLAS += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMATLAS += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHbbSMATLAS += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFbbSMATLAS += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFbbSMATLAS += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMATLAS += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHbbSMATLAS += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHbbSMATLAS += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFtataSMATLAS += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFtataSMATLAS += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMATLAS += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMATLAS += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHtataSMATLAS += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFtataSMATLAS += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFtataSMATLAS += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMATLAS += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHtataSMATLAS += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHtataSMATLAS += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZZSMATLAS += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZZSMATLAS += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMATLAS += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMATLAS += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZZSMATLAS += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZZSMATLAS += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZZSMATLAS += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMATLAS += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZZSMATLAS += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZZSMATLAS += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFWWSMATLAS += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFWWSMATLAS += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMATLAS += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMATLAS += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHWWSMATLAS += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFWWSMATLAS += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFWWSMATLAS += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMATLAS += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHWWSMATLAS += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHWWSMATLAS += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZySMATLAS += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZySMATLAS += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMATLAS += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMATLAS += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZySMATLAS += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFZySMATLAS += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFZySMATLAS += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMATLAS += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHZySMATLAS += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHZySMATLAS += ATLASZy8*atof(line.substr(start,end-start).c_str());
			}
			if ( line.find("Params",0) != string::npos ) {
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				fl.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				fh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mh0.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ma0.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mhh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				tB.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				tG.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				t12.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				t13.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				sigfact.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				kG.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				kY.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				kL.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Cyy.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Cgg.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CZZ.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CWW.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Cbb.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Ctata.push_back(pow(atof(line.substr(start,end-start).c_str()),2.0));
				
				
				
				
				//cout << fl.size() << " ";
				//cout << fh.size() << " ";
				//cout << mh0.size() << " ";
				//cout << ma0.size() << " ";
				//cout << mhh.size() << " ";
				//cout << tB.size() << " ";
				//cout << tG.size() << " ";
				//cout << t12.size() << " ";
				//cout << t13.size() << " ";
				//cout << sigfact.size() << " ";
				//cout << kG.size() << " ";
				//cout << kY.size() << " ";
				//cout << kL.size() << endl;
			}
			if ( line.find("lighth",0) != string::npos ) {
				//cout << "Entering lighth" << endl;
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				hggFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.at(ggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				hggFhyyATLAS.at(hggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				hVBFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyATLAS.at(VBFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				hVBFhyyATLAS.at(hVBFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				hVHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				hVHhyyATLAS.at(hVHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				hVHhyyATLAS.at(hVHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				hVHhyyATLAS.at(hVHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				httHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.at(ttHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				httHhyyATLAS.at(httHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				hggFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.at(ggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				hggFhbbATLAS.at(hggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				hVBFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbATLAS.at(VBFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				hVBFhbbATLAS.at(hVBFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				hVHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				hVHhbbATLAS.at(hVHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				hVHhbbATLAS.at(hVHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				hVHhbbATLAS.at(hVHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				httHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.at(ttHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				httHhbbATLAS.at(httHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				hggFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.at(ggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				hggFhtataATLAS.at(hggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				hVBFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataATLAS.at(VBFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				hVBFhtataATLAS.at(hVBFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				hVHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				hVHhtataATLAS.at(hVHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				hVHhtataATLAS.at(hVHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				hVHhtataATLAS.at(hVHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				httHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.at(ttHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				httHhtataATLAS.at(httHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				hggFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZATLAS.at(ggFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				hggFhZZATLAS.at(hggFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				hVBFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZATLAS.at(VBFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				hVBFhZZATLAS.at(hVBFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				hVHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				hVHhZZATLAS.at(hVHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				hVHhZZATLAS.at(hVHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				hVHhZZATLAS.at(hVHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				httHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZATLAS.at(ttHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				httHhZZATLAS.at(httHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				hggFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWATLAS.at(ggFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				hggFhWWATLAS.at(hggFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				hVBFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWATLAS.at(VBFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				hVBFhWWATLAS.at(hVBFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				hVHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				hVHhWWATLAS.at(hVHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				hVHhWWATLAS.at(hVHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				hVHhWWATLAS.at(hVHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				httHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWATLAS.at(ttHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				httHhWWATLAS.at(httHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				hggFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.at(ggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				hggFhZyATLAS.at(hggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				hVBFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyATLAS.at(VBFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				hVBFhZyATLAS.at(hVBFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				hVHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				hVHhZyATLAS.at(hVHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				hVHhZyATLAS.at(hVHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				hVHhZyATLAS.at(hVHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				httHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.at(ttHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				httHhZyATLAS.at(httHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.push_back(CMSyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.at(ggFhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyCMS.push_back(CMSyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyCMS.at(VBFhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.push_back(CMSyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.push_back(CMSyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.at(ttHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.push_back(CMSbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.at(ggFhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbCMS.push_back(CMSbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbCMS.at(VBFhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.push_back(CMSbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.push_back(CMSbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.at(ttHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.push_back(CMStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.at(ggFhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataCMS.push_back(CMStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataCMS.at(VBFhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.push_back(CMStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.push_back(CMStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.at(ttHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZCMS.push_back(CMSZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZCMS.at(ggFhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZCMS.push_back(CMSZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZCMS.at(VBFhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.push_back(CMSZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZCMS.push_back(CMSZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZCMS.at(ttHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWCMS.push_back(CMSWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWCMS.at(ggFhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWCMS.push_back(CMSWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWCMS.at(VBFhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.push_back(CMSWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWCMS.push_back(CMSWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWCMS.at(ttHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.push_back(CMSZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.at(ggFhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyCMS.push_back(CMSZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyCMS.at(VBFhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.push_back(CMSZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.push_back(CMSZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.at(ttHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				
				//cout << ggFhyyATLAS.size() << " ";
				//cout << VBFhyyATLAS.size() << " ";
				//cout << VHhyyATLAS.size() << " ";
				//cout << ttHhyyATLAS.size() << " ";
				//cout << ggFhbbATLAS.size() << " ";
				//cout << VBFhbbATLAS.size() << " ";
				//cout << VHhbbATLAS.size() << " ";
				//cout << ttHhbbATLAS.size() << " ";
				//cout << ggFhtataATLAS.size() << " ";
				//cout << VBFhtataATLAS.size() << " ";
				//cout << VHhtataATLAS.size() << " ";
				//cout << ttHhtataATLAS.size() << " ";
				//cout << ggFhZZATLAS.size() << " ";
				//cout << VBFhZZATLAS.size() << " ";
				//cout << VHhZZATLAS.size() << " ";
				//cout << ttHhZZATLAS.size() << " ";
				//cout << ggFhWWATLAS.size() << " ";
				//cout << VBFhWWATLAS.size() << " ";
				//cout << VHhWWATLAS.size() << " ";
				//cout << ttHhWWATLAS.size() << endl;
				//cout << ggFhyyCMS.size() << " ";
				//cout << VBFhyyCMS.size() << " ";
				//cout << VHhyyCMS.size() << " ";
				//cout << ttHhyyCMS.size() << " ";
				//cout << ggFhbbCMS.size() << " ";
				//cout << VBFhbbCMS.size() << " ";
				//cout << VHhbbCMS.size() << " ";
				//cout << ttHhbbCMS.size() << " ";
				//cout << ggFhtataCMS.size() << " ";
				//cout << VBFhtataCMS.size() << " ";
				//cout << VHhtataCMS.size() << " ";
				//cout << ttHhtataCMS.size() << " ";
				//cout << ggFhZZCMS.size() << " ";
				//cout << VBFhZZCMS.size() << " ";
				//cout << VHhZZCMS.size() << " ";
				//cout << ttHhZZCMS.size() << " ";
				//cout << ggFhWWCMS.size() << " ";
				//cout << VBFhWWCMS.size() << " ";
				//cout << VHhWWCMS.size() << " ";
				//cout << ttHhWWCMS.size() << endl;
			}
			
			if ( line.find("pseudo",0) != string::npos && ma0.at(ma0.size()-1) < 900 ) {
				//cout << "Entering Pseudo" << endl;
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.at(ggFhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				AggFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.at(ggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				AggFhyyATLAS.at(AggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.at(ttHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				AttHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.at(ttHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				AttHhyyATLAS.at(AttHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.at(ggFhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				AggFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.at(ggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				AggFhbbATLAS.at(AggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.at(ttHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				AttHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.at(ttHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				AttHhbbATLAS.at(AttHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.at(ggFhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				AggFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.at(ggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				AggFhtataATLAS.at(AggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.at(ttHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				AttHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.at(ttHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				AttHhtataATLAS.at(AttHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.at(ggFhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				AggFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.at(ggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				AggFhZyATLAS.at(AggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.at(ttHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				AttHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.at(ttHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				AttHhZyATLAS.at(AttHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.at(ggFhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.at(ggFhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.at(ttHhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.at(ttHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.at(ggFhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.at(ggFhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.at(ttHhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.at(ttHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.at(ggFhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.at(ggFhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.at(ttHhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.at(ttHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.at(ggFhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.at(ggFhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.at(ttHhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.at(ttHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				
			}
			else if ( line.find("pseudo",0) != string::npos && ma0.at(ma0.size()-1) > 900 ) {
				AggFhyyATLAS.push_back(0);
				AttHhyyATLAS.push_back(0);
				AggFhbbATLAS.push_back(0);
				AttHhbbATLAS.push_back(0);
				AggFhtataATLAS.push_back(0);
				AttHhtataATLAS.push_back(0);
				AggFhZyATLAS.push_back(0);
				AttHhZyATLAS.push_back(0);
			}
			
			if ( line.find("heavyh",0) != string::npos && mhh.at(mhh.size()-1) < 900 ) {
				//cout << "Entering heavyh" << endl;
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.at(ggFhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				HggFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyATLAS.at(ggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				HggFhyyATLAS.at(HggFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyATLAS.at(VBFhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				HVBFhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyATLAS.at(VBFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				HVBFhyyATLAS.at(HVBFhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				HVHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				HVHhyyATLAS.at(HVHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				HVHhyyATLAS.at(HVHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyATLAS.at(VHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				HVHhyyATLAS.at(HVHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.at(ttHhyyATLAS.size()-1) += ATLASyy7*atof(line.substr(start,end-start).c_str());
				HttHhyyATLAS.push_back(ATLASyy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyATLAS.at(ttHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				HttHhyyATLAS.at(HttHhyyATLAS.size()-1) += ATLASyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.at(ggFhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				HggFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbATLAS.at(ggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				HggFhbbATLAS.at(HggFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbATLAS.at(VBFhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				HVBFhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbATLAS.at(VBFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				HVBFhbbATLAS.at(HVBFhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				HVHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				HVHhbbATLAS.at(HVHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				HVHhbbATLAS.at(HVHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbATLAS.at(VHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				HVHhbbATLAS.at(HVHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.at(ttHhbbATLAS.size()-1) += ATLASbb7*atof(line.substr(start,end-start).c_str());
				HttHhbbATLAS.push_back(ATLASbb7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbATLAS.at(ttHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				HttHhbbATLAS.at(HttHhbbATLAS.size()-1) += ATLASbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.at(ggFhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				HggFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataATLAS.at(ggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				HggFhtataATLAS.at(HggFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataATLAS.at(VBFhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				HVBFhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataATLAS.at(VBFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				HVBFhtataATLAS.at(HVBFhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				HVHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				HVHhtataATLAS.at(HVHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				HVHhtataATLAS.at(HVHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataATLAS.at(VHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				HVHhtataATLAS.at(HVHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.at(ttHhtataATLAS.size()-1) += ATLAStata7*atof(line.substr(start,end-start).c_str());
				HttHhtataATLAS.push_back(ATLAStata7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataATLAS.at(ttHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				HttHhtataATLAS.at(HttHhtataATLAS.size()-1) += ATLAStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZATLAS.at(ggFhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				HggFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZATLAS.at(ggFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				HggFhZZATLAS.at(HggFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZATLAS.at(VBFhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				HVBFhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZATLAS.at(VBFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				HVBFhZZATLAS.at(HVBFhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				HVHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				HVHhZZATLAS.at(HVHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				HVHhZZATLAS.at(HVHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZATLAS.at(VHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				HVHhZZATLAS.at(HVHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZATLAS.at(ttHhZZATLAS.size()-1) += ATLASZZ7*atof(line.substr(start,end-start).c_str());
				HttHhZZATLAS.push_back(ATLASZZ7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZATLAS.at(ttHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				HttHhZZATLAS.at(HttHhZZATLAS.size()-1) += ATLASZZ8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWATLAS.at(ggFhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				HggFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWATLAS.at(ggFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				HggFhWWATLAS.at(HggFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				Htest.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWATLAS.at(VBFhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				HVBFhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWATLAS.at(VBFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				HVBFhWWATLAS.at(HVBFhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				HVHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				HVHhWWATLAS.at(HVHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				HVHhWWATLAS.at(HVHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWATLAS.at(VHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				HVHhWWATLAS.at(HVHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWATLAS.at(ttHhWWATLAS.size()-1) += ATLASWW7*atof(line.substr(start,end-start).c_str());
				HttHhWWATLAS.push_back(ATLASWW7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWATLAS.at(ttHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				HttHhWWATLAS.at(HttHhWWATLAS.size()-1) += ATLASWW8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.at(ggFhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				HggFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyATLAS.at(ggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				HggFhZyATLAS.at(HggFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyATLAS.at(VBFhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				HVBFhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyATLAS.at(VBFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				HVBFhZyATLAS.at(HVBFhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				HVHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				HVHhZyATLAS.at(HVHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				HVHhZyATLAS.at(HVHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyATLAS.at(VHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				HVHhZyATLAS.at(HVHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.at(ttHhZyATLAS.size()-1) += ATLASZy7*atof(line.substr(start,end-start).c_str());
				HttHhZyATLAS.push_back(ATLASZy7*atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyATLAS.at(ttHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				HttHhZyATLAS.at(HttHhZyATLAS.size()-1) += ATLASZy8*atof(line.substr(start,end-start).c_str());
				
				start = 7;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.at(ggFhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhyyCMS.at(ggFhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyCMS.at(VBFhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhyyCMS.at(VBFhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhyyCMS.at(VHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.at(ttHhyyCMS.size()-1) += CMSyy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhyyCMS.at(ttHhyyCMS.size()-1) += CMSyy8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.at(ggFhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhbbCMS.at(ggFhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbCMS.at(VBFhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhbbCMS.at(VBFhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhbbCMS.at(VHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.at(ttHhbbCMS.size()-1) += CMSbb7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhbbCMS.at(ttHhbbCMS.size()-1) += CMSbb8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.at(ggFhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhtataCMS.at(ggFhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataCMS.at(VBFhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhtataCMS.at(VBFhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhtataCMS.at(VHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.at(ttHhtataCMS.size()-1) += CMStata7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhtataCMS.at(ttHhtataCMS.size()-1) += CMStata8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZCMS.at(ggFhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZZCMS.at(ggFhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZCMS.at(VBFhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZZCMS.at(VBFhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZZCMS.at(VHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZCMS.at(ttHhZZCMS.size()-1) += CMSZZ7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZZCMS.at(ttHhZZCMS.size()-1) += CMSZZ8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWCMS.at(ggFhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhWWCMS.at(ggFhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWCMS.at(VBFhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhWWCMS.at(VBFhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhWWCMS.at(VHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWCMS.at(ttHhWWCMS.size()-1) += CMSWW7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhWWCMS.at(ttHhWWCMS.size()-1) += CMSWW8*atof(line.substr(start,end-start).c_str());
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.at(ggFhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ggFhZyCMS.at(ggFhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyCMS.at(VBFhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VBFhZyCMS.at(VBFhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				VHhZyCMS.at(VHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.at(ttHhZyCMS.size()-1) += CMSZy7*atof(line.substr(start,end-start).c_str());
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ttHhZyCMS.at(ttHhZyCMS.size()-1) += CMSZy8*atof(line.substr(start,end-start).c_str());
			}
			else if (  line.find("heavyh",0) != string::npos && mhh.at(mhh.size()-1) >= 900 ) {
				Htest.push_back(0);
				HggFhyyATLAS.push_back(0);
				HVBFhyyATLAS.push_back(0);
				HVHhyyATLAS.push_back(0);
				HttHhyyATLAS.push_back(0);
				HggFhbbATLAS.push_back(0);
				HVBFhbbATLAS.push_back(0);
				HVHhbbATLAS.push_back(0);
				HttHhbbATLAS.push_back(0);
				HggFhtataATLAS.push_back(0);
				HVBFhtataATLAS.push_back(0);
				HVHhtataATLAS.push_back(0);
				HttHhtataATLAS.push_back(0);
				HggFhZZATLAS.push_back(0);
				HVBFhZZATLAS.push_back(0);
				HVHhZZATLAS.push_back(0);
				HttHhZZATLAS.push_back(0);
				HggFhWWATLAS.push_back(0);
				HVBFhWWATLAS.push_back(0);
				HVHhWWATLAS.push_back(0);
				HttHhWWATLAS.push_back(0);
				HggFhZyATLAS.push_back(0);
				HVBFhZyATLAS.push_back(0);
				HVHhZyATLAS.push_back(0);
				HttHhZyATLAS.push_back(0);
			}
			iline++;
			if ( iline%1000 == 0 ) {
				cout << iline << endl;
			}
		}
	}
	
	//cout << "Done reading in data" << endl;
	
	
	cout << "Vector size ggFhyyATLAS: " << ggFhyyATLAS.size() << endl;
	cout << "Vector size VBFhyyATLAS: " << VBFhyyATLAS.size() << endl;
	cout << "Vector size VHhyyATLAS: " << VHhyyATLAS.size() << endl;
	cout << "Vector size ttHhyyATLAS: " << ttHhyyATLAS.size() << endl;
	cout << "Vector size ggFhbbATLAS: " << ggFhbbATLAS.size() << endl;
	cout << "Vector size VBFhbbATLAS: " << VBFhbbATLAS.size() << endl;
	cout << "Vector size VHhbbATLAS: " << VHhbbATLAS.size() << endl;
	cout << "Vector size ttHhbbATLAS: " << ttHhbbATLAS.size() << endl;
	cout << "Vector size ggFhtataATLAS: " << ggFhtataATLAS.size() << endl;
	cout << "Vector size VBFhtataATLAS: " << VBFhtataATLAS.size() << endl;
	cout << "Vector size VHhtataATLAS: " << VHhtataATLAS.size() << endl;
	cout << "Vector size ttHhtataATLAS: " << ttHhtataATLAS.size() << endl;
	cout << "Vector size ggFhZZATLAS: " << ggFhZZATLAS.size() << endl;
	cout << "Vector size VBFhZZATLAS: " << VBFhZZATLAS.size() << endl;
	cout << "Vector size VHhZZATLAS: " << VHhZZATLAS.size() << endl;
	cout << "Vector size ttHhZZATLAS: " << ttHhZZATLAS.size() << endl;
	cout << "Vector size ggFhWWATLAS: " << ggFhWWATLAS.size() << endl;
	cout << "Vector size VBFhWWATLAS: " << VBFhWWATLAS.size() << endl;
	cout << "Vector size VHhWWATLAS: " << VHhWWATLAS.size() << endl;
	cout << "Vector size ttHhWWATLAS: " << ttHhWWATLAS.size() << endl;
	cout << "Vector size ggFhZyATLAS: " << ggFhZyATLAS.size() << endl;
	cout << "Vector size VBFhZyATLAS: " << VBFhZyATLAS.size() << endl;
	cout << "Vector size VHhZyATLAS: " << VHhZyATLAS.size() << endl;
	cout << "Vector size ttHhZyATLAS: " << ttHhZyATLAS.size() << endl;
	
	cout << "Vector size hggFhyyATLAS: " << hggFhyyATLAS.size() << endl;
	cout << "Vector size hVBFhyyATLAS: " << hVBFhyyATLAS.size() << endl;
	cout << "Vector size hVHhyyATLAS: " << hVHhyyATLAS.size() << endl;
	cout << "Vector size httHhyyATLAS: " << httHhyyATLAS.size() << endl;
	cout << "Vector size hggFhbbATLAS: " << hggFhbbATLAS.size() << endl;
	cout << "Vector size hVBFhbbATLAS: " << hVBFhbbATLAS.size() << endl;
	cout << "Vector size hVHhbbATLAS: " << hVHhbbATLAS.size() << endl;
	cout << "Vector size httHhbbATLAS: " << httHhbbATLAS.size() << endl;
	cout << "Vector size hggFhtataATLAS: " << hggFhtataATLAS.size() << endl;
	cout << "Vector size hVBFhtataATLAS: " << hVBFhtataATLAS.size() << endl;
	cout << "Vector size hVHhtataATLAS: " << hVHhtataATLAS.size() << endl;
	cout << "Vector size httHhtataATLAS: " << httHhtataATLAS.size() << endl;
	cout << "Vector size hggFhZZATLAS: " << hggFhZZATLAS.size() << endl;
	cout << "Vector size hVBFhZZATLAS: " << hVBFhZZATLAS.size() << endl;
	cout << "Vector size hVHhZZATLAS: " << hVHhZZATLAS.size() << endl;
	cout << "Vector size httHhZZATLAS: " << httHhZZATLAS.size() << endl;
	cout << "Vector size hggFhWWATLAS: " << hggFhWWATLAS.size() << endl;
	cout << "Vector size hVBFhWWATLAS: " << hVBFhWWATLAS.size() << endl;
	cout << "Vector size hVHhWWATLAS: " << hVHhWWATLAS.size() << endl;
	cout << "Vector size httHhWWATLAS: " << httHhWWATLAS.size() << endl;
	cout << "Vector size hggFhZyATLAS: " << hggFhZyATLAS.size() << endl;
	cout << "Vector size hVBFhZyATLAS: " << hVBFhZyATLAS.size() << endl;
	cout << "Vector size hVHhZyATLAS: " << hVHhZyATLAS.size() << endl;
	cout << "Vector size httHhZyATLAS: " << httHhZyATLAS.size() << endl;
	
	cout << "Vector size AggFhyyATLAS: " << AggFhyyATLAS.size() << endl;
	cout << "Vector size AttHhyyATLAS: " << AttHhyyATLAS.size() << endl;
	cout << "Vector size AggFhbbATLAS: " << AggFhbbATLAS.size() << endl;
	cout << "Vector size AttHhbbATLAS: " << AttHhbbATLAS.size() << endl;
	cout << "Vector size AggFhtataATLAS: " << AggFhtataATLAS.size() << endl;
	cout << "Vector size AttHhtataATLAS: " << AttHhtataATLAS.size() << endl;
	cout << "Vector size AggFhZyATLAS: " << AggFhZyATLAS.size() << endl;
	cout << "Vector size AttHhZyATLAS: " << AttHhZyATLAS.size() << endl;
	
	cout << "Vector size HggFhyyATLAS: " << HggFhyyATLAS.size() << endl;
	cout << "Vector size HVBFhyyATLAS: " << HVBFhyyATLAS.size() << endl;
	cout << "Vector size HVHhyyATLAS: " << HVHhyyATLAS.size() << endl;
	cout << "Vector size HttHhyyATLAS: " << HttHhyyATLAS.size() << endl;
	cout << "Vector size HggFhbbATLAS: " << HggFhbbATLAS.size() << endl;
	cout << "Vector size HVBFhbbATLAS: " << HVBFhbbATLAS.size() << endl;
	cout << "Vector size HVHhbbATLAS: " << HVHhbbATLAS.size() << endl;
	cout << "Vector size HttHhbbATLAS: " << HttHhbbATLAS.size() << endl;
	cout << "Vector size HggFhtataATLAS: " << HggFhtataATLAS.size() << endl;
	cout << "Vector size HVBFhtataATLAS: " << HVBFhtataATLAS.size() << endl;
	cout << "Vector size HVHhtataATLAS: " << HVHhtataATLAS.size() << endl;
	cout << "Vector size HttHhtataATLAS: " << HttHhtataATLAS.size() << endl;
	cout << "Vector size HggFhZZATLAS: " << HggFhZZATLAS.size() << endl;
	cout << "Vector size HVBFhZZATLAS: " << HVBFhZZATLAS.size() << endl;
	cout << "Vector size HVHhZZATLAS: " << HVHhZZATLAS.size() << endl;
	cout << "Vector size HttHhZZATLAS: " << HttHhZZATLAS.size() << endl;
	cout << "Vector size HggFhWWATLAS: " << HggFhWWATLAS.size() << endl;
	cout << "Vector size HVBFhWWATLAS: " << HVBFhWWATLAS.size() << endl;
	cout << "Vector size HVHhWWATLAS: " << HVHhWWATLAS.size() << endl;
	cout << "Vector size HttHhWWATLAS: " << HttHhWWATLAS.size() << endl;
	cout << "Vector size HggFhZyATLAS: " << HggFhZyATLAS.size() << endl;
	cout << "Vector size HVBFhZyATLAS: " << HVBFhZyATLAS.size() << endl;
	cout << "Vector size HVHhZyATLAS: " << HVHhZyATLAS.size() << endl;
	cout << "Vector size HttHhZyATLAS: " << HttHhZyATLAS.size() << endl;
	
	cout << "Vector size ggFhyyCMS: " << ggFhyyCMS.size() << endl;
	cout << "Vector size VBFhyyCMS: " << VBFhyyCMS.size() << endl;
	cout << "Vector size VHhyyCMS: " << VHhyyCMS.size() << endl;
	cout << "Vector size ttHhyyCMS: " << ttHhyyCMS.size() << endl;
	cout << "Vector size ggFhbbCMS: " << ggFhbbCMS.size() << endl;
	cout << "Vector size VBFhbbCMS: " << VBFhbbCMS.size() << endl;
	cout << "Vector size VHhbbCMS: " << VHhbbCMS.size() << endl;
	cout << "Vector size ttHhbbCMS: " << ttHhbbCMS.size() << endl;
	cout << "Vector size ggFhtataCMS: " << ggFhtataCMS.size() << endl;
	cout << "Vector size VBFhtataCMS: " << VBFhtataCMS.size() << endl;
	cout << "Vector size VHhtataCMS: " << VHhtataCMS.size() << endl;
	cout << "Vector size ttHhtataCMS: " << ttHhtataCMS.size() << endl;
	cout << "Vector size ggFhZZCMS: " << ggFhZZCMS.size() << endl;
	cout << "Vector size VBFhZZCMS: " << VBFhZZCMS.size() << endl;
	cout << "Vector size VHhZZCMS: " << VHhZZCMS.size() << endl;
	cout << "Vector size ttHhZZCMS: " << ttHhZZCMS.size() << endl;
	cout << "Vector size ggFhWWCMS: " << ggFhWWCMS.size() << endl;
	cout << "Vector size VBFhWWCMS: " << VBFhWWCMS.size() << endl;
	cout << "Vector size VHhWWCMS: " << VHhWWCMS.size() << endl;
	cout << "Vector size ttHhWWCMS: " << ttHhWWCMS.size() << endl;
	cout << "Vector size ggFhZyCMS: " << ggFhZyCMS.size() << endl;
	cout << "Vector size VBFhZyCMS: " << VBFhZyCMS.size() << endl;
	cout << "Vector size VHhZyCMS: " << VHhZyCMS.size() << endl;
	cout << "Vector size ttHhZyCMS: " << ttHhZyCMS.size() << endl;
	
	for ( int i=0; i<fl.size(); i++) {
		deltam.push_back(ma0.at(i)-mh0.at(i));
		muCMSyy1blh=(ggFhyyCMS.at(i)+VBFhyyCMS.at(i)+VHhyyCMS.at(i)+ttHhyyCMS.at(i))/(ggFyySMCMS+VBFyySMCMS+VHyySMCMS+ttHyySMCMS);
		muCMSZZ1blh=(ggFhZZCMS.at(i)+VBFhZZCMS.at(i)+VHhZZCMS.at(i)+ttHhZZCMS.at(i))/(ggFZZSMCMS+VBFZZSMCMS+VHZZSMCMS+ttHZZSMCMS);
		muCMSZZ2blh=(VBFhZZCMS.at(i))/(VBFZZSMCMS);
		muCMSWW1blh=(ggFhWWCMS.at(i))/(ggFWWSMCMS);
		muCMSbb1blh=(VHhbbCMS.at(i))/(VHbbSMCMS);
		muCMStata1blh=(ggFhtataCMS.at(i)+VBFhtataCMS.at(i)+VHhtataCMS.at(i)+ttHhtataCMS.at(i))/(ggFtataSMCMS+VBFtataSMCMS+VHtataSMCMS+ttHtataSMCMS);
		muCMStata2blh=(VBFhtataCMS.at(i))/(VBFtataSMCMS);
		muCMSZy1blh=(ggFhZyCMS.at(i)+VBFhZyCMS.at(i)+VHhZyCMS.at(i)+ttHhZyCMS.at(i))/(ggFZySMCMS+VBFZySMCMS+VHZySMCMS+ttHZySMCMS);
		
		muATLASyy1blh=(ggFhyyATLAS.at(i)+VBFhyyATLAS.at(i)+VHhyyATLAS.at(i)+ttHhyyATLAS.at(i))/(ggFyySMATLAS+VBFyySMATLAS+VHyySMATLAS+ttHyySMATLAS);
		muATLASZZ1blh=(ggFhZZATLAS.at(i)+VBFhZZATLAS.at(i)+VHhZZATLAS.at(i)+ttHhZZATLAS.at(i))/(ggFZZSMATLAS+VBFZZSMATLAS+VHZZSMATLAS+ttHZZSMATLAS);
		muATLASWW1blh=(ggFhWWATLAS.at(i)+VBFhWWATLAS.at(i))/(ggFWWSMATLAS+VBFWWSMATLAS);
		muATLASWW2blh=(VBFhWWATLAS.at(i))/(VBFWWSMATLAS);
		muATLASWW3blh=(ggFhWWATLAS.at(i))/(ggFWWSMATLAS);
		muATLASbb1blh=(VHhbbATLAS.at(i))/(VHbbSMATLAS);
		muATLAStata1blh=(ggFhtataATLAS.at(i)+VBFhtataATLAS.at(i)+VHhtataATLAS.at(i)+ttHhtataATLAS.at(i))/(ggFtataSMATLAS+VBFtataSMATLAS+VHtataSMATLAS+ttHtataSMATLAS);
		muATLASZy1blh=(ggFhZyATLAS.at(i)+VBFhZyATLAS.at(i)+VHhZyATLAS.at(i)+ttHhZyATLAS.at(i))/(ggFZySMATLAS+VBFZySMATLAS+VHZySMATLAS+ttHZySMATLAS);
		
		hmuATLASyy1blh.push_back((hggFhyyATLAS.at(i)+hVBFhyyATLAS.at(i)+hVHhyyATLAS.at(i)+httHhyyATLAS.at(i))/(ggFyySMATLAS+VBFyySMATLAS+VHyySMATLAS+ttHyySMATLAS));
		hmuATLASZZ1blh.push_back((hggFhZZATLAS.at(i)+hVBFhZZATLAS.at(i)+hVHhZZATLAS.at(i)+httHhZZATLAS.at(i))/(ggFZZSMATLAS+VBFZZSMATLAS+VHZZSMATLAS+ttHZZSMATLAS));
		hmuATLASWW1blh.push_back((hggFhWWATLAS.at(i)+hVBFhWWATLAS.at(i))/(ggFWWSMATLAS+VBFWWSMATLAS));
		hmuATLASWW2blh.push_back((hVBFhWWATLAS.at(i))/(VBFWWSMATLAS));
		hmuATLASWW3blh.push_back((hggFhWWATLAS.at(i))/(ggFWWSMATLAS));
		hmuATLASbb1blh.push_back((hVHhbbATLAS.at(i))/(VHbbSMATLAS));
		hmuATLAStata1blh.push_back((hggFhtataATLAS.at(i)+hVBFhtataATLAS.at(i)+hVHhtataATLAS.at(i)+httHhtataATLAS.at(i))/(ggFtataSMATLAS+VBFtataSMATLAS+VHtataSMATLAS+ttHtataSMATLAS));
		hmuATLASZy1blh.push_back((hggFhZyATLAS.at(i)+hVBFhZyATLAS.at(i)+hVHhZyATLAS.at(i)+httHhZyATLAS.at(i))/(ggFZySMATLAS+VBFZySMATLAS+VHZySMATLAS+ttHZySMATLAS));
		
		AmuATLASyy1blh.push_back((AggFhyyATLAS.at(i)+AttHhyyATLAS.at(i))/(ggFyySMATLAS+VBFyySMATLAS+VHyySMATLAS+ttHyySMATLAS));
		AmuATLAStata1blh.push_back((AggFhtataATLAS.at(i)+AttHhtataATLAS.at(i))/(ggFtataSMATLAS+VBFtataSMATLAS+VHtataSMATLAS+ttHtataSMATLAS));
		AmuATLASZy1blh.push_back((AggFhZyATLAS.at(i)+AttHhZyATLAS.at(i))/(ggFZySMATLAS+VBFZySMATLAS+VHZySMATLAS+ttHZySMATLAS));
		
		HmuATLASyy1blh.push_back((HggFhyyATLAS.at(i)+HVBFhyyATLAS.at(i)+HVHhyyATLAS.at(i)+HttHhyyATLAS.at(i))/(ggFyySMATLAS+VBFyySMATLAS+VHyySMATLAS+ttHyySMATLAS));
		HmuATLASZZ1blh.push_back((HggFhZZATLAS.at(i)+HVBFhZZATLAS.at(i)+HVHhZZATLAS.at(i)+HttHhZZATLAS.at(i))/(ggFZZSMATLAS+VBFZZSMATLAS+VHZZSMATLAS+ttHZZSMATLAS));
		HmuATLASWW1blh.push_back((HggFhWWATLAS.at(i)+HVBFhWWATLAS.at(i))/(ggFWWSMATLAS+VBFWWSMATLAS));
		HmuATLASWW2blh.push_back((HVBFhWWATLAS.at(i))/(VBFWWSMATLAS));
		HmuATLASWW3blh.push_back((HggFhWWATLAS.at(i))/(ggFWWSMATLAS));
		HmuATLASbb1blh.push_back((HVHhbbATLAS.at(i))/(VHbbSMATLAS));
		HmuATLAStata1blh.push_back((HggFhtataATLAS.at(i)+HVBFhtataATLAS.at(i)+HVHhtataATLAS.at(i)+HttHhtataATLAS.at(i))/(ggFtataSMATLAS+VBFtataSMATLAS+VHtataSMATLAS+ttHtataSMATLAS));
		HmuATLASZy1blh.push_back((HggFhZyATLAS.at(i)+HVBFhZyATLAS.at(i)+HVHhZyATLAS.at(i)+HttHhZyATLAS.at(i))/(ggFZySMATLAS+VBFZySMATLAS+VHZySMATLAS+ttHZySMATLAS));
		
		
		vmuCMSyy1blh.push_back((ggFhyyCMS.at(i)+VBFhyyCMS.at(i)+VHhyyCMS.at(i)+ttHhyyCMS.at(i))/(ggFyySMCMS+VBFyySMCMS+VHyySMCMS+ttHyySMCMS));
		vmuCMSZZ1blh.push_back((ggFhZZCMS.at(i)+VBFhZZCMS.at(i)+VHhZZCMS.at(i)+ttHhZZCMS.at(i))/(ggFZZSMCMS+VBFZZSMCMS+VHZZSMCMS+ttHZZSMCMS));
		vmuCMSZZ2blh.push_back((VBFhZZCMS.at(i))/(VBFZZSMCMS));
		vmuCMSWW1blh.push_back((ggFhWWCMS.at(i))/(ggFWWSMCMS));
		vmuCMSbb1blh.push_back((VHhbbCMS.at(i))/(VHbbSMCMS));
		vmuCMStata1blh.push_back((ggFhtataCMS.at(i)+VBFhtataCMS.at(i)+VHhtataCMS.at(i)+ttHhtataCMS.at(i))/(ggFtataSMCMS+VBFtataSMCMS+VHtataSMCMS+ttHtataSMCMS));
		vmuCMStata2blh.push_back((VBFhtataCMS.at(i))/(VBFtataSMCMS));
		vmuCMSZy1blh.push_back((ggFhZyCMS.at(i)+VBFhZyCMS.at(i)+VHhZyCMS.at(i)+ttHhZyCMS.at(i))/(ggFZySMCMS+VBFZySMCMS+VHZySMCMS+ttHZySMCMS));
		
		vmuATLASyy1blh.push_back((ggFhyyATLAS.at(i)+VBFhyyATLAS.at(i)+VHhyyATLAS.at(i)+ttHhyyATLAS.at(i))/(ggFyySMATLAS+VBFyySMATLAS+VHyySMATLAS+ttHyySMATLAS));
		vmuATLASZZ1blh.push_back((ggFhZZATLAS.at(i)+VBFhZZATLAS.at(i)+VHhZZATLAS.at(i)+ttHhZZATLAS.at(i))/(ggFZZSMATLAS+VBFZZSMATLAS+VHZZSMATLAS+ttHZZSMATLAS));
		vmuATLASWW1blh.push_back((ggFhWWATLAS.at(i)+VBFhWWATLAS.at(i))/(ggFWWSMATLAS+VBFWWSMATLAS));
		vmuATLASWW2blh.push_back((VBFhWWATLAS.at(i))/(VBFWWSMATLAS));
		vmuATLASWW3blh.push_back((ggFhWWATLAS.at(i))/(ggFWWSMATLAS));
		vmuATLASbb1blh.push_back((VHhbbATLAS.at(i))/(VHbbSMATLAS));
		vmuATLAStata1blh.push_back((ggFhtataATLAS.at(i)+VBFhtataATLAS.at(i)+VHhtataATLAS.at(i)+ttHhtataATLAS.at(i))/(ggFtataSMATLAS+VBFtataSMATLAS+VHtataSMATLAS+ttHtataSMATLAS));
		vmuATLASZy1blh.push_back((ggFhZyATLAS.at(i)+VBFhZyATLAS.at(i)+VHhZyATLAS.at(i)+ttHhZyATLAS.at(i))/(ggFZySMATLAS+VBFZySMATLAS+VHZySMATLAS+ttHZySMATLAS));
		
		chisqCMS.push_back(pow((muCMSyy1blh-muCMSyy1)/emuCMSyy1,2.0));
		chisqCMS.at(i)+=pow((muCMSZZ1blh-muCMSZZ1)/emuCMSZZ1,2.0);
		chisqCMS.at(i)+=pow((muCMSZZ2blh-muCMSZZ2)/emuCMSZZ2,2.0);
		chisqCMS.at(i)+=pow((muCMSWW1blh-muCMSWW1)/emuCMSWW1,2.0);
		chisqCMS.at(i)+=pow((muCMSbb1blh-muCMSbb1)/emuCMSbb1,2.0);
		chisqCMS.at(i)+=pow((muCMStata1blh-muCMStata1)/emuCMStata1,2.0);
		chisqCMS.at(i)+=pow((muCMStata2blh-muCMStata2)/emuCMStata2,2.0);
		
		chisqATLAS.push_back(pow((muATLASyy1blh-muATLASyy1)/emuATLASyy1,2.0));
		chisqATLAS.at(i)+=pow((muATLASZZ1blh-muATLASZZ1)/emuATLASZZ1,2.0);
		chisqATLAS.at(i)+=pow((muATLASWW1blh-muATLASWW1)/emuATLASWW1,2.0);
		chisqATLAS.at(i)+=pow((muATLASWW2blh-muATLASWW2)/emuATLASWW2,2.0);
		chisqATLAS.at(i)+=pow((muATLASWW3blh-muATLASWW3)/emuATLASWW3,2.0);
		chisqATLAS.at(i)+=pow((muATLASbb1blh-muATLASbb1)/emuATLASbb1,2.0);
		chisqATLAS.at(i)+=pow((muATLAStata1blh-muATLAStata1)/emuATLAStata1,2.0);
		
		chisq.push_back(pow((muCMSyy1blh-muCMSyy1)/emuCMSyy1,2.0));
		chisq.at(i)+=pow((muCMSZZ1blh-muCMSZZ1)/emuCMSZZ1,2.0);
		chisq.at(i)+=pow((muCMSZZ2blh-muCMSZZ2)/emuCMSZZ2,2.0);
		chisq.at(i)+=pow((muCMSWW1blh-muCMSWW1)/emuCMSWW1,2.0);
		chisq.at(i)+=pow((muCMSbb1blh-muCMSbb1)/emuCMSbb1,2.0);
		chisq.at(i)+=pow((muCMStata1blh-muCMStata1)/emuCMStata1,2.0);
		chisq.at(i)+=pow((muCMStata2blh-muCMStata2)/emuCMStata2,2.0);
		chisq.at(i)+=pow((muATLASyy1blh-muATLASyy1)/emuATLASyy1,2.0);
		chisq.at(i)+=pow((muATLASZZ1blh-muATLASZZ1)/emuATLASZZ1,2.0);
		chisq.at(i)+=pow((muATLASWW1blh-muATLASWW1)/emuATLASWW1,2.0);
		chisq.at(i)+=pow((muATLASWW2blh-muATLASWW2)/emuATLASWW2,2.0);
		chisq.at(i)+=pow((muATLASWW3blh-muATLASWW3)/emuATLASWW3,2.0);
		chisq.at(i)+=pow((muATLASbb1blh-muATLASbb1)/emuATLASbb1,2.0);
		chisq.at(i)+=pow((muATLAStata1blh-muATLAStata1)/emuATLAStata1,2.0);
		
		
		rchisqCMS.push_back(pow((muCMSyy1blh-muCMSyy1)/emuCMSyy1,2.0));
		rchisqCMS.at(i)+=pow((muCMSZZ1blh-muCMSZZ1)/emuCMSZZ1,2.0);
		rchisqCMS.at(i)+=pow((muCMSWW1blh-muCMSWW1)/emuCMSWW1,2.0);
		
		rchisqATLAS.push_back(pow((muATLASyy1blh-muATLASyy1)/emuATLASyy1,2.0));
		rchisqATLAS.at(i)+=pow((muATLASWW1blh-muATLASWW1)/emuATLASWW1,2.0);
		rchisqATLAS.at(i)+=pow((muATLASWW3blh-muATLASWW3)/emuATLASWW3,2.0);
		
		rchisq.push_back(pow((muCMSyy1blh-muCMSyy1)/emuCMSyy1,2.0));
		rchisq.at(i)+=pow((muCMSZZ1blh-muCMSZZ1)/emuCMSZZ1,2.0);
		rchisq.at(i)+=pow((muCMSWW1blh-muCMSWW1)/emuCMSWW1,2.0);
		rchisq.at(i)+=pow((muATLASyy1blh-muATLASyy1)/emuATLASyy1,2.0);
		rchisq.at(i)+=pow((muATLASWW1blh-muATLASWW1)/emuATLASWW1,2.0);
		rchisq.at(i)+=pow((muATLASWW3blh-muATLASWW3)/emuATLASWW3,2.0);
		
		if ( i%1000 == 0 ) {
			cout << "Chi sq calc: " << i << endl;
		}
		//cout << "Chisq calculated for point " << i << endl;
	}
	
	double mhlimit;
	vector<double> b_ratiotest;
	
	
	for ( int i=0; i<fl.size(); i++ ) {
		if ( mhh.at(i) >= 128.2 && mhh.at(i) < 160 ) {
			mhlimit=-18813.1 - 5.91403E7/(mhh.at(i)*mhh.at(i)) + 1.72403E6/(mhh.at(i)) + 91.227 *(mhh.at(i)) - 0.165836 *(mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 160 && mhh.at(i) < 200 ) {
			mhlimit=49437.5 + 2.69503E8/(mhh.at(i)*mhh.at(i)) - 5.96676E6/(mhh.at(i)) - 181.663 *(mhh.at(i)) + 0.249854 *(mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 200 && mhh.at(i) < 300 ) {
			mhlimit=8.10975 - 0.0266974 *(mhh.at(i));
		}
		else if ( mhh.at(i) >= 300 && mhh.at(i) < 600 ) {
			mhlimit=0.084537 + 0.0000532794 *(mhh.at(i));
		}
		else if ( mhh.at(i) >= 600 ) {
			mhlimit=-1.14292 + 0.00209904 *(mhh.at(i));
		}
		else {
			mhlimit = 100;
		}
		b_ratiotest.push_back(true);
		if ( Htest.at(i) > mhlimit ) {
			b_ratiotest.at(i) = false;
		}
	}
	
	cout << "There are " << fl.size() << " data points." << endl;
	
	minATLAS=99999999;
	minCMS=99999999;
	min=99999999;
	for ( int i=0; i<fl.size(); i++) {
		if ( chisqCMS.at(i) < minCMS ) {
			minCMS = chisqCMS.at(i);
		}
		if ( chisqATLAS.at(i) < minATLAS ) {
			minATLAS = chisqATLAS.at(i);
		}
		if ( chisq.at(i) < min ) {
			min = chisq.at(i);
		}
	}
	
	rminATLAS=99999999;
	rminCMS=99999999;
	rmin=99999999;
	for ( int i=0; i<fl.size(); i++) {
		if ( rchisqCMS.at(i) < rminCMS ) {
			rminCMS = rchisqCMS.at(i);
		}
		if ( rchisqATLAS.at(i) < rminATLAS ) {
			rminATLAS = rchisqATLAS.at(i);
		}
		if ( rchisq.at(i) < rmin ) {
			rmin = rchisq.at(i);
		}
	}
	
	//cout << "Minimums found" << endl;
	
	cout << "The chisq minimum values are: " << endl;
	cout << "CMS \t" << minCMS << endl;
	cout << "ATLAS \t" << minATLAS << endl;
	cout << "Comb \t" << min << endl;
	
	cout << "The rchisq minimum values are: " << endl;
	cout << "CMS \t" << rminCMS << endl;
	cout << "ATLAS \t" << rminATLAS << endl;
	cout << "Comb \t" << rmin << endl;
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	
	vsm = 1.0/sqrt(sqrt(2)*0.0000116637);
	ofstream data;
	 data.open("dataofinterest2.txt");
	 for ( int i=0; i<fl.size(); i++) {
		 sg = tG.at(i)/sqrt(1.0+tG.at(i)*tG.at(i));
		 sB = tB.at(i)/sqrt(1.0+tB.at(k)*tB.at(i));
		 cB = sqrt(1.0-sB*sB);
		 vblh=vsm*sqrt(1. + ( 1. /6. + ( 5. /4. - 2.*sg*sg)*( (fl.at(i)*fl.at(i)) /(fl.at(i)*fl.at(i) + fh.at(i)*fh.at(i))))  vsm*vsm/(fl.at(i)*fl.at(i)));
		 lambda0 = (mh0.at(i)*mh0.at(i)/(vblh*vblh))*((mh0.at(i)*mh0.at(i) - ma0.at(i)*ma0.at(i))/(mh0.at(i)*mh0.at(i) - 2.0*sB*cB* ma0.at(i)*ma0.at(i)));
		 Bmu =  (lambda0 vblh*vblh + ma0.at(i)*ma0.at(i))*(2.0*sB*cB)/2.0;
		 ta = (1.0/(Bmu - lambda0 *vblh*vblh*2.0*sB*cB))*((Bmu*(cB/(2.0*sB) - sB/(2.0*cB)))+sqrt( Bmu*Bmu/(4.0*sB*sB*cB*cB)   -2.0*lambda0*Bmu*vblh*vblh*2.0*sB*cB  +lambda0*lambda0*vblh*vblh*vblh*vblh*4.0*sB*sB*cB*cB     ));
	 //if ( vmuATLAStata1blh.at(i) > 4.2 && vmuATLAStata1blh.at(i) < 4.6 && vmuATLASZy1blh.at(i) > 0.3 && vmuATLASZy1blh.at(i) < 0.6 ) {
	 data << fl.at(i) << " \t ";
	 data << fh.at(i) << " \t ";
	 data << mh0.at(i) << " \t ";
	 data << ma0.at(i) << " \t ";
	 data << mhh.at(i) << " \t ";
	 data << tB.at(i) << " \t ";
	 data << tG.at(i) << " \t ";
	 data << t12.at(i) << " \t ";
	 data << t13.at(i) << " \t ";
	 data << sigfact.at(i) << " \t ";
	 data << kG.at(i) << " \t ";
	 data << kY.at(i) << " \t ";
	 data << kL.at(i) << " \t ";
	 data << Cyy.at(i) << " \t ";
	 data << Cgg.at(i) << " \t ";
	 data << CWW.at(i) << " \t ";
	 data << CZZ.at(i) << " \t ";
	 data << Cbb.at(i) << " \t ";
	 data << Ctata.at(i) << " \t ";
	 data << hmuATLASyy1blh.at(i) << " \t ";
	 data << AmuATLASyy1blh.at(i) << " \t ";
	 data << HmuATLASyy1blh.at(i) << " \t ";
	 data << hmuATLAStata1blh.at(i) << " \t ";
	 data << AmuATLAStata1blh.at(i) << " \t ";
	 data << HmuATLAStata1blh.at(i) << " \t ";
	 data << hmuATLASWW1blh.at(i) << " \t ";
	 data << HmuATLASWW1blh.at(i) << " \t ";
	 data << endl;
	// }
	 }
	 data.close();
	
	
	for ( int i=0; i<fl.size(); i++) {
		if ( chisqCMS.at(i)-minCMS < 1 ) {
			CMSlt1++;
		}
		else if ( chisqCMS.at(i)-minCMS < 4 ) {
			CMSlt2++;
		}
		else if ( chisqCMS.at(i)-minCMS < 9 ) {
			CMSlt3++;
		}
		else {
			CMSgt3++;
		}
		if ( chisqATLAS.at(i)-minATLAS < 1 ) {
			ATLASlt1++;
		}
		else if ( chisqATLAS.at(i)-minATLAS < 4 ) {
			ATLASlt2++;
		}
		else if ( chisqATLAS.at(i)-minATLAS < 9 ) {
			ATLASlt3++;
		}
		else {
			ATLASgt3++;
		}
		if ( chisq.at(i)-min < 1 ) {
			lt1++;
		}
		else if ( chisq.at(i)-min < 4 ) {
			lt2++;
		}
		else if ( chisq.at(i)-min < 9 ) {
			lt3++;
		}
		else {
			gt3++;
		}
	}
	
	rCMSlt1=0;
	rCMSlt2=0;
	rCMSlt3=0;
	rCMSgt3=0;
	rATLASlt1=0;
	rATLASlt2=0;
	rATLASlt3=0;
	rATLASgt3=0;
	rlt1=0;
	rlt2=0;
	rlt3=0;
	rgt3=0;
	for ( int i=0; i<fl.size(); i++) {
		if ( rchisqCMS.at(i)-rminCMS < 1 ) {
			rCMSlt1++;
		}
		else if ( rchisqCMS.at(i)-rminCMS < 4 ) {
			rCMSlt2++;
		}
		else if ( rchisqCMS.at(i)-rminCMS < 9 ) {
			rCMSlt3++;
		}
		else {
			rCMSgt3++;
		}
		if ( rchisqATLAS.at(i)-rminATLAS < 1 ) {
			rATLASlt1++;
		}
		else if ( rchisqATLAS.at(i)-rminATLAS < 4 ) {
			rATLASlt2++;
		}
		else if ( rchisqATLAS.at(i)-rminATLAS < 9 ) {
			rATLASlt3++;
		}
		else {
			rATLASgt3++;
		}
		if ( rchisq.at(i)-rmin < 1 ) {
			rlt1++;
		}
		else if ( rchisq.at(i)-rmin < 4 ) {
			rlt2++;
		}
		else if ( rchisq.at(i)-rmin < 9 ) {
			rlt3++;
		}
		else {
			rgt3++;
		}
	}
	
	cout << "CMS numbers: " << rCMSlt1 << " " << rCMSlt2 << " " << rCMSlt3 << " " << rCMSgt3 << endl;
	cout << "ATLAS numbers: " << rATLASlt1 << " " << rATLASlt2 << " " << rATLASlt3 << " " << rATLASgt3 << endl;
	
	/*cout << fl.size() << endl;
	 cout << fh.size() << endl;
	 cout << mh0.size() << endl;
	 cout << ma0.size() << endl;
	 cout << tB.size() << endl;
	 cout << tG.size() << endl;
	 cout << rDSgg.size() << endl;
	 cout << rSgg.size() << endl;
	 cout << rDSyy.size() << endl;
	 cout << rSyy.size() << endl;
	 cout << Cyy.size() << endl;
	 cout << Cgg.size() << endl;
	 cout << CWW.size() << endl;
	 cout << CZZ.size() << endl;
	 cout << Cbb.size() << endl;
	 cout << Ctata.size() << endl;
	 cout << vmuCMSyy1blh.size() << endl;
	 cout << vmuCMSZZ1blh.size() << endl;
	 cout << vmuCMSWW1blh.size() << endl;
	 cout << vmuCMSbb1blh.size() << endl;
	 cout << vmuCMStata1blh.size() << endl;
	 cout << vmuCMSZy1blh.size() << endl;
	 cout << vmuATLASyy1blh.size() << endl;
	 cout << vmuATLASZZ1blh.size() << endl;
	 cout << vmuATLASWW1blh.size() << endl;
	 cout << vmuATLASbb1blh.size() << endl;
	 cout << vmuATLAStata1blh.size() << endl;
	 cout << vmuATLASZy1blh.size() << endl;*/
	
	
	
	for ( int i=0; i<fl.size(); i++) {
		//cout << i << endl;
		if ( fl.at(i) < flmin ) {
			flmin = fl.at(i);
		}
		else if (fl.at(i) > flmax ) {
			flmax = fl.at(i);
		}
		if ( fh.at(i) < fhmin ) {
			fhmin = fh.at(i);
		}
		else if (fh.at(i) > fhmax ) {
			fhmax = fh.at(i);
		}
		if ( mh0.at(i) < mh0min ) {
			mh0min = mh0.at(i);
		}
		else if (mh0.at(i) > mh0max ) {
			mh0max = mh0.at(i);
		}
		if ( ma0.at(i) < ma0min ) {
			ma0min = ma0.at(i);
		}
		else if (ma0.at(i) > ma0max ) {
			ma0max = ma0.at(i);
		}
		if ( mhh.at(i) > mhhmax ) {
			mhhmax = mhh.at(i);
		}
		else if ( mhh.at(i) < mhhmin ) {
			mhhmin = mhh.at(i);
		}
		if ( deltam.at(i) < deltammin ) {
			deltammin = deltam.at(i);
		}
		else if (deltam.at(i) > deltammax ) {
			deltammax = deltam.at(i);
		}
		if ( tB.at(i) < tBmin ) {
			tBmin = tB.at(i);
		}
		else if (tB.at(i) > tBmax ) {
			tBmax = tB.at(i);
		}
		if ( tG.at(i) < tGmin ) {
			tGmin = tG.at(i);
		}
		else if (tG.at(i) > tGmax ) {
			tGmax = tG.at(i);
		}
		if ( rDSgg.at(i) < DSggmin ) {
			DSggmin = rDSgg.at(i);
		}
		else if (rDSgg.at(i) > DSggmax ) {
			DSggmax = rDSgg.at(i);
		}
		if ( rSgg.at(i) < Sggmin ) {
			Sggmin = rSgg.at(i);
		}
		else if (rSgg.at(i) > Sggmax ) {
			Sggmax = rSgg.at(i);
		}
		if ( rDSyy.at(i) < DSyymin ) {
			DSyymin = rDSyy.at(i);
		}
		else if (rDSyy.at(i) > DSyymax ) {
			DSyymax = rDSyy.at(i);
		}
		if ( rSyy.at(i) < Syymin ) {
			Syymin = rSyy.at(i);
		}
		else if (rSyy.at(i) > Syymax ) {
			Syymax = rSyy.at(i);
		}
		if ( Cyy.at(i) < Cyymin ) {
			Cyymin = Cyy.at(i);
		}
		else if (Cyy.at(i) > Cyymax ) {
			Cyymax = Cyy.at(i);
		}
		if ( Cgg.at(i) < Cggmin ) {
			Cggmin = Cgg.at(i);
		}
		else if (Cgg.at(i) > Cggmax ) {
			Cggmax = Cgg.at(i);
		}
		if ( CWW.at(i) < CWWmin ) {
			CWWmin = CWW.at(i);
		}
		else if (CWW.at(i) > CWWmax ) {
			CWWmax = CWW.at(i);
		}
		if ( CZZ.at(i) < CZZmin ) {
			CZZmin = CZZ.at(i);
		}
		else if (CZZ.at(i) > CZZmax ) {
			CZZmax = CZZ.at(i);
		}
		if ( Cbb.at(i) < Cbbmin ) {
			Cbbmin = Cbb.at(i);
		}
		else if (Cbb.at(i) > Cbbmax ) {
			Cbbmax = Cbb.at(i);
		}
		if ( Ctata.at(i) < Ctatamin ) {
			Ctatamin = Ctata.at(i);
		}
		else if (Ctata.at(i) > Ctatamax ) {
			Ctatamax = Ctata.at(i);
		}
		if ( vmuCMSyy1blh.at(i) < muyymin ) {
			muyymin = vmuCMSyy1blh.at(i);
		}
		else if (vmuCMSyy1blh.at(i) > muyymax ) {
			muyymax = vmuCMSyy1blh.at(i);
		}
		if ( vmuCMSZZ1blh.at(i) < muZZmin ) {
			muZZmin = vmuCMSZZ1blh.at(i);
		}
		else if (vmuCMSZZ1blh.at(i) > muZZmax ) {
			muZZmax = vmuCMSZZ1blh.at(i);
		}
		if ( vmuCMSWW1blh.at(i) < muWWmin ) {
			muWWmin = vmuCMSWW1blh.at(i);
		}
		else if (vmuCMSWW1blh.at(i) > muWWmax ) {
			muWWmax = vmuCMSWW1blh.at(i);
		}
		if ( vmuCMStata1blh.at(i) < mutatamin ) {
			mutatamin = vmuCMStata1blh.at(i);
		}
		else if (vmuCMStata1blh.at(i) > mutatamax ) {
			mutatamax = vmuCMStata1blh.at(i);
		}
		if ( vmuATLASyy1blh.at(i) < muyymin ) {
			muyymin = vmuATLASyy1blh.at(i);
		}
		else if (vmuATLASyy1blh.at(i) > muyymax ) {
			muyymax = vmuATLASyy1blh.at(i);
		}
		if ( vmuATLASZZ1blh.at(i) < muZZmin ) {
			muZZmin = vmuATLASZZ1blh.at(i);
		}
		else if (vmuATLASZZ1blh.at(i) > muZZmax ) {
			muZZmax = vmuATLASZZ1blh.at(i);
		}
		if ( vmuATLASWW1blh.at(i) < muWWmin ) {
			muWWmin = vmuATLASWW1blh.at(i);
		}
		else if (vmuATLASWW1blh.at(i) > muWWmax ) {
			muWWmax = vmuATLASWW1blh.at(i);
		}
		if ( vmuATLASbb1blh.at(i) < mubbmin ) {
			mubbmin = vmuATLASbb1blh.at(i);
		}
		else if (vmuATLASbb1blh.at(i) > mubbmax ) {
			mubbmax = vmuATLASbb1blh.at(i);
		}
		if ( vmuATLAStata1blh.at(i) < mutatamin ) {
			mutatamin = vmuATLAStata1blh.at(i);
		}
		else if (vmuATLAStata1blh.at(i) > mutatamax ) {
			mutatamax = vmuATLAStata1blh.at(i);
		}
		if ( vmuCMSbb1blh.at(i) < mubbmin ) {
			mubbmin = vmuCMSbb1blh.at(i);
		}
		else if (vmuCMSbb1blh.at(i) > mubbmax ) {
			mubbmax = vmuCMSbb1blh.at(i);
		}
		if ( vmuCMSZy1blh.at(i) < muZymin ) {
			muZymin = vmuCMSZy1blh.at(i);
		}
		else if (vmuCMSZy1blh.at(i) > muZymax ) {
			muZymax = vmuCMSZy1blh.at(i);
		}
		if ( vmuATLASZy1blh.at(i) < muZymin ) {
			muZymin = vmuATLASZy1blh.at(i);
		}
		else if (vmuATLASZy1blh.at(i) > muZymax ) {
			muZymax = vmuATLASZy1blh.at(i);
		}
	}
	if ( !degen ) {
		cout << "Setting non-degenerate maximums." << endl;
		Cyymax = 4;
		muZZmax = 1.5;
		muWWmax = 1.5;
		muyymax = 2.5;
		muZymax = 1.5;
		mubbmax = 1.5;
		mutatamax = 6.0;
		Syymin = 0;
		Syymax = 10;
		DSyymax = 4;
		DSyymin = 0;
		Sggmin = 0;
		Sggmax = 1;
		CWWmax = 1.2;
		CZZmax = 1.2;
	}
	else {
		cout << "Setting degenerate maximums." << endl;
		Cyymax = 2.0;
		muZZmax = 1.5;
		muWWmax = 1.5;
		muyymax = 2.5;
		muZymax = 1.5;
		mubbmax = 1.5;
		mutatamax = 6.0;
		Syymin = 0;
		Syymax = 10;
		DSyymax = 4;
		DSyymin = 0;
		Sggmin = 0;
		Sggmax = 1;
		Cggmax = 2.2;
		Cbbmax = 2.2;
		CWWmax = 1.2;
		CZZmax = 1.2;
		Ctatamax = 2.2;
		mhhmax = 300;
	}
	
	cout << "MH maximum is " << mhhmax << endl;
	
	
	
	//Keep these lines the same. There should be no need to change them.
	plot figures;
	figures.ratiotest = b_ratiotest;
	figures.chisqCMS = chisqCMS;
	figures.chisqATLAS = chisqATLAS;
	figures.chisq = chisq;
	figures.minCMS = minCMS;
	figures.minATLAS = minATLAS;
	figures.min = min;
	figures.CMSlt1=CMSlt1;
	figures.CMSlt2=CMSlt2;
	figures.CMSlt3=CMSlt3;
	figures.CMSgt3=CMSgt3;
	figures.ATLASlt1=ATLASlt1;
	figures.ATLASlt2=ATLASlt2;
	figures.ATLASlt3=ATLASlt3;
	figures.ATLASgt3=ATLASgt3;
	figures.lt1=lt1;
	figures.lt2=lt2;
	figures.lt3=lt3;
	figures.gt3=gt3;
	
	for ( int i=0; i<rDSyy.size(); i++) {
		DSyy.push_back(sqrt(rDSyy.at(i)*rDSyy.at(i)+iDSyy.at(i)*iDSyy.at(i)));
		Syy.push_back(sqrt(rSyy.at(i)*rSyy.at(i)+iSyy.at(i)*iSyy.at(i)));
		DSgg.push_back(sqrt(rDSgg.at(i)*rDSgg.at(i)+iDSgg.at(i)*iDSgg.at(i)));
		Sgg.push_back(sqrt(rSgg.at(i)*rSgg.at(i)+iSgg.at(i)*iSgg.at(i)));
	}
	
	vector<double> muZZcombined;
	for ( int i=0; i<vmuCMSZZ1blh.size(); i++) {
		muZZcombined.push_back(vmuCMSZZ1blh.at(i)/2.0+vmuATLASZZ1blh.at(i)/2.0);
	}
	vector<double> muWWcombined;
	for ( int i=0; i<vmuCMSWW1blh.size(); i++) {
		muWWcombined.push_back(vmuCMSWW1blh.at(i)/2.0+vmuATLASWW1blh.at(i)/2.0);
	}
	vector<double> mubbcombined;
	for ( int i=0; i<vmuCMSbb1blh.size(); i++) {
		mubbcombined.push_back(vmuCMSbb1blh.at(i)/2.0+vmuATLASbb1blh.at(i)/2.0);
	}
	vector<double> mutatacombined;
	for ( int i=0; i<vmuCMStata1blh.size(); i++) {
		mutatacombined.push_back(vmuCMStata1blh.at(i)/2.0+vmuATLAStata1blh.at(i)/2.0);
	}
	vector<double> muZycombined;
	for ( int i=0; i<vmuCMSZy1blh.size(); i++) {
		muZycombined.push_back(vmuCMSZy1blh.at(i)/2.0+vmuATLASZy1blh.at(i)/2.0);
	}
	vector<double> muyycombined;
	for ( int i=0; i<vmuCMSyy1blh.size(); i++) {
		muyycombined.push_back(vmuCMSyy1blh.at(i)/2.0+vmuATLASyy1blh.at(i)/2.0);
	}
	
	//To make a new plot, this is the parts that need to be altered.
	//If a value you want is not in a vector, you can create a new vector at any point and calculate the values.
	//The combined value in the third plot is taken as (CMS+ATLAS)/2.0 (the basic average).
	//The format for generate_plot is:
	//2igures.generate_plot("filename","xlabel","ylabel",xmin,xmax,ymin,ymax);
	
	
	//Plots involving CZZ
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = CZZ;
	figures.ATLASx = fl;
	figures.ATLASy = CZZ;
	figures.COMBx = fl;
	figures.COMBy = CZZ;
	if ( degen ) {
		figures.generate_plot("./deg/xfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		figures.generate_plot("./norm/xfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	cout << "Testing" << endl;
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = CZZ;
	figures.ATLASx = mh0;
	figures.ATLASy = CZZ;
	figures.COMBx = mh0;
	figures.COMBy = CZZ;
	
	cout << "Testing" << endl;
	
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		figures.generate_plot("./norm/xmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = CZZ;
	figures.ATLASx = deltam;
	figures.ATLASy = CZZ;
	figures.COMBx = deltam;
	figures.COMBy = CZZ;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}*/
	
	//Plots involving CWW
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = CWW;
	figures.ATLASx = fl;
	figures.ATLASy = CWW;
	figures.COMBx = fl;
	figures.COMBy = CWW;
	if ( degen ) {
		figures.generate_plot("./deg/xfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		figures.generate_plot("./norm/xfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = CWW;
	figures.ATLASx = mh0;
	figures.ATLASy = CWW;
	figures.COMBx = mh0;
	figures.COMBy = CWW;
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		figures.generate_plot("./norm/xmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = CWW;
	figures.ATLASx = deltam;
	figures.ATLASy = CWW;
	figures.COMBx = deltam;
	figures.COMBy = CWW;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}*/
	
	
	//Plots involving Cyy
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = Cyy;
	figures.ATLASx = fl;
	figures.ATLASy = Cyy;
	figures.COMBx = fl;
	figures.COMBy = Cyy;
	if ( degen ) {
		figures.generate_plot("./deg/xfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		figures.generate_plot("./norm/xfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = Cyy;
	figures.ATLASx = mh0;
	figures.ATLASy = Cyy;
	figures.COMBx = mh0;
	figures.COMBy = Cyy;
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		figures.generate_plot("./norm/xmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = Cyy;
	figures.ATLASx = deltam;
	figures.ATLASy = Cyy;
	figures.COMBx = deltam;
	figures.COMBy = Cyy;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}*/
	
	
	//Plots involving Cgg
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = Cgg;
	figures.ATLASx = fl;
	figures.ATLASy = Cgg;
	figures.COMBx = fl;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = Cgg;
	figures.ATLASx = mh0;
	figures.ATLASy = Cgg;
	figures.COMBx = mh0;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = Cgg;
	figures.ATLASx = deltam;
	figures.ATLASy = Cgg;
	figures.COMBx = deltam;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy = Cgg;
	figures.ATLASx = tB;
	figures.ATLASy = Cgg;
	figures.COMBx = tB;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = Cbb;
	figures.ATLASx = fl;
	figures.ATLASy = Cbb;
	figures.COMBx = fl;
	figures.COMBy = Cbb;
	if ( degen ) {
		figures.generate_plot("./deg/xfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		figures.generate_plot("./norm/xfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = Cbb;
	figures.ATLASx = mh0;
	figures.ATLASy = Cbb;
	figures.COMBx = mh0;
	figures.COMBy = Cbb;
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		figures.generate_plot("./norm/xmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = Cbb;
	figures.ATLASx = deltam;
	figures.ATLASy = Cbb;
	figures.COMBx = deltam;
	figures.COMBy = Cbb;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}*/
	
	
	//Plots involving Ctata
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = Ctata;
	figures.ATLASx = fl;
	figures.ATLASy = Ctata;
	figures.COMBx = fl;
	figures.COMBy = Ctata;
	if ( degen ) {
		figures.generate_plot("./deg/xfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		figures.generate_plot("./norm/xfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = mh0;
	figures.CMSy = Ctata;
	figures.ATLASx = mh0;
	figures.ATLASy = Ctata;
	figures.COMBx = mh0;
	figures.COMBy = Ctata;
	if ( degen ) {
		figures.generate_plot("./deg/xmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		figures.generate_plot("./norm/xmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = Ctata;
	figures.ATLASx = deltam;
	figures.ATLASy = Ctata;
	figures.COMBx = deltam;
	figures.COMBy = Ctata;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}*/
	
	//parameter vs parameter plots
	
	/*figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = tB;
	figures.ATLASx = fl;
	figures.ATLASy = tB;
	figures.COMBx = fl;
	figures.COMBy = tB;
	if ( degen ) {
		figures.generate_plot("./deg/xfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		figures.generate_plot("./norm/xfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = deltam;
	figures.ATLASx = fl;
	figures.ATLASy = deltam;
	figures.COMBx = fl;
	figures.COMBy = deltam;
	if ( degen ) {
		figures.generate_plot("./deg/xfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		figures.generate_plot("./norm/xfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = tB;
	figures.ATLASx = deltam;
	figures.ATLASy = tB;
	figures.COMBx = deltam;
	figures.COMBy = tB;
	if ( degen ) {
		figures.generate_plot("./deg/xdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		figures.generate_plot("./norm/xdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = DSgg;
	figures.CMSy = DSyy;
	figures.ATLASx = DSgg;
	figures.ATLASy = DSyy;
	figures.COMBx = DSgg;
	figures.COMBy = DSyy;
	if ( degen ) {
		figures.generate_plot("./deg/xdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		figures.generate_plot("./norm/xdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Syy;
	figures.CMSy = DSyy;
	figures.ATLASx = Syy;
	figures.ATLASy = DSyy;
	figures.COMBx = Syy;
	figures.COMBy = DSyy;
	if ( degen ) {
		figures.generate_plot("./deg/xsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		figures.generate_plot("./norm/xsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Sgg;
	figures.CMSy = DSgg;
	figures.ATLASx = Sgg;
	figures.ATLASy = DSgg;
	figures.COMBx = Sgg;
	figures.COMBy = DSgg;
	if ( degen ) {
		figures.generate_plot("./deg/xsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		figures.generate_plot("./norm/xsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}*/
	
	/*figures.clear_plot_data();
	figures.CMSx = Cyy;
	figures.CMSy = Cgg;
	figures.ATLASx = Cyy;
	figures.ATLASy = Cgg;
	figures.COMBx = Cyy;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = CZZ;
	figures.CMSy = Cyy;
	figures.ATLASx = CZZ;
	figures.ATLASy = Cyy;
	figures.COMBx = CZZ;
	figures.COMBy = Cyy;
	if ( degen ) {
		figures.generate_plot("./deg/xczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		figures.generate_plot("./norm/xczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = CZZ;
	figures.CMSy = Cgg;
	figures.ATLASx = CZZ;
	figures.ATLASy = Cgg;
	figures.COMBx = CZZ;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = CZZ;
	figures.CMSy = CWW;
	figures.ATLASx = CZZ;
	figures.ATLASy = CWW;
	figures.COMBx = CZZ;
	figures.COMBy = CWW;
	if ( degen ) {
		figures.generate_plot("./deg/xczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		figures.generate_plot("./norm/xczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Cbb;
	figures.CMSy = Cgg;
	figures.ATLASx = Cbb;
	figures.ATLASy = Cgg;
	figures.COMBx = Cbb;
	figures.COMBy = Cgg;
	if ( degen ) {
		figures.generate_plot("./deg/xcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		figures.generate_plot("./norm/xcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	figures.clear_plot_data();
	figures.CMSx = Cgg;
	figures.CMSy = vmuCMSyy1blh;
	figures.ATLASx = Cgg;
	figures.ATLASy = vmuATLASyy1blh;
	figures.COMBx = Cgg;
	figures.COMBy = muyycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		figures.generate_plot("./norm/xcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Cyy;
	figures.CMSy = vmuCMSyy1blh;
	figures.ATLASx = Cyy;
	figures.ATLASy = vmuATLASyy1blh;
	figures.COMBx = Cyy;
	figures.COMBy = muyycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		figures.generate_plot("./norm/xcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Cgg;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = Cgg;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = Cgg;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/xcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = CWW;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = CWW;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = CWW;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/xcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = Cgg;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = Cgg;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = Cgg;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/xcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = CZZ;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = CZZ;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = CZZ;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/xczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//mu vs mu plots
	
	/*figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/xmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/xmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/xmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/xmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/xmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/xmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/xmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/xmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/xmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/xmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/xmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/xmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSbb1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASbb1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = mubbcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/xmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSbb1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASbb1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = mubbcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/xmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMStata1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLAStata1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = mutatacombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/xmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/xmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}*/
	
	
	
	
	
	//#####################################################################
	//#####################################################################
	//#####################################################################
	//#####################################################################
	//#####################################################################
	//#####################################################################
	//#####################################################################
	//Below are the plots generated by only using the measured mu values
	//with errors less than 0.4. The reason for looking at these values is
	//that the CMS and ATLAS results display two very distinct and
	//important possibilities: CMS results suggest something very SM-like,
	//while ATLAS results suggest something distinctly non-SM-like, with
	//Cyy != 1, that is understandable in the near-degenerate scenario.
	
	
	
	plot rfigures;
	rfigures.ratiotest = b_ratiotest;
	rfigures.chisqCMS = rchisqCMS;
	rfigures.chisqATLAS = rchisqATLAS;
	rfigures.chisq = rchisq;
	rfigures.minCMS = rminCMS;
	rfigures.minATLAS = rminATLAS;
	rfigures.min = rmin;
	rfigures.CMSlt1=rCMSlt1;
	rfigures.CMSlt2=rCMSlt2;
	rfigures.CMSlt3=rCMSlt3;
	rfigures.CMSgt3=rCMSgt3;
	rfigures.ATLASlt1=rATLASlt1;
	rfigures.ATLASlt2=rATLASlt2;
	rfigures.ATLASlt3=rATLASlt3;
	rfigures.ATLASgt3=rATLASgt3;
	rfigures.lt1=rlt1;
	rfigures.lt2=rlt2;
	rfigures.lt3=rlt3;
	rfigures.gt3=rgt3;
	
	
	
	//Plots involving CZZ
	/*rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = CZZ;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = CZZ;
	rfigures.COMBx = fl;
	rfigures.COMBy = CZZ;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = CZZ;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = CZZ;
	rfigures.COMBx = mh0;
	rfigures.COMBy = CZZ;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = CZZ;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = CZZ;
	rfigures.COMBx = deltam;
	rfigures.COMBy = CZZ;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = CZZ;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = CZZ;
	rfigures.COMBx = tB;
	rfigures.COMBy = CZZ;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//Plots involving CWW
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = CWW;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = CWW;
	rfigures.COMBx = fl;
	rfigures.COMBy = CWW;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = CWW;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = CWW;
	rfigures.COMBx = mh0;
	rfigures.COMBy = CWW;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = CWW;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = CWW;
	rfigures.COMBx = deltam;
	rfigures.COMBy = CWW;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = CWW;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = CWW;
	rfigures.COMBx = tB;
	rfigures.COMBy = CWW;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = Cyy;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = Cyy;
	rfigures.COMBx = fl;
	rfigures.COMBy = Cyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = Cyy;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = Cyy;
	rfigures.COMBx = mh0;
	rfigures.COMBy = Cyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = Cyy;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = Cyy;
	rfigures.COMBx = deltam;
	rfigures.COMBy = Cyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = Cyy;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = Cyy;
	rfigures.COMBx = tB;
	rfigures.COMBy = Cyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = fl;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = mh0;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = deltam;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = tB;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = Cbb;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = Cbb;
	rfigures.COMBx = fl;
	rfigures.COMBy = Cbb;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = Cbb;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = Cbb;
	rfigures.COMBx = mh0;
	rfigures.COMBy = Cbb;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = Cbb;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = Cbb;
	rfigures.COMBx = deltam;
	rfigures.COMBy = Cbb;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = Cbb;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = Cbb;
	rfigures.COMBx = tB;
	rfigures.COMBy = Cbb;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//Plots involving Ctata
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = Ctata;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = Ctata;
	rfigures.COMBx = fl;
	rfigures.COMBy = Ctata;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = mh0;
	rfigures.CMSy = Ctata;
	rfigures.ATLASx = mh0;
	rfigures.ATLASy = Ctata;
	rfigures.COMBx = mh0;
	rfigures.COMBy = Ctata;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = Ctata;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = Ctata;
	rfigures.COMBx = deltam;
	rfigures.COMBy = Ctata;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = Ctata;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = Ctata;
	rfigures.COMBx = tB;
	rfigures.COMBy = Ctata;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	
	//f vs mu plots
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMSZZ1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLASZZ1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = muZZcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMSWW1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMSyy1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLASyy1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = muyycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMSbb1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLASbb1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = mubbcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = fl;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	
	//parameter vs parameter plots
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = tB;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = tB;
	rfigures.COMBx = fl;
	rfigures.COMBy = tB;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = fl;
	rfigures.CMSy = deltam;
	rfigures.ATLASx = fl;
	rfigures.ATLASy = deltam;
	rfigures.COMBx = fl;
	rfigures.COMBy = deltam;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		rfigures.generate_plot("./rnorm/rfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = deltam;
	rfigures.CMSy = tB;
	rfigures.ATLASx = deltam;
	rfigures.ATLASy = tB;
	rfigures.COMBx = deltam;
	rfigures.COMBy = tB;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = DSgg;
	rfigures.CMSy = DSyy;
	rfigures.ATLASx = DSgg;
	rfigures.ATLASy = DSyy;
	rfigures.COMBx = DSgg;
	rfigures.COMBy = DSyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Syy;
	rfigures.CMSy = DSyy;
	rfigures.ATLASx = Syy;
	rfigures.ATLASy = DSyy;
	rfigures.COMBx = Syy;
	rfigures.COMBy = DSyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Sgg;
	rfigures.CMSy = DSgg;
	rfigures.ATLASx = Sgg;
	rfigures.ATLASy = DSgg;
	rfigures.COMBx = Sgg;
	rfigures.COMBy = DSgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*rfigures.clear_plot_data();
	rfigures.CMSx = Cyy;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = Cyy;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = Cyy;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = CZZ;
	rfigures.CMSy = Cyy;
	rfigures.ATLASx = CZZ;
	rfigures.ATLASy = Cyy;
	rfigures.COMBx = CZZ;
	rfigures.COMBy = Cyy;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = CZZ;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = CZZ;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = CZZ;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = CZZ;
	rfigures.CMSy = CWW;
	rfigures.ATLASx = CZZ;
	rfigures.ATLASy = CWW;
	rfigures.COMBx = CZZ;
	rfigures.COMBy = CWW;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Cbb;
	rfigures.CMSy = Cgg;
	rfigures.ATLASx = Cbb;
	rfigures.ATLASy = Cgg;
	rfigures.COMBx = Cbb;
	rfigures.COMBy = Cgg;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Cgg;
	rfigures.CMSy = vmuCMSyy1blh;
	rfigures.ATLASx = Cgg;
	rfigures.ATLASy = vmuATLASyy1blh;
	rfigures.COMBx = Cgg;
	rfigures.COMBy = muyycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Cyy;
	rfigures.CMSy = vmuCMSyy1blh;
	rfigures.ATLASx = Cyy;
	rfigures.ATLASy = vmuATLASyy1blh;
	rfigures.COMBx = Cyy;
	rfigures.COMBy = muyycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Cgg;
	rfigures.CMSy = vmuCMSWW1blh;
	rfigures.ATLASx = Cgg;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = Cgg;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = CWW;
	rfigures.CMSy = vmuCMSWW1blh;
	rfigures.ATLASx = CWW;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = CWW;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = Cgg;
	rfigures.CMSy = vmuCMSZZ1blh;
	rfigures.ATLASx = Cgg;
	rfigures.ATLASy = vmuATLASZZ1blh;
	rfigures.COMBx = Cgg;
	rfigures.COMBy = muZZcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = CZZ;
	rfigures.CMSy = vmuCMSZZ1blh;
	rfigures.ATLASx = CZZ;
	rfigures.ATLASy = vmuATLASZZ1blh;
	rfigures.COMBx = CZZ;
	rfigures.COMBy = muZZcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//mu vs mu plots
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSyy1blh;
	rfigures.CMSy = vmuCMSZZ1blh;
	rfigures.ATLASx = vmuATLASyy1blh;
	rfigures.ATLASy = vmuATLASZZ1blh;
	rfigures.COMBx = muyycombined;
	rfigures.COMBy = muZZcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSyy1blh;
	rfigures.CMSy = vmuCMSWW1blh;
	rfigures.ATLASx = vmuATLASyy1blh;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = muyycombined;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSyy1blh;
	rfigures.CMSy = vmuCMSbb1blh;
	rfigures.ATLASx = vmuATLASyy1blh;
	rfigures.ATLASy = vmuATLASbb1blh;
	rfigures.COMBx = muyycombined;
	rfigures.COMBy = mubbcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSyy1blh;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = vmuATLASyy1blh;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = muyycombined;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSyy1blh;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = vmuATLASyy1blh;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = muyycombined;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSZZ1blh;
	rfigures.CMSy = vmuCMSWW1blh;
	rfigures.ATLASx = vmuATLASZZ1blh;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = muZZcombined;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSZZ1blh;
	rfigures.CMSy = vmuCMSbb1blh;
	rfigures.ATLASx = vmuATLASZZ1blh;
	rfigures.ATLASy = vmuATLASbb1blh;
	rfigures.COMBx = muZZcombined;
	rfigures.COMBy = mubbcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSZZ1blh;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = vmuATLASZZ1blh;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = muZZcombined;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSZZ1blh;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = vmuATLASZZ1blh;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = muZZcombined;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSWW1blh;
	rfigures.CMSy = vmuCMSbb1blh;
	rfigures.ATLASx = vmuATLASWW1blh;
	rfigures.ATLASy = vmuATLASbb1blh;
	rfigures.COMBx = muWWcombined;
	rfigures.COMBy = mubbcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSWW1blh;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = vmuATLASWW1blh;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = muWWcombined;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSWW1blh;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = vmuATLASWW1blh;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = muWWcombined;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSbb1blh;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = vmuATLASbb1blh;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = mubbcombined;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMSbb1blh;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = vmuATLASbb1blh;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = mubbcombined;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = vmuCMStata1blh;
	rfigures.CMSy = vmuCMSZy1blh;
	rfigures.ATLASx = vmuATLAStata1blh;
	rfigures.ATLASy = vmuATLASZy1blh;
	rfigures.COMBx = mutatacombined;
	rfigures.COMBy = muZycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = vmuCMSZZ1blh;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = vmuATLASZZ1blh;
	rfigures.COMBx = tB;
	rfigures.COMBy = muZZcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy =vmuCMSWW1blh;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = vmuATLASWW1blh;
	rfigures.COMBx = tB;
	rfigures.COMBy = muWWcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy =vmuCMSyy1blh;
	rfigures.ATLASx = tB;
	rfigures.ATLASy =vmuATLASyy1blh;
	rfigures.COMBx = tB;
	rfigures.COMBy = muyycombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = vmuCMSbb1blh;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = vmuATLASbb1blh;
	rfigures.COMBx = tB;
	rfigures.COMBy = mubbcombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = vmuCMStata1blh;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = vmuATLAStata1blh;
	rfigures.COMBx = tB;
	rfigures.COMBy = mutatacombined;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	//This section produces plots in which the colouration is due to different values of f, rather than the chisq.
	
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	for ( int i=0; i<fl.size(); i++) {
		if ( fl.at(i) < 1000 ) {
			CMSlt1++;
			ATLASlt1++;
			lt1++;
		}
		else if ( fl.at(i) < 1500 ) {
			CMSlt2++;
			ATLASlt2++;
			lt2++;
		}
		else if ( fl.at(i) < 2000 ) {
			CMSlt3++;
			ATLASlt3++;
			lt3++;
		}
		else {
			CMSgt3++;
			ATLASgt3++;
			gt3++;
		}
	}
	
	plot ffigures;
	ffigures.ratiotest = b_ratiotest;
	ffigures.bool_chisq = false;
	ffigures.bool_f = true;
	if ( degen ) {
		ffigures.bool_deg = true;
	}
	ffigures.chisqCMS = fl;
	ffigures.chisqATLAS = fl;
	ffigures.chisq = fl;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//Plots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//Plots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//parameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ffvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./fnorm/ffvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//mu vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/fmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./fnorm/fmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./fdeg/ftbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./fnorm/ftbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	
	
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	//#########################################################################
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	for ( int i=0; i<fl.size(); i++) {
		if ( !degen ) {
			if ( tB.at(i) < 1.25 ) {
				CMSlt1++;
				ATLASlt1++;
				lt1++;
			}
			else if ( tB.at(i) < 2.5 ) {
				CMSlt2++;
				ATLASlt2++;
				lt2++;
			}
			else if ( tB.at(i) < 3.75 ) {
				CMSlt3++;
				ATLASlt3++;
				lt3++;
			}
			else {
				CMSgt3++;
				ATLASgt3++;
				gt3++;
			}
		}
		else {
			if ( tB.at(i) < 1.01 ) {
				CMSlt1++;
				ATLASlt1++;
				lt1++;
			}
			else if ( tB.at(i) < 1.05 ) {
				CMSlt2++;
				ATLASlt2++;
				lt2++;
			}
			else if ( tB.at(i) < 1.10 ) {
				CMSlt3++;
				ATLASlt3++;
				lt3++;
			}
			else {
				CMSgt3++;
				ATLASgt3++;
				gt3++;
			}
		}
		
	}
	
	//plot ffigures;
	ffigures.bool_f = false;
	ffigures.bool_tB = true;
	ffigures.chisqCMS = tB;
	ffigures.chisqATLAS = tB;
	ffigures.chisq = tB;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//Plots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//Plots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//parameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./tnorm/tfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//mu vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/tmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./tnorm/tmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./tdeg/ttbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./tnorm/ttbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	for ( int i=0; i<fl.size(); i++) {
		if ( !degen ) {
			if ( deltam.at(i) < 25 ) {
				CMSlt1++;
				ATLASlt1++;
				lt1++;
			}
			else if ( deltam.at(i) < 100 ) {
				CMSlt2++;
				ATLASlt2++;
				lt2++;
			}
			else if ( deltam.at(i) < 200 ) {
				CMSlt3++;
				ATLASlt3++;
				lt3++;
			}
			else {
				CMSgt3++;
				ATLASgt3++;
				gt3++;
			}
		}
		else {
			if ( deltam.at(i) < 0.5 ) {
				CMSlt1++;
				ATLASlt1++;
				lt1++;
			}
			else if ( deltam.at(i) < 1 ) {
				CMSlt2++;
				ATLASlt2++;
				lt2++;
			}
			else if ( deltam.at(i) < 1.5 ) {
				CMSlt3++;
				ATLASlt3++;
				lt3++;
			}
			else {
				CMSgt3++;
				ATLASgt3++;
				gt3++;
			}
		}
		
	}
	
	//plot ffigures;
	ffigures.bool_tB = false;
	ffigures.bool_ma0 = true;
	ffigures.chisqCMS = deltam;
	ffigures.chisqATLAS = deltam;
	ffigures.chisq = deltam;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//Plots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//Plots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//parameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./mnorm/mfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/msyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/msyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/msggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/msggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//mu vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./mdeg/mtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./mnorm/mtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	
	//plot ffigures;
	//ffigures.bool_chisq = false;
	//ffigures.bool_deg = degen;
	
	
	vector<double> primary;
	
	for ( int i=0; i<fl.size(); i++) {
		if ( hmuATLASyy1blh.at(i) > AmuATLASyy1blh.at(i) && hmuATLASyy1blh.at(i) > HmuATLASyy1blh.at(i) ) {
			primary.push_back(1.0);
		}
		else if ( AmuATLASyy1blh.at(i) > hmuATLASyy1blh.at(i) && AmuATLASyy1blh.at(i) > HmuATLASyy1blh.at(i) ) {
			primary.push_back(3.0);
		}
		else if ( HmuATLASyy1blh.at(i) > AmuATLASyy1blh.at(i) && HmuATLASyy1blh.at(i) > hmuATLASyy1blh.at(i) ) {
			primary.push_back(2.0);
		}
		else {
			primary.push_back(4.0);
		}
	}
	
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	
	
	
	for ( int i=0; i<fl.size(); i++) {
		if ( primary.at(i) == 1.0  ) {
			CMSlt1++;
			ATLASlt1++;
			lt1++;
		}
		else if ( primary.at(i) == 2.0 ) {
			CMSlt2++;
			ATLASlt2++;
			lt2++;
		}
		else if ( primary.at(i) == 3.0 ) {
			CMSlt3++;
			ATLASlt3++;
			lt3++;
		}
		else {
			CMSgt3++;
			ATLASgt3++;
			gt3++;
		}
	}
	
	//plot ffigures;
	ffigures.bool_ma0 = false;
	ffigures.bool_prim = true;
	ffigures.chisqCMS = primary;
	ffigures.chisqATLAS = primary;
	ffigures.chisq = primary;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//Plots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//Plots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//parameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./pnorm/pfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pdsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/psyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/psyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/psggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/psggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//pu vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./pnorm/pmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLAStata1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLAStata1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLAStata1blh;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./pnorm/pmhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/pmhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./pnorm/pmhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,muyymin,1.3);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = hmuATLASyy1blh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = hmuATLASyy1blh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = hmuATLASyy1blh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/phmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}
	else {
		ffigures.generate_plot("./pnorm/phmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./pdeg/ptbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./pnorm/ptbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	
	
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	
	
	
	
	
	
	
	
	
	
	
	vector<double> tbdm;
	
	
	for ( int i=0; i<fl.size(); i++) {
		if ( !degen ) {
			if ( tB.at(i) < 1.25 ) {
				tbdm.push_back(0.0);
			}
			else if ( tB.at(i) < 2.5 ) {
				tbdm.push_back(1.0);
			}
			else if ( tB.at(i) < 3.75 ) {
				tbdm.push_back(2.0);
			}
			else {
				tbdm.push_back(3.0);
			}
		}
		else {
			if ( tB.at(i) < 1.01 ) {
				tbdm.push_back(0.0);
			}
			else if ( tB.at(i) < 1.05 ) {
				tbdm.push_back(1.0);
			}
			else if ( tB.at(i) < 1.10 ) {
				tbdm.push_back(2.0);
			}
			else {
				tbdm.push_back(3.0);
			}
		}
		
	}
	
	for ( int i=0; i<fl.size(); i++) {
		if ( !degen ) {
			if ( deltam.at(i) < 25 ) {
				tbdm.at(i) += 0.0;
			}
			else if ( deltam.at(i) < 100 ) {
				tbdm.at(i) += 4.0;
			}
			else if ( deltam.at(i) < 200 ) {
				tbdm.at(i) += 8.0;
			}
			else {
				tbdm.at(i) += 12.0;
			}
		}
		else {
			if ( deltam.at(i) < 0.5 ) {
				tbdm.at(i) += 0.0;
			}
			else if ( deltam.at(i) < 1 ) {
				tbdm.at(i) += 4.0;
			}
			else if ( deltam.at(i) < 1.5 ) {
				tbdm.at(i) += 8.0;
			}
			else {
				tbdm.at(i) += 12.0;
			}
		}
		
	}
	
	
	//plot ffigures;
	ffigures.bool_prim = false;
	ffigures.bool_double = true;
	ffigures.chisqCMS = tbdm;
	ffigures.chisqATLAS = tbdm;
	ffigures.chisq = tbdm;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	
	
	//dlots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//dlots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//dlots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//dlots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	
	
	//dlots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//darameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./dnorm/dfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddeltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/ddsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/ddsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dsyyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dsggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dczzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dczzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dczzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dcggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dczzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//du vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmuWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dmutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLAStata1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLAStata1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLAStata1blh;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./dnorm/dmhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dmhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./dnorm/dmhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,muyymin,1.3);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = hmuATLASyy1blh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = hmuATLASyy1blh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = hmuATLASyy1blh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dhmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}
	else {
		ffigures.generate_plot("./dnorm/dhmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./ddeg/dtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./dnorm/dtbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	
	
	
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	//#########################################################################################################
	
	
	
	
	
	
	vector<double> prim2;
	
	for ( int i=0; i<fl.size(); i++) {
		if ( hmuATLASyy1blh.at(i) > AmuATLASyy1blh.at(i) && hmuATLASyy1blh.at(i) > HmuATLASyy1blh.at(i) ) {
			prim2.push_back(0.0);
		}
		else if ( AmuATLASyy1blh.at(i) > hmuATLASyy1blh.at(i) && AmuATLASyy1blh.at(i) > HmuATLASyy1blh.at(i) ) {
			prim2.push_back(2.0);
		}
		else if ( HmuATLASyy1blh.at(i) > AmuATLASyy1blh.at(i) && HmuATLASyy1blh.at(i) > hmuATLASyy1blh.at(i) ) {
			prim2.push_back(1.0);
		}
		
		if ( hmuATLAStata1blh.at(i) > HmuATLAStata1blh.at(i) && hmuATLAStata1blh.at(i) > AmuATLAStata1blh.at(i) ) {
			prim2.at(i) += 0.0;
		}
		else if ( HmuATLAStata1blh.at(i) > hmuATLAStata1blh.at(i) && HmuATLAStata1blh.at(i) > AmuATLAStata1blh.at(i) ) {
			prim2.at(i) += 3.0;
		}
		else if ( AmuATLAStata1blh.at(i) > hmuATLAStata1blh.at(i) && AmuATLAStata1blh.at(i) > HmuATLAStata1blh.at(i) ) {
			prim2.at(i) += 6.0;
		}
	}
	
	
	
	CMSlt1=0;
	CMSlt2=0;
	CMSlt3=0;
	CMSgt3=0;
	ATLASlt1=0;
	ATLASlt2=0;
	ATLASlt3=0;
	ATLASgt3=0;
	lt1=0;
	lt2=0;
	lt3=0;
	gt3=0;
	
	cout << "Test1" << endl;
	
	
	//plot ffigures;
	ffigures.bool_double = false;
	ffigures.bool_prim2 = true;
	ffigures.chisqCMS = prim2;
	ffigures.chisqATLAS = prim2;
	ffigures.chisq = prim2;
	ffigures.minCMS = minCMS;
	ffigures.minATLAS = minATLAS;
	ffigures.min = min;
	ffigures.CMSlt1=CMSlt1;
	ffigures.CMSlt2=CMSlt2;
	ffigures.CMSlt3=CMSlt3;
	ffigures.CMSgt3=CMSgt3;
	ffigures.ATLASlt1=ATLASlt1;
	ffigures.ATLASlt2=ATLASlt2;
	ffigures.ATLASlt3=ATLASlt3;
	ffigures.ATLASgt3=ATLASgt3;
	ffigures.lt1=lt1;
	ffigures.lt2=lt2;
	ffigures.lt3=lt3;
	ffigures.gt3=gt3;
	
	cout << "Test2" << endl;

	
	//Plots involving CZZ
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = fl;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsczz.png","f (GeV)","r_{ZZ}",flmin,flmax,CZZmin,CZZmax);
	}
	
	cout << "Test3" << endl;

	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vsczz.png","mh0 (GeV)","r_{ZZ}",mh0min,mh0max,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvsczz.png","#Delta M_{h,A} (GeV)","r_{ZZ}",deltammin,deltammax,CZZmin,CZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CZZ;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CZZ;
	ffigures.COMBx = tB;
	ffigures.COMBy = CZZ;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsczz.png","tan#beta","r_{ZZ}",tBmin,tBmax,CZZmin,CZZmax);
	}
	
	cout << "Test4" << endl;

	
	//Plots involving CWW
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = fl;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvscww.png","f (GeV)","r_{WW}",flmin,flmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = mh0;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vscww.png","mh0 (GeV)","r_{WW}",mh0min,mh0max,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = deltam;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvscww.png","#Delta M_{h,A} (GeV)","r_{WW}",deltammin,deltammax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = tB;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvscww.png","tan#beta","r_{WW}",tBmin,tBmax,CWWmin,CWWmax);
	}
	
	
	
	//Plots involving Cyy
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvscyy.png","f (GeV)","r_{#gamma#gamma}",flmin,flmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vscyy.png","mh0 (GeV)","r_{#gamma#gamma}",mh0min,mh0max,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvscyy.png","#Delta M_{h,A} (GeV)","r_{#gamma#gamma}",deltammin,deltammax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvscyy.png","tan#beta","r_{#gamma#gamma}",tBmin,tBmax,Cyymin,Cyymax);
	}
	
	
	
	//Plots involving Cgg
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvscgg.png","f (GeV)","r_{hgg}",flmin,flmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vscgg.png","mh0 (GeV)","r_{hgg}",mh0min,mh0max,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvscgg.png","#Delta M_{h,A} (GeV)","r_{hgg}",deltammin,deltammax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvscgg.png","tan#beta","r_{hgg}",tBmin,tBmax,Cggmin,Cggmax);
	}
	
	
	
	//Plots involving Cbb
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = fl;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvscbb.png","f (GeV)","r_{bb}",flmin,flmax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vscbb.png","mh0 (GeV)","r_{bb}",mh0min,mh0max,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvscbb.png","#Delta M_{h,A} (GeV)","r_{bb}",deltammin,deltammax,Cbbmin,Cbbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Cbb;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Cbb;
	ffigures.COMBx = tB;
	ffigures.COMBy = Cbb;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvscbb.png","tan#beta","r_{bb}",tBmin,tBmax,Cbbmin,Cbbmax);
	}
	
	cout << "Test5" << endl;

	
	
	//Plots involving Ctata
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = fl;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsctata.png","f (GeV)","r_{#tau#tau}",flmin,flmax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mh0;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = mh0;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = mh0;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2mh0vsctata.png","mh0 (GeV)","r_{#tau#tau}",mh0min,mh0max,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = deltam;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvsctata.png","#Delta M_{h,A} (GeV)","r_{#tau#tau}",deltammin,deltammax,Ctatamin,Ctatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = Ctata;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = Ctata;
	ffigures.COMBx = tB;
	ffigures.COMBy = Ctata;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsctata.png","tan#beta","r_{#tau#tau}",tBmin,tBmax,Ctatamin,Ctatamax);
	}*/
	
	//f vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = fl;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	//Parameter vs parameter plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = tB;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = tB;
	ffigures.COMBx = fl;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = fl;
	ffigures.CMSy = deltam;
	ffigures.ATLASx = fl;
	ffigures.ATLASy = deltam;
	ffigures.COMBx = fl;
	ffigures.COMBy = deltam;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		ffigures.generate_plot("./2norm/2fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = deltam;
	ffigures.CMSy = tB;
	ffigures.ATLASx = deltam;
	ffigures.ATLASy = tB;
	ffigures.COMBx = deltam;
	ffigures.COMBy = tB;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2deltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		ffigures.generate_plot("./2norm/2deltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = DSgg;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = DSgg;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = DSgg;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2dsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2dsggvsdsyy.png","#Delta S_{gg}","#Delta S_{#gamma#gamma}",DSggmin,DSggmax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Syy;
	ffigures.CMSy = DSyy;
	ffigures.ATLASx = Syy;
	ffigures.ATLASy = DSyy;
	ffigures.COMBx = Syy;
	ffigures.COMBy = DSyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2syyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2syyvsdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Sgg;
	ffigures.CMSy = DSgg;
	ffigures.ATLASx = Sgg;
	ffigures.ATLASy = DSgg;
	ffigures.COMBx = Sgg;
	ffigures.COMBy = DSgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2sggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2sggvsdsgg.png","S_{gg}","#Delta S_{gg}",Sggmin,Sggmax,DSggmin,DSggmax);
	}
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2cyyvscgg.png","r_{#gamma#gamma}","r_{hgg}",Cyymin,Cyymax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cyy;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cyy;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cyy;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2czzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2czzvscyy.png","r_{ZZ}","r_{#gamma#gamma}",CZZmin,CZZmax,Cyymin,Cyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2czzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2czzvscgg.png","r_{ZZ}","r_{hgg}",CZZmin,CZZmax,Cggmin,Cggmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = CWW;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = CWW;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = CWW;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2czzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2czzvscww.png","r_{ZZ}","r_{WW}",CZZmin,CZZmax,CWWmin,CWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cbb;
	ffigures.CMSy = Cgg;
	ffigures.ATLASx = Cbb;
	ffigures.ATLASy = Cgg;
	ffigures.COMBx = Cbb;
	ffigures.COMBy = Cgg;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	else {
		ffigures.generate_plot("./2norm/2cbbvscgg.png","r_{bb}","r_{hgg}",Cbbmin,Cbbmax,Cggmin,Cggmax);
	}
	
	
	//C vs mu plots
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2cggvsmuyy.png","r_{hgg}","#mu_{#gamma#gamma}",Cggmin,Cggmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cyy;
	ffigures.CMSy = vmuCMSyy1blh;
	ffigures.ATLASx = Cyy;
	ffigures.ATLASy = vmuATLASyy1blh;
	ffigures.COMBx = Cyy;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2cyyvsmuyy.png","r_{#gamma#gamma}","#mu_{#gamma#gamma}",Cyymin,Cyymax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2cggvsmuww.png","r_{hgg}","#mu_{WW}",Cggmin,Cggmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CWW;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = CWW;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = CWW;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2cwwvsmuww.png","r_{WW}","#mu_{WW}",CWWmin,CWWmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = Cgg;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = Cgg;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = Cgg;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2cggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2cggvsmuzz.png","r_{hgg}","#mu_{ZZ}",Cggmin,Cggmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = CZZ;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = CZZ;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = CZZ;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2czzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2czzvsmuzz.png","r_{ZZ}","#mu_{ZZ}",CZZmin,CZZmax,muZZmin,muZZmax);
	}*/
	
	//Pu vs mu plots
	
	/*ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2muyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSyy1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASyy1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muyycombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2muyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSWW1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2muZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSZZ1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASZZ1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muZZcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2muZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2muWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2muWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSWW1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASWW1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = muWWcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2muWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2muWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2mubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMSbb1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLASbb1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mubbcombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2mubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = vmuCMStata1blh;
	ffigures.CMSy = vmuCMSZy1blh;
	ffigures.ATLASx = vmuATLAStata1blh;
	ffigures.ATLASy = vmuATLASZy1blh;
	ffigures.COMBx = mutatacombined;
	ffigures.COMBy = muZycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		ffigures.generate_plot("./2norm/2mutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLAStata1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLAStata1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLAStata1blh;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./2norm/2mhhvsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",mhhmin,mhhmax,mutatamin,1.3);
	}
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = mhh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = mhh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = mhh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2mhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,mutatamin,1.3);
	}
	else {
		ffigures.generate_plot("./2norm/2mhhvsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",mhhmin,mhhmax,muyymin,1.3);
	}
	
	
	
	
	ffigures.clear_plot_data();
	ffigures.CMSx = hmuATLASyy1blh;
	ffigures.CMSy = HmuATLASyy1blh;
	ffigures.ATLASx = hmuATLASyy1blh;
	ffigures.ATLASy = HmuATLASyy1blh;
	ffigures.COMBx = hmuATLASyy1blh;
	ffigures.COMBy = HmuATLASyy1blh;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2hmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}
	else {
		ffigures.generate_plot("./2norm/2hmuyyvsHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	}

	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSZZ1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASZZ1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muZZcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSWW1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASWW1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muWWcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy =vmuCMSyy1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy =vmuATLASyy1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = muyycombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMSbb1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLASbb1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mubbcombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsmubb.png","tan#beta","#Mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	ffigures.clear_plot_data();
	ffigures.CMSx = tB;
	ffigures.CMSy = vmuCMStata1blh;
	ffigures.ATLASx = tB;
	ffigures.ATLASy = vmuATLAStata1blh;
	ffigures.COMBx = tB;
	ffigures.COMBy = mutatacombined;
	if ( degen ) {
		ffigures.generate_plot("./2deg/2tbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		ffigures.generate_plot("./2norm/2tbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}

	
	return 0;
}

#endif