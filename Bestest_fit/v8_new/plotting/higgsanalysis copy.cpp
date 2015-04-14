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
		datafile.open("bestest_varyall_deg.txt");
	}
	else {
		datafile.open("bestest_varyall_ndeg.txt");
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
	double emuCMSyy1=sqrt(0.28*0.28+0.26*0.26)/sqrt(2.0);
	double emuCMSZZ1=sqrt(0.30*0.30+0.24*0.24)/sqrt(2.0), emuCMSZZ2=sqrt(0.84*0.84+0.57*0.57)/sqrt(2.0);
	double emuCMSWW1=0.21;
	double emuCMSbb1=sqrt(0.7*0.7+0.6*0.6)/sqrt(2.0);
	double emuCMStata1=0.4, emuCMStata2=0.6;
	
	double muATLASyy1=1.65;
	double muATLASZZ1=1.7, muATLASZZ2=1.5;
	double muATLASWW1=1.01, muATLASWW2=1.66, muATLASWW3=0.82;
	double muATLASbb1=-0.4;
	double muATLAStata1=0.7;
	double emuATLASyy1=sqrt(0.34*0.34+0.30*0.30)/sqrt(2.0);
	double emuATLASZZ1=sqrt(0.5*0.5+0.4*0.4)/sqrt(2.0), emuATLASZZ2=0.4;
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
	vector<double> tA;
	vector<double> tG;
	vector<double> t12;
	vector<double> t13;
	vector<double> sigfact;
	vector<double> kG;
	vector<double> kY;
	vector<double> kL;
//	vector<double> rDSgg;
//	vector<double> rSgg;
//	vector<double> rDSyy;
//	vector<double> rSyy;
//	vector<double> iDSgg;
//	vector<double> iSgg;
//	vector<double> iDSyy;
//	vector<double> iSyy;
//	vector<double> DSgg;
//	vector<double> Sgg;
//	vector<double> DSyy;
//	vector<double> Syy;
//	vector<double> Cyy;
//	vector<double> Cgg;
//	vector<double> CZZ;
//	vector<double> CWW;
//	vector<double> Cbb;
//	vector<double> Ctata;
	vector<double> mhc;
	vector<double> mphi0;
	vector<double> mphic;
	vector<double> meta0;
	vector<double> metac;
	vector<double> msigma;
	vector<double> mtau;
	vector<double> muq;
	vector<double> mdq;
	vector<double> msq;
	vector<double> mcq;
	vector<double> mbq;
	vector<double> mtq;
	vector<double> mTu;
	vector<double> mTd;
	vector<double> mTb5;
	vector<double> mTb2;
	vector<double> mT5;
	vector<double> mT6;
	vector<double> mWblh;
	vector<double> mZblh;
	vector<double> mwpblh;
	vector<double> mzpblh;
	vector<double> talpha;
	vector<double> sapb;
	vector<double> rDShgg;
	vector<double> iDShgg;
	vector<double> rShgg;
	vector<double> iShgg;
	vector<double> rDShyy;
	vector<double> iDShyy;
	vector<double> rShyy;
	vector<double> iShyy;
	vector<double> rDSHgg;
	vector<double> iDSHgg;
	vector<double> rSHgg;
	vector<double> iSHgg;
	vector<double> rDSHyy;
	vector<double> iDSHyy;
	vector<double> rSHyy;
	vector<double> iSHyy;
	vector<double> rDSAgg;
	vector<double> iDSAgg;
	vector<double> rSAgg;
	vector<double> iSAgg;
	vector<double> rDSAyy;
	vector<double> iDSAyy;
	vector<double> rSAyy;
	vector<double> iSAyy;
	vector<double> yv;
	vector<double> ratiowmass;
	vector<double> ywff;
	vector<double> yzff;
	vector<double> yhHH;
	vector<double> yhphiphi;
	vector<double> yhetaeta;
	vector<double> yHHH;
	vector<double> yHphiphi;
	vector<double> yHetaeta;
	vector<double> yAHH;
	vector<double> yAphiphi;
	vector<double> yAetaeta;
	vector<double> yhtata;
	vector<double> yHtata;
	vector<double> yAtata;
	vector<double> yhuu;
	vector<double> yhdd;
	vector<double> yhcc;
	vector<double> yhss;
	vector<double> yhbb;
	vector<double> yhtt;
	vector<double> yhTuTu;
	vector<double> yhTdTd;
	vector<double> yhTb5Tb5;
	vector<double> yhTb2Tb2;
	vector<double> yhT5T5;
	vector<double> yhT6T6;
	vector<double> yHuu;
	vector<double> yHdd;
	vector<double> yHcc;
	vector<double> yHss;
	vector<double> yHbb;
	vector<double> yHtt;
	vector<double> yHTuTu;
	vector<double> yHTdTd;
	vector<double> yHTb5Tb5;
	vector<double> yHTb2Tb2;
	vector<double> yHT5T5;
	vector<double> yHT6T6;
	vector<double> yAuu;
	vector<double> yAdd;
	vector<double> yAcc;
	vector<double> yAss;
	vector<double> yAbb;
	vector<double> yAtt;
	vector<double> yATuTu;
	vector<double> yATdTd;
	vector<double> yATb5Tb5;
	vector<double> yATb2Tb2;
	vector<double> yAT5T5;
	vector<double> yAT6T6;
	vector<double> yhww;
	vector<double> yHww;
	vector<double> yAww;
	vector<double> yhwpwp;
	vector<double> yHwpwp;
	vector<double> yAwpwp;
	vector<double> Chtt;
	vector<double> Chbb;
	vector<double> Chcc;
	vector<double> Chtata;
	vector<double> ChWW;
	vector<double> ChZZ;
	vector<double> Chyy;
	vector<double> ChyZ;
	vector<double> Chgg;
	vector<double> CHtt;
	vector<double> CHbb;
	vector<double> CHcc;
	vector<double> CHtata;
	vector<double> CHWW;
	vector<double> CHZZ;
	vector<double> CHyy;
	vector<double> CHyZ;
	vector<double> CHgg;
	vector<double> CAtt;
	vector<double> CAbb;
	vector<double> CAcc;
	vector<double> CAtata;
	vector<double> CAWW;
	vector<double> CAZZ;
	vector<double> CAyy;
	vector<double> CAyZ;
	vector<double> CAgg;
	vector<double> rShyyu;
	vector<double> iShyyu;
	vector<double> rShyyd;
	vector<double> iShyyd;
	vector<double> rShyys;
	vector<double> iShyys;
	vector<double> rShyyc;
	vector<double> iShyyc;
	vector<double> rShyyb;
	vector<double> iShyyb;
	vector<double> rShyyt;
	vector<double> iShyyt;
	vector<double> rShyyTu;
	vector<double> iShyyTu;
	vector<double> rShyyTd;
	vector<double> iShyyTd;
	vector<double> rShyyTb5;
	vector<double> iShyyTb5;
	vector<double> rShyyTb2;
	vector<double> iShyyTb2;
	vector<double> rShyyT5;
	vector<double> iShyyT5;
	vector<double> rShyyT6;
	vector<double> iShyyT6;
	vector<double> rShggu;
	vector<double> iShggu;
	vector<double> rShggd;
	vector<double> iShggd;
	vector<double> rShggs;
	vector<double> iShggs;
	vector<double> rShggc;
	vector<double> iShggc;
	vector<double> rShggb;
	vector<double> iShggb;
	vector<double> rShggt;
	vector<double> iShggt;
	vector<double> rShggTu;
	vector<double> iShggTu;
	vector<double> rShggTd;
	vector<double> iShggTd;
	vector<double> rShggTb5;
	vector<double> iShggTb5;
	vector<double> rShggTb2;
	vector<double> iShggTb2;
	vector<double> rShggT5;
	vector<double> iShggT5;
	vector<double> rShggT6;
	vector<double> iShggT6;
	vector<double> rSHyyu;
	vector<double> iSHyyu;
	vector<double> rSHyyd;
	vector<double> iSHyyd;
	vector<double> rSHyys;
	vector<double> iSHyys;
	vector<double> rSHyyc;
	vector<double> iSHyyc;
	vector<double> rSHyyb;
	vector<double> iSHyyb;
	vector<double> rSHyyt;
	vector<double> iSHyyt;
	vector<double> rSHyyTu;
	vector<double> iSHyyTu;
	vector<double> rSHyyTd;
	vector<double> iSHyyTd;
	vector<double> rSHyyTb5;
	vector<double> iSHyyTb5;
	vector<double> rSHyyTb2;
	vector<double> iSHyyTb2;
	vector<double> rSHyyT5;
	vector<double> iSHyyT5;
	vector<double> rSHyyT6;
	vector<double> iSHyyT6;
	vector<double> rSHggu;
	vector<double> iSHggu;
	vector<double> rSHggd;
	vector<double> iSHggd;
	vector<double> rSHggs;
	vector<double> iSHggs;
	vector<double> rSHggc;
	vector<double> iSHggc;
	vector<double> rSHggb;
	vector<double> iSHggb;
	vector<double> rSHggt;
	vector<double> iSHggt;
	vector<double> rSHggTu;
	vector<double> iSHggTu;
	vector<double> rSHggTd;
	vector<double> iSHggTd;
	vector<double> rSHggTb5;
	vector<double> iSHggTb5;
	vector<double> rSHggTb2;
	vector<double> iSHggTb2;
	vector<double> rSHggT5;
	vector<double> iSHggT5;
	vector<double> rSHggT6;
	vector<double> iSHggT6;
	vector<double> rSAyyu;
	vector<double> iSAyyu;
	vector<double> rSAyyd;
	vector<double> iSAyyd;
	vector<double> rSAyys;
	vector<double> iSAyys;
	vector<double> rSAyyc;
	vector<double> iSAyyc;
	vector<double> rSAyyb;
	vector<double> iSAyyb;
	vector<double> rSAyyt;
	vector<double> iSAyyt;
	vector<double> rSAyyTu;
	vector<double> iSAyyTu;
	vector<double> rSAyyTd;
	vector<double> iSAyyTd;
	vector<double> rSAyyTb5;
	vector<double> iSAyyTb5;
	vector<double> rSAyyTb2;
	vector<double> iSAyyTb2;
	vector<double> rSAyyT5;
	vector<double> iSAyyT5;
	vector<double> rSAyyT6;
	vector<double> iSAyyT6;
	vector<double> rSAggu;
	vector<double> iSAggu;
	vector<double> rSAggd;
	vector<double> iSAggd;
	vector<double> rSAggs;
	vector<double> iSAggs;
	vector<double> rSAggc;
	vector<double> iSAggc;
	vector<double> rSAggb;
	vector<double> iSAggb;
	vector<double> rSAggt;
	vector<double> iSAggt;
	vector<double> rSAggTu;
	vector<double> iSAggTu;
	vector<double> rSAggTd;
	vector<double> iSAggTd;
	vector<double> rSAggTb5;
	vector<double> iSAggTb5;
	vector<double> rSAggTb2;
	vector<double> iSAggTb2;
	vector<double> rSAggT5;
	vector<double> iSAggT5;
	vector<double> rSAggT6;
	vector<double> iSAggT6;
	vector<double> rShyyta;
	vector<double> iShyyta;
	vector<double> rSHyyta;
	vector<double> iSHyyta;
	vector<double> rSAyyta;
	vector<double> iSAyyta;
	vector<double> rShyyw;
	vector<double> iShyyw;
	vector<double> rShyywp;
	vector<double> iShyywp;
	vector<double> rSHyyw;
	vector<double> iSHyyw;
	vector<double> rSHyywp;
	vector<double> iSHyywp;
	vector<double> rShyyH;
	vector<double> iShyyH;
	vector<double> rShyyeta;
	vector<double> iShyyeta;
	vector<double> rShyyphi;
	vector<double> iShyyphi;
	vector<double> rSHyyH;
	vector<double> iSHyyH;
	vector<double> rSHyyeta;
	vector<double> iSHyyeta;
	vector<double> rSHyyphi;
	vector<double> iSHyyphi;
	
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
	double tAmin=9999, tAmax=0;
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
	
	vector<double> HtestW;
	vector<double> Htesty;
	vector<double> Atesty;
	
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
				mhc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mphi0.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mphic.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				meta0.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				metac.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				msigma.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mtau.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				muq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mdq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				msq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mcq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mbq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mtq.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mWblh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mZblh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mwpblh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				mzpblh.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				tB.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				talpha.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				sapb.push_back(atof(line.substr(start,end-start).c_str()));
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
				kL.push_back(atof(line.substr(start,end-start).c_str()));//38
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDShgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDShgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDShyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDShyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSHgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSHgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSHyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSHyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSAgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSAgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rDSAyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iDSAyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyy.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yv.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ratiowmass.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ywff.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yzff.push_back(atof(line.substr(start,end-start).c_str()));
			}
			if ( line.find("yvalues",0) != string::npos ) {
				start = 8;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhHH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhphiphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhetaeta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHHH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHphiphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHetaeta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAHH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAphiphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAetaeta.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhtata.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHtata.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAtata.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhuu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHuu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAuu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhdd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHdd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAdd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhss.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHss.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAss.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhTuTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHTuTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yATuTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhTdTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHTdTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yATdTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhTb5Tb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHTb5Tb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yATb5Tb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhTb2Tb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHTb2Tb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yATb2Tb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhT5T5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHT5T5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAT5T5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhT6T6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHT6T6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAT6T6.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhww.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHww.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAww.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yhwpwp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yHwpwp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				yAwpwp.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chtata.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ChWW.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ChZZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				ChyZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				Chgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHtata.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHWW.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHZZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHyZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CHgg.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAtt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAbb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAcc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAtata.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAWW.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAZZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAyy.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAyZ.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				CAgg.push_back(atof(line.substr(start,end-start).c_str()));
				
				
				
				
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
			if ( line.find("loopfactsh",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyw.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyw.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyywp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyywp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShyyeta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShyyeta.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rShggT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iShggT6.push_back(atof(line.substr(start,end-start).c_str()));
				
			}
			if ( line.find("loopfactsH",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyw.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyw.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyywp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyywp.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyH.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyphi.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHyyeta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHyyeta.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSHggT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSHggT6.push_back(atof(line.substr(start,end-start).c_str()));
				
			}
			if ( line.find("loopfactsA",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyys.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyT6.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAyyta.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAyyta.push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggc.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggs.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggt.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggb.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggTu.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggTd.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggTb5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggTb2.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggT5.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				rSAggT6.push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				iSAggT6.push_back(atof(line.substr(start,end-start).c_str()));
				
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
				Atesty.push_back(atof(line.substr(start,end-start).c_str()));
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
				Atesty.push_back(0);
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
				Htesty.push_back(atof(line.substr(start,end-start).c_str()));
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
				HtestW.push_back(atof(line.substr(start,end-start).c_str()));
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
				HtestW.push_back(0);
				Htesty.push_back(0);
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
	
	cout << "Vector size Htestw: " << HtestW.size() << endl;
	cout << "Vector size Htesty: " << Htesty.size() << endl;
	cout << "Vector size Atesty: " << Atesty.size() << endl;
	
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
		
		/*		if ( chisqATLAS.at(i) < 16 ) {
		 cout << "ATLAS ";
		 cout << fl.at(i) << " ";
		 cout << fh.at(i) << " ";
		 cout << mh0.at(i) << " ";
		 cout << ma0.at(i) << " ";
		 cout << mhh.at(i) << " ";
		 cout << tB.at(i) << " ";
		 cout << tG.at(i) << " ";
		 cout << t12.at(i) << " ";
		 cout << t13.at(i) << " ";
		 cout << sigfact.at(i) << " ";
		 cout << kG.at(i) << " ";
		 cout << kY.at(i) << " ";
		 cout << kL.at(i) << " ";
		 cout << Cyy.at(i) << " ";
		 cout << Cgg.at(i) << " ";
		 cout << CWW.at(i) << " ";
		 cout << CZZ.at(i) << " ";
		 cout << Cbb.at(i) << " ";
		 cout << Ctata.at(i) << " ";
		 cout << hmuATLASyy1blh.at(i) << " ";
		 cout << AmuATLASyy1blh.at(i) << " ";
		 cout << HmuATLASyy1blh.at(i) << " ";
		 cout << hmuATLAStata1blh.at(i) << " ";
		 cout << AmuATLAStata1blh.at(i) << " ";
		 cout << HmuATLAStata1blh.at(i) << " ";
		 cout << hmuATLASWW1blh.at(i) << " ";
		 cout << HmuATLASWW1blh.at(i) << " ";
		 cout << endl;
		 }
		 
		 if ( chisqCMS.at(i) < 3.97 ) {
		 cout << "CMS ";
		 cout << fl.at(i) << " ";
		 cout << fh.at(i) << " ";
		 cout << mh0.at(i) << " ";
		 cout << ma0.at(i) << " ";
		 cout << mhh.at(i) << " ";
		 cout << tB.at(i) << " ";
		 cout << tG.at(i) << " ";
		 cout << t12.at(i) << " ";
		 cout << t13.at(i) << " ";
		 cout << sigfact.at(i) << " ";
		 cout << kG.at(i) << " ";
		 cout << kY.at(i) << " ";
		 cout << kL.at(i) << " ";
		 cout << Cyy.at(i) << " ";
		 cout << Cgg.at(i) << " ";
		 cout << CWW.at(i) << " ";
		 cout << CZZ.at(i) << " ";
		 cout << Cbb.at(i) << " ";
		 cout << Ctata.at(i) << " ";
		 cout << hmuATLASyy1blh.at(i) << " ";
		 cout << AmuATLASyy1blh.at(i) << " ";
		 cout << HmuATLASyy1blh.at(i) << " ";
		 cout << hmuATLAStata1blh.at(i) << " ";
		 cout << AmuATLAStata1blh.at(i) << " ";
		 cout << HmuATLAStata1blh.at(i) << " ";
		 cout << hmuATLASWW1blh.at(i) << " ";
		 cout << HmuATLASWW1blh.at(i) << " ";
		 cout << endl;
		 }*/
		
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
	
	cout << "The chisq of the CMS SM is: " << pow((1.0-muCMSyy1)/emuCMSyy1,2.0)+pow((1.0-muCMSZZ1)/emuCMSZZ1,2.0)+pow((1.0-muCMSWW1)/emuCMSWW1,2.0)+pow((1.0-muCMStata2)/emuCMStata2,2.0)+pow((1.0-muCMSZZ2)/emuCMSZZ2,2.0)+pow((1.0-muCMSbb1)/emuCMSbb1,2.0)+pow((1.0-muCMStata1)/emuCMStata1,2.0) << endl;
	cout << muCMSyy1 << " " << emuCMSyy1 << " " << pow((1.0-muCMSyy1)/emuCMSyy1,2.0) << endl;
	cout << muCMSZZ1 << " " << emuCMSZZ1 << " " << pow((1.0-muCMSZZ1)/emuCMSZZ1,2.0) << endl;
	cout << muCMSZZ2 << " " << emuCMSZZ2 << " " << pow((1.0-muCMSZZ2)/emuCMSZZ2,2.0) << endl;
	cout << muCMSWW1 << " " << emuCMSWW1 << " " << pow((1.0-muCMSWW1)/emuCMSWW1,2.0) << endl;
	cout << muCMSbb1 << " " << emuCMSbb1 << " " << pow((1.0-muCMSbb1)/emuCMSbb1,2.0) << endl;
	cout << muCMStata1 << " " << emuCMStata1 << " " << pow((1.0-muCMStata1)/emuCMStata1,2.0) << endl;
	cout << muCMStata2 << " " << emuCMStata2 << " " << pow((1.0-muCMStata2)/emuCMStata2,2.0) << endl;
	cout << "The chisq of the ATLAS SM is: " << pow((1.0-muATLASyy1)/emuATLASyy1,2.0)+pow((1.0-muATLASZZ1)/emuATLASZZ1,2.0)/*+pow((1.0-muATLASZZ2)/emuATLASZZ2,2.0)*/+pow((1.0-muATLASWW1)/emuATLASWW1,2.0)+pow((1.0-muATLASWW2)/emuATLASWW2,2.0)+pow((1.0-muATLASWW3)/emuATLASWW3,2.0)+pow((1.0-muATLASbb1)/emuATLASbb1,2.0)+pow((1.0-muATLAStata1)/emuATLAStata1,2.0) << endl;
	cout << muATLASyy1 << " " << emuATLASyy1 << " " << pow((1.0-muATLASyy1)/emuATLASyy1,2.0) << endl;
	cout << muATLASZZ1 << " " << emuATLASZZ1 << " " << pow((1.0-muATLASZZ1)/emuATLASZZ1,2.0) << endl;
	//cout << muATLASZZ2 << " " << emuATLASZZ2 << " " << pow((1.0-muATLASZZ2)/emuATLASZZ2,2.0) << endl;
	cout << muATLASWW1 << " " << emuATLASWW1 << " " << pow((1.0-muATLASWW1)/emuATLASWW1,2.0) << endl;
	cout << muATLASWW2 << " " << emuATLASWW2 << " " << pow((1.0-muATLASWW2)/emuATLASWW2,2.0) << endl;
	cout << muATLASWW3 << " " << emuATLASWW3 << " " << pow((1.0-muATLASWW3)/emuATLASWW3,2.0) << endl;
	cout << muATLASbb1 << " " << emuATLASbb1 << " " << pow((1.0-muATLASbb1)/emuATLASbb1,2.0) << endl;
	cout << muATLAStata1 << " " << emuATLAStata1 << " " << pow((1.0-muATLAStata1)/emuATLAStata1,2.0) << endl;
	
	
	double mhwwlimit,mhyylimit,mayylimit;
	vector<double> b_ratiotest;
	
	double vsm, sg,sB,cB,vblh,lambda0,Bmu,ta,alphaem,MZ,MW,MWp,MWpconstraint,gA,gB,sw,g;
	
	
	vsm = 1.0/sqrt(sqrt(2)*0.0000116637);
	alphaem=1/127.9;
	MZ=91.1876;
	
	for ( int i=0; i<fl.size(); i++ ) {
		vblh = vsm*sqrt(1. + (1./6. + (5./4. - 2.*sg*sg)*fl.at(i)*fl.at(i)/(fl.at(i)*fl.at(i) + fh.at(i)*fh.at(i)))*vsm*vsm/(fl.at(i)*fl.at(i)));
		sw = sqrt(0.5 - 0.5*sqrt(1.0 - (4.0*3.14159*alphaem*vblh*vblh )/(MZ*MZ)));
		MW = sqrt(1.0 - sw*sw)*MZ;
		sg = tG.at(i)/sqrt(1.0+tG.at(i)*tG.at(i));
		g = sqrt(4.0*3.14159*alphaem)/sw;
		gA = g/sqrt(1.0-sg*sg);
		gB = g/sg;
		MWp=sqrt(0.25* (gA*gA + gB*gB)*(fl.at(i)*fl.at(i) + fh.at(i)*fh.at(i)) - MW*MW);
		MWpconstraint=14.182931524204026 - 23.378448106386056*sg - 0.042544864480156076*sg*sg + 17.296487025431006*sg*sg*sg*sg;
		
		
		if ( mhh.at(i) >= 128.2 && mhh.at(i) < 150 ) {
			mhwwlimit=pow(10.0,-1.38188 + 85958./(mhh.at(i)*mhh.at(i)) - 1219.03/mhh.at(i) + 0.353545*sqrt(mhh.at(i)) + 0.0523065*mhh.at(i) - 0.000258021*mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 150 && mhh.at(i) < 170 ) {
			mhwwlimit=pow(10.0,-2.70814 + 104083./(mhh.at(i)*mhh.at(i)) - 147.355/mhh.at(i) - 0.186522*sqrt(mhh.at(i)) - 0.00739727*mhh.at(i) + 0.000139588*mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 170 && mhh.at(i) < 200 ) {
			mhwwlimit=pow(10.0,19.2623 - 1746340/(mhh.at(i)*mhh.at(i)) + 18418.7/mhh.at(i) - 2.90129*sqrt(mhh.at(i)) - 0.393658*mhh.at(i) + 0.00131283*mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 200 && mhh.at(i) < 300 ) {
			mhwwlimit=pow(10.0,0.0135371 + 28937.7/(mhh.at(i)*mhh.at(i)) + 73.4727/mhh.at(i) - 0.00883728*sqrt(mhh.at(i)) - 0.00112041*mhh.at(i) - 0.0000076754*mhh.at(i)*mhh.at(i));
		}
		else if ( mhh.at(i) >= 300 ) {
			mhwwlimit=pow(10.0,269.072 + 1569510/(mhh.at(i)*mhh.at(i)) - 25921.9/mhh.at(i) - 17.797*sqrt(mhh.at(i)) + 0.376611*mhh.at(i) - 0.000060398*mhh.at(i)*mhh.at(i));
		}
		else {
			mhwwlimit = 100;
		}
		
		if ( mhh.at(i) >=128.2 && mhh.at(i) < 130 ) {
			mhyylimit=0.999769 - 42763.9/(mhh.at(i)*mhh.at(i)) + 718.678/mhh.at(i) - 0.222259*sqrt(mhh.at(i)) - 0.0358014*mhh.at(i) + 0.000190164*mhh.at(i)*mhh.at(i);
		}
		else if ( mhh.at(i) >= 130 && mhh.at(i) < 150) {
			mhyylimit=-5.01711 - 30608.8/(mhh.at(i)*mhh.at(i)) + 716.536/mhh.at(i) - 0.191598*sqrt(mhh.at(i)) + 0.0315056*mhh.at(i) - 0.0000338228*mhh.at(i)*mhh.at(i);
		}
		else if	( mhh.at(i) >= 150 ) {
			mhyylimit=0.017654;
		}
		else {
			mhyylimit=100;
		}
		
		if ( ma0.at(i) >=128.2 && ma0.at(i) < 130 ) {
			mayylimit=0.999769 - 42763.9/(ma0.at(i)*ma0.at(i)) + 718.678/ma0.at(i) - 0.222259*sqrt(ma0.at(i)) - 0.0358014*ma0.at(i) + 0.000190164*ma0.at(i)*ma0.at(i);
		}
		else if ( ma0.at(i) >= 130 && ma0.at(i) < 150) {
			mayylimit=-5.01711 - 30608.8/(ma0.at(i)*ma0.at(i)) + 716.536/ma0.at(i) - 0.191598*sqrt(ma0.at(i)) + 0.0315056*ma0.at(i) - 0.0000338228*ma0.at(i)*ma0.at(i);
		}
		else if	( ma0.at(i) >= 150 ) {
			mayylimit=0.017654;
		}
		else {
			mayylimit=100;
		}
		
		
		b_ratiotest.push_back(true);
		if ( HtestW.at(i) > mhwwlimit || MWp < MWpconstraint || Htesty.at(i) > mhyylimit || Atesty.at(i) > mayylimit ) {
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
	
	cout << "Starting data of interest" << endl;
	ofstream data;
	data.open("dataofinterest2.txt");
	data << "fl,fh,tB,tG,t12,t13,sigfact,kG,kY,kL,tA,sin(a+b),yv,mWblh/mWSM,mh,mA,mH,,mH+,mphi,mphi+,meta,meta+,msigma,mTu,mTd,";
	data << "mTb5,mTb2,mT5,mT6,rDShgg,iDShgg,rShgg,iShgg,rDShyy,iDShyy,rShyy,iShyy,rDSHgg,iDSHgg,rSHgg,iSHgg,rDSHyy,iDSHyy,rSHyy,";
	data << "iSHyy,rDSAgg,iDSAgg,rSAgg,iSAgg,rDSAyy,iDSAyy,rSAyy,iSAyy,Cyy,Cgg,CWW,CZZ,Cbb,Ctata,muhyy,muAyy,muHyy,muhtata,muAtata,muHtata,muhWW,muHWW" << endl;
	for ( int i=0; i<fl.size(); i++) {
		sg = tG.at(i)/sqrt(1.0+tG.at(i)*tG.at(i));
		sB = tB.at(i)/sqrt(1.0+tB.at(i)*tB.at(i));
		cB = sqrt(1.0-sB*sB);
		vblh=vsm*sqrt(1. + (1./6. + (5./4. - 2.*sg*sg)*fl.at(i)*fl.at(i)/(fl.at(i)*fl.at(i) + fh.at(i)*fh.at(i)))*vsm*vsm/(fl.at(i)*fl.at(i)));
		lambda0 = (mh0.at(i)*mh0.at(i)/(vblh*vblh))*((mh0.at(i)*mh0.at(i) - ma0.at(i)*ma0.at(i))/(mh0.at(i)*mh0.at(i) - 2.0*sB*cB* ma0.at(i)*ma0.at(i)));
		Bmu =  (lambda0*vblh*vblh + ma0.at(i)*ma0.at(i))*(2.0*sB*cB)/2.0;
		ta = (1.0/(Bmu - lambda0 *vblh*vblh*2.0*sB*cB))*((Bmu*(cB/(2.0*sB) - sB/(2.0*cB)))+sqrt( Bmu*Bmu/(4.0*sB*sB*cB*cB)   -2.0*lambda0*Bmu*vblh*vblh*2.0*sB*cB  +lambda0*lambda0*vblh*vblh*vblh*vblh*4.0*sB*sB*cB*cB     ));
		tA.push_back(ta);
		if ( rchisqATLAS.at(i) < 1 && b_ratiotest.at(i) == true ) {
			data << fl.at(i) << ",";
			data << fh.at(i) << ",";
			data << tB.at(i) << ",";
			data << tG.at(i) << ",";
			data << t12.at(i) << ",";
			data << t13.at(i) << ",";
			data << sigfact.at(i) << ",";
			data << kG.at(i) << ",";
			data << kY.at(i) << ",";
			data << kL.at(i) << ",";

			data << tA.at(i) << ",";
			data << sapb.at(i) << ",";
			data << yv.at(i) << ",";
			data << ratiowmass.at(i) << ",";
			
			data << mh0.at(i) << ",";
			data << ma0.at(i) << ",";
			data << mhh.at(i) << ",";
			data << mhc.at(i) << ",";
			data << mphi0.at(i) << ",";
			data << mphic.at(i) << ",";
			data << meta0.at(i) << ",";
			data << metac.at(i) << ",";
			data << msigma.at(i) << ",";
			data << mTu.at(i) << ",";
			data << mTd.at(i) << ",";
			data << mTb5.at(i) << ",";
			data << mTb2.at(i) << ",";
			data << mT5.at(i) << ",";
			data << mT6.at(i) << ",";
			data << mwpblh.at(i) << ",";
			data << mzpblh.at(i) << ",";
			
			
			data << rDShgg.at(i) << ",";
			data << iDShgg.at(i) << ",";
			data << rShgg.at(i) << ",";
			data << iShgg.at(i) << ",";
			data << rDShyy.at(i) << ",";
			data << iDShyy.at(i) << ",";
			data << rShyy.at(i) << ",";
			data << iShyy.at(i) << ",";
			data << rDSHgg.at(i) << ",";
			data << iDSHgg.at(i) << ",";
			data << rSHgg.at(i) << ",";
			data << iSHgg.at(i) << ",";
			data << rDSHyy.at(i) << ",";
			data << iDSHyy.at(i) << ",";
			data << rSHyy.at(i) << ",";
			data << iSHyy.at(i) << ",";
			data << rDSAgg.at(i) << ",";
			data << iDSAgg.at(i) << ",";
			data << rSAgg.at(i) << ",";
			data << iSAgg.at(i) << ",";
			data << rDSAyy.at(i) << ",";
			data << iDSAyy.at(i) << ",";
			data << rSAyy.at(i) << ",";
			data << iSAyy.at(i) << ",";
			
			data << Chyy.at(i) << ",";
			data << Chgg.at(i) << ",";
			data << ChWW.at(i) << ",";
			data << ChZZ.at(i) << ",";
			data << Chbb.at(i) << ",";
			data << Chtata.at(i) << ",";
			data << hmuATLASyy1blh.at(i) << ",";
			data << AmuATLASyy1blh.at(i) << ",";
			data << HmuATLASyy1blh.at(i) << ",";
			data << hmuATLAStata1blh.at(i) << ",";
			data << AmuATLAStata1blh.at(i) << ",";
			data << HmuATLAStata1blh.at(i) << ",";
			data << hmuATLASWW1blh.at(i) << ",";
			data << HmuATLASWW1blh.at(i) << ",";
			data << endl;
		}
	}
	data.close();
	
	cout << "Done data of interest, starting min/max determinations." << endl;
	
	
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
		if ( tA.at(i) < tAmin ) {
			tAmin = tA.at(i);
		}
		else if (tA.at(i) > tAmax ) {
			tAmax = tA.at(i);
		}
		if ( tG.at(i) < tGmin ) {
			tGmin = tG.at(i);
		}
		else if (tG.at(i) > tGmax ) {
			tGmax = tG.at(i);
		}
		if ( rDShgg.at(i) < DSggmin ) {
			DSggmin = rDShgg.at(i);
		}
		else if (rDShgg.at(i) > DSggmax ) {
			DSggmax = rDShgg.at(i);
		}
		if ( rShgg.at(i) < Sggmin ) {
			Sggmin = rShgg.at(i);
		}
		else if (rShgg.at(i) > Sggmax ) {
			Sggmax = rShgg.at(i);
		}
		if ( rDShyy.at(i) < DSyymin ) {
			DSyymin = rDShyy.at(i);
		}
		else if (rDShyy.at(i) > DSyymax ) {
			DSyymax = rDShyy.at(i);
		}
		if ( rShyy.at(i) < Syymin ) {
			Syymin = rShyy.at(i);
		}
		else if (rShyy.at(i) > Syymax ) {
			Syymax = rShyy.at(i);
		}
		if ( Chyy.at(i) < Cyymin ) {
			Cyymin = Chyy.at(i);
		}
		else if (Chyy.at(i) > Cyymax ) {
			Cyymax = Chyy.at(i);
		}
		if ( Chgg.at(i) < Cggmin ) {
			Cggmin = Chgg.at(i);
		}
		else if (Chgg.at(i) > Cggmax ) {
			Cggmax = Chgg.at(i);
		}
		if ( ChWW.at(i) < CWWmin ) {
			CWWmin = ChWW.at(i);
		}
		else if (ChWW.at(i) > CWWmax ) {
			CWWmax = ChWW.at(i);
		}
		if ( ChZZ.at(i) < CZZmin ) {
			CZZmin = ChZZ.at(i);
		}
		else if (ChZZ.at(i) > CZZmax ) {
			CZZmax = ChZZ.at(i);
		}
		if ( Chbb.at(i) < Cbbmin ) {
			Cbbmin = Chbb.at(i);
		}
		else if (Chbb.at(i) > Cbbmax ) {
			Cbbmax = Chbb.at(i);
		}
		if ( Chtata.at(i) < Ctatamin ) {
			Ctatamin = Chtata.at(i);
		}
		else if (Chtata.at(i) > Ctatamax ) {
			Ctatamax = Chtata.at(i);
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
		muZymax = 2.5;
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
	
	
	
	for ( int i=0; i<fl.size(); i++) {
		if ( chisqCMS.at(i) < 2.8915 ) {
			CMSlt1++;
		}
		else if ( chisqCMS.at(i) < 3.8915 ) {
			CMSlt2++;
		}
		else if ( chisqCMS.at(i) < 6.8915 ) {
			CMSlt3++;
		}
		else {
			CMSgt3++;
		}
		if ( chisqATLAS.at(i) < 11.1553 ) {
			ATLASlt1++;
		}
		else if ( chisqATLAS.at(i)< 12.1553 ) {
			ATLASlt2++;
		}
		else if ( chisqATLAS.at(i) < 15.1553 ) {
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
	
	
	//Keep these lines the same. There should be no need to change them.
	plot figures;
	figures.ratiotest = b_ratiotest;
	figures.bool_delta = true;
	if ( degen ) {
		figures.bool_degen = true;
	}
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
	
/*
	
	
	//f vs mu plots
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = fl;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/fvsmuzz.png","f (GeV)","#mu_{ZZ}",flmin,flmax,muZZmin,muZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = fl;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/fvsmuww.png","f (GeV)","#mu_{WW}",flmin,flmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMSyy1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLASyy1blh;
	figures.COMBx = fl;
	figures.COMBy = muyycombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	else {
		figures.generate_plot("./norm/fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",flmin,flmax,muyymin,muyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = fl;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",flmin,flmax,muZymin,muZymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = fl;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/fvsmubb.png","f (GeV)","#mu_{bb}",flmin,flmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = fl;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = fl;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/fvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/fvsmutata.png","f (GeV)","#mu_{#tau#tau}",flmin,flmax,mutatamin,mutatamax);
	}
	
	
	
	//parameter vs parameter plots
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = tB;
	figures.ATLASx = fl;
	figures.ATLASy = tB;
	figures.COMBx = fl;
	figures.COMBy = tB;
	if ( degen ) {
		figures.generate_plot("./deg/fvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	else {
		figures.generate_plot("./norm/fvstB.png","f (GeV)","tan#beta",flmin,flmax,tBmin,tBmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = fl;
	figures.CMSy = deltam;
	figures.ATLASx = fl;
	figures.ATLASy = deltam;
	figures.COMBx = fl;
	figures.COMBy = deltam;
	if ( degen ) {
		figures.generate_plot("./deg/fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	else {
		figures.generate_plot("./norm/fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",flmin,flmax,deltammin,deltammax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = deltam;
	figures.CMSy = tB;
	figures.ATLASx = deltam;
	figures.ATLASy = tB;
	figures.COMBx = deltam;
	figures.COMBy = tB;
	if ( degen ) {
		figures.generate_plot("./deg/deltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	else {
		figures.generate_plot("./norm/deltamvstB.png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,tBmin,tBmax);
	}
	

	
	//mu vs mu plots
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/muyyvsmuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/muyyvsmuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/muyyvsmubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/muyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/muyyvsmutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSyy1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASyy1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muyycombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/muyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/muyyvsmuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSWW1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/muZZvsmuWW.png","#mu_{ZZ}","#mu_{WW}",muZZmin,muZZmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/muZZvsmubb.png","#mu_{ZZ}","#mu_{bb}",muZZmin,muZZmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/muZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/muZZvsmutata.png","#mu_{ZZ}","#mu_{#tau#tau}",muZZmin,muZZmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSZZ1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASZZ1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muZZcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/muZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/muZZvsmuZy.png","#mu_{ZZ}","#mu_{Z#gamma}",muZZmin,muZZmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/muWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/muWWvsmubb.png","#mu_{WW}","#mu_{bb}",muWWmin,muWWmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/muWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/muWWvsmutata.png","#mu_{WW}","#mu_{#tau#tau}",muWWmin,muWWmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSWW1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASWW1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = muWWcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/muWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/muWWvsmuZy.png","#mu_{WW}","#mu_{Z#gamma}",muWWmin,muWWmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSbb1blh;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = vmuATLASbb1blh;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = mubbcombined;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/mubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/mubbvsmutata.png","#mu_{bb}","#mu_{#tau#tau}",mubbmin,mubbmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMSbb1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLASbb1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = mubbcombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/mubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/mubbvsmuZy.png","#mu_{bb}","#mu_{Z#gamma}",mubbmin,mubbmax,muZymin,muZymax);
	}
	
	
	
	
	
	figures.clear_plot_data();
	figures.CMSx = vmuCMStata1blh;
	figures.CMSy = vmuCMSZy1blh;
	figures.ATLASx = vmuATLAStata1blh;
	figures.ATLASy = vmuATLASZy1blh;
	figures.COMBx = mutatacombined;
	figures.COMBy = muZycombined;
	if ( degen ) {
		figures.generate_plot("./deg/mutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	else {
		figures.generate_plot("./norm/mutatavsmuZy.png","#mu_{#tau#tau}","#mu_{Z#gamma}",mutatamin,mutatamax,muZymin,muZymax);
	}
	
	
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = tB;
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = tB;
	figures.COMBy = muZZcombined;
	if ( degen ) {
		figures.generate_plot("./deg/tbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/tbvsmuzz.png","tan#beta","#mu_{ZZ}",tBmin,tBmax,muZZmin,muZZmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy =vmuCMSWW1blh;
	figures.ATLASx = tB;
	figures.ATLASy = vmuATLASWW1blh;
	figures.COMBx = tB;
	figures.COMBy = muWWcombined;
	if ( degen ) {
		figures.generate_plot("./deg/tbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	else {
		figures.generate_plot("./norm/tbvsmuww.png","tan#beta","#mu_{WW}",tBmin,tBmax,muWWmin,muWWmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy =vmuCMSyy1blh;
	figures.ATLASx = tB;
	figures.ATLASy =vmuATLASyy1blh;
	figures.COMBx = tB;
	figures.COMBy = muyycombined;
	if ( degen ) {
		figures.generate_plot("./deg/tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	else {
		figures.generate_plot("./norm/tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",tBmin,tBmax,muyymin,muyymax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy = vmuCMSbb1blh;
	figures.ATLASx = tB;
	figures.ATLASy = vmuATLASbb1blh;
	figures.COMBx = tB;
	figures.COMBy = mubbcombined;
	if ( degen ) {
		figures.generate_plot("./deg/tbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	else {
		figures.generate_plot("./norm/tbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy = vmuCMStata1blh;
	figures.ATLASx = tB;
	figures.ATLASy = vmuATLAStata1blh;
	figures.COMBx = tB;
	figures.COMBy = mutatacombined;
	if ( degen ) {
		figures.generate_plot("./deg/tbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	else {
		figures.generate_plot("./norm/tbvsmutata.png","tan#beta","#mu_{#tau#tau}",tBmin,tBmax,mutatamin,mutatamax);
	}
	
	figures.clear_plot_data();
	figures.CMSx = tB;
	figures.CMSy = tA;
	figures.ATLASx = tB;
	figures.ATLASy = tA;
	figures.COMBx = tB;
	figures.COMBy = tA;
	if ( degen ) {
		figures.generate_plot("./rdeg/rtbvsta.png","tan(#beta)","tan(#alpha)",tBmin,tBmax,tAmin,tAmax);
	}
	else {
		figures.generate_plot("./rnorm/rtbvsta.png","tan(#beta)","tan(#alpha)",tBmin,tBmax,tAmin,tAmax);
	}
	
	
	
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
	
	
	rfigures.clear_plot_data();
	rfigures.CMSx = tB;
	rfigures.CMSy = tA;
	rfigures.ATLASx = tB;
	rfigures.ATLASy = tA;
	rfigures.COMBx = tB;
	rfigures.COMBy = tA;
	if ( degen ) {
		rfigures.generate_plot("./rdeg/rtbvsta.png","tan(#beta)","tan(#alpha)",tBmin,tBmax,tAmin,tAmax);
	}
	else {
		rfigures.generate_plot("./rnorm/rtbvsta.png","tan(#beta)","tan(#alpha)",tBmin,tBmax,tAmin,tAmax);
	}
	
	
	
	
	
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
		rfigures.generate_plot("./rnorm/rtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	
	
	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
	
	//mu vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./fnorm/ftbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
	
	//mu vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./tnorm/ttbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	
	
	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
	
	//mu vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./mnorm/mtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	
	
	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
		
	//pu vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./pnorm/ptbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	

	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
	
	//du vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./dnorm/dtbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	
	
	
	//f vs mu plots
	
	ffigures.clear_plot_data();
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
	
	
	//Pu vs mu plots
	
	ffigures.clear_plot_data();
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
		ffigures.generate_plot("./2norm/2tbvsmubb.png","tan#beta","#mu_{bb}",tBmin,tBmax,mubbmin,mubbmax);
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
	*/
	
	return 0;
}

#endif