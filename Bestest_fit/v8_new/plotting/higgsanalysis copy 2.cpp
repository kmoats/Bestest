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
	
	vector<double> param[318];
	double parammax[318];
	double parammin[318];
	
	//Param relations
	/*
	 param[0]<->fl;
	 param[1]<->fh;
	 param[2]<->mh0;
	 param[3]<->ma0;
	 param[4]<->mhh;
	 param[5]<->tB;
	 param[6]<->tA;
	 param[7]<->tG;
	 param[8]<->t12;
	 param[9]<->t13;
	 param[10]<->sigfact;
	 param[11]<->kG;
	 param[12]<->kY;
	 param[13]<->kL;
	 param[14]<->mhc;
	 param[15]<->mphi0;
	 param[16]<->mphic;
	 param[17]<->meta0;
	 param[18]<->metac;
	 param[19]<->msigma;
	 param[20]<->mtau;
	 param[21]<->muq;
	 param[22]<->mdq;
	 param[23]<->msq;
	 param[24]<->mcq;
	 param[25]<->mbq;
	 param[26]<->mtq;
	 param[27]<->mTu;
	 param[28]<->mTd;
	 param[29]<->mTb5;
	 param[30]<->mTb2;
	 param[31]<->mT5;
	 param[32]<->mT6;
	 param[33]<->mWblh;
	 param[34]<->mZblh;
	 param[35]<->mwpblh;
	 param[36]<->mzpblh;
	 param[37]<->talpha;
	 param[38]<->sapb;
	 param[39]<->rDShgg;
	 param[40]<->iDShgg;
	 param[41]<->rShgg;
	 param[42]<->iShgg;
	 param[43]<->rDShyy;
	 param[44]<->iDShyy;
	 param[45]<->rShyy;
	 param[46]<->iShyy;
	 param[47]<->rDSHgg;
	 param[48]<->iDSHgg;
	 param[49]<->rSHgg;
	 param[50]<->iSHgg;
	 param[51]<->rDSHyy;
	 param[52]<->iDSHyy;
	 param[53]<->rSHyy;
	 param[54]<->iSHyy;
	 param[55]<->rDSAgg;
	 param[56]<->iDSAgg;
	 param[57]<->rSAgg;
	 param[58]<->iSAgg;
	 param[59]<->rDSAyy;
	 param[60]<->iDSAyy;
	 param[61]<->rSAyy;
	 param[62]<->iSAyy;
	 param[63]<->yv;
	 param[64]<->ratiowmass;
	 param[65]<->ywff;
	 param[66]<->yzff;
	 param[67]<->yhHH;
	 param[68]<->yhphiphi;
	 param[69]<->yhetaeta;
	 param[70]<->yHHH;
	 param[71]<->yHphiphi;
	 param[72]<->yHetaeta;
	 param[73]<->yAHH;
	 param[74]<->yAphiphi;
	 param[75]<->yAetaeta;
	 param[76]<->yhtata;
	 param[77]<->yHtata;
	 param[78]<->yAtata;
	 param[79]<->yhuu;
	 param[80]<->yhdd;
	 param[81]<->yhcc;
	 param[82]<->yhss;
	 param[83]<->yhbb;
	 param[84]<->yhtt;
	 param[85]<->yhTuTu;
	 param[86]<->yhTdTd;
	 param[87]<->yhTb5Tb5;
	 param[88]<->yhTb2Tb2;
	 param[89]<->yhT5T5;
	 param[90]<->yhT6T6;
	 param[91]<->yHuu;
	 param[92]<->yHdd;
	 param[93]<->yHcc;
	 param[94]<->yHss;
	 param[95]<->yHbb;
	 param[96]<->yHtt;
	 param[97]<->yHTuTu;
	 param[98]<->yHTdTd;
	 param[99]<->yHTb5Tb5;
	 param[100]<->yHTb2Tb2;
	 param[101]<->yHT5T5;
	 param[102]<->yHT6T6;
	 param[103]<->yAuu;
	 param[104]<->yAdd;
	 param[105]<->yAcc;
	 param[106]<->yAss;
	 param[107]<->yAbb;
	 param[108]<->yAtt;
	 param[109]<->yATuTu;
	 param[110]<->yATdTd;
	 param[111]<->yATb5Tb5;
	 param[112]<->yATb2Tb2;
	 param[113]<->yAT5T5;
	 param[114]<->yAT6T6;
	 param[115]<->yhww;
	 param[116]<->yHww;
	 param[117]<->yAww;
	 param[118]<->yhwpwp;
	 param[119]<->yHwpwp;
	 param[120]<->yAwpwp;
	 param[121]<->Chtt;
	 param[122]<->Chbb;
	 param[123]<->Chcc;
	 param[124]<->Chtata;
	 param[125]<->ChWW;
	 param[126]<->ChZZ;
	 param[127]<->Chyy;
	 param[128]<->ChyZ;
	 param[129]<->Chgg;
	 param[130]<->CHtt;
	 param[131]<->CHbb;
	 param[132]<->CHcc;
	 param[133]<->CHtata;
	 param[134]<->CHWW;
	 param[135]<->CHZZ;
	 param[136]<->CHyy;
	 param[137]<->CHyZ;
	 param[138]<->CHgg;
	 param[139]<->CAtt;
	 param[140]<->CAbb;
	 param[141]<->CAcc;
	 param[142]<->CAtata;
	 param[143]<->CAWW;
	 param[144]<->CAZZ;
	 param[145]<->CAyy;
	 param[146]<->CAyZ;
	 param[147]<->CAgg;
	 param[148]<->rShyyu;
	 param[149]<->iShyyu;
	 param[150]<->rShyyd;
	 param[151]<->iShyyd;
	 param[152]<->rShyys;
	 param[153]<->iShyys;
	 param[154]<->rShyyc;
	 param[155]<->iShyyc;
	 param[156]<->rShyyb;
	 param[157]<->iShyyb;
	 param[158]<->rShyyt;
	 param[159]<->iShyyt;
	 param[160]<->rShyyTu;
	 param[161]<->iShyyTu;
	 param[162]<->rShyyTd;
	 param[163]<->iShyyTd;
	 param[164]<->rShyyTb5;
	 param[165]<->iShyyTb5;
	 param[166]<->rShyyTb2;
	 param[167]<->iShyyTb2;
	 param[168]<->rShyyT5;
	 param[169]<->iShyyT5;
	 param[170]<->rShyyT6;
	 param[171]<->iShyyT6;
	 param[172]<->rShggu;
	 param[173]<->iShggu;
	 param[174]<->rShggd;
	 param[175]<->iShggd;
	 param[176]<->rShggs;
	 param[177]<->iShggs;
	 param[178]<->rShggc;
	 param[179]<->iShggc;
	 param[180]<->rShggb;
	 param[181]<->iShggb;
	 param[182]<->rShggt;
	 param[183]<->iShggt;
	 param[184]<->rShggTu;
	 param[185]<->iShggTu;
	 param[186]<->rShggTd;
	 param[187]<->iShggTd;
	 param[188]<->rShggTb5;
	 param[189]<->iShggTb5;
	 param[190]<->rShggTb2;
	 param[191]<->iShggTb2;
	 param[192]<->rShggT5;
	 param[193]<->iShggT5;
	 param[194]<->rShggT6;
	 param[195]<->iShggT6;
	 param[196]<->rSHyyu;
	 param[197]<->iSHyyu;
	 param[198]<->rSHyyd;
	 param[199]<->iSHyyd;
	 param[200]<->rSHyys;
	 param[201]<->iSHyys;
	 param[202]<->rSHyyc;
	 param[203]<->iSHyyc;
	 param[204]<->rSHyyb;
	 param[205]<->iSHyyb;
	 param[206]<->rSHyyt;
	 param[207]<->iSHyyt;
	 param[208]<->rSHyyTu;
	 param[209]<->iSHyyTu;
	 param[210]<->rSHyyTd;
	 param[211]<->iSHyyTd;
	 param[212]<->rSHyyTb5;
	 param[213]<->iSHyyTb5;
	 param[214]<->rSHyyTb2;
	 param[215]<->iSHyyTb2;
	 param[216]<->rSHyyT5;
	 param[217]<->iSHyyT5;
	 param[218]<->rSHyyT6;
	 param[219]<->iSHyyT6;
	 param[220]<->rSHggu;
	 param[221]<->iSHggu;
	 param[222]<->rSHggd;
	 param[223]<->iSHggd;
	 param[224]<->rSHggs;
	 param[225]<->iSHggs;
	 param[226]<->rSHggc;
	 param[227]<->iSHggc;
	 param[228]<->rSHggb;
	 param[229]<->iSHggb;
	 param[230]<->rSHggt;
	 param[231]<->iSHggt;
	 param[232]<->rSHggTu;
	 param[233]<->iSHggTu;
	 param[234]<->rSHggTd;
	 param[235]<->iSHggTd;
	 param[236]<->rSHggTb5;
	 param[237]<->iSHggTb5;
	 param[238]<->rSHggTb2;
	 param[239]<->iSHggTb2;
	 param[240]<->rSHggT5;
	 param[241]<->iSHggT5;
	 param[242]<->rSHggT6;
	 param[243]<->iSHggT6;
	 param[244]<->rSAyyu;
	 param[245]<->iSAyyu;
	 param[246]<->rSAyyd;
	 param[247]<->iSAyyd;
	 param[248]<->rSAyys;
	 param[249]<->iSAyys;
	 param[250]<->rSAyyc;
	 param[251]<->iSAyyc;
	 param[252]<->rSAyyb;
	 param[253]<->iSAyyb;
	 param[254]<->rSAyyt;
	 param[255]<->iSAyyt;
	 param[256]<->rSAyyTu;
	 param[257]<->iSAyyTu;
	 param[258]<->rSAyyTd;
	 param[259]<->iSAyyTd;
	 param[260]<->rSAyyTb5;
	 param[261]<->iSAyyTb5;
	 param[262]<->rSAyyTb2;
	 param[263]<->iSAyyTb2;
	 param[264]<->rSAyyT5;
	 param[265]<->iSAyyT5;
	 param[266]<->rSAyyT6;
	 param[267]<->iSAyyT6;
	 param[268]<->rSAggu;
	 param[269]<->iSAggu;
	 param[270]<->rSAggd;
	 param[271]<->iSAggd;
	 param[272]<->rSAggs;
	 param[273]<->iSAggs;
	 param[274]<->rSAggc;
	 param[275]<->iSAggc;
	 param[276]<->rSAggb;
	 param[277]<->iSAggb;
	 param[278]<->rSAggt;
	 param[279]<->iSAggt;
	 param[280]<->rSAggTu;
	 param[281]<->iSAggTu;
	 param[282]<->rSAggTd;
	 param[283]<->iSAggTd;
	 param[284]<->rSAggTb5;
	 param[285]<->iSAggTb5;
	 param[286]<->rSAggTb2;
	 param[287]<->iSAggTb2;
	 param[288]<->rSAggT5;
	 param[289]<->iSAggT5;
	 param[290]<->rSAggT6;
	 param[291]<->iSAggT6;
	 param[292]<->rShyyta;
	 param[293]<->iShyyta;
	 param[294]<->rSHyyta;
	 param[295]<->iSHyyta;
	 param[296]<->rSAyyta;
	 param[297]<->iSAyyta;
	 param[298]<->rShyyw;
	 param[299]<->iShyyw;
	 param[300]<->rShyywp;
	 param[301]<->iShyywp;
	 param[302]<->rSHyyw;
	 param[303]<->iSHyyw;
	 param[304]<->rSHyywp;
	 param[305]<->iSHyywp;
	 param[306]<->rShyyH;
	 param[307]<->iShyyH;
	 param[308]<->rShyyeta;
	 param[309]<->iShyyeta;
	 param[310]<->rShyyphi;
	 param[311]<->iShyyphi;
	 param[312]<->rSHyyH;
	 param[313]<->iSHyyH;
	 param[314]<->rSHyyeta;
	 param[315]<->iSHyyeta;
	 param[316]<->rSHyyphi;
	 param[317]<->iSHyyphi;
	 */
	string nametext[318];
	nametext[0] = "fl";
	nametext[1] = "fh";
	nametext[2] = "mh0";
	nametext[3] = "ma0";
	nametext[4] = "mhh";
	nametext[5] = "tB";
	nametext[6] = "tA";
	nametext[7] = "tG";
	nametext[8] = "t12";
	nametext[9] = "t13";
	nametext[10] = "sigfact";
	nametext[11] = "kG";
	nametext[12] = "kY";
	nametext[13] = "kL";
	nametext[14] = "mhc";
	nametext[15] = "mphi0";
	nametext[16] = "mphic";
	nametext[17] = "meta0";
	nametext[18] = "metac";
	nametext[19] = "msigma";
	nametext[20] = "mtau";
	nametext[21] = "muq";
	nametext[22] = "mdq";
	nametext[23] = "msq";
	nametext[24] = "mcq";
	nametext[25] = "mbq";
	nametext[26] = "mtq";
	nametext[27] = "mTu";
	nametext[28] = "mTd";
	nametext[29] = "mTb5";
	nametext[30] = "mTb2";
	nametext[31] = "mT5";
	nametext[32] = "mT6";
	nametext[33] = "mWblh";
	nametext[34] = "mZblh";
	nametext[35] = "mwpblh";
	nametext[36] = "mzpblh";
	nametext[37] = "talpha";
	nametext[38] = "sapb";
	nametext[39] = "rDShgg";
	nametext[40] = "iDShgg";
	nametext[41] = "rShgg";
	nametext[42] = "iShgg";
	nametext[43] = "rDShyy";
	nametext[44] = "iDShyy";
	nametext[45] = "rShyy";
	nametext[46] = "iShyy";
	nametext[47] = "rDSHgg";
	nametext[48] = "iDSHgg";
	nametext[49] = "rSHgg";
	nametext[50] = "iSHgg";
	nametext[51] = "rDSHyy";
	nametext[52] = "iDSHyy";
	nametext[53] = "rSHyy";
	nametext[54] = "iSHyy";
	nametext[55] = "rDSAgg";
	nametext[56] = "iDSAgg";
	nametext[57] = "rSAgg";
	nametext[58] = "iSAgg";
	nametext[59] = "rDSAyy";
	nametext[60] = "iDSAyy";
	nametext[61] = "rSAyy";
	nametext[62] = "iSAyy";
	nametext[63] = "yv";
	nametext[64] = "ratiowmass";
	nametext[65] = "ywff";
	nametext[66] = "yzff";
	nametext[67] = "yhHH";
	nametext[68] = "yhphiphi";
	nametext[69] = "yhetaeta";
	nametext[70] = "yHHH";
	nametext[71] = "yHphiphi";
	nametext[72] = "yHetaeta";
	nametext[73] = "yAHH";
	nametext[74] = "yAphiphi";
	nametext[75] = "yAetaeta";
	nametext[76] = "yhtata";
	nametext[77] = "yHtata";
	nametext[78] = "yAtata";
	nametext[79] = "yhuu";
	nametext[80] = "yhdd";
	nametext[81] = "yhcc";
	nametext[82] = "yhss";
	nametext[83] = "yhbb";
	nametext[84] = "yhtt";
	nametext[85] = "yhTuTu";
	nametext[86] = "yhTdTd";
	nametext[87] = "yhTb5Tb5";
	nametext[88] = "yhTb2Tb2";
	nametext[89] = "yhT5T5";
	nametext[90] = "yhT6T6";
	nametext[91] = "yHuu";
	nametext[92] = "yHdd";
	nametext[93] = "yHcc";
	nametext[94] = "yHss";
	nametext[95] = "yHbb";
	nametext[96] = "yHtt";
	nametext[97] = "yHTuTu";
	nametext[98] = "yHTdTd";
	nametext[99] = "yHTb5Tb5";
	nametext[100] = "yHTb2Tb2";
	nametext[101] = "yHT5T5";
	nametext[102] = "yHT6T6";
	nametext[103] = "yAuu";
	nametext[104] = "yAdd";
	nametext[105] = "yAcc";
	nametext[106] = "yAss";
	nametext[107] = "yAbb";
	nametext[108] = "yAtt";
	nametext[109] = "yATuTu";
	nametext[110] = "yATdTd";
	nametext[111] = "yATb5Tb5";
	nametext[112] = "yATb2Tb2";
	nametext[113] = "yAT5T5";
	nametext[114] = "yAT6T6";
	nametext[115] = "yhww";
	nametext[116] = "yHww";
	nametext[117] = "yAww";
	nametext[118] = "yhwpwp";
	nametext[119] = "yHwpwp";
	nametext[120] = "yAwpwp";
	nametext[121] = "Chtt";
	nametext[122] = "Chbb";
	nametext[123] = "Chcc";
	nametext[124] = "Chtata";
	nametext[125] = "ChWW";
	nametext[126] = "ChZZ";
	nametext[127] = "Chyy";
	nametext[128] = "ChyZ";
	nametext[129] = "Chgg";
	nametext[130] = "CHtt";
	nametext[131] = "CHbb";
	nametext[132] = "CHcc";
	nametext[133] = "CHtata";
	nametext[134] = "CHWW";
	nametext[135] = "CHZZ";
	nametext[136] = "CHyy";
	nametext[137] = "CHyZ";
	nametext[138] = "CHgg";
	nametext[139] = "CAtt";
	nametext[140] = "CAbb";
	nametext[141] = "CAcc";
	nametext[142] = "CAtata";
	nametext[143] = "CAWW";
	nametext[144] = "CAZZ";
	nametext[145] = "CAyy";
	nametext[146] = "CAyZ";
	nametext[147] = "CAgg";
	nametext[148] = "rShyyu";
	nametext[149] = "iShyyu";
	nametext[150] = "rShyyd";
	nametext[151] = "iShyyd";
	nametext[152] = "rShyys";
	nametext[153] = "iShyys";
	nametext[154] = "rShyyc";
	nametext[155] = "iShyyc";
	nametext[156] = "rShyyb";
	nametext[157] = "iShyyb";
	nametext[158] = "rShyyt";
	nametext[159] = "iShyyt";
	nametext[160] = "rShyyTu";
	nametext[161] = "iShyyTu";
	nametext[162] = "rShyyTd";
	nametext[163] = "iShyyTd";
	nametext[164] = "rShyyTb5";
	nametext[165] = "iShyyTb5";
	nametext[166] = "rShyyTb2";
	nametext[167] = "iShyyTb2";
	nametext[168] = "rShyyT5";
	nametext[169] = "iShyyT5";
	nametext[170] = "rShyyT6";
	nametext[171] = "iShyyT6";
	nametext[172] = "rShggu";
	nametext[173] = "iShggu";
	nametext[174] = "rShggd";
	nametext[175] = "iShggd";
	nametext[176] = "rShggs";
	nametext[177] = "iShggs";
	nametext[178] = "rShggc";
	nametext[179] = "iShggc";
	nametext[180] = "rShggb";
	nametext[181] = "iShggb";
	nametext[182] = "rShggt";
	nametext[183] = "iShggt";
	nametext[184] = "rShggTu";
	nametext[185] = "iShggTu";
	nametext[186] = "rShggTd";
	nametext[187] = "iShggTd";
	nametext[188] = "rShggTb5";
	nametext[189] = "iShggTb5";
	nametext[190] = "rShggTb2";
	nametext[191] = "iShggTb2";
	nametext[192] = "rShggT5";
	nametext[193] = "iShggT5";
	nametext[194] = "rShggT6";
	nametext[195] = "iShggT6";
	nametext[196] = "rSHyyu";
	nametext[197] = "iSHyyu";
	nametext[198] = "rSHyyd";
	nametext[199] = "iSHyyd";
	nametext[200] = "rSHyys";
	nametext[201] = "iSHyys";
	nametext[202] = "rSHyyc";
	nametext[203] = "iSHyyc";
	nametext[204] = "rSHyyb";
	nametext[205] = "iSHyyb";
	nametext[206] = "rSHyyt";
	nametext[207] = "iSHyyt";
	nametext[208] = "rSHyyTu";
	nametext[209] = "iSHyyTu";
	nametext[210] = "rSHyyTd";
	nametext[211] = "iSHyyTd";
	nametext[212] = "rSHyyTb5";
	nametext[213] = "iSHyyTb5";
	nametext[214] = "rSHyyTb2";
	nametext[215] = "iSHyyTb2";
	nametext[216] = "rSHyyT5";
	nametext[217] = "iSHyyT5";
	nametext[218] = "rSHyyT6";
	nametext[219] = "iSHyyT6";
	nametext[220] = "rSHggu";
	nametext[221] = "iSHggu";
	nametext[222] = "rSHggd";
	nametext[223] = "iSHggd";
	nametext[224] = "rSHggs";
	nametext[225] = "iSHggs";
	nametext[226] = "rSHggc";
	nametext[227] = "iSHggc";
	nametext[228] = "rSHggb";
	nametext[229] = "iSHggb";
	nametext[230] = "rSHggt";
	nametext[231] = "iSHggt";
	nametext[232] = "rSHggTu";
	nametext[233] = "iSHggTu";
	nametext[234] = "rSHggTd";
	nametext[235] = "iSHggTd";
	nametext[236] = "rSHggTb5";
	nametext[237] = "iSHggTb5";
	nametext[238] = "rSHggTb2";
	nametext[239] = "iSHggTb2";
	nametext[240] = "rSHggT5";
	nametext[241] = "iSHggT5";
	nametext[242] = "rSHggT6";
	nametext[243] = "iSHggT6";
	nametext[244] = "rSAyyu";
	nametext[245] = "iSAyyu";
	nametext[246] = "rSAyyd";
	nametext[247] = "iSAyyd";
	nametext[248] = "rSAyys";
	nametext[249] = "iSAyys";
	nametext[250] = "rSAyyc";
	nametext[251] = "iSAyyc";
	nametext[252] = "rSAyyb";
	nametext[253] = "iSAyyb";
	nametext[254] = "rSAyyt";
	nametext[255] = "iSAyyt";
	nametext[256] = "rSAyyTu";
	nametext[257] = "iSAyyTu";
	nametext[258] = "rSAyyTd";
	nametext[259] = "iSAyyTd";
	nametext[260] = "rSAyyTb5";
	nametext[261] = "iSAyyTb5";
	nametext[262] = "rSAyyTb2";
	nametext[263] = "iSAyyTb2";
	nametext[264] = "rSAyyT5";
	nametext[265] = "iSAyyT5";
	nametext[266] = "rSAyyT6";
	nametext[267] = "iSAyyT6";
	nametext[268] = "rSAggu";
	nametext[269] = "iSAggu";
	nametext[270] = "rSAggd";
	nametext[271] = "iSAggd";
	nametext[272] = "rSAggs";
	nametext[273] = "iSAggs";
	nametext[274] = "rSAggc";
	nametext[275] = "iSAggc";
	nametext[276] = "rSAggb";
	nametext[277] = "iSAggb";
	nametext[278] = "rSAggt";
	nametext[279] = "iSAggt";
	nametext[280] = "rSAggTu";
	nametext[281] = "iSAggTu";
	nametext[282] = "rSAggTd";
	nametext[283] = "iSAggTd";
	nametext[284] = "rSAggTb5";
	nametext[285] = "iSAggTb5";
	nametext[286] = "rSAggTb2";
	nametext[287] = "iSAggTb2";
	nametext[288] = "rSAggT5";
	nametext[289] = "iSAggT5";
	nametext[290] = "rSAggT6";
	nametext[291] = "iSAggT6";
	nametext[292] = "rShyyta";
	nametext[293] = "iShyyta";
	nametext[294] = "rSHyyta";
	nametext[295] = "iSHyyta";
	nametext[296] = "rSAyyta";
	nametext[297] = "iSAyyta";
	nametext[298] = "rShyyw";
	nametext[299] = "iShyyw";
	nametext[300] = "rShyywp";
	nametext[301] = "iShyywp";
	nametext[302] = "rSHyyw";
	nametext[303] = "iSHyyw";
	nametext[304] = "rSHyywp";
	nametext[305] = "iSHyywp";
	nametext[306] = "rShyyH";
	nametext[307] = "iShyyH";
	nametext[308] = "rShyyeta";
	nametext[309] = "iShyyeta";
	nametext[310] = "rShyyphi";
	nametext[311] = "iShyyphi";
	nametext[312] = "rSHyyH";
	nametext[313] = "iSHyyH";
	nametext[314] = "rSHyyeta";
	nametext[315] = "iSHyyeta";
	nametext[316] = "rSHyyphi";
	nametext[317] = "iSHyyphi";
	
	string labeltext[318];
	labeltext[0] = "f (GeV)";
	labeltext[1] = "F (GeV)";
	labeltext[2] = "m_{h} (GeV)";
	labeltext[3] = "m_{A} (GeV)";
	labeltext[4] = "m_{H} (GeV) ";
	labeltext[5] = "tan#beta";
	labeltext[6] = "tan#alpha";
	labeltext[7] = "tan#theta_{g}";
	labeltext[8] = "tan#theta_{12}";
	labeltext[9] = "tan#theta_{13}";
	labeltext[10] = "K_{#sigma}";
	labeltext[11] = "K_{G}";
	labeltext[12] = "K_{Y}";
	labeltext[13] = "K_{L}";
	labeltext[14] = "m_{H^{+}} (GeV)";
	labeltext[15] = "m_{#phi} (GeV)";
	labeltext[16] = "m_{#phi^{+}} (GeV)";
	labeltext[17] = "m_{#eta} (GeV)";
	labeltext[18] = "m_{#eta^{+}} (GeV)";
	labeltext[19] = "m_{#sigma} (GeV)";
	labeltext[20] = "m_{#tau} (GeV)";
	labeltext[21] = "m_{u} (GeV)";
	labeltext[22] = "m_{d} (GeV)";
	labeltext[23] = "m_{s} (GeV)";
	labeltext[24] = "m_{c} (GeV)";
	labeltext[25] = "m_{b} (GeV)";
	labeltext[26] = "m_{t} (GeV)";
	labeltext[27] = "m_{T}  (GeV)";
	labeltext[28] = "m_{B}  (GeV)";
	labeltext[29] = "m_{T_{b5}} (GeV)";
	labeltext[30] = "m_{T_{b2}} (GeV)";
	labeltext[31] = "m_{T_{5}} (GeV)";
	labeltext[32] = "m_{T_{6}} (GeV)";
	labeltext[33] = "m_{W} (GeV)";
	labeltext[34] = "m_{Z} (GeV)";
	labeltext[35] = "m_{W\'} (GeV)";
	labeltext[36] = "m_{Z\'} (GeV)";
	labeltext[37] = "tan#alpha";
	labeltext[38] = "sin(#alpha+#beta)";
	labeltext[39] = "Re(#Delta S_{hgg})";
	labeltext[40] = "Im(#Delta S_{hgg})";
	labeltext[41] = "Re(S_{hgg})";
	labeltext[42] = "Im(S_{hgg})";
	labeltext[43] = "Re(#Delta S_{h#gamma#gamma})";
	labeltext[44] = "Im(#Delta S_{h#gamma#gamma})";
	labeltext[45] = "Re(S_{h#gamma#gamma})";
	labeltext[46] = "Im(S_{h#gamma#gamma})";
	labeltext[47] = "Re(#Delta S_{Hgg})";
	labeltext[48] = "Im(#Delta S_{Hgg})";
	labeltext[49] = "Re(S_{Hgg})";
	labeltext[50] = "Im(S_{Hgg})";
	labeltext[51] = "Re(#Delta S_{H#gamma#gamma})";
	labeltext[52] = "Im(#Delta S_{H#gamma#gamma})";
	labeltext[53] = "Re(S_{H#gamma#gamma})";
	labeltext[54] = "Im(S_{H#gamma#gamma})";
	labeltext[55] = "Re(#Delta S_{Agg})";
	labeltext[56] = "Im(#Delta S_{Agg})";
	labeltext[57] = "Re(S_{Agg})";
	labeltext[58] = "Im(S_{Agg})";
	labeltext[59] = "Re(#Delta S_{A#gamma#gamma})";
	labeltext[60] = "Im(#Delta S_{A#gamma#gamma})";
	labeltext[61] = "Re(S_{A#gamma#gamma})";
	labeltext[62] = "Im(S_{A#gamma#gamma})";
	labeltext[63] = "y_{v}";
	labeltext[64] = "m_{W}^{BLH}/m_{W}^{SM}";
	labeltext[65] = "y_{Wff}";
	labeltext[66] = "y_{Zff}";
	labeltext[67] = "y_{hHH}";
	labeltext[68] = "y_{h#phi#phi}";
	labeltext[69] = "y_{h#eta#eta}";
	labeltext[70] = "y_{HHH}";
	labeltext[71] = "y_{H#phi#phi}";
	labeltext[72] = "y_{H#eta#eta}";
	labeltext[73] = "y_{AHH}";
	labeltext[74] = "y_{Ah#phi#phi}";
	labeltext[75] = "y_{A#eta#eta}";
	labeltext[76] = "y_{h#tau#tau}";
	labeltext[77] = "y_{H#tau#tau}";
	labeltext[78] = "y_{A#tau#tau}";
	labeltext[79] = "y_{huu}";
	labeltext[80] = "y_{hdd}";
	labeltext[81] = "y_{hcc}";
	labeltext[82] = "y_{hss}";
	labeltext[83] = "y_{hbb}";
	labeltext[84] = "y_{htt}";
	labeltext[85] = "y_{hTT}";
	labeltext[86] = "y_{hBB}";
	labeltext[87] = "y_{hT_{b5}T_{b5}}";
	labeltext[88] = "y_{hT_{b2}T_{b2}}";
	labeltext[89] = "y_{hT_{5}T_{5}}";
	labeltext[90] = "y_{hT_{6}T_{6}}";
	labeltext[91] = "y_{Huu}";
	labeltext[92] = "y_{Hdd}";
	labeltext[93] = "y_{Hcc}";
	labeltext[94] = "y_{Hss}";
	labeltext[95] = "y_{Hbb}";
	labeltext[96] = "y_{Htt}";
	labeltext[97] = "y_{HTT}";
	labeltext[98] = "y_{HBB}";
	labeltext[99] = "y_{HT_{b5}T_{b5}}";
	labeltext[100] = "y_{HT_{b2}T_{b2}}";
	labeltext[101] = "y_{HT_{5}T_{5}}";
	labeltext[102] = "y_{HT_{6}T_{6}}";
	labeltext[103] = "y_{Auu}";
	labeltext[104] = "y_{Add}";
	labeltext[105] = "y_{Acc}";
	labeltext[106] = "y_{Ass}";
	labeltext[107] = "y_{Abb}";
	labeltext[108] = "y_{Att}";
	labeltext[109] = "y_{ATT}";
	labeltext[110] = "y_{ABB}";
	labeltext[111] = "y_{AT_{b5}T_{b5}}";
	labeltext[112] = "y_{AT_{b2}T_{b2}}";
	labeltext[113] = "y_{AT_{5}T_{5}}";
	labeltext[114] = "y_{AT_{6}T_{6}}";
	labeltext[115] = "y_{hWW}";
	labeltext[116] = "y_{HWW}";
	labeltext[117] = "y_{AWW}";
	labeltext[118] = "y_{hW\'W\'}";
	labeltext[119] = "y_{HW\'W\'}";
	labeltext[120] = "y_{AW\'W\'}";
	labeltext[121] = "C_{htt}";
	labeltext[122] = "C_{hbb}";
	labeltext[123] = "C_{hcc}";
	labeltext[124] = "C_{h#tau#tau}";
	labeltext[125] = "C_{hWW}";
	labeltext[126] = "C_{hZZ}";
	labeltext[127] = "C_{h#gamma#gamma}";
	labeltext[128] = "C_{h#gammaZ}";
	labeltext[129] = "C_{hgg}";
	labeltext[130] = "C_{Htt}";
	labeltext[131] = "C_{Hbb}";
	labeltext[132] = "C_{Hcc}";
	labeltext[133] = "C_{H#tau#tau}";
	labeltext[134] = "C_{HWW}";
	labeltext[135] = "C_{HZZ}";
	labeltext[136] = "C_{H#gamma#gamma}";
	labeltext[137] = "C_{H#gammaZ}";
	labeltext[138] = "C_{Hgg}";
	labeltext[139] = "C_{Att}";
	labeltext[140] = "C_{Abb}";
	labeltext[141] = "C_{Acc}";
	labeltext[142] = "C_{A#tau#tau}";
	labeltext[143] = "C_{AWW}";
	labeltext[144] = "C_{AZZ}";
	labeltext[145] = "C_{A#gamma#gamma}";
	labeltext[146] = "C_{A#gammaZ}";
	labeltext[147] = "C_{Agg}";
	labeltext[148] = "Re(S_{h#gamma#gamma-u})";
	labeltext[149] = "Im(S_{h#gamma#gamma-u})";
	labeltext[150] = "Re(S_{h#gamma#gamma-d})";
	labeltext[151] = "Im(S_{h#gamma#gamma-d})";
	labeltext[152] = "Re(S_{h#gamma#gamma-s})";
	labeltext[153] = "Im(S_{h#gamma#gamma-s})";
	labeltext[154] = "Re(S_{h#gamma#gamma-c})";
	labeltext[155] = "Im(S_{h#gamma#gamma-c})";
	labeltext[156] = "Re(S_{h#gamma#gamma-b})";
	labeltext[157] = "Im(S_{h#gamma#gamma-b})";
	labeltext[158] = "Re(S_{h#gamma#gamma-t})";
	labeltext[159] = "Im(S_{h#gamma#gamma-t})";
	labeltext[160] = "Re(S_{h#gamma#gamma-T})";
	labeltext[161] = "Im(S_{h#gamma#gamma-T})";
	labeltext[162] = "Re(S_{h#gamma#gamma-B})";
	labeltext[163] = "Im(S_{h#gamma#gamma-B})";
	labeltext[164] = "Re(S_{h#gamma#gamma-Tb5})";
	labeltext[165] = "Im(S_{h#gamma#gamma-Tb5})";
	labeltext[166] = "Re(S_{h#gamma#gamma-Tb2})";
	labeltext[167] = "Im(S_{h#gamma#gamma-Tb2})";
	labeltext[168] = "Re(S_{h#gamma#gamma-T5})";
	labeltext[169] = "Im(S_{h#gamma#gamma-T5})";
	labeltext[170] = "Re(S_{h#gamma#gamma-T6})";
	labeltext[171] = "Im(S_{h#gamma#gamma-T6})";
	labeltext[172] = "Re(S_{hgg-u})";
	labeltext[173] = "Im(S_{hgg-u})";
	labeltext[174] = "Re(S_{hgg-d})";
	labeltext[175] = "Im(S_{hgg-d})";
	labeltext[176] = "Re(S_{hgg-s})";
	labeltext[177] = "Im(S_{hgg-s})";
	labeltext[178] = "Re(S_{hgg-c})";
	labeltext[179] = "Im(S_{hgg-c})";
	labeltext[180] = "Re(S_{hgg-b})";
	labeltext[181] = "Im(S_{hgg-b})";
	labeltext[182] = "Re(S_{hgg-t})";
	labeltext[183] = "Im(S_{hgg-t})";
	labeltext[184] = "Re(S_{hgg-T})";
	labeltext[185] = "Im(S_{hgg-T})";
	labeltext[186] = "Re(S_{hgg-B})";
	labeltext[187] = "Im(S_{hgg-B})";
	labeltext[188] = "Re(S_{hgg-Tb5})";
	labeltext[189] = "Im(S_{hgg-Tb5})";
	labeltext[190] = "Re(S_{hgg-Tb2})";
	labeltext[191] = "Im(S_{hgg-Tb2})";
	labeltext[192] = "Re(S_{hgg-T5})";
	labeltext[193] = "Im(S_{hgg-T5})";
	labeltext[194] = "Re(S_{hgg-T6})";
	labeltext[195] = "Im(S_{hgg-T6})";
	labeltext[196] = "Re(S_{H#gamma#gamma-u})";
	labeltext[197] = "Im(S_{H#gamma#gamma-u})";
	labeltext[198] = "Re(S_{H#gamma#gamma-d})";
	labeltext[199] = "Im(S_{H#gamma#gamma-d})";
	labeltext[200] = "Re(S_{H#gamma#gamma-s})";
	labeltext[201] = "Im(S_{H#gamma#gamma-s})";
	labeltext[202] = "Re(S_{H#gamma#gamma-c})";
	labeltext[203] = "Im(S_{H#gamma#gamma-c})";
	labeltext[204] = "Re(S_{H#gamma#gamma-b})";
	labeltext[205] = "Im(S_{H#gamma#gamma-b})";
	labeltext[206] = "Re(S_{H#gamma#gamma-t})";
	labeltext[207] = "Im(S_{H#gamma#gamma-t})";
	labeltext[208] = "Re(S_{H#gamma#gamma-T})";
	labeltext[209] = "Im(S_{H#gamma#gamma-T})";
	labeltext[210] = "Re(S_{H#gamma#gamma-B})";
	labeltext[211] = "Im(S_{H#gamma#gamma-B})";
	labeltext[212] = "Re(S_{H#gamma#gamma-Tb5})";
	labeltext[213] = "Im(S_{H#gamma#gamma-Tb5})";
	labeltext[214] = "Re(S_{H#gamma#gamma-Tb2})";
	labeltext[215] = "Im(S_{H#gamma#gamma-Tb2})";
	labeltext[216] = "Re(S_{H#gamma#gamma-T5})";
	labeltext[217] = "Im(S_{H#gamma#gamma-T5})";
	labeltext[218] = "Re(S_{H#gamma#gamma-T6})";
	labeltext[219] = "Im(S_{H#gamma#gamma-T6})";
	labeltext[220] = "Re(S_{Hgg-u})";
	labeltext[221] = "Im(S_{Hgg-u})";
	labeltext[222] = "Re(S_{Hgg-d})";
	labeltext[223] = "Im(S_{Hgg-d})";
	labeltext[224] = "Re(S_{Hgg-s})";
	labeltext[225] = "Im(S_{Hgg-s})";
	labeltext[226] = "Re(S_{Hgg-c})";
	labeltext[227] = "Im(S_{Hgg-c})";
	labeltext[228] = "Re(S_{Hgg-b})";
	labeltext[229] = "Im(S_{Hgg-b})";
	labeltext[230] = "Re(S_{Hgg-t})";
	labeltext[231] = "Im(S_{Hgg-t})";
	labeltext[232] = "Re(S_{Hgg-T})";
	labeltext[233] = "Im(S_{Hgg-T})";
	labeltext[234] = "Re(S_{Hgg-B})";
	labeltext[235] = "Im(S_{Hgg-B})";
	labeltext[236] = "Re(S_{Hgg-Tb5})";
	labeltext[237] = "Im(S_{Hgg-Tb5})";
	labeltext[238] = "Re(S_{Hgg-Tb2})";
	labeltext[239] = "Im(S_{Hgg-Tb2})";
	labeltext[240] = "Re(S_{Hgg-T5})";
	labeltext[241] = "Im(S_{Hgg-T5})";
	labeltext[242] = "Re(S_{Hgg-T6})";
	labeltext[243] = "Im(S_{Hgg-T6})";
	labeltext[244] = "Re(S_{A#gamma#gamma-u})";
	labeltext[245] = "Im(S_{A#gamma#gamma-u})";
	labeltext[246] = "Re(S_{A#gamma#gamma-d})";
	labeltext[247] = "Im(S_{A#gamma#gamma-d})";
	labeltext[248] = "Re(S_{A#gamma#gamma-s})";
	labeltext[249] = "Im(S_{A#gamma#gamma-s})";
	labeltext[250] = "Re(S_{A#gamma#gamma-c})";
	labeltext[251] = "Im(S_{A#gamma#gamma-c})";
	labeltext[252] = "Re(S_{A#gamma#gamma-b})";
	labeltext[253] = "Im(S_{A#gamma#gamma-b})";
	labeltext[254] = "Re(S_{A#gamma#gamma-t})";
	labeltext[255] = "Im(S_{A#gamma#gamma-t})";
	labeltext[256] = "Re(S_{A#gamma#gamma-T})";
	labeltext[257] = "Im(S_{A#gamma#gamma-T})";
	labeltext[258] = "Re(S_{A#gamma#gamma-B})";
	labeltext[259] = "Im(S_{A#gamma#gamma-B})";
	labeltext[260] = "Re(S_{A#gamma#gamma-Tb5})";
	labeltext[261] = "Im(S_{A#gamma#gamma-Tb5})";
	labeltext[262] = "Re(S_{A#gamma#gamma-Tb2})";
	labeltext[263] = "Im(S_{A#gamma#gamma-Tb2})";
	labeltext[264] = "Re(S_{A#gamma#gamma-T5})";
	labeltext[265] = "Im(S_{A#gamma#gamma-T5})";
	labeltext[266] = "Re(S_{A#gamma#gamma-T6})";
	labeltext[267] = "Im(S_{A#gamma#gamma-T6})";
	labeltext[268] = "Re(S_{Agg-u})";
	labeltext[269] = "Im(S_{Agg-u})";
	labeltext[270] = "Re(S_{Agg-d})";
	labeltext[271] = "Im(S_{Agg-d})";
	labeltext[272] = "Re(S_{Agg-s})";
	labeltext[273] = "Im(S_{Agg-s})";
	labeltext[274] = "Re(S_{Agg-c})";
	labeltext[275] = "Im(S_{Agg-c})";
	labeltext[276] = "Re(S_{Agg-b})";
	labeltext[277] = "Im(S_{Agg-b})";
	labeltext[278] = "Re(S_{Agg-t})";
	labeltext[279] = "Im(S_{Agg-t})";
	labeltext[280] = "Re(S_{Agg-T})";
	labeltext[281] = "Im(S_{Agg-T})";
	labeltext[282] = "Re(S_{Agg-B})";
	labeltext[283] = "Im(S_{Agg-B})";
	labeltext[284] = "Re(S_{Agg-Tb5})";
	labeltext[285] = "Im(S_{Agg-Tb5})";
	labeltext[286] = "Re(S_{Agg-Tb2})";
	labeltext[287] = "Im(S_{Agg-Tb2})";
	labeltext[288] = "Re(S_{Agg-T5})";
	labeltext[289] = "Im(S_{Agg-T5})";
	labeltext[290] = "Re(S_{Agg-T6})";
	labeltext[291] = "Im(S_{Agg-T6})";
	labeltext[292] = "Re(S_{h#gamma#gamma-#tau})";
	labeltext[293] = "Im(S_{h#gamma#gamma-#tau})";
	labeltext[294] = "Re(S_{H#gamma#gamma-#tau})";
	labeltext[295] = "Im(S_{H#gamma#gamma-#tau})";
	labeltext[296] = "Re(S_{A#gamma#gamma-#tau})";
	labeltext[297] = "Im(S_{A#gamma#gamma-#tau})";
	labeltext[298] = "Re(S_{h#gamma#gamma-W})";
	labeltext[299] = "Im(S_{h#gamma#gamma-W})";
	labeltext[300] = "Re(S_{h#gamma#gamma-W\'})";
	labeltext[301] = "Im(S_{h#gamma#gamma-W\'})";
	labeltext[302] = "Re(S_{H#gamma#gamma-W})";
	labeltext[303] = "Im(S_{H#gamma#gamma-W})";
	labeltext[304] = "Re(S_{H#gamma#gamma-W\'})";
	labeltext[305] = "Im(S_{H#gamma#gamma-W\'})";
	labeltext[306] = "Re(S_{h#gamma#gamma-H})";
	labeltext[307] = "Im(S_{h#gamma#gamma-H})";
	labeltext[308] = "Re(S_{h#gamma#gamma-#eta})";
	labeltext[309] = "Im(S_{h#gamma#gamma-#eta})";
	labeltext[310] = "Re(S_{h#gamma#gamma-#phi})";
	labeltext[311] = "Im(S_{h#gamma#gamma-#phi})";
	labeltext[312] = "Re(S_{H#gamma#gamma-H})";
	labeltext[313] = "Im(S_{H#gamma#gamma-H})";
	labeltext[314] = "Re(S_{H#gamma#gamma-#eta})";
	labeltext[315] = "Im(S_{H#gamma#gamma-#eta})";
	labeltext[316] = "Re(S_{H#gamma#gamma-#phi})";
	labeltext[317] = "Im(S_{H#gamma#gamma-#phi})";
	
	
	
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
	
	/*double param[0]min=9999, param[0]max=0;
	double param[1]min=9999, param[1]max=0;
	double param[2]min=9999, param[2]max=0;
	double param[3]min=9999, param[3]max=0;
	double param[5]min=9999, param[5]max=0;
	double param[6]min=9999, param[6]max=0;
	double param[7]min=9999, param[7]max=0;
	double DSggmin=9999, DSggmax=0;
	double Sggmin=9999, Sggmax=0;
	double DSyymin=9999, DSyymax=0;
	double Syymin=9999, Syymax=0;
	double Cyymin=9999, Cyymax=0;
	double Cggmin=9999, Cggmax=0;
	double CWWmin=9999, CWWmax=0;
	double CZZmin=9999, CZZmax=0;
	double Cbbmin=9999, Cbbmax=0;
	double Ctatamin=9999, Ctatamax=0;*/
	double deltammin=9999, deltammax=0;
	double muyymin=9999, muyymax=0;
	double muZZmin=9999, muZZmax=0;
	double muWWmin=9999, muWWmax=0;
	double mubbmin=9999, mubbmax=0;
	double mutatamin=9999, mutatamax=0;
	double muZymin=9999, muZymax=0;
	//double param[4]min=9999, param[4]max=0;
	
	
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
				param[0].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[1].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[2].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[3].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[4].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[14].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[15].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[16].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[17].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[18].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[19].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[20].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[21].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[22].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[23].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[24].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[25].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[26].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[27].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[28].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[29].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[30].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[31].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[32].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[33].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[34].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[35].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[36].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[5].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[37].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[38].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[7].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[8].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[9].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[10].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[11].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[12].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[13].push_back(atof(line.substr(start,end-start).c_str()));//38
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[39].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[40].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[41].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[42].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[43].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[44].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[45].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[46].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[47].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[48].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[49].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[50].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[51].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[52].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[53].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[54].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[55].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[56].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[57].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[58].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[59].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[60].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[61].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[62].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[63].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[64].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[65].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[66].push_back(atof(line.substr(start,end-start).c_str()));
			}
			if ( line.find("yvalues",0) != string::npos ) {
				start = 8;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[67].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[68].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[69].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[70].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[71].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[72].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[73].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[74].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[75].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[76].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[77].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[78].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[79].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[91].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[103].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[80].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[92].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[104].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[81].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[93].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[105].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[82].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[94].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[106].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[84].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[96].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[108].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[83].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[95].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[107].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[85].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[97].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[109].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[86].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[98].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[110].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[87].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[99].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[111].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[88].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[100].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[112].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[89].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[101].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[113].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[90].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[102].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[114].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[115].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[116].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[117].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[118].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[119].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[120].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[121].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[122].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[123].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[124].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[125].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[126].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[127].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[128].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[129].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[130].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[131].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[132].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[133].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[134].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[135].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[136].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[137].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[138].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[139].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[140].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[141].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[142].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[143].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[144].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[145].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[146].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[147].push_back(atof(line.substr(start,end-start).c_str()));
				
				
				
				
				//cout << param[0].size() << " ";
				//cout << param[1].size() << " ";
				//cout << param[2].size() << " ";
				//cout << param[3].size() << " ";
				//cout << param[4].size() << " ";
				//cout << param[5].size() << " ";
				//cout << param[7].size() << " ";
				//cout << param[8].size() << " ";
				//cout << param[9].size() << " ";
				//cout << param[10].size() << " ";
				//cout << param[11].size() << " ";
				//cout << param[12].size() << " ";
				//cout << param[13].size() << endl;
			}
			if ( line.find("loopfactsh",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[148].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[149].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[150].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[151].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[154].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[155].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[152].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[153].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[158].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[159].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[156].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[157].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[160].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[161].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[162].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[163].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[164].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[165].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[166].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[167].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[168].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[169].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[170].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[171].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[292].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[293].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[298].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[299].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[300].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[301].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[306].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[307].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[310].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[311].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[308].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[309].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[172].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[173].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[174].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[175].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[178].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[179].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[176].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[177].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[182].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[183].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[180].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[181].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[184].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[185].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[186].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[187].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[188].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[189].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[190].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[191].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[192].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[193].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[194].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[195].push_back(atof(line.substr(start,end-start).c_str()));
				
			}
			if ( line.find("loopfactsH",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[196].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[197].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[198].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[199].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[202].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[203].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[200].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[201].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[206].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[207].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[204].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[205].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[208].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[209].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[210].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[211].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[212].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[213].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[214].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[215].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[216].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[217].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[218].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[219].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[294].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[295].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[302].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[303].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[304].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[305].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[312].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[313].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[316].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[317].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[314].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[315].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[220].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[221].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[222].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[223].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[226].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[227].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[224].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[225].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[230].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[231].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[228].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[229].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[232].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[233].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[234].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[235].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[236].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[237].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[238].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[239].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[240].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[241].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[242].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[243].push_back(atof(line.substr(start,end-start).c_str()));
				
			}
			if ( line.find("loopfactsA",0) != string::npos ) {
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[244].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[245].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[246].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[247].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[250].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[251].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[248].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[249].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[254].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[255].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[252].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[253].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[256].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[257].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[258].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[259].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[260].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[261].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[262].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[263].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[264].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[265].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[266].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[267].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[296].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[297].push_back(atof(line.substr(start,end-start).c_str()));
				
				start = 11;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[268].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[269].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[270].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[271].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[274].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[275].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[272].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[273].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[278].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[279].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[276].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[277].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[280].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[281].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[282].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[283].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[284].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[285].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[286].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[287].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[288].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[289].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[290].push_back(atof(line.substr(start,end-start).c_str()));
				start = end+1;
				end = line.find_first_not_of("0123456789Ee.-+",start);
				param[291].push_back(atof(line.substr(start,end-start).c_str()));
				
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
			
			if ( line.find("pseudo",0) != string::npos && param[3].at(param[3].size()-1) < 900 ) {
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
			else if ( line.find("pseudo",0) != string::npos && param[3].at(param[3].size()-1) > 900 ) {
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
			
			if ( line.find("heavyh",0) != string::npos && param[4].at(param[4].size()-1) < 900 ) {
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
			else if (  line.find("heavyh",0) != string::npos && param[4].at(param[4].size()-1) >= 900 ) {
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
	
	for ( int i=0; i<param[0].size(); i++) {
		deltam.push_back(param[3].at(i)-param[2].at(i));
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
		 cout << param[0].at(i) << " ";
		 cout << param[1].at(i) << " ";
		 cout << param[2].at(i) << " ";
		 cout << param[3].at(i) << " ";
		 cout << param[4].at(i) << " ";
		 cout << param[5].at(i) << " ";
		 cout << param[7].at(i) << " ";
		 cout << param[8].at(i) << " ";
		 cout << param[9].at(i) << " ";
		 cout << param[10].at(i) << " ";
		 cout << param[11].at(i) << " ";
		 cout << param[12].at(i) << " ";
		 cout << param[13].at(i) << " ";
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
		 cout << param[0].at(i) << " ";
		 cout << param[1].at(i) << " ";
		 cout << param[2].at(i) << " ";
		 cout << param[3].at(i) << " ";
		 cout << param[4].at(i) << " ";
		 cout << param[5].at(i) << " ";
		 cout << param[7].at(i) << " ";
		 cout << param[8].at(i) << " ";
		 cout << param[9].at(i) << " ";
		 cout << param[10].at(i) << " ";
		 cout << param[11].at(i) << " ";
		 cout << param[12].at(i) << " ";
		 cout << param[13].at(i) << " ";
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
	
	for ( int i=0; i<param[0].size(); i++ ) {
		vblh = vsm*sqrt(1. + (1./6. + (5./4. - 2.*sg*sg)*param[0].at(i)*param[0].at(i)/(param[0].at(i)*param[0].at(i) + param[1].at(i)*param[1].at(i)))*vsm*vsm/(param[0].at(i)*param[0].at(i)));
		sw = sqrt(0.5 - 0.5*sqrt(1.0 - (4.0*3.14159*alphaem*vblh*vblh )/(MZ*MZ)));
		MW = sqrt(1.0 - sw*sw)*MZ;
		sg = param[7].at(i)/sqrt(1.0+param[7].at(i)*param[7].at(i));
		g = sqrt(4.0*3.14159*alphaem)/sw;
		gA = g/sqrt(1.0-sg*sg);
		gB = g/sg;
		MWp=sqrt(0.25* (gA*gA + gB*gB)*(param[0].at(i)*param[0].at(i) + param[1].at(i)*param[1].at(i)) - MW*MW);
		MWpconstraint=14.182931524204026 - 23.378448106386056*sg - 0.042544864480156076*sg*sg + 17.296487025431006*sg*sg*sg*sg;
		
		
		if ( param[4].at(i) >= 128.2 && param[4].at(i) < 150 ) {
			mhwwlimit=pow(10.0,-1.38188 + 85958./(param[4].at(i)*param[4].at(i)) - 1219.03/param[4].at(i) + 0.353545*sqrt(param[4].at(i)) + 0.0523065*param[4].at(i) - 0.000258021*param[4].at(i)*param[4].at(i));
		}
		else if ( param[4].at(i) >= 150 && param[4].at(i) < 170 ) {
			mhwwlimit=pow(10.0,-2.70814 + 104083./(param[4].at(i)*param[4].at(i)) - 147.355/param[4].at(i) - 0.186522*sqrt(param[4].at(i)) - 0.00739727*param[4].at(i) + 0.000139588*param[4].at(i)*param[4].at(i));
		}
		else if ( param[4].at(i) >= 170 && param[4].at(i) < 200 ) {
			mhwwlimit=pow(10.0,19.2623 - 1746340/(param[4].at(i)*param[4].at(i)) + 18418.7/param[4].at(i) - 2.90129*sqrt(param[4].at(i)) - 0.393658*param[4].at(i) + 0.00131283*param[4].at(i)*param[4].at(i));
		}
		else if ( param[4].at(i) >= 200 && param[4].at(i) < 300 ) {
			mhwwlimit=pow(10.0,0.0135371 + 28937.7/(param[4].at(i)*param[4].at(i)) + 73.4727/param[4].at(i) - 0.00883728*sqrt(param[4].at(i)) - 0.00112041*param[4].at(i) - 0.0000076754*param[4].at(i)*param[4].at(i));
		}
		else if ( param[4].at(i) >= 300 ) {
			mhwwlimit=pow(10.0,269.072 + 1569510/(param[4].at(i)*param[4].at(i)) - 25921.9/param[4].at(i) - 17.797*sqrt(param[4].at(i)) + 0.376611*param[4].at(i) - 0.000060398*param[4].at(i)*param[4].at(i));
		}
		else {
			mhwwlimit = 100;
		}
		
		if ( param[4].at(i) >=128.2 && param[4].at(i) < 130 ) {
			mhyylimit=0.999769 - 42763.9/(param[4].at(i)*param[4].at(i)) + 718.678/param[4].at(i) - 0.222259*sqrt(param[4].at(i)) - 0.0358014*param[4].at(i) + 0.000190164*param[4].at(i)*param[4].at(i);
		}
		else if ( param[4].at(i) >= 130 && param[4].at(i) < 150) {
			mhyylimit=-5.01711 - 30608.8/(param[4].at(i)*param[4].at(i)) + 716.536/param[4].at(i) - 0.191598*sqrt(param[4].at(i)) + 0.0315056*param[4].at(i) - 0.0000338228*param[4].at(i)*param[4].at(i);
		}
		else if	( param[4].at(i) >= 150 ) {
			mhyylimit=0.017654;
		}
		else {
			mhyylimit=100;
		}
		
		if ( param[3].at(i) >=128.2 && param[3].at(i) < 130 ) {
			mayylimit=0.999769 - 42763.9/(param[3].at(i)*param[3].at(i)) + 718.678/param[3].at(i) - 0.222259*sqrt(param[3].at(i)) - 0.0358014*param[3].at(i) + 0.000190164*param[3].at(i)*param[3].at(i);
		}
		else if ( param[3].at(i) >= 130 && param[3].at(i) < 150) {
			mayylimit=-5.01711 - 30608.8/(param[3].at(i)*param[3].at(i)) + 716.536/param[3].at(i) - 0.191598*sqrt(param[3].at(i)) + 0.0315056*param[3].at(i) - 0.0000338228*param[3].at(i)*param[3].at(i);
		}
		else if	( param[3].at(i) >= 150 ) {
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
	
	cout << "There are " << param[0].size() << " data points." << endl;
	
	minATLAS=99999999;
	minCMS=99999999;
	min=99999999;
	for ( int i=0; i<param[0].size(); i++) {
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
	for ( int i=0; i<param[0].size(); i++) {
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
	
	
	for ( int j=0; j<318; j++) {
		cout << "param " << j << " has size " << param[j].size() << endl;
	}
	
	
	cout << "Starting data of interest" << endl;
	ofstream data;
	data.open("dataofinterest2.txt");
	data << "f,F,tB,tG,t12,t13,sigfact,KG,KY,KL,tA,sina+b,yv,mWblh/mWSM,mh,mA,mH,,mH+,mphi,mphi+,meta,meta+,msigma,mT,mB,";
	data << "mTb5,mTb2,mT5,mT6,rDShgg,iDShgg,rShgg,iShgg,rDShyy,iDShyy,rShyy,iShyy,rDSHgg,iDSHgg,rSHgg,iSHgg,rDSHyy,iDSHyy,rSHyy,iSHyy,";
	data << "rDSAgg,iDSAgg,rSAgg,iSAgg,rDSAyy,iDSAyy,rSAyy,iSAyy,Cyy,Cgg,CWW,CZZ,Cbb,Ctata,muhyy,muAyy,muHyy,muhtata,muAtata,muHtata,muhWW,muHWW" << endl;
	for ( int i=0; i<param[0].size(); i++) {
		sg = param[7].at(i)/sqrt(1.0+param[7].at(i)*param[7].at(i));
		sB = param[5].at(i)/sqrt(1.0+param[5].at(i)*param[5].at(i));
		cB = sqrt(1.0-sB*sB);
		vblh=vsm*sqrt(1. + (1./6. + (5./4. - 2.*sg*sg)*param[0].at(i)*param[0].at(i)/(param[0].at(i)*param[0].at(i) + param[1].at(i)*param[1].at(i)))*vsm*vsm/(param[0].at(i)*param[0].at(i)));
		lambda0 = (param[2].at(i)*param[2].at(i)/(vblh*vblh))*((param[2].at(i)*param[2].at(i) - param[3].at(i)*param[3].at(i))/(param[2].at(i)*param[2].at(i) - 2.0*sB*cB* param[3].at(i)*param[3].at(i)));
		Bmu =  (lambda0*vblh*vblh + param[3].at(i)*param[3].at(i))*(2.0*sB*cB)/2.0;
		ta = (1.0/(Bmu - lambda0 *vblh*vblh*2.0*sB*cB))*((Bmu*(cB/(2.0*sB) - sB/(2.0*cB)))+sqrt( Bmu*Bmu/(4.0*sB*sB*cB*cB)   -2.0*lambda0*Bmu*vblh*vblh*2.0*sB*cB  +lambda0*lambda0*vblh*vblh*vblh*vblh*4.0*sB*sB*cB*cB     ));
		param[6].push_back(ta);
		if ( rchisqATLAS.at(i) < 1 && b_ratiotest.at(i) == true ) {
			data << param[0].at(i) << ",";
			data << param[1].at(i) << ",";
			data << param[5].at(i) << ",";
			data << param[7].at(i) << ",";
			data << param[8].at(i) << ",";
			data << param[9].at(i) << ",";
			data << param[10].at(i) << ",";
			data << param[11].at(i) << ",";
			data << param[12].at(i) << ",";
			data << param[13].at(i) << ",";
			
			data << param[6].at(i) << ",";
			data << param[38].at(i) << ",";
			data << param[63].at(i) << ",";
			data << param[64].at(i) << ",";
			
			data << param[2].at(i) << ",";
			data << param[3].at(i) << ",";
			data << param[4].at(i) << ",";
			data << param[14].at(i) << ",";
			data << param[15].at(i) << ",";
			data << param[16].at(i) << ",";
			data << param[17].at(i) << ",";
			data << param[18].at(i) << ",";
			data << param[19].at(i) << ",";
			data << param[27].at(i) << ",";
			data << param[28].at(i) << ",";
			data << param[29].at(i) << ",";
			data << param[30].at(i) << ",";
			data << param[31].at(i) << ",";
			data << param[32].at(i) << ",";
			data << param[35].at(i) << ",";
			data << param[36].at(i) << ",";
			
			
			data << param[39].at(i) << ",";
			data << param[40].at(i) << ",";
			data << param[41].at(i) << ",";
			data << param[42].at(i) << ",";
			data << param[43].at(i) << ",";
			data << param[44].at(i) << ",";
			data << param[45].at(i) << ",";
			data << param[46].at(i) << ",";
			data << param[47].at(i) << ",";
			data << param[48].at(i) << ",";
			data << param[49].at(i) << ",";
			data << param[50].at(i) << ",";
			data << param[51].at(i) << ",";
			data << param[52].at(i) << ",";
			data << param[53].at(i) << ",";
			data << param[54].at(i) << ",";
			data << param[55].at(i) << ",";
			data << param[56].at(i) << ",";
			data << param[57].at(i) << ",";
			data << param[58].at(i) << ",";
			data << param[59].at(i) << ",";
			data << param[60].at(i) << ",";
			data << param[61].at(i) << ",";
			data << param[62].at(i) << ",";
			
			data << param[127].at(i) << ",";
			data << param[129].at(i) << ",";
			data << param[125].at(i) << ",";
			data << param[126].at(i) << ",";
			data << param[122].at(i) << ",";
			data << param[124].at(i) << ",";
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
	
	
	/*cout << param[0].size() << endl;
	 cout << param[1].size() << endl;
	 cout << param[2].size() << endl;
	 cout << param[3].size() << endl;
	 cout << param[5].size() << endl;
	 cout << param[7].size() << endl;
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
	
	
	
	for ( int j=0; j<318; j++) {
		parammin[j]=9999;
		parammax[j]=-9999;
		for ( int i=0; i<param[j].size(); i++) {
			if ( param[j].at(i) < parammin[j] ) {
				parammin[j] = param[j].at(i);
			}
			else if (param[j].at(i) > parammax[j] ) {
				parammax[j] = param[j].at(i);
			}
		}
	}
	
	for ( int i=0; i<param[0].size(); i++) {
		if ( deltam.at(i) < deltammin ) {
			deltammin = deltam.at(i);
		}
		else if (deltam.at(i) > deltammax ) {
			deltammax = deltam.at(i);
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
		muZZmax = 1.5;
		muWWmax = 1.5;
		muyymax = 2.5;
		muZymax = 1.5;
		mubbmax = 1.5;
		mutatamax = 6.0;
	}
	else {
		cout << "Setting degenerate maximums." << endl;
		muZZmax = 1.5;
		muWWmax = 1.5;
		muyymax = 2.5;
		muZymax = 2.5;
		mubbmax = 1.5;
		mutatamax = 6.0;
		parammax[4] = 300;
	}
	
	cout << "MH maximum is " << parammax[4] << endl;
	
	/*for ( int i=0; i<rDSyy.size(); i++) {
	 DSyy.push_back(sqrt(rDSyy.at(i)*rDSyy.at(i)+iDSyy.at(i)*iDSyy.at(i)));
	 Syy.push_back(sqrt(rSyy.at(i)*rSyy.at(i)+iSyy.at(i)*iSyy.at(i)));
	 DSgg.push_back(sqrt(rDSgg.at(i)*rDSgg.at(i)+iDSgg.at(i)*iDSgg.at(i)));
	 Sgg.push_back(sqrt(rSgg.at(i)*rSgg.at(i)+iSgg.at(i)*iSgg.at(i)));
	 }*/
	
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

	
	
	
	for ( int i=0; i<param[0].size(); i++) {
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
	
	
	//To make a new plot, this is the parts that need to be altered.
	//If a value you want is not in a vector, you can create a new vector at any point and calculate the values.
	//The combined value in the third plot is taken as (CMS+ATLAS)/2.0 (the basic average).
	//The format for generate_plot is:
	//2igures.generate_plot("filename","xlabel","ylabel",xmin,xmax,ymin,ymax);
		
	figures.clear_plot_data();
	figures.CMSx = param[xset];
	figures.CMSy = vmuCMSZZ1blh;
	figures.ATLASx = param[xset];
	figures.ATLASy = vmuATLASZZ1blh;
	figures.COMBx = param[xset];
	figures.COMBy = muZZcombined;
	ss.str("");
	ss << "./deg/" << nametext[xset] << "vs" << nametext[yset] << ".png";
	
	if ( degen ) {
		figures.generate_plot("./deg/fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	}
	else {
		figures.generate_plot("./norm/fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	}

	
	/*
	 
	 
	 //f vs mu plots
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMSZZ1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLASZZ1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = muZZcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMSWW1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLASWW1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = muWWcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMSyy1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLASyy1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = muyycombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMSZy1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLASZy1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = muZycombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMSbb1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLASbb1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = mubbcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = vmuCMStata1blh;
	 figures.ATLASx = param[0];
	 figures.ATLASy = vmuATLAStata1blh;
	 figures.COMBx = param[0];
	 figures.COMBy = mutatacombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 
	 //parameter vs parameter plots
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = param[5];
	 figures.ATLASx = param[0];
	 figures.ATLASy = param[5];
	 figures.COMBx = param[0];
	 figures.COMBy = param[5];
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 figures.generate_plot("./norm/fvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[0];
	 figures.CMSy = deltam;
	 figures.ATLASx = param[0];
	 figures.ATLASy = deltam;
	 figures.COMBx = param[0];
	 figures.COMBy = deltam;
	 if ( degen ) {
	 figures.generate_plot("./deg/fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 figures.generate_plot("./norm/fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = deltam;
	 figures.CMSy = param[5];
	 figures.ATLASx = deltam;
	 figures.ATLASy = param[5];
	 figures.COMBx = deltam;
	 figures.COMBy = param[5];
	 if ( degen ) {
	 figures.generate_plot("./deg/deltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 figures.generate_plot("./norm/deltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 figures.generate_plot("./deg/muyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 figures.generate_plot("./norm/muyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = vmuCMSyy1blh;
	 figures.CMSy = vmuCMSWW1blh;
	 figures.ATLASx = vmuATLASyy1blh;
	 figures.ATLASy = vmuATLASWW1blh;
	 figures.COMBx = muyycombined;
	 figures.COMBy = muWWcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/muyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 figures.generate_plot("./norm/muyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = vmuCMSyy1blh;
	 figures.CMSy = vmuCMSbb1blh;
	 figures.ATLASx = vmuATLASyy1blh;
	 figures.ATLASy = vmuATLASbb1blh;
	 figures.COMBx = muyycombined;
	 figures.COMBy = mubbcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/muyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 figures.generate_plot("./norm/muyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = vmuCMSyy1blh;
	 figures.CMSy = vmuCMStata1blh;
	 figures.ATLASx = vmuATLASyy1blh;
	 figures.ATLASy = vmuATLAStata1blh;
	 figures.COMBx = muyycombined;
	 figures.COMBy = mutatacombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/muyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 figures.generate_plot("./norm/muyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = vmuCMSyy1blh;
	 figures.CMSy = vmuCMSZy1blh;
	 figures.ATLASx = vmuATLASyy1blh;
	 figures.ATLASy = vmuATLASZy1blh;
	 figures.COMBx = muyycombined;
	 figures.COMBy = muZycombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/muyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 figures.generate_plot("./norm/muyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 figures.CMSx = param[5];
	 figures.CMSy = vmuCMSZZ1blh;
	 figures.ATLASx = param[5];
	 figures.ATLASy = vmuATLASZZ1blh;
	 figures.COMBx = param[5];
	 figures.COMBy = muZZcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/tbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 figures.generate_plot("./norm/tbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[5];
	 figures.CMSy =vmuCMSWW1blh;
	 figures.ATLASx = param[5];
	 figures.ATLASy = vmuATLASWW1blh;
	 figures.COMBx = param[5];
	 figures.COMBy = muWWcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/tbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 figures.generate_plot("./norm/tbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[5];
	 figures.CMSy =vmuCMSyy1blh;
	 figures.ATLASx = param[5];
	 figures.ATLASy =vmuATLASyy1blh;
	 figures.COMBx = param[5];
	 figures.COMBy = muyycombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 figures.generate_plot("./norm/tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[5];
	 figures.CMSy = vmuCMSbb1blh;
	 figures.ATLASx = param[5];
	 figures.ATLASy = vmuATLASbb1blh;
	 figures.COMBx = param[5];
	 figures.COMBy = mubbcombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/tbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 figures.generate_plot("./norm/tbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[5];
	 figures.CMSy = vmuCMStata1blh;
	 figures.ATLASx = param[5];
	 figures.ATLASy = vmuATLAStata1blh;
	 figures.COMBx = param[5];
	 figures.COMBy = mutatacombined;
	 if ( degen ) {
	 figures.generate_plot("./deg/tbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 figures.generate_plot("./norm/tbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 
	 figures.clear_plot_data();
	 figures.CMSx = param[5];
	 figures.CMSy = param[6];
	 figures.ATLASx = param[5];
	 figures.ATLASy = param[6];
	 figures.COMBx = param[5];
	 figures.COMBy = param[6];
	 if ( degen ) {
	 figures.generate_plot("./rdeg/rtbvsta.png","tan(#beta)","tan(#alpha)",param[5]min,param[5]max,param[6]min,param[6]max);
	 }
	 else {
	 figures.generate_plot("./rnorm/rtbvsta.png","tan(#beta)","tan(#alpha)",param[5]min,param[5]max,param[6]min,param[6]max);
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
	 for ( int i=0; i<param[0].size(); i++) {
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
	 rfigures.CMSx = param[5];
	 rfigures.CMSy = param[6];
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy = param[6];
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = param[6];
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsta.png","tan(#beta)","tan(#alpha)",param[5]min,param[5]max,param[6]min,param[6]max);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsta.png","tan(#beta)","tan(#alpha)",param[5]min,param[5]max,param[6]min,param[6]max);
	 }
	 
	 
	 
	 
	 
	 //f vs mu plots
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMSZZ1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLASZZ1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = muZZcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMSWW1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLASWW1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = muWWcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMSyy1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLASyy1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = muyycombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMSZy1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLASZy1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = muZycombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMSbb1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLASbb1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = mubbcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = vmuCMStata1blh;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = vmuATLAStata1blh;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = mutatacombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 
	 //parameter vs parameter plots
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = param[5];
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = param[5];
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = param[5];
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[0];
	 rfigures.CMSy = deltam;
	 rfigures.ATLASx = param[0];
	 rfigures.ATLASy = deltam;
	 rfigures.COMBx = param[0];
	 rfigures.COMBy = deltam;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = deltam;
	 rfigures.CMSy = param[5];
	 rfigures.ATLASx = deltam;
	 rfigures.ATLASy = param[5];
	 rfigures.COMBx = deltam;
	 rfigures.COMBy = param[5];
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 rfigures.generate_plot("./rdeg/rsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 rfigures.generate_plot("./rdeg/rmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = vmuCMSyy1blh;
	 rfigures.CMSy = vmuCMSWW1blh;
	 rfigures.ATLASx = vmuATLASyy1blh;
	 rfigures.ATLASy = vmuATLASWW1blh;
	 rfigures.COMBx = muyycombined;
	 rfigures.COMBy = muWWcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = vmuCMSyy1blh;
	 rfigures.CMSy = vmuCMSbb1blh;
	 rfigures.ATLASx = vmuATLASyy1blh;
	 rfigures.ATLASy = vmuATLASbb1blh;
	 rfigures.COMBx = muyycombined;
	 rfigures.COMBy = mubbcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = vmuCMSyy1blh;
	 rfigures.CMSy = vmuCMStata1blh;
	 rfigures.ATLASx = vmuATLASyy1blh;
	 rfigures.ATLASy = vmuATLAStata1blh;
	 rfigures.COMBx = muyycombined;
	 rfigures.COMBy = mutatacombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = vmuCMSyy1blh;
	 rfigures.CMSy = vmuCMSZy1blh;
	 rfigures.ATLASx = vmuATLASyy1blh;
	 rfigures.ATLASy = vmuATLASZy1blh;
	 rfigures.COMBx = muyycombined;
	 rfigures.COMBy = muZycombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 rfigures.CMSx = param[5];
	 rfigures.CMSy = vmuCMSZZ1blh;
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy = vmuATLASZZ1blh;
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = muZZcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[5];
	 rfigures.CMSy =vmuCMSWW1blh;
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy = vmuATLASWW1blh;
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = muWWcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[5];
	 rfigures.CMSy =vmuCMSyy1blh;
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy =vmuATLASyy1blh;
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = muyycombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[5];
	 rfigures.CMSy = vmuCMSbb1blh;
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy = vmuATLASbb1blh;
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = mubbcombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 rfigures.clear_plot_data();
	 rfigures.CMSx = param[5];
	 rfigures.CMSy = vmuCMStata1blh;
	 rfigures.ATLASx = param[5];
	 rfigures.ATLASy = vmuATLAStata1blh;
	 rfigures.COMBx = param[5];
	 rfigures.COMBy = mutatacombined;
	 if ( degen ) {
	 rfigures.generate_plot("./rdeg/rtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 rfigures.generate_plot("./rnorm/rtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
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
	 for ( int i=0; i<param[0].size(); i++) {
	 if ( param[0].at(i) < 1000 ) {
	 CMSlt1++;
	 ATLASlt1++;
	 lt1++;
	 }
	 else if ( param[0].at(i) < 1500 ) {
	 CMSlt2++;
	 ATLASlt2++;
	 lt2++;
	 }
	 else if ( param[0].at(i) < 2000 ) {
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
	 ffigures.chisqCMS = param[0];
	 ffigures.chisqATLAS = param[0];
	 ffigures.chisq = param[0];
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //parameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ffvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ffvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/fdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./fdeg/fsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./fdeg/fmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/fmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/fmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/fmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/fmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/fmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ftbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ftbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ftbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ftbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ftbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ftbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ftbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ftbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./fdeg/ftbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./fnorm/ftbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
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
	 for ( int i=0; i<param[0].size(); i++) {
	 if ( !degen ) {
	 if ( param[5].at(i) < 1.25 ) {
	 CMSlt1++;
	 ATLASlt1++;
	 lt1++;
	 }
	 else if ( param[5].at(i) < 2.5 ) {
	 CMSlt2++;
	 ATLASlt2++;
	 lt2++;
	 }
	 else if ( param[5].at(i) < 3.75 ) {
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
	 if ( param[5].at(i) < 1.01 ) {
	 CMSlt1++;
	 ATLASlt1++;
	 lt1++;
	 }
	 else if ( param[5].at(i) < 1.05 ) {
	 CMSlt2++;
	 ATLASlt2++;
	 lt2++;
	 }
	 else if ( param[5].at(i) < 1.10 ) {
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
	 ffigures.bool_param[5] = true;
	 ffigures.chisqCMS = param[5];
	 ffigures.chisqATLAS = param[5];
	 ffigures.chisq = param[5];
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //parameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./tdeg/tsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./tdeg/tmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/tmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/tmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/ttbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/ttbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/ttbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/ttbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/ttbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/ttbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/ttbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/ttbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./tdeg/ttbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./tnorm/ttbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
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
	 for ( int i=0; i<param[0].size(); i++) {
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
	 ffigures.bool_param[5] = false;
	 ffigures.bool_param[3] = true;
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //parameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./mdeg/msyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/msyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./mdeg/mmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./mdeg/mtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./mnorm/mtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 
	 
	 
	 //plot ffigures;
	 //ffigures.bool_chisq = false;
	 //ffigures.bool_deg = degen;
	 
	 
	 vector<double> primary;
	 
	 for ( int i=0; i<param[0].size(); i++) {
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
	 
	 
	 
	 for ( int i=0; i<param[0].size(); i++) {
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
	 ffigures.bool_param[3] = false;
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //parameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pdeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./pdeg/psyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/psyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./pdeg/pmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLAStata1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLAStata1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLAStata1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pparam[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pparam[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/pparam[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/pparam[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,muyymin,1.3);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = hmuATLASyy1blh;
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = hmuATLASyy1blh;
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = hmuATLASyy1blh;
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/phmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/phmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/ptbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/ptbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/ptbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/ptbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/ptbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/ptbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/ptbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/ptbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./pdeg/ptbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./pnorm/ptbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
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
	 
	 
	 for ( int i=0; i<param[0].size(); i++) {
	 if ( !degen ) {
	 if ( param[5].at(i) < 1.25 ) {
	 tbdm.push_back(0.0);
	 }
	 else if ( param[5].at(i) < 2.5 ) {
	 tbdm.push_back(1.0);
	 }
	 else if ( param[5].at(i) < 3.75 ) {
	 tbdm.push_back(2.0);
	 }
	 else {
	 tbdm.push_back(3.0);
	 }
	 }
	 else {
	 if ( param[5].at(i) < 1.01 ) {
	 tbdm.push_back(0.0);
	 }
	 else if ( param[5].at(i) < 1.05 ) {
	 tbdm.push_back(1.0);
	 }
	 else if ( param[5].at(i) < 1.10 ) {
	 tbdm.push_back(2.0);
	 }
	 else {
	 tbdm.push_back(3.0);
	 }
	 }
	 
	 }
	 
	 for ( int i=0; i<param[0].size(); i++) {
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //darameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dfvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/ddeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/ddeltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./ddeg/dsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dsyparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./ddeg/dmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dmuyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dmuyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dmuyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dmuyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dmuyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLAStata1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLAStata1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLAStata1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dparam[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dparam[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dparam[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dparam[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,muyymin,1.3);
	 }
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = hmuATLASyy1blh;
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = hmuATLASyy1blh;
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = hmuATLASyy1blh;
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dhmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dhmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dtbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dtbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dtbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dtbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./ddeg/dtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./dnorm/dtbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
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
	 
	 for ( int i=0; i<param[0].size(); i++) {
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
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmuzz.png","f (GeV)","#mu_{ZZ}",param[0]min,param[0]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmuww.png","f (GeV)","#mu_{WW}",param[0]min,param[0]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSyy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASyy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmuyy.png","f (GeV)","#mu_{#gamma#gamma}",param[0]min,param[0]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmuzy.png","f (GeV)","#mu_{Z#gamma}",param[0]min,param[0]max,muZymin,muZymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmubb.png","f (GeV)","#mu_{bb}",param[0]min,param[0]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsmutata.png","f (GeV)","#mu_{#tau#tau}",param[0]min,param[0]max,mutatamin,mutatamax);
	 }
	 
	 
	 //Parameter vs parameter plots
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsparam[5].png","f (GeV)","tan#beta",param[0]min,param[0]max,param[5]min,param[5]max);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[0];
	 ffigures.CMSy = deltam;
	 ffigures.ATLASx = param[0];
	 ffigures.ATLASy = deltam;
	 ffigures.COMBx = param[0];
	 ffigures.COMBy = deltam;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2fvsdeltam.png","f (GeV)","#Delta M_{h,A} (GeV)",param[0]min,param[0]max,deltammin,deltammax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = deltam;
	 ffigures.CMSy = param[5];
	 ffigures.ATLASx = deltam;
	 ffigures.ATLASy = param[5];
	 ffigures.COMBx = deltam;
	 ffigures.COMBy = param[5];
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2deltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2deltamvsparam[5].png","#Delta M_{h,A} (GeV)","tan#beta",deltammin,deltammax,param[5]min,param[5]max);
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
	 ffigures.generate_plot("./2deg/2syparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2syparam[63]sdsyy.png","S_{#gamma#gamma}","#Delta S_{#gamma#gamma}",Syymin,Syymax,DSyymin,DSyymax);
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
	 ffigures.generate_plot("./2deg/2muyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2muyparam[63]smuZZ.png","#mu_{#gamma#gamma}","#mu_{ZZ}",muyymin,muyymax,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSWW1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2muyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2muyparam[63]smuWW.png","#mu_{#gamma#gamma}","#mu_{WW}",muyymin,muyymax,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2muyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2muyparam[63]smubb.png","#mu_{#gamma#gamma}","#mu_{bb}",muyymin,muyymax,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2muyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2muyparam[63]smutata.png","#mu_{#gamma#gamma}","#mu_{#tau#tau}",muyymin,muyymax,mutatamin,mutatamax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = vmuCMSyy1blh;
	 ffigures.CMSy = vmuCMSZy1blh;
	 ffigures.ATLASx = vmuATLASyy1blh;
	 ffigures.ATLASy = vmuATLASZy1blh;
	 ffigures.COMBx = muyycombined;
	 ffigures.COMBy = muZycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2muyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2muyparam[63]smuZy.png","#mu_{#gamma#gamma}","#mu_{Z#gamma}",muyymin,muyymax,muZymin,muZymax);
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
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLAStata1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLAStata1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLAStata1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2param[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2param[4]vsmutata.png","m_{H}","#mu^{H}_{#tau#tau}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[4];
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = param[4];
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = param[4];
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2param[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,mutatamin,1.3);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2param[4]vsmuyy.png","m_{H}","#mu^{H}_{#gamma#gamma}",param[4]min,param[4]max,muyymin,1.3);
	 }
	 
	 
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = hmuATLASyy1blh;
	 ffigures.CMSy = HmuATLASyy1blh;
	 ffigures.ATLASx = hmuATLASyy1blh;
	 ffigures.ATLASy = HmuATLASyy1blh;
	 ffigures.COMBx = hmuATLASyy1blh;
	 ffigures.COMBy = HmuATLASyy1blh;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2hmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2hmuyparam[63]sHmuyy.png","#mu^{h}_{#gamma#gamma}","#mu^{H}_{#gamma#gamma}",0,1.6,0,1.6);
	 }
	 
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSZZ1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASZZ1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muZZcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2tbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2tbvsmuzz.png","tan#beta","#mu_{ZZ}",param[5]min,param[5]max,muZZmin,muZZmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSWW1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASWW1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muWWcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2tbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2tbvsmuww.png","tan#beta","#mu_{WW}",param[5]min,param[5]max,muWWmin,muWWmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy =vmuCMSyy1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy =vmuATLASyy1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = muyycombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2tbvsmuyy.png","tan#beta","#mu_{#gamma#gamma}",param[5]min,param[5]max,muyymin,muyymax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMSbb1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLASbb1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mubbcombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2tbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2tbvsmubb.png","tan#beta","#mu_{bb}",param[5]min,param[5]max,mubbmin,mubbmax);
	 }
	 
	 ffigures.clear_plot_data();
	 ffigures.CMSx = param[5];
	 ffigures.CMSy = vmuCMStata1blh;
	 ffigures.ATLASx = param[5];
	 ffigures.ATLASy = vmuATLAStata1blh;
	 ffigures.COMBx = param[5];
	 ffigures.COMBy = mutatacombined;
	 if ( degen ) {
	 ffigures.generate_plot("./2deg/2tbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 else {
	 ffigures.generate_plot("./2norm/2tbvsmutata.png","tan#beta","#mu_{#tau#tau}",param[5]min,param[5]max,mutatamin,mutatamax);
	 }
	 */
	
	return 0;
}

#endif
