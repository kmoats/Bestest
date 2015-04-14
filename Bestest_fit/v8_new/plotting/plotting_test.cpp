#ifndef plotting_cpp
#define plotting_cpp

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "plotting.h"
#include "style.h"

using namespace std;

plot::~plot() {
}


plot::plot() {
	style_Theory myStyle;
	myStyle.SetStyle();
	bool_chisq=true;
	bool_f=false;
	bool_tB=false;
	bool_ma0=false;
	bool_deg=false;
	bool_prim=false;
	bool_double=false;
};

void plot::clear_plot_data() {
	CMSx.clear();
	CMSy.clear();
	ATLASx.clear();
	ATLASy.clear();
	COMBx.clear();
	COMBy.clear();
}

void plot::generate_plot(std::string plotname, std::string xtext, std::string ytext, double xmin, double xmax, double ymin, double ymax) {
	
	for ( int c=0; c<1; c++) {
		if ( bool_chisq ) {
			canv = new TCanvas("canv","Canvas",100,100,2400,1200);
			canv->SetBorderSize(0);
			canv->SetFillColor(10);
			canv->SetBorderMode(0);
			pad1 = new TPad("pad1","lhepad1",0,0,0.5,1.0,10);
			pad2 = new TPad("pad2","lhepad2",0.5,0,1.0,1.0,10);
			//	pad1 = new TPad("pad1","lhepad1",0,0,0.333,1.0,10);
			//	pad2 = new TPad("pad2","lhepad2",0.333,0,0.667,1.0,10);
			//pad3 = new TPad("pad3","lhepad3",0.667,0,1.0,1.0,10);
			pad1->Draw();
			pad2->Draw();
			//pad3->Draw();
		}
		else if ( !bool_double )  {
			canv = new TCanvas("canv","Canvas",100,100,1200,1200);
			canv->SetBorderSize(0);
			canv->SetFillColor(10);
			canv->SetBorderMode(0);
			pad1 = new TPad("pad1","lhepad1",0,0,1.0,1.0,10);
			pad1->Draw();
		}
		else {
			canv = new TCanvas("canv","Canvas",100,100,1500,1200);
			canv->SetBorderSize(0);
			canv->SetFillColor(10);
			canv->SetBorderMode(0);
			pad1 = new TPad("pad1","lhepad1",0,0,0.8,1.0,10);
			pad1->Draw();
			
		}
		
		
		h_axis = new TH1D("axis","axis",1,xmin,xmax);
		h_axis->GetXaxis()->SetLabelSize(0.04);
		h_axis->GetXaxis()->SetNdivisions(508);
		h_axis->GetYaxis()->SetLabelSize(0.04);
		h_axis->GetXaxis()->SetTitleSize(0.07);
		h_axis->GetYaxis()->SetTitleSize(0.07);
		h_axis->GetXaxis()->SetTitleOffset(1.1);
		h_axis->GetYaxis()->SetTitleOffset(1.1);
		h_axis->GetXaxis()->SetLabelColor(1);
		h_axis->GetXaxis()->SetAxisColor(1);
		h_axis->GetXaxis()->SetTitleColor(1);
		h_axis->GetYaxis()->SetLabelColor(1);
		h_axis->GetYaxis()->SetAxisColor(1);
		h_axis->GetYaxis()->SetTitleColor(1);
		h_axis->GetXaxis()->SetTitle(xtext.c_str());
		h_axis->GetYaxis()->SetTitle(ytext.c_str());
		h_axis->SetBinContent(1,ymax);
		h_axis->SetLineColor(10);
		h_axis->GetYaxis()->SetRangeUser(ymin,ymax);
		
		
		gCMS1 = new TGraph(CMSlt1);
		gCMS2 = new TGraph(CMSlt2);
		gCMS3 = new TGraph(CMSlt3);
		gCMS4 = new TGraph(CMSgt3);
		gATLAS1 = new TGraph(ATLASlt1);
		gATLAS2 = new TGraph(ATLASlt2);
		gATLAS3 = new TGraph(ATLASlt3);
		gATLAS4 = new TGraph(ATLASgt3);
		g1 = new TGraph(lt1);
		g2 = new TGraph(lt2);
		g3 = new TGraph(lt3);
		g4 = new TGraph(gt3);
		
		gc0 = new TGraph(20000);
		gc1 = new TGraph(20000);
		gc2 = new TGraph(20000);
		gc3 = new TGraph(20000);
		gc4 = new TGraph(20000);
		gc5 = new TGraph(20000);
		gc6 = new TGraph(20000);
		gc7 = new TGraph(20000);
		gc8 = new TGraph(20000);
		gc9 = new TGraph(20000);
		gc10 = new TGraph(20000);
		gc11 = new TGraph(20000);
		gc12 = new TGraph(20000);
		gc13 = new TGraph(20000);
		gc14 = new TGraph(20000);
		gc15 = new TGraph(20000);
		
		reject = new TGraph(20000);
		reject->SetMarkerColor(kRed-10);
		reject->SetFillColor(kRed-10);
		reject->SetMarkerStyle(7);

		
		if ( c==0 ) {
			if ( bool_chisq ) {
				gCMS1->SetMarkerColor(kGreen+1);
				gCMS1->SetFillColor(kGreen+1);
				gCMS1->SetMarkerStyle(6);
				gCMS2->SetMarkerColor(kOrange);
				gCMS2->SetFillColor(kOrange);
				gCMS2->SetMarkerStyle(6);
				gCMS3->SetMarkerColor(kOrange+10);
				gCMS3->SetFillColor(kOrange+10);
				gCMS3->SetMarkerStyle(6);
				gCMS4->SetMarkerColor(kGray);
				gCMS4->SetFillColor(kGray);
				gCMS4->SetMarkerStyle(6);
				gATLAS1->SetMarkerColor(kGreen+1);
				gATLAS1->SetFillColor(kGreen+1);
				gATLAS1->SetMarkerStyle(6);
				gATLAS2->SetMarkerColor(kOrange);
				gATLAS2->SetFillColor(kOrange);
				gATLAS2->SetMarkerStyle(6);
				gATLAS3->SetMarkerColor(kOrange+10);
				gATLAS3->SetFillColor(kOrange+10);
				gATLAS3->SetMarkerStyle(6);
				gATLAS4->SetMarkerColor(kGray);
				gATLAS4->SetFillColor(kGray);
				gATLAS4->SetMarkerStyle(6);
				g1->SetMarkerColor(kGreen+1);
				g1->SetFillColor(kGreen+1);
				g1->SetMarkerStyle(6);
				g2->SetMarkerColor(kOrange);
				g2->SetFillColor(kOrange);
				g2->SetMarkerStyle(6);
				g3->SetMarkerColor(kOrange+10);
				g3->SetFillColor(kOrange+10);
				g3->SetMarkerStyle(6);
				g4->SetMarkerColor(kGray);
				g4->SetFillColor(kGray);
				g4->SetMarkerStyle(6);
			}
			else if ( !bool_double && !bool_prim2 ) {
				gCMS4->SetMarkerColor(kOrange-2);
				gCMS4->SetFillColor(kOrange-2);
				gCMS4->SetMarkerStyle(6);
				gCMS3->SetMarkerColor(kOrange+8);
				gCMS3->SetFillColor(kOrange+8);
				gCMS3->SetMarkerStyle(6);
				gCMS2->SetMarkerColor(kViolet+8);
				gCMS2->SetFillColor(kViolet+8);
				gCMS2->SetMarkerStyle(6);
				gCMS1->SetMarkerColor(kAzure+8);
				gCMS1->SetFillColor(kAzure+8);
				gCMS1->SetMarkerStyle(6);
				gATLAS4->SetMarkerColor(kOrange-2);
				gATLAS4->SetFillColor(kOrange-2);
				gATLAS4->SetMarkerStyle(6);
				gATLAS3->SetMarkerColor(kOrange+8);
				gATLAS3->SetFillColor(kOrange+8);
				gATLAS3->SetMarkerStyle(6);
				gATLAS2->SetMarkerColor(kViolet+8);
				gATLAS2->SetFillColor(kViolet+8);
				gATLAS2->SetMarkerStyle(6);
				gATLAS1->SetMarkerColor(kAzure+8);
				gATLAS1->SetFillColor(kAzure+8);
				gATLAS1->SetMarkerStyle(6);
				g4->SetMarkerColor(kOrange-2);
				g4->SetFillColor(kOrange-2);
				g4->SetMarkerStyle(6);
				g3->SetMarkerColor(kOrange+8);
				g3->SetFillColor(kOrange+8);
				g3->SetMarkerStyle(6);
				g2->SetMarkerColor(kViolet+8);
				g2->SetFillColor(kViolet+8);
				g2->SetMarkerStyle(6);
				g1->SetMarkerColor(kAzure+8);
				g1->SetFillColor(kAzure+8);
				g1->SetMarkerStyle(6);
			}
			else if ( bool_double ) {
				gc0->SetMarkerColor(kOrange+3);
				gc0->SetFillColor(kOrange+3);
				gc0->SetMarkerStyle(6);
				gc1->SetMarkerColor(kOrange-6);
				gc1->SetFillColor(kOrange-6);
				gc1->SetMarkerStyle(6);
				gc2->SetMarkerColor(kOrange-3);
				gc2->SetFillColor(kOrange-3);
				gc2->SetMarkerStyle(6);
				gc3->SetMarkerColor(kOrange);
				gc3->SetFillColor(kOrange);
				gc3->SetMarkerStyle(6);
				
				gc4->SetMarkerColor(kRed+3);
				gc4->SetFillColor(kRed+3);
				gc4->SetMarkerStyle(6);
				gc5->SetMarkerColor(kRed+1);
				gc5->SetFillColor(kRed+1);
				gc5->SetMarkerStyle(6);
				gc6->SetMarkerColor(kRed-4);
				gc6->SetFillColor(kRed-4);
				gc6->SetMarkerStyle(6);
				gc7->SetMarkerColor(kRed-9);
				gc7->SetFillColor(kRed-9);
				gc7->SetMarkerStyle(6);
				
				gc8->SetMarkerColor(kViolet+3);
				gc8->SetFillColor(kViolet+3);
				gc8->SetMarkerStyle(6);
				gc9->SetMarkerColor(kViolet-6);
				gc9->SetFillColor(kViolet-6);
				gc9->SetMarkerStyle(6);
				gc10->SetMarkerColor(kViolet-3);
				gc10->SetFillColor(kViolet-3);
				gc10->SetMarkerStyle(6);
				gc11->SetMarkerColor(kViolet);
				gc11->SetFillColor(kViolet);
				gc11->SetMarkerStyle(6);
				
				gc12->SetMarkerColor(kAzure+3);
				gc12->SetFillColor(kAzure+3);
				gc12->SetMarkerStyle(6);
				gc13->SetMarkerColor(kAzure-6);
				gc13->SetFillColor(kAzure-6);
				gc13->SetMarkerStyle(6);
				gc14->SetMarkerColor(kAzure-3);
				gc14->SetFillColor(kAzure-3);
				gc14->SetMarkerStyle(6);
				gc15->SetMarkerColor(kCyan-4);
				gc15->SetFillColor(kCyan-4);
				gc15->SetMarkerStyle(6);
			}
			else if ( bool_prim2 ) {
				gc0->SetMarkerColor(kOrange+3);
				gc0->SetFillColor(kOrange+3);
				gc0->SetMarkerStyle(6);
				gc1->SetMarkerColor(kOrange-3);
				gc1->SetFillColor(kOrange-3);
				gc1->SetMarkerStyle(6);
				gc2->SetMarkerColor(kOrange+10);
				gc2->SetFillColor(kOrange+10);
				gc2->SetMarkerStyle(6);
				
				gc3->SetMarkerColor(kSpring-7);
				gc3->SetFillColor(kSpring+-7);
				gc3->SetMarkerStyle(6);
				gc4->SetMarkerColor(kSpring-5);
				gc4->SetFillColor(kSpring-5);
				gc4->SetMarkerStyle(6);
				gc5->SetMarkerColor(kSpring+10);
				gc5->SetFillColor(kSpring+10);
				gc5->SetMarkerStyle(6);

				gc6->SetMarkerColor(kAzure+3);
				gc6->SetFillColor(kAzure+3);
				gc6->SetMarkerStyle(6);
				gc7->SetMarkerColor(kAzure-3);
				gc7->SetFillColor(kAzure-3);
				gc7->SetMarkerStyle(6);
				gc8->SetMarkerColor(kAzure+10);
				gc8->SetFillColor(kAzure+10);
				gc8->SetMarkerStyle(6);
			}
		}
		else {
			plotname_start = plotname.substr(0,plotname.find(".",1));
			plotname_end = plotname.substr(plotname.find(".",1),plotname.size());
			plotname = plotname_start + "_bw" + plotname_end;
			if ( bool_chisq ) {
				gCMS1->SetMarkerColor(kBlack);
				gCMS1->SetFillColor(kBlack);
				gCMS1->SetMarkerStyle(6);
				gCMS2->SetMarkerColor(kGray+2);
				gCMS2->SetFillColor(kGray+2);
				gCMS2->SetMarkerStyle(6);
				gCMS3->SetMarkerColor(kGray+1);
				gCMS3->SetFillColor(kGray+1);
				gCMS3->SetMarkerStyle(6);
				gCMS4->SetMarkerColor(kGray);
				gCMS4->SetFillColor(kGray);
				gCMS4->SetMarkerStyle(6);
				gATLAS1->SetMarkerColor(kBlack);
				gATLAS1->SetFillColor(kBlack);
				gATLAS1->SetMarkerStyle(6);
				gATLAS2->SetMarkerColor(kGray+2);
				gATLAS2->SetFillColor(kGray+2);
				gATLAS2->SetMarkerStyle(6);
				gATLAS3->SetMarkerColor(kGray+1);
				gATLAS3->SetFillColor(kGray+1);
				gATLAS3->SetMarkerStyle(6);
				gATLAS4->SetMarkerColor(kGray);
				gATLAS4->SetFillColor(kGray);
				gATLAS4->SetMarkerStyle(6);
				g1->SetMarkerColor(kBlack);
				g1->SetFillColor(kBlack);
				g1->SetMarkerStyle(6);
				g2->SetMarkerColor(kGray+2);
				g2->SetFillColor(kGray+2);
				g2->SetMarkerStyle(6);
				g3->SetMarkerColor(kGray+1);
				g3->SetFillColor(kGray+1);
				g3->SetMarkerStyle(6);
				g4->SetMarkerColor(kGray);
				g4->SetFillColor(kGray);
				g4->SetMarkerStyle(6);
			}
			else {
				gCMS4->SetMarkerColor(kGray);
				gCMS4->SetFillColor(kGray);
				gCMS4->SetMarkerStyle(6);
				gCMS3->SetMarkerColor(kGray+1);
				gCMS3->SetFillColor(kGray+1);
				gCMS3->SetMarkerStyle(6);
				gCMS2->SetMarkerColor(kGray+2);
				gCMS2->SetFillColor(kGray+2);
				gCMS2->SetMarkerStyle(6);
				gCMS1->SetMarkerColor(kGray+3);
				gCMS1->SetFillColor(kGray+3);
				gCMS1->SetMarkerStyle(6);
				gATLAS4->SetMarkerColor(kGray);
				gATLAS4->SetFillColor(kGray);
				gATLAS4->SetMarkerStyle(6);
				gATLAS3->SetMarkerColor(kGray+1);
				gATLAS3->SetFillColor(kGray+1);
				gATLAS3->SetMarkerStyle(6);
				gATLAS2->SetMarkerColor(kGray+2);
				gATLAS2->SetFillColor(kGray+2);
				gATLAS2->SetMarkerStyle(6);
				gATLAS1->SetMarkerColor(kGray+3);
				gATLAS1->SetFillColor(kGray+3);
				gATLAS1->SetMarkerStyle(6);
				g4->SetMarkerColor(kGray);
				g4->SetFillColor(kGray);
				g4->SetMarkerStyle(6);
				g3->SetMarkerColor(kGray+1);
				g3->SetFillColor(kGray+1);
				g3->SetMarkerStyle(6);
				g2->SetMarkerColor(kGray+2);
				g2->SetFillColor(kGray+2);
				g2->SetMarkerStyle(6);
				g1->SetMarkerColor(kGray+3);
				g1->SetFillColor(kGray+3);
				g1->SetMarkerStyle(6);
			}
		}
		
		
		
		//cout << "Plots created and set up" << endl;
		
		iCMSlt1=0;
		iCMSlt2=0;
		iCMSlt3=0;
		iCMSgt3=0;
		iATLASlt1=0;
		iATLASlt2=0;
		iATLASlt3=0;
		iATLASgt3=0;
		ilt1=0;
		ilt2=0;
		ilt3=0;
		igt3=0;
		ic0=0;
		ic1=0;
		ic2=0;
		ic3=0;
		ic4=0;
		ic5=0;
		ic6=0;
		ic7=0;
		ic8=0;
		ic9=0;
		ic10=0;
		ic11=0;
		ic12=0;
		ic13=0;
		ic14=0;
		ic15=0;
		
		for ( int i=0; i<CMSx.size(); i++) {
			if ( ratiotest.at(i) ) {
				if ( bool_chisq ) {
					if ( chisqCMS.at(i)-minCMS < 1 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i)-minCMS < 4 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i)-minCMS < 9 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i)-minATLAS < 1 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i)-minATLAS < 4 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i)-minATLAS < 9 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i)-min < 1 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i)-min < 4 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i)-min < 9 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_f ) {
					if ( chisqCMS.at(i) < 1000 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 1500 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 2000 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) < 1000 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 1500 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 2000 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) < 1000 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 1500 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 2000 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_tB && !bool_deg ) {
					if ( chisqCMS.at(i) < 1.25 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 2.5 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 3.75 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) < 1.25 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 2.5 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 3.75 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) < 1.25 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 2.5 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 3.75 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_tB && bool_deg ) {
					if ( chisqCMS.at(i) < 1.01 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 1.05 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 1.1 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) < 1.01 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 1.05 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 1.1 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) < 1.01 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 1.05 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 1.1 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_ma0 && !bool_deg ) {
					if ( chisqCMS.at(i) < 25 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 100 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 200 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) < 25 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 100 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 200 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) < 25 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 100 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 200 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_ma0 && bool_deg ) {
					if ( chisqCMS.at(i) < 0.5 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 1 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) < 1.5 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) < 0.5 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 1 ) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) < 1.5 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) < 0.5 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 1 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) < 1.5 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_prim ) {
					if ( chisqCMS.at(i) == 1.0 ) {
						iCMSlt1++;
						gCMS1->SetPoint(iCMSlt1,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) == 2.0 ) {
						iCMSlt2++;
						gCMS2->SetPoint(iCMSlt2,CMSx.at(i),CMSy.at(i));
					}
					else if ( chisqCMS.at(i) == 3.0 ) {
						iCMSlt3++;
						gCMS3->SetPoint(iCMSlt3,CMSx.at(i),CMSy.at(i));
					}
					else {
						iCMSgt3++;
						gCMS4->SetPoint(iCMSgt3,CMSx.at(i),CMSy.at(i));
					}
					if ( chisqATLAS.at(i) == 1.0 ) {
						iATLASlt1++;
						gATLAS1->SetPoint(iATLASlt1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) == 2.0) {
						iATLASlt2++;
						gATLAS2->SetPoint(iATLASlt2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisqATLAS.at(i) == 3.0 ) {
						iATLASlt3++;
						gATLAS3->SetPoint(iATLASlt3,ATLASx.at(i),ATLASy.at(i));
					}
					else {
						iATLASgt3++;
						gATLAS4->SetPoint(iATLASgt3,ATLASx.at(i),ATLASy.at(i));
					}
					if ( chisq.at(i) == 1.0 ) {
						ilt1++;
						g1->SetPoint(ilt1,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) == 2.0 ) {
						ilt2++;
						g2->SetPoint(ilt2,COMBx.at(i),COMBy.at(i));
					}
					else if ( chisq.at(i) == 3.0 ) {
						ilt3++;
						g3->SetPoint(ilt3,COMBx.at(i),COMBy.at(i));
					}
					else {
						igt3++;
						g4->SetPoint(igt3,COMBx.at(i),COMBy.at(i));
					}
				}
				else if ( bool_double || bool_prim2 ) {
					if ( chisq.at(i) == 0.0 ) {
						ic0++;
						gc0->SetPoint(ic0,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 1.0 ) {
						ic1++;
						gc1->SetPoint(ic1,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 2.0 ) {
						ic2++;
						gc2->SetPoint(ic2,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 3.0 ) {
						ic3++;
						gc3->SetPoint(ic3,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 4.0 ) {
						ic4++;
						gc4->SetPoint(ic4,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 5.0 ) {
						ic5++;
						gc5->SetPoint(ic5,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 6.0 ) {
						ic6++;
						gc6->SetPoint(ic6,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 7.0 ) {
						ic7++;
						gc7->SetPoint(ic7,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 8.0 ) {
						ic8++;
						gc8->SetPoint(ic8,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 9.0 ) {
						ic9++;
						gc9->SetPoint(ic9,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 10.0 ) {
						ic10++;
						gc10->SetPoint(ic10,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 11.0 ) {
						ic11++;
						gc11->SetPoint(ic11,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 12.0 ) {
						ic12++;
						gc12->SetPoint(ic12,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 13.0 ) {
						ic13++;
						gc13->SetPoint(ic13,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 14.0 ) {
						ic14++;
						gc14->SetPoint(ic14,ATLASx.at(i),ATLASy.at(i));
					}
					else if ( chisq.at(i) == 15.0 ) {
						ic15++;
						gc15->SetPoint(ic15,ATLASx.at(i),ATLASy.at(i));
					}
				}
			}
			else {
				ireject++;
				reject->SetPoint(ireject,ATLASx.at(i),ATLASy.at(i));
			}
		}
		
		pt1 = new TPaveText(0.75,0.2,0.90,0.25, "NDC");
		pt1->SetFillColor(4000);
		pt1->SetTextSize(0.05);
		pt1->SetTextAlign(12);
		pt1->SetBorderSize(0);
		pt1->AddText("CMS");
		pt2 = new TPaveText(0.75,0.2,0.90,0.25, "NDC");
		pt2->SetFillColor(4000);
		pt2->SetTextSize(0.05);
		pt2->SetTextAlign(12);
		pt2->SetBorderSize(0);
		pt2->AddText("ATLAS");
		//pt3 = new TPaveText(0.75,0.2,0.90,0.25, "NDC");
		//pt3->SetFillColor(4000);
		//pt3->SetTextSize(0.05);
		//pt3->SetTextAlign(12);
		//pt3->SetBorderSize(0);
		//pt3->AddText("COMB");
		
		if ( !bool_double ) {
			leg1 = new TLegend(0.72,0.8,0.92,0.95);
			leg2 = new TLegend(0.72,0.8,0.92,0.95);
			leg3 = new TLegend(0.72,0.8,0.92,0.95);
		}
		else {
			leg1 = new TLegend(0.95,0.5,1.2,0.95);
			leg2 = new TLegend(0.72,0.5,0.92,0.95);
			leg3 = new TLegend(0.72,0.5,0.92,0.95);
		}
		
		if ( bool_chisq ) {
			leg1->AddEntry(gCMS1,"#Delta#chi^{2} < 1","f");
			leg1->AddEntry(gCMS2,"#Delta#chi^{2} < 4","f");
			leg1->AddEntry(gCMS3,"#Delta#chi^{2} < 9","f");
			leg1->AddEntry(gCMS4,"#Delta#chi^{2} > 9","f");
			leg2->AddEntry(gATLAS1,"#Delta#chi^{2} < 1","f");
			leg2->AddEntry(gATLAS2,"#Delta#chi^{2} < 4","f");
			leg2->AddEntry(gATLAS3,"#Delta#chi^{2} < 9","f");
			leg2->AddEntry(gATLAS4,"#Delta#chi^{2} > 9","f");
			leg3->AddEntry(g1,"#Delta#chi^{2} < 1","f");
			leg3->AddEntry(g2,"#Delta#chi^{2} < 4","f");
			leg3->AddEntry(g3,"#Delta#chi^{2} < 9","f");
			leg3->AddEntry(g4,"#Delta#chi^{2} > 9","f");
		}
		else if ( bool_f ) {
			leg1->AddEntry(gCMS1,"f < 1000","f");
			leg1->AddEntry(gCMS2,"f < 1500","f");
			leg1->AddEntry(gCMS3,"f < 2000","f");
			leg1->AddEntry(gCMS4,"f > 2000","f");
			leg2->AddEntry(gATLAS1,"f < 1000","f");
			leg2->AddEntry(gATLAS2,"f < 1500","f");
			leg2->AddEntry(gATLAS3,"f < 2000","f");
			leg2->AddEntry(gATLAS4,"f > 2000","f");
			leg3->AddEntry(g1,"f < 1000","f");
			leg3->AddEntry(g2,"f < 1500","f");
			leg3->AddEntry(g3,"f < 2000","f");
			leg3->AddEntry(g4,"f > 2000","f");
		}
		else if ( bool_tB && !bool_deg ) {
			leg1->AddEntry(gCMS1,"tan#beta < 1.25","f");
			leg1->AddEntry(gCMS2,"tan#beta < 2.5","f");
			leg1->AddEntry(gCMS3,"tan#beta < 3.75","f");
			leg1->AddEntry(gCMS4,"tan#beta > 3.75","f");
			leg2->AddEntry(gATLAS1,"tan#beta < 1.25","f");
			leg2->AddEntry(gATLAS2,"tan#beta < 2.5","f");
			leg2->AddEntry(gATLAS3,"tan#beta < 3.75","f");
			leg2->AddEntry(gATLAS4,"tan#beta > 3.75","f");
			leg3->AddEntry(g1,"tan#beta < 1.25","f");
			leg3->AddEntry(g2,"tan#beta < 2.5","f");
			leg3->AddEntry(g3,"tan#beta < 3.75","f");
			leg3->AddEntry(g4,"tan#beta > 3.75","f");
		}
		else if ( bool_tB && bool_deg ) {
			leg1->AddEntry(gCMS1,"tan#beta < 1.01","f");
			leg1->AddEntry(gCMS2,"tan#beta < 1.05","f");
			leg1->AddEntry(gCMS3,"tan#beta < 1.10","f");
			leg1->AddEntry(gCMS4,"tan#beta > 1.10","f");
			leg2->AddEntry(gATLAS1,"tan#beta < 1.01","f");
			leg2->AddEntry(gATLAS2,"tan#beta < 1.05","f");
			leg2->AddEntry(gATLAS3,"tan#beta < 1.10","f");
			leg2->AddEntry(gATLAS4,"tan#beta > 1.10","f");
			leg3->AddEntry(g1,"tan#beta < 1.01","f");
			leg3->AddEntry(g2,"tan#beta < 1.05","f");
			leg3->AddEntry(g3,"tan#beta < 1.10","f");
			leg3->AddEntry(g4,"tan#beta > 1.10","f");
		}
		else if ( bool_ma0 && !bool_deg ) {
			leg1->AddEntry(gCMS1,"#Delta M < 25","f");
			leg1->AddEntry(gCMS2,"#Delta M < 100","f");
			leg1->AddEntry(gCMS3,"#Delta M < 200","f");
			leg1->AddEntry(gCMS4,"#Delta M > 200","f");
			leg2->AddEntry(gATLAS1,"#Delta M < 25","f");
			leg2->AddEntry(gATLAS2,"#Delta M < 100","f");
			leg2->AddEntry(gATLAS3,"#Delta M < 200","f");
			leg2->AddEntry(gATLAS4,"#Delta M > 200","f");
			leg3->AddEntry(g1,"#Delta M < 25","f");
			leg3->AddEntry(g2,"#Delta M < 100","f");
			leg3->AddEntry(g3,"#Delta M < 200","f");
			leg3->AddEntry(g4,"#Delta M > 200","f");
		}
		else if ( bool_ma0 && bool_deg ) {
			leg1->AddEntry(gCMS1,"#Delta M < 0.5","f");
			leg1->AddEntry(gCMS2,"#Delta M < 1","f");
			leg1->AddEntry(gCMS3,"#Delta M < 1.5","f");
			leg1->AddEntry(gCMS4,"#Delta M > 1.5","f");
			leg2->AddEntry(gATLAS1,"#Delta M < 0.5","f");
			leg2->AddEntry(gATLAS2,"#Delta M < 1","f");
			leg2->AddEntry(gATLAS3,"#Delta M < 1.5","f");
			leg2->AddEntry(gATLAS4,"#Delta M > 1.5","f");
			leg3->AddEntry(g1,"#Delta M < 0.5","f");
			leg3->AddEntry(g2,"#Delta M < 1","f");
			leg3->AddEntry(g3,"#Delta M < 1.5","f");
			leg3->AddEntry(g4,"#Delta M > 1.5","f");
		}
		else if ( bool_prim ) {
			leg1->AddEntry(gCMS1,"h > #gamma#gamma","f");
			leg1->AddEntry(gCMS2,"H > #gamma#gamma","f");
			leg1->AddEntry(gCMS3,"A > #gamma#gamma","f");
			leg2->AddEntry(gATLAS1,"h > #gamma#gamma","f");
			leg2->AddEntry(gATLAS2,"H > #gamma#gamma","f");
			leg2->AddEntry(gATLAS3,"A > #gamma#gamma","f");
			leg3->AddEntry(g1,"h > #gamma#gamma","f");
			leg3->AddEntry(g2,"H > #gamma#gamma","f");
			leg3->AddEntry(g3,"A > #gamma#gamma","f");
		}
		else if ( bool_double && bool_deg ) {
			leg1->AddEntry(gc0,"tan#beta < 1.01 & #Delta M < 0.5","f");
			leg1->AddEntry(gc1,"tan#beta < 1.05 & #Delta M < 0.5","f");
			leg1->AddEntry(gc2,"tan#beta < 1.10 & #Delta M < 0.5","f");
			leg1->AddEntry(gc3,"tan#beta > 1.10 & #Delta M < 0.5","f");
			leg1->AddEntry(gc4,"tan#beta < 1.01 & #Delta M < 1.0","f");
			leg1->AddEntry(gc5,"tan#beta < 1.05 & #Delta M < 1.0","f");
			leg1->AddEntry(gc6,"tan#beta < 1.10 & #Delta M < 1.0","f");
			leg1->AddEntry(gc7,"tan#beta > 1.10 & #Delta M < 1.0","f");
			leg1->AddEntry(gc8,"tan#beta < 1.01 & #Delta M < 1.5","f");
			leg1->AddEntry(gc9,"tan#beta < 1.05 & #Delta M < 1.5","f");
			leg1->AddEntry(gc10,"tan#beta < 1.10 & #Delta M < 1.5","f");
			leg1->AddEntry(gc11,"tan#beta > 1.10 & #Delta M < 1.5","f");
			leg1->AddEntry(gc12,"tan#beta < 1.01 & #Delta M > 1.5","f");
			leg1->AddEntry(gc13,"tan#beta < 1.05 & #Delta M > 1.5","f");
			leg1->AddEntry(gc14,"tan#beta < 1.10 & #Delta M > 1.5","f");
			leg1->AddEntry(gc15,"tan#beta > 1.10 & #Delta M > 1.5","f");
		}
		else if ( bool_double && !bool_deg ) {
			leg1->AddEntry(gc0,"tan#beta < 1.25 & #Delta M < 25","f");
			leg1->AddEntry(gc1,"tan#beta < 2.5 & #Delta M < 25","f");
			leg1->AddEntry(gc2,"tan#beta < 3.75 & #Delta M < 25","f");
			leg1->AddEntry(gc3,"tan#beta > 3.75 & #Delta M < 25","f");
			leg1->AddEntry(gc4,"tan#beta < 1.25 & #Delta M < 100","f");
			leg1->AddEntry(gc5,"tan#beta < 2.5 & #Delta M < 100","f");
			leg1->AddEntry(gc6,"tan#beta < 3.75 & #Delta M < 100","f");
			leg1->AddEntry(gc7,"tan#beta > 3.75 & #Delta M < 100","f");
			leg1->AddEntry(gc8,"tan#beta < 1.25 & #Delta M < 200","f");
			leg1->AddEntry(gc9,"tan#beta < 2.5 & #Delta M < 200","f");
			leg1->AddEntry(gc10,"tan#beta < 3.75 & #Delta M < 200","f");
			leg1->AddEntry(gc11,"tan#beta > 3.75 & #Delta M < 200","f");
			leg1->AddEntry(gc12,"tan#beta < 1.25 & #Delta M > 200","f");
			leg1->AddEntry(gc13,"tan#beta < 2.5 & #Delta M > 200","f");
			leg1->AddEntry(gc14,"tan#beta < 3.75 & #Delta M > 200","f");
			leg1->AddEntry(gc15,"tan#beta > 3.75 & #Delta M > 200","f");
		}
		else if ( bool_prim2 && bool_deg ) {
			if ( ic0 > 0 ) {
				leg1->AddEntry(gc0,"h > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic1 > 0 ) {
				leg1->AddEntry(gc1,"H > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic2 > 0 ) {
				leg1->AddEntry(gc2,"A > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic3 > 0 ) {
				leg1->AddEntry(gc3,"h > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic4 > 0 ) {
				leg1->AddEntry(gc4,"H > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic5 > 0 ) {
				leg1->AddEntry(gc5,"A > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic6 > 0 ) {
				leg1->AddEntry(gc6,"h > #gamma#gamma, A > #tau#tau","f");
			}
			if ( ic7 > 0 ) {
				leg1->AddEntry(gc7,"H > #gamma#gamma, A > #tau#tau","f");
			}
			if ( ic8 > 0 ) {
				leg1->AddEntry(gc8,"A > #gamma#gamma, A > #tau#tau","f");
			}
		}
		else if ( bool_prim2 && !bool_deg ) {
			if ( ic0 > 0 ) {
				leg1->AddEntry(gc0,"h > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic1 > 0 ) {
				leg1->AddEntry(gc1,"H > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic2 > 0 ) {
				leg1->AddEntry(gc2,"A > #gamma#gamma, h > #tau#tau","f");
			}
			if ( ic3 > 0 ) {
				leg1->AddEntry(gc3,"h > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic4 > 0 ) {
				leg1->AddEntry(gc4,"H > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic5 > 0 ) {
				leg1->AddEntry(gc5,"A > #gamma#gamma, H > #tau#tau","f");
			}
			if ( ic6 > 0 ) {
				leg1->AddEntry(gc6,"h > #gamma#gamma, A > #tau#tau","f");
			}
			if ( ic7 > 0 ) {
				leg1->AddEntry(gc7,"H > #gamma#gamma, A > #tau#tau","f");
			}
			if ( ic8 > 0 ) {
				leg1->AddEntry(gc8,"A > #gamma#gamma, A > #tau#tau","f");
			}
		}
		
		
		
		if ( bool_chisq ) {
			pad1->cd();
			h_axis->Draw("AXIS");
			reject->Draw("P same");
			gCMS4->Draw("P same");
			gCMS3->Draw("P same");
			gCMS2->Draw("P same");
			gCMS1->Draw("P same");
			pt1->Draw();
			leg1->Draw();
			pad2->cd();
			h_axis->Draw("AXIS");
			reject->Draw("P same");
			gATLAS4->Draw("P same");
			gATLAS3->Draw("P same");
			gATLAS2->Draw("P same");
			gATLAS1->Draw("P same");
			pt2->Draw();
			leg2->Draw();
			//	pad3->cd();
			//	h_axis->Draw("AXIS");
			//	g4->Draw("P same");
			//	g3->Draw("P same");
			//	g2->Draw("P same");
			//	g1->Draw("P same");
			//	pt3->Draw();
			//	leg3->Draw();
		}
		else if ( !bool_double && !bool_prim2 ) {
			pad1->cd();
			h_axis->Draw("AXIS");
			reject->Draw("P same");
			gATLAS4->Draw("P same");
			gATLAS3->Draw("P same");
			gATLAS2->Draw("P same");
			gATLAS1->Draw("P same");
			//pt2->Draw();
			leg2->Draw();
		}
		else if ( bool_double ) {
			pad1->cd();
			h_axis->Draw("AXIS");
			reject->Draw("P same");
			gc0->Draw("P same");
			gc1->Draw("P same");
			gc2->Draw("P same");
			gc3->Draw("P same");
			gc4->Draw("P same");
			gc5->Draw("P same");
			gc6->Draw("P same");
			gc7->Draw("P same");
			gc8->Draw("P same");
			gc9->Draw("P same");
			gc10->Draw("P same");
			gc11->Draw("P same");
			gc12->Draw("P same");
			gc13->Draw("P same");
			gc14->Draw("P same");
			gc15->Draw("P same");
			leg1->Draw();
		}
		else if ( bool_prim2 ) {
			pad1->cd();
			h_axis->Draw("AXIS");
			reject->Draw("P same");
			gc8->Draw("P same");
			gc7->Draw("P same");
			gc6->Draw("P same");
			gc0->Draw("P same");
			gc1->Draw("P same");
			gc2->Draw("P same");
			gc3->Draw("P same");
			gc4->Draw("P same");
			gc5->Draw("P same");
			leg1->Draw();
		}
		
		
		
		canv->Update();
		canv->SaveAs(plotname.c_str());
		canv->Clear();
		
		
	}
}

#endif