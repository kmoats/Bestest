/********************************************************************
 * style_Theory.h; Created Aug, 2006 by JPA
 * Created to make a consistent style for graphs for papers, etc
 * The style is taken from a Theory style
 ********************************************************************/

#ifndef style_Theory_h
#define style_Theory_h

#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "TROOT.h"
#include "TStyle.h"

class style_Theory : TStyle
	{
	public:
		style_Theory();
		~style_Theory();
		
		void SetStyle();
		
	};

style_Theory::style_Theory()
{
}

style_Theory::~style_Theory()
{
}

void style_Theory::SetStyle()
{
	std::cout << "Setting up Theory Styles" << std::endl;
	//
	// based on a style file from BaBar
	//
	
	//..BABAR style from RooLogon.C in workdir
	//  TStyle *TheoryStyle= new TStyle("Theory","Theory style");
	
	// use plain black on white colors
	Int_t icol=0;
	gStyle->SetLineWidth(5);
	gStyle->SetFrameLineWidth(5);
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetStatColor(icol);
	gStyle->SetFillColor(icol);
	
	// set the paper & margin sizes
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.05);
	//  gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadRightMargin(0.075);
	gStyle->SetPadBottomMargin(0.16);
	// gStyle->SetPadLeftMargin(0.12); //default
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetEndErrorSize(10);
	
	// use large fonts
	//Int_t font=72;
	Int_t font=42;
	Double_t tsize=0.05;
	gStyle->SetTextFont(font);
	
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font,"xyz");
	gStyle->SetTitleFont(font,"xyz");
	
	gStyle->SetLabelSize(tsize,"xyz");
	gStyle->SetTitleSize(tsize,"xyz");
	
	//use bold lines and markers
	//  gStyle->SetMarkerStyle(20); // This one hides the error bars
	gStyle->SetMarkerStyle(6);
	gStyle->SetMarkerSize(0.4);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
	
	// Added by JPA
	// get rid of borders 
	Int_t isize = 0;
	gStyle->SetTitleBorderSize(isize);
	gStyle->SetLegendBorderSize(isize+1);
	gStyle->SetStatBorderSize(isize+1);
	
	//set Axis style
	gStyle->SetNdivisions(410,"xyz");
	
	//get rid of X error bars and y error bar caps
	gStyle->SetErrorX(0.001);
	gStyle->SetEndErrorSize(1);
	// End JPA
	
	
	//do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
	
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	
	//gROOT->SetStyle("Plain");
	
}

#endif

