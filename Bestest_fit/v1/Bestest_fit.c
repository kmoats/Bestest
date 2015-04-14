{
	gROOT->Reset();

	gStyle->SetOptStat(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPadTickX(1);	
	gStyle->SetPadTickY(1);
	gStyle->SetPadTopMargin(0.13);
	//gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.1);
	//gStyle->SetTitleOffset(1.2,"X");
	gStyle->SetTitleOffset(1.8,"Y");
	gStyle->SetTitleOffset(1.5,"Z");
	gStyle->SetTitleFillColor(kWhite);
	//gStyle->SetTitleBorderSize(0);
	//gStyle->SetTitleX(0.16);
	//gStyle->SetTitleY(0.87);

	#include <cmath>
	#include <fstream>

	Int_t pi = 3.14159;
	
	Int_t Nmax = 99999;
	
	Double_t alphaem[Nmax], GF[Nmax], yG2[Nmax], thetaWSM[Nmax], thetaWBLH[Nmax], gSM[Nmax], gBLH[Nmax], vSM[Nmax], vBLH[Nmax], gA[Nmax], gB[Nmax], gY[Nmax], xs[Nmax], v1[Nmax], v2[Nmax], m1[Nmax], m2[Nmax], m4[Nmax], m5[Nmax], m6[Nmax], aS[Nmax], aG[Nmax], aY[Nmax], Lambda[Nmax];
	Double_t y1[Nmax], y2[Nmax], y3[Nmax], yu[Nmax], yd[Nmax], ys[Nmax], yc[Nmax], yb[Nmax], yt[Nmax], lambda0[Nmax], Bmu[Nmax], sigfact[Nmax], m56[Nmax], m65[Nmax], lambda56[Nmax], lambda65[Nmax];
	Double_t tanbeta[Nmax], tanalpha[Nmax], tanthetag[Nmax], tantheta12[Nmax], tantheta13[Nmax], f[Nmax], F[Nmax];
	Double_t Mtau[Nmax], Mu[Nmax], Md[Nmax], Ms[Nmax], Mc[Nmax], Mb[Nmax], Mt[Nmax], McPOLE[Nmax], MbPOLE[Nmax]; 
	Double_t MTb5[Nmax], MTb2[Nmax], MT6[Nmax], MT5[Nmax], MTad[Nmax], MTau[Nmax];
	Double_t Mh[Nmax], MA0[Nmax], MHc[Nmax], MH0[Nmax], Mphi0[Nmax], Mphic[Nmax], Meta0[Nmax], Metac[Nmax], Msigma[Nmax];
	Double_t Mgamma[Nmax], MZ[Nmax], MWSM[Nmax], MWBLH[Nmax], MZprime[Nmax], MWprime[Nmax];
	Double_t ReShyySM[Nmax], ImShyySM[Nmax], DeltaShyy[Nmax], ReShyy[Nmax], ImShyy[Nmax], Chyy[Nmax];
	Double_t ReShyZSM[Nmax], ImShyZSM[Nmax], DeltaShyZ[Nmax], ReShyZ[Nmax], ImShyZ[Nmax], ChyZ[Nmax];
	Double_t ReShggSM[Nmax], ImShggSM[Nmax], DeltaShgg[Nmax], ReShgg[Nmax], ImShgg[Nmax], Chgg[Nmax];
	Double_t chisq_yy[Nmax], chisq_ZZ[Nmax], chisq_WW[Nmax], chisq_bb[Nmax], chisq_tautau[Nmax], chisq_tot[Nmax]; 
	
	Double_t deltachisq[Nmax];
	
	ifstream infile;
	ofstream outfile;
	
	//open input file
	infile.open("Bestest_fit_nb_out.txt");
	if (infile.fail())
	{
		cout << "\n Could not open input file" << endl;
		exit(1);
	}
	
	//open output file
	outfile.open("Bestest_fit_c_out.txt");
	if (outfile.fail())
	{
		cout << "\n Could not open output file" << endl;
		exit(1);
	}
	
	
	//read data from data file
	Int_t i = 0;
	while (!infile.eof())	
	{
		infile >> alphaem[i] >> GF[i] >> yG2[i] >> thetaWSM[i] >> thetaWBLH[i] >> gSM[i] >> gBLH[i] >> vSM[i] >> vBLH[i] >> gA[i] >> gB[i] >> gY[i] >> xs[i] >> v1[i] >> v2[i] >> m1[i] >> m2[i] >> m4[i] >> m5[i] >> m6[i] >> aS[i] >> aG[i] >> aY[i] >> Lambda[i];
		infile >> y1[i] >> y2[i] >> y3[i] >> yu[i] >> yd[i] >> ys[i] >> yc[i] >> yb[i] >> yt[i] >> lambda0[i] >> Bmu[i] >> sigfact[i] >> m56[i] >> m65[i] >> lambda56[i] >> lambda65[i] >> tanbeta[i] >> tanalpha[i] >> tanthetag[i] >> tantheta12[i] >> tantheta13[i] >> f[i] >> F[i]; 
		infile >> Mtau[i] >> Mu[i] >> Md[i] >> Ms[i] >> Mc[i] >> Mb[i] >> Mt[i] >> McPOLE[i] >> MbPOLE[i] >> MTb5[i] >> MTb2[i] >> MT6[i] >> MT5[i] >> MTad[i] >> MTau[i] >> Mh[i] >> MA0[i] >> MHc[i] >> MH0[i] >> Mphi0[i] >> Mphic[i] >> Meta0[i] >> Metac[i] >> Msigma[i] >> Mgamma[i] >> MZ[i] >> MWSM[i] >> MWBLH[i] >> MZprime[i] >> MWprime[i]; 
		infile >> ReShyySM[i] >> ImShyySM[i] >> DeltaShyy[i] >> ReShyy[i] >> ImShyy[i] >> Chyy[i];
		infile >> ReShyZSM[i] >> ImShyZSM[i] >> DeltaShyZ[i] >> ReShyZ[i] >> ImShyZ[i] >> ChyZ[i];
		infile >> ReShggSM[i] >> ImShggSM[i] >> DeltaShgg[i] >> ReShgg[i] >> ImShgg[i] >> Chgg[i];
		infile >> chisq_yy[i] >> chisq_ZZ[i] >> chisq_WW[i] >> chisq_bb[i] >> chisq_tautau[i] >> chisq_tot[i];

		i++;
	
	}
	
	Int_t N = i;
	cout << endl << "Number of data points = " << N << endl << endl;
	outfile << endl << "Number of data points = " << N << endl << endl;
	infile.close();
	

	Double_t chisqmin = 100000;
	
	Double_t Chyyfit, Chggfit;

	Double_t Chyy1sigmalow = 100000;
	Double_t Chyy1sigmahigh = 0;
	Double_t Chgg1sigmalow = 100000;
	Double_t Chgg1sigmahigh = 0;
	
	Double_t Chyy2sigmalow = 100000;
	Double_t Chyy2sigmahigh = 0;
	Double_t Chgg2sigmalow = 100000;
	Double_t Chgg2sigmahigh = 0;
	
	Double_t Chyy3sigmalow = 100000;
	Double_t Chyy3sigmahigh = 0;
	Double_t Chgg3sigmalow = 100000;
	Double_t Chgg3sigmahigh = 0;
	
	Double_t Chyy1sigmaminus, Chyy1sigmaplus, Chyy2sigmaminus, Chyy2sigmaplus, Chyy3sigmaminus, Chyy3sigmaplus;
	Double_t Chgg1sigmaminus, Chgg1sigmaplus, Chgg2sigmaminus, Chgg2sigmaplus, Chgg3sigmaminus, Chgg3sigmaplus;	

	Double_t deltachisq1sigma = 1.00;
	Double_t deltachisq2sigma = 3.84;
	Double_t deltachisq3sigma = 9.00;
	
//	Double_t deltachisq1sigma = 2.30;
//	Double_t deltachisq2sigma = 5.99;
//	Double_t deltachisq3sigma = 11.83;
	
	
	for(Int_t i = 0; i < N; i++)
	{		
		if(chisq_tot[i] < chisqmin)
		{
			chisqmin = chisq_tot[i];
			Chyyfit = Chyy[i];
			Chggfit = Chgg[i];
		}
	}	
	
	for(Int_t i = 0; i < N; i++)
	{		
		deltachisq[i] = chisq_tot[i]-chisqmin;
		
		if(deltachisq[i] <= deltachisq1sigma)
		{
			
			if(Chyy[i] < Chyy1sigmalow)
			{
				Chyy1sigmalow = Chyy[i];
			}	
			
			if(Chyy[i] > Chyy1sigmahigh)
			{
				Chyy1sigmahigh = Chyy[i];
			}

			if(Chgg[i] < Chgg1sigmalow)
			{
				Chgg1sigmalow = Chgg[i];
			}	
			
			if(Chgg[i] > Chgg1sigmahigh)
			{
				Chgg1sigmahigh = Chgg[i];
			}
		}
		
		if(deltachisq[i] <= deltachisq2sigma)
		{
			
			if(Chyy[i] < Chyy2sigmalow)
			{
				Chyy2sigmalow = Chyy[i];
			}	
			
			if(Chyy[i] > Chyy2sigmahigh)
			{
				Chyy2sigmahigh = Chyy[i];
			}
			
			if(Chgg[i] < Chgg2sigmalow)
			{
				Chgg2sigmalow = Chgg[i];
			}	
			
			if(Chgg[i] > Chgg2sigmahigh)
			{
				Chgg2sigmahigh = Chgg[i];
			}
		}

		if(deltachisq[i] <= deltachisq3sigma)
		{
			
			if(Chyy[i] < Chyy3sigmalow)
			{
				Chyy3sigmalow = Chyy[i];
			}	
			
			if(Chyy[i] > Chyy3sigmahigh)
			{
				Chyy3sigmahigh = Chyy[i];
			}
			
			if(Chgg[i] < Chgg3sigmalow)
			{
				Chgg3sigmalow = Chgg[i];
			}	
			
			if(Chgg[i] > Chgg3sigmahigh)
			{
				Chgg3sigmahigh = Chgg[i];
			}
		}		
	}	
	
	Chyy1sigmaminus = Chyyfit - Chyy1sigmalow;
	Chgg1sigmaminus = Chggfit - Chgg1sigmalow;

	Chyy1sigmaplus = Chyy1sigmahigh - Chyyfit;
	Chgg1sigmaplus = Chgg1sigmahigh - Chggfit;
	
	Chyy2sigmaminus = Chyyfit - Chyy2sigmalow;
	Chgg2sigmaminus = Chggfit - Chgg2sigmalow;
	
	Chyy2sigmaplus = Chyy2sigmahigh - Chyyfit;
	Chgg2sigmaplus = Chgg2sigmahigh - Chggfit;
	
	Chyy3sigmaminus = Chyyfit - Chyy3sigmalow;
	Chgg3sigmaminus = Chggfit - Chgg3sigmalow;
	
	Chyy3sigmaplus = Chyy3sigmahigh - Chyyfit;
	Chgg3sigmaplus = Chgg3sigmahigh - Chggfit;
	
	cout << "Minimum chi-squared = " << chisqmin  << endl << endl;
	cout << "Chyy = " << Chyyfit <<endl;
	cout << "+1 sigma = " << Chyy1sigmaplus <<  ", -1 sigma = " << Chyy1sigmaminus << endl;
	cout << "+2 sigma = " << Chyy2sigmaplus <<  ", -2 sigma = " << Chyy2sigmaminus << endl;
	cout << "+3 sigma = " << Chyy3sigmaplus <<  ", -3 sigma = " << Chyy3sigmaminus << endl << endl;
	
	cout << "Chgg = " << Chggfit <<endl;
	cout << "+1 sigma = " << Chgg1sigmaplus <<  ", -1 sigma = " << Chgg1sigmaminus << endl;
	cout << "+2 sigma = " << Chgg2sigmaplus <<  ", -2 sigma = " << Chgg2sigmaminus << endl;
	cout << "+3 sigma = " << Chgg3sigmaplus <<  ", -3 sigma = " << Chgg3sigmaminus << endl << endl;
	
	outfile << "Minimum chi-squared = " << chisqmin  << endl << endl;
	outfile << "Chyy = " << Chyyfit <<endl;
	outfile << "+1 sigma = " << Chyy1sigmaplus <<  ", -1 sigma = " << Chyy1sigmaminus << endl;
	outfile << "+2 sigma = " << Chyy2sigmaplus <<  ", -2 sigma = " << Chyy2sigmaminus << endl;
	outfile << "+3 sigma = " << Chyy3sigmaplus <<  ", -3 sigma = " << Chyy3sigmaminus << endl << endl;
	
	outfile << "Chgg = " << Chggfit <<endl;
	outfile << "+1 sigma = " << Chgg1sigmaplus <<  ", -1 sigma = " << Chgg1sigmaminus << endl;
	outfile << "+2 sigma = " << Chgg2sigmaplus <<  ", -2 sigma = " << Chgg2sigmaminus << endl;
	outfile << "+3 sigma = " << Chgg3sigmaplus <<  ", -3 sigma = " << Chgg3sigmaminus << endl << endl;
	
	
	Int_t n_canvas = 3;
	Char_t *canvas_title[] = {"parameters","chisq","DeltaS & C"};
	
	Int_t n_pad_x[] = {4,4,4,4,4,4};
	Int_t n_pad_y[] = {2,2,2,2,2,2};
	Int_t n_pads[] = {8,8,8,8,8,8};
//	Int_t n_pads[n_canvas];
	Int_t canvas_width[n_canvas];
	Int_t canvas_height[n_canvas];
	
	Int_t n_pads_max = 0;
	
	for(Int_t canvas = 0; canvas < n_canvas; canvas++)
	{
//		n_pads[canvas] = n_pad_x[canvas]*n_pad_y[canvas];
		canvas_width[canvas]  = 425*n_pad_x[canvas];
		canvas_height[canvas] = 400*n_pad_y[canvas];
		
		if(n_pads[canvas] > n_pads_max)
		{
			n_pads_max = n_pads[canvas];
		}
	}
		
	Char_t *pad_title[][n_pads_max] = { {"tan#beta vs. f","tan#beta vs. M_{A}","M_{A} vs. f","tan#theta_{13} vs. tan#theta_{12}","tan#alpha vs. tan#beta","","",""},
										{"#Delta#chi^{2} vs. f","#Delta#chi^{2} vs. M_{A}","#Delta#chi^{2} vs. tan#beta","#Delta#chi^{2} vs. F","#Delta#chi^{2} vs. #DeltaS_{h#gamma#gamma}","#Delta#chi^{2} vs. #DeltaS_{hgg}","#Delta#chi^{2} vs. C_{h#gamma#gamma}","#Delta#chi^{2} vs. C_{hgg}"},
										{"#DeltaS_{hgg} vs. #DeltaS_{h#gamma#gamma}","C_{hgg} vs. C_{h#gamma#gamma}","C_{h#gamma#gamma} vs. f","C_{hgg} vs. f","","","",""}
									  };
	Char_t *x_title[][n_pads_max] = {	{"f (GeV)","M_{A} (GeV)","f (GeV)", "tan#theta_{12}", "tan#beta","","",""},
										{"f (GeV)","M_{A} (GeV)","tan#beta","F (GeV)", "#DeltaS_{h#gamma#gamma}","#DeltaS_{hgg}","C_{h#gamma#gamma}","C_{hgg}"},
										{"#DeltaS_{h#gamma#gamma}","C_{h#gamma#gamma}","f (GeV)","f (GeV)","","","",""}
									};
	Char_t *y_title[][n_pads_max] = {	{"tan#beta","tan#beta","M_{A} (GeV)", "tan#theta_{13}", "tan#alpha","","",""},
										{"#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}","#Delta#chi^{2}"},
										{"#DeltaS_{hgg}","C_{hgg}","C_{h#gamma#gamma}","C_{hgg}","","","",""}
									};
	
	Int_t n_plots = 5;
	TGraph *graph[n_canvas][n_pads_max][n_plots+1];
	
	
	Double_t x_var[n_canvas][n_pads_max][N];
	Double_t y_var[n_canvas][n_pads_max][N];

	
	for(Int_t canvas = 0; canvas < n_canvas; canvas++)
	{
		for(Int_t pad = 0; pad < n_pads[canvas]; pad++)
		{
			for(Int_t plot = 0; plot < n_plots; plot++)
			{
				graph[canvas][pad][plot] = new TGraph(N);
				graph[canvas][pad][plot]->SetMarkerStyle(6);
				graph[canvas][pad][plot]->SetMarkerColor(plot+2);
			}
		}
		
		for(Int_t i = 0; i < N; i++)
		{
			x_var[0][0][i] = f[i];			
			y_var[0][0][i] = tanbeta[i];

			x_var[0][1][i] = MA0[i];
			y_var[0][1][i] = tanbeta[i];
			
			x_var[0][2][i] = f[i];
			y_var[0][2][i] = MA0[i];
			
			x_var[0][3][i] = tantheta12[i];
			y_var[0][3][i] = tantheta13[i];
			
			x_var[0][4][i] = tanbeta[i];
			y_var[0][4][i] = tanalpha[i];
			
			x_var[1][0][i] = f[i];			
			y_var[1][0][i] = deltachisq[i];
			
			x_var[1][1][i] = MA0[i];			
			y_var[1][1][i] = deltachisq[i];

			x_var[1][2][i] = tanbeta[i];			
			y_var[1][2][i] = deltachisq[i];	
			
			x_var[1][3][i] = F[i];			
			y_var[1][3][i] = deltachisq[i];	
			
			x_var[1][4][i] = DeltaShyy[i];			
			y_var[1][4][i] = deltachisq[i];	
			
			x_var[1][5][i] = DeltaShgg[i];			
			y_var[1][5][i] = deltachisq[i];	
			
			x_var[1][6][i] = Chyy[i];			
			y_var[1][6][i] = deltachisq[i];	

			x_var[1][7][i] = Chgg[i];			
			y_var[1][7][i] = deltachisq[i];	
			
			x_var[2][0][i] = DeltaShyy[i];			
			y_var[2][0][i] = DeltaShgg[i];
			
			x_var[2][1][i] = Chyy[i];			
			y_var[2][1][i] = Chgg[i];
			
			x_var[2][2][i] = f[i];			
			y_var[2][2][i] = Chyy[i];
			
			x_var[2][3][i] = f[i];			
			y_var[2][3][i] = Chgg[i];

			
		}
		
	}
	
	Double_t x_min[][n_pads_max] = {{500,125,500,0,1},
									{500,125,1,500,-1,-1,0.5,0.5},
									{-1,0.5,500,500}
								  };
	
	Double_t x_max[][n_pads_max] = {{3000,500,3000,2,7},
									{3000,500,7,5000,1,1,1.5,1.5},
									{1,1.5,3000,3000}
								  };
	
	Double_t y_min[][n_pads_max] = {{1,1,125,0,-2},
									{0,0,0,0,0,0,0,0},
									{-0.1,0.5,0.5,0.5}
								  };
	
	Double_t y_max[][n_pads_max] = {{7,7,500,2,2},
									{15,15,15,15,15,15,15,15},
									{0.1,1.5,1.5,1.5}
								  };
		
	for(Int_t canvas = 0; canvas < n_canvas; canvas++)
	{
		TCanvas *C = new TCanvas(canvas_title[canvas],canvas_title[canvas],canvas_width[canvas],canvas_height[canvas]);
		C->Divide(n_pad_x[canvas],n_pad_y[canvas]);
		
		for(Int_t pad = 0; pad < n_pads[canvas]; pad++)
		{
			C->cd(pad+1);
		
			for(Int_t i = 0; i < N; i++)
			{				
				if(lambda0[i] < 4*pi && lambda56[i] < 4*pi && lambda65[i] < 4*pi && y1[i] < 4*pi && y2[i] < 4*pi && y3[i] < 4*pi && gA[i] < 4*pi && gB[i] < 4*pi && gY[i] < 4*pi)
				{
					graph[canvas][pad][0]->SetPoint(i,x_var[canvas][pad][i],y_var[canvas][pad][i]);

					if(chisq_tot[i] < chisqmin + deltachisq3sigma)
					{
						graph[canvas][pad][1]->SetPoint(i,x_var[canvas][pad][i],y_var[canvas][pad][i]);
					
						if(chisq_tot[i] < chisqmin + deltachisq2sigma)
						{
							graph[canvas][pad][2]->SetPoint(i,x_var[canvas][pad][i],y_var[canvas][pad][i]);
					
							if(chisq_tot[i] < chisqmin + deltachisq1sigma)
							{
								graph[canvas][pad][3]->SetPoint(i,x_var[canvas][pad][i],y_var[canvas][pad][i]);
							}
						}
					}
				}
			}
			
			if(canvas == 2 && pad == 1)			
			{
				graph[canvas][pad][4]->SetPoint(1,Chyyfit,Chggfit);
				graph[canvas][pad][4]->SetMarkerStyle(20);
				graph[canvas][pad][4]->SetMarkerColor(kBlack);
			}
		
			graph[canvas][pad][0]->SetTitle(Form("%s",pad_title[canvas][pad]));
			graph[canvas][pad][0]->GetXaxis()->SetTitle(Form("%s",x_title[canvas][pad]));
			graph[canvas][pad][0]->GetYaxis()->SetTitle(Form("%s",y_title[canvas][pad]));
			
			graph[canvas][pad][0]->GetXaxis()->SetRangeUser(x_min[canvas][pad],x_max[canvas][pad]);	
			graph[canvas][pad][0]->GetYaxis()->SetRangeUser(y_min[canvas][pad],y_max[canvas][pad]);	
			
			graph[canvas][pad][0]->Draw("AP");
			graph[canvas][pad][1]->Draw("P");
			graph[canvas][pad][2]->Draw("P");
			graph[canvas][pad][3]->Draw("P");
			graph[canvas][pad][4]->Draw("P");
			
			gPad->RedrawAxis();
			
			legend = new TLegend(0.75,0.67,0.86,0.84);
			legend->SetFillColor(kWhite);
		//	legend->SetBorderSize(0);
		//	legend->SetEntrySeparation(0.24);
		//	legend->Clear();
			
			legend->AddEntry(graph[canvas][pad][3],"< 1#sigma","p");
			legend->AddEntry(graph[canvas][pad][2],"< 2#sigma","p");
			legend->AddEntry(graph[canvas][pad][1],"< 3#sigma","p");
			legend->AddEntry(graph[canvas][pad][0],"> 3#sigma","p");
			
			legend->Draw();
			

			
		}
	}
	
}