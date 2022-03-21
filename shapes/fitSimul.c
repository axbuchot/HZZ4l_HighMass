void fitSimul(int year=2018){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TString chanName[3]={"4mu","4e","2e2mu"};
	//int mass [5]={125,124,120,126,130};
  int mass [18]={200,250,300,350,400,450,500,550,600,700,750,800,900,1000,1500,2000,2500,3000}

    const int nstage=4; //24;
//      TString stageName [nstage]={"VBF_Rest", "VBF_2j_mjj_GT700_2j", "VBF_GT200", "VBF_2j_mjj_350_700_2j", "VBF_2j_mjj_GT350_3j", "VH_Had", "VH_Lep_0_150", "VH_Lep_GT150", "BBH", "ggH_1j_0_60", "ggH_0j_10_200", "ggH_0j_0_10", "ggH_VBF", "ggH_1j_60_120", "ggH_1j_120_200", "ggH_2j_60_120", "ggH_GT200", "ggH_2j_120_200", "ggH_2j_0_60", "TH", "TTH","VBF_GT200_VHori","VBF_Rest_VHori","VH_Had_VBFori"};
//      int stage [nstage]= {201, 209, 226, 207, 208, 212, 301, 303, 701, 107, 106, 105, 113, 108, 109, 111, 150, 112, 110, 801, 601, 231, 206, 217};   

    // Split VH in WH, ZH
    TString stageName [nstage]={"ggH","VBFH_1jet","VBFH_2jet","inclusive"};

    int stage [nstage]= {201, 202, 203, 204};

	RooRealVar *ZZMass =new RooRealVar("ZZMass","",150,3500);
	RooRealVar *stagev=new RooRealVar("highmass_categories","",0,1000);
	RooRealVar *chan=new RooRealVar("chan","",0,1000);
	RooRealVar *weight=new RooRealVar("weight","",1.e-9,1.e9);

	//TString inputTrees = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/ReducedTrees/";
  TString inputTrees = '/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/2018/VBFH1000/VBFH1000_reducedTree_MC_2018.root'
	TString path = "./";

	// std::vector<string> sigName;
	// sigName.push_back("VBFH");  sigName.push_back("WminusH"); 
	// sigName.push_back("WplusH"); sigName.push_back("ZH"); 
	// sigName.push_back("bbH"); sigName.push_back("ggH"); sigName.push_back("ttH");

	for (int j =0;j<3;j++){
		ofstream outmass ;
		outmass.open(Form(path + "shape/sim_massParam_highmass_categoriest%s%d.txt",chanName[j].Data(),year));
		RooDataSet *data[5];
		for (int i =0;i<5;i++){
			TChain *t = new TChain ("SelectedTree");
			t->Add(Form(inputTrees + "*H%d_redTree_*%d*.root",mass[i],year));
			// Uncomment to use only TTrees for which we have all the samples (mass variated and tuneup/down)
    		// for(int l = 0; l < sigName.size(); l++) {
      //       		TString name = sigName.at(l);
      //       		t->Add(Form(inputTrees+name+"%d_redTree_*%d*.root", mass[i],year));
      //       		cout << sigName.at(l) << " added to TChain" << endl;
    		// }

			data[i]=new RooDataSet(Form("data%d_%d",i,mass[i]),"",t,RooArgSet(*ZZMass,*weight,*chan,*stagev),Form("chan==%d",j+1),"weight");
		}

		for (int k =0;k<nstage;k++){
			//	MH->setVal(mass[i]);
			//	MH->setConstant(kTRUE);
			RooCategory massrc("massrc","");
			RooSimultaneous simPdf("simPdf","simultaneous pdf",massrc) ;
			RooDataSet* data_stage[5]; 
			RooDataSet* data_stage_cat[5]; 
			RooDataSet* data_all;
			RooDoubleCBFast* cpdf[5]; 
			RooAddPdf* addpdf[5]; 

			for (int i=0;i<5;i++){
				TString cat_type =Form("mh%d",mass[i]);
				massrc.defineType(cat_type,mass[i]);
			}
			RooRealVar *a1_1=new RooRealVar("a1_1","",0.00111,-0.2,0.2);
			RooRealVar *a1_2=new RooRealVar("a1_2","",0.00146,-0.1,0.1);
			RooRealVar *a1_0=new RooRealVar("a1_0","",1.2917,1.,5.);
			RooRealVar *a2_1=new RooRealVar("a2_1","",0.00130,-0.2,0.2);
			RooRealVar *a2_2=new RooRealVar("a2_2","",0.00127,-0.1,0.1);
			RooRealVar *a2_0=new RooRealVar("a2_0","",1.8125,1.,5.);
			RooRealVar *n1_1=new RooRealVar("n1_1","",0.00138,-0.2,0.2);
			RooRealVar *n1_2=new RooRealVar("n1_2","",0.00185,-0.1,0.1);
			RooRealVar *n1_0=new RooRealVar("n1_0","",2.0312,1.,7.);
			RooRealVar *n2_1=new RooRealVar("n2_1","",0.001,-0.2,0.2);
			RooRealVar *n2_2=new RooRealVar("n2_2","",-0.001,-0.1,0.1);
			RooRealVar *n2_0=new RooRealVar("n2_0","",3.0369,1.,5.);

			RooRealVar *mean_1=new RooRealVar("mean_1","",1,0,1.5);
			RooRealVar *mean_2=new RooRealVar("mean_2","",0.99945,-1.,1.);
			RooRealVar *mean_0=new RooRealVar("mean_0","",124.8424,120.,130);

			RooRealVar *sigma_1=new RooRealVar("sigma_1","",0.00860,-0.5,0.5);
			RooRealVar *sigma_2=new RooRealVar("sigma_2","",0.00026,-1.,1.);
			RooRealVar *sigma_0=new RooRealVar("sigma_0","",1.5,0.5,3.);

			RooRealVar* mean_l_0 = new RooRealVar("landau_mean_0","",130,110,140);
			RooRealVar* mean_l_1 = new RooRealVar("landau_mean_1","",0,-1.5,1.5);
			RooRealVar* mean_l_2 = new RooRealVar("landau_mean_2","",0,-1,1);
			RooRealVar* sigma_l_0 = new RooRealVar("landau_sigma_0","",15,2,20);
			RooRealVar* sigma_l_1 = new RooRealVar("landau_sigma_1","",0.,-1,1);
			RooRealVar* sigma_l_2 = new RooRealVar("landau_sigma_2","",0,-1,1);
			RooRealVar* frac_0 = new RooRealVar("frac_0","",0.65,0,1);
			RooRealVar* frac_1 = new RooRealVar("frac_1","",-0.1,0.1);
			RooRealVar* frac_2 = new RooRealVar("frac_2","",-0.1,0.1);
			for (int i=0;i<5;i++){
				RooConstVar *MH=new RooConstVar("MH","",mass[i]);
				// As of AN: the CB parameters are extracted by a fit p = a1+a11(mh-125)+a12(mh-125)**2
				// RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*a1_0,*a1_1,*a1_2,*MH));
				// RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*a2_0,*a2_1,*a2_2,*MH));
				// RooFormulaVar* n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*n1_0,*n1_1,*n1_2,*MH));
				// RooFormulaVar* n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*n2_0,*n2_1,*n2_2,*MH));
				RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*a1_0,*a1_1,*MH));
				RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*a2_0,*a2_1,*MH));
				RooFormulaVar* n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*n1_0,*n1_1,*MH));
				RooFormulaVar* n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*n2_0,*n2_1,*MH));
				RooFormulaVar* mean=new RooFormulaVar(Form("mean_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*mean_0,*mean_1,*MH));
				RooFormulaVar* sigma=new RooFormulaVar(Form("sigma_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*sigma_0,*sigma_1,*MH));
				RooFormulaVar* sigma_l=new RooFormulaVar(Form("sigma_l_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*sigma_l_0,*sigma_l_1,*MH));
				RooFormulaVar* mean_l=new RooFormulaVar(Form("mean_l_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*mean_l_0,*mean_l_1,*MH));
				RooFormulaVar* frac=new RooFormulaVar(Form("frac_%d",mass[i]),"","@0",RooArgList(*frac_0));
				}

				TString cat_type =Form("mh%d",mass[i]);
				data_stage[i] = (RooDataSet*)data[i]->reduce(Form("highmass_categories==%d",stage[k])); 
				data_stage_cat[i] = new RooDataSet(Form("data_%d",i),"",RooArgSet(*ZZMass,*weight),Index(massrc),Import(cat_type,*data_stage[i]),WeightVar("weight")); 
				cpdf[i] = new RooDoubleCBFast(Form("DCBall_%d",i),"",*ZZMass,*mean,*sigma,*a1,*n1,*a2,*n2);
				RooLandau *pdf_landau = new RooLandau("landau","",*ZZMass, *mean_l, *sigma_l);
				addpdf[i] = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdf[i],*pdf_landau),*frac);
				cout<< data_stage_cat[i]->sumEntries()<<endl;
				if (stageName[k].Contains("TTH")|| stageName[k].Contains("VH_Lep") || stageName[k].Contains("WH_Lep") || stageName[k].Contains("ZH_Lep")) // WH, ZH commented
					simPdf.addPdf(*addpdf[i],cat_type); // Crystal Ball + Landau in TTH and VH to cope with the rising tail
				else
					simPdf.addPdf(*cpdf[i],cat_type); // Crystal Ball in all categories
				if(i==0)
					data_all = data_stage_cat[i];
				else
					data_all->append(*data_stage_cat[i]);
				if(i==0){ // mH = 125 GeV
					// Quoting from analysis note: "The initial value for the parmsCB0 is obtained by fitting 125 GeV mass sample alone,
					// and refitted in the simultaneous fits. " 
					RooDoubleCBFast *cpdftmp = new RooDoubleCBFast("DCBall","",*ZZMass,*mean_0,*sigma_0,*a1_0,*n1_0,*a2_0,*n2_0);
					RooLandau *pdf_landautmp = new RooLandau("landautmp","",*ZZMass, *mean_l_0, *sigma_l_0);
					RooAddPdf *addpdftmp = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdftmp,*pdf_landautmp),*frac);
					if (stageName[k].Contains("TTH")|| stageName[k].Contains("VH_Lep") || stageName[k].Contains("WH_Lep") || stageName[k].Contains("ZH_Lep"))
						addpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					else
						cpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					//a1_0->setConstant(1);
					//a2_0->setConstant(1);
					//n1_0->setConstant(1);
					//n2_0->setConstant(1);
					// The CB and landau means, as well as CB sigma are set to constants for the simultaneous fit, 
					//using the parameters obtained from 125GeV fit.
					sigma_0->setConstant(1);
					mean_0->setConstant(1);
					mean_l_0->setConstant(1);
					//sigma_l_0->setConstant(1);
				}
			}
			simPdf.Print("v");
			data_all->Print("v");
			//			MH->setConstant(1);
			simPdf.fitTo(*data_all,InitialHesse(true),Strategy(2)) ;
			TString a1outs= Form ("a1_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,a1_0->getVal(),a1_1->getVal());
			TString a2outs= Form ("a2_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,a2_0->getVal(),a2_1->getVal());
			TString n1outs= Form ("n1_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,n1_0->getVal(),n1_1->getVal());
			TString n2outs= Form ("n2_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,n2_0->getVal(),n2_1->getVal());
			TString meanouts= Form ("mean_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,mean_0->getVal(),mean_1->getVal());
			TString sigmaouts= Form ("sigma_%s_%s%d \t %.4f+%.5f*(MH-125)",stageName[k].Data(),chanName[j].Data(),year,sigma_0->getVal(),sigma_1->getVal());

			outmass<< a1outs<<endl;
			outmass<< a2outs<<endl;
			outmass<< meanouts<<endl;
			outmass<< sigmaouts<<endl;
			outmass<< n1outs<<endl;
			outmass<< n2outs<<endl;

			cout<<"where1"<<endl;
			for (int i =0;i<5;i++){
				RooPlot *frame = ZZMass->frame();
				TString cat_type =Form("mh%d",mass[i]);
				data_all->plotOn(frame,Cut(Form("massrc==massrc::mh%d",mass[i]))) ;
				simPdf.plotOn(frame,Slice(massrc,cat_type),ProjWData(massrc,*data_all)) ;
				frame->Draw();
				gPad->Print(Form(path + "shape/fig_%d/simFit_%s_%d_%s_%d.png", year, stageName[k].Data(),mass[i],chanName[j].Data(),year));
			}
		}
		outmass.close();
	}
}

