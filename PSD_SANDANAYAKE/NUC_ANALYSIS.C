#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TH1.h"
#include "Riostream.h"
#include "TMath.h"
#include "TSystem.h"
#include <stdlib.h>
#include "TCanvas.h"
#include <ctype.h>
using namespace std;
//function to fix files and count the no. of lines
int fixer(TString filename)
{
    //initial nametag creation
    TString filename_fixed    = filename + TString("_fixed.txt");
    TString filename_new      = filename + TString(".txt");
    ifstream read_filename;     read_filename.open(filename); 
    rename(filename,filename_new);
    ofstream ff (filename_fixed);
    //fix to fill the gap at 2*J=<<_>>X and N=  <<_>>X, and make new files
    FILE* file = fopen(filename_new, "r");
    char line[256];
    char space = ' ';
    char zero = '0';
    int count =0; //count no. of lines
    while (fgets(line, sizeof(line), file))//make necessary replacements
    {
       if (line[16] ==  space){line[16] = zero;}
           line[15] =   space;
           line[14] =   '=';
           line[13] =   'j';
           line[12] =   '2';
       if (line[27] == space){line[27] = zero;} 
       ff << line;
       ++count;
    }
    ff.close();
    fclose(file);
    return count;
}

//histogram creator for the level density plots
TH1F* plotter(TString filename, double binwidth, double Edif, int Jstat, int Jval)
{
    TString filename_fixed    = filename + TString("_fixed.txt");
    ifstream read;
    read.open(filename_fixed);
    if (!read){cout << "the file "<< filename <<" does not exist" << endl;}   
    int J,N_eig,C, N,Z;
    double EX,E,E_MAX=0,SUM;
    TString NZ, TJ, P, CEQ, NEQ, EXC, EQ, EEQ;
    //reading from file
    while (read >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {if(EX>E_MAX){E_MAX=EX;}}
    read.close();
    E_MAX = ceil(E_MAX+Edif);
    //The histogram
    TH1F* rho = new TH1F(filename + TString(Form("%d",Jstat)) + TString(Form("%d",Jval)), filename, (E_MAX)/binwidth, 0, E_MAX);
    gStyle->SetOptStat(kFALSE); 
    //reading from file
    read.open(filename_fixed);
    while (read >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {
        if(Jstat==0)        {rho-> Fill(EX+Edif);}
        else if (Jstat==1)  {if(J/2 == Jval){rho-> Fill(EX+Edif);}} 
    }
    read.close();
    //write data to new text files
    ofstream ss; 
    if(Jstat==0)        {ss.open( filename + TString("_hist.txt") ); }
    else if(Jstat==1)   {ss.open( filename + TString("_J_") + TString( Form("%d",Jval) ) + TString("_hist.txt") ); }
    ss << "Bin no.  "<<"   Bin center value   "<<"  no. of levels   "<<"value"<<endl;
    for(int r=0; r<=rho->GetXaxis()->GetNbins(); ++r) 
    {   
        double Bin_cont = rho->GetBinContent( r )/rho->GetBinWidth( 1 );
        ss << r<<"             "<< rho->GetBinCenter( r )<<"      "<<rho->GetBinContent( r ) <<"             "<< Bin_cont <<endl;
        rho-> SetBinContent(r, Bin_cont);
    }
    ss.close();
    return rho;
}

//histogram creator for B
TH1F* Bplotter(TString name, double binwidth,TString sgn, int cutstat,  double LL, double UL, int Ji_stat, int J_i)
{
    //initializations
    TString s = ToLower(name);
    TString filename_new,filename1,filename2, filename_new2;
    TString tag, tail;
    //rename data files and also have a .txt extension
    if(name == TString("TE128"))
    {
        if(sgn==TString("pos")){filename1 = s + TString("pos_m1");  filename_new = s + TString("pos_m1.txt");   filename2 = s + TString("neg_m1");  filename_new2 = s + TString("neg_m1.txt");  tail=TString(" for (+) parity");}
        if(sgn==TString("neg")){filename1 = s + TString("neg_m1");  filename_new = s + TString("neg_m1.txt");   filename1 = s + TString("pos_m1");  filename_new2 = s + TString("pos_m1.txt");  tail=TString(" for (-) parity");}
        rename(filename1,filename_new);
        rename(filename2,filename_new2);
    }
    if(name == TString("SC44"))
    {
        filename1 = s + TString("e1");
        filename_new = s + TString("e1.txt");
        rename(filename1,filename_new);
    }
    ifstream read;
    read.open(filename_new);  
    int Ji,Jf,c=0;
    double Ei,Egamma,Bif, Bfi,B_sum,Egamma_max=0;
    //read spectrum file to set the range for histograms
    while(read >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
    {
        if(Egamma_max < abs(Egamma)) {Egamma_max = abs(Egamma);}
    }
    read.close();
    Egamma_max = ceil(Egamma_max); //upper limit for histogram
    int nbins;  
    nbins = Egamma_max/binwidth; //no.of bins
    if(name == TString("TE128")){tag = TString("(M1)");}//strings to set histogram names
    if(name == TString("SC44")){tag = TString("(E1)");  tail = TString("");}
    TH1F* bm1 = new TH1F(filename_new + sgn + TString(Form("%d",J_i)) + TString(Form("%d",Ji_stat)), TString("Graph of B") + tag + TString(" vs. Gamma de-excitation energy for ") + name + tail, nbins, 0 , Egamma_max ); // the histogram
    gStyle->SetOptStat(kFALSE); 
    bm1->GetXaxis()->SetTitle("Excitation energy (MeV)"); 
	bm1->GetYaxis()->SetTitle("Density of States (MeV^{-1})");
    bm1->SetLineWidth(1);
    ofstream ff;//write histogram to file
    if(cutstat==0)//choose what to do depending on whether user asks cuts in energy/J_i
    {
        if(Ji_stat==0){ff.open(name + sgn +TString("_B" +tag + ".txt"));}
        if(Ji_stat==1){ff.open(name + sgn +TString("_B") + tag + "J_i=" + Form("%d",J_i) + ".txt");}
    }
    if(cutstat==1)
    {   if(Ji_stat==0){ff.open(name + sgn + TString("_") + Form("%.0f",LL) + TString("-") + Form("%.0f",UL) + TString("MeV") + TString("_B" +tag + ".txt"));}
        if(Ji_stat==1){ff.open(name + sgn + TString("_") + Form("%.0f",LL) + TString("-") + Form("%.0f",UL) + TString("MeV") + TString("_B") + tag + "J_i=" + Form("%d",J_i) + ".txt");}
    }
    ff << "Bin no.  "<<    "Bin center value    "<<"  no. of entries       B"<<tag<<endl;
    //reading from file to find the minimum excitation energy
    int i=0;
    double emin=0, efmin=0, refp=0, refn=0;
    read.open(filename_new);
            while(read >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
            {
                refn = Ei-Egamma;
                if(Ei<emin){emin=Ei;}
                if(refn<efmin){efmin = refn;}
            }
            emin = min(emin,efmin);
    read.close();
    int counts[120];
    for(int z=0; z<120; ++z){counts[z]=0;}
    if(cutstat==0){UL=1000; LL=-1000;} //set extreme limits if no cuts in energy are done
    read.open(filename_new);
    int t = 0;
    while(read >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)//fill histogram depending on user preference
    {
        if( (Egamma > 0) && ((Ei-emin) < UL) && ((Ei-emin) > LL) && (Ji_stat==0) )   
            { 
                t = bm1-> Fill( Egamma , Bif );
                ++counts[t];
            }
        if( (Egamma < 0) && ((Ei-emin-Egamma) < UL) && ((Ei-emin-Egamma) > LL) && (Ji_stat==0) )   
            { 
                t = bm1-> Fill( abs(Egamma) , Bfi );
                ++counts[t];   
            }
        if( (Egamma > 0) && ((Ei-emin) < UL) && ((Ei-emin) > LL) && (Ji_stat==1) && (J_i==Ji/2) )   
            { 
                t = bm1-> Fill( Egamma , Bif );
                ++counts[t];
            }
        if( (Egamma < 0) && ((Ei-emin-Egamma) < UL) && ((Ei-emin-Egamma) > LL) && (Ji_stat==1) && (J_i==Ji/2) )   
            { 
                t = bm1-> Fill( abs(Egamma) , Bfi );
                ++counts[t];   
            }    
        t = 0;
    }
    for(int f=0; f< ( bm1->GetXaxis()->GetNbins() ); ++f )
    {
            if(counts[f]==0){bm1->SetBinContent(f,0);}//set bin count to zero if no events are in a bin
            else{bm1->SetBinContent(f, (bm1->GetBinContent(f)/counts[f]) );}//calculate <B(XL)>
            ff << f <<"    "<<bm1->GetBinCenter(f)<<"     "<<counts[f]<<"     "<<bm1->GetBinContent(f)<<endl;//output to data file
    }   
    ff<<"****The -nan states were replaced by a value 0 when filling the histogram.****"<<endl;    
    ff.close();
    if(name == TString("SC44"))//setting the plot environment
    {
        bm1 ->GetXaxis()->SetTitle("E_{#gamma} (MeV)"); 
	    bm1 ->GetYaxis()->SetTitle("<B(E1)> (e^{2}.fm^{2})");
        bm1 ->SetLineWidth(1); 
    }
    TString sym;
    if(name == TString("TE128"))
    {
        if(sgn == TString("pos")){sym = TString("+");}
        if(sgn == TString("neg")){sym = TString("-");}
        bm1 ->GetXaxis()->SetTitle("E_{#gamma}  (MeV)"); 
	    bm1 ->GetYaxis()->SetTitle(TString("<B(M1_{(") + sym + TString(")})> (#mu_{N}^{2})"));
        bm1 ->SetLineWidth(1); 
    } 
    return bm1;
}

//The main one
void NUC_ANALYSIS()
{
    //define Canvas event	
    TCanvas *c1 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("c1"); 
    delete c1; 
    c1 = new TCanvas("c1");
    c1->SetLogy();
    //define variables
    TString name, filename1, filename2;
    double BW = 0;
    int Jstat, Jnum, colors[] = {1,616,634,3,619,803,419,416};//array storing colors (ex: kWhite + an integer in the array gives a specific color)
    //ask user input
    cout <<"******* Start of phase I: LEVEL DENSITY plots *******"<<endl;
    cout <<"Input the STARTING STRING of the files you need to analyze (IN UPPERCASE. For this project, please enter either TE128 or SC44):"<<endl;
    cin >> name;
    cout << "Input the bin width (in MeV s):"<<endl;
    cin >> BW;
    cout << "Do you wish to make plots specific in J for level density (1 if YES, 0 if NO)?"<<endl;
    cin >> Jstat;
    //make nametags depending on the parity, for level scheme files 
    filename1 = name + TString("POS");
    filename2 = name + TString("NEG");
    //run fixer
    int count_1 = fixer(filename1);
    int count_2 = fixer(filename2);
    //read the reference energy values from 2 files
    double bin_max = 0;
    int J,N_eig,C, N,Z, h=0;
    double EX,E,E_MAX=0, E_refn,E_refp,E_dif;
    TString NZ, TJ, P, CEQ, NEQ, EXC, EQ, EEQ;
   ifstream r1;//open file of + parity to read the reference energy E_refp (in line 1)
    r1.open(filename1 + TString("_fixed.txt"));
    while (r1 >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {
        if(h==0){E_refp = E; ++h;}
        if(h>0){break;}
    }
    r1.close();
    h=0;
    ifstream r2;//open file of - parity to read the reference energy E_refn (in line 1)
    r2.open(filename2 + TString("_fixed.txt"));
    while (r2 >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {
        if(h==0){E_refn = E; ++h;}
        if(h>0){break;}
    }
    r2.close();
    E_dif = E_refn - E_refp; //E_refp is the GND state reference (0+)
    //run plotter and create 2 histograms for +/- parities
    TH1F *rho_plot_1 = plotter(filename1,BW,0,0,0);
    TH1F *rho_plot_2 = plotter(filename2,BW,E_dif,0,0);
    rho_plot_1->SetLineColor(kRed);
    rho_plot_2->SetLineColor(kBlue);
    auto legend = new TLegend(0.2,0.1,0.2,0.1);
    //create the total level density histogram
    double bin_max1 = rho_plot_1->GetBinCenter(rho_plot_1->GetXaxis()->GetNbins()); 
    double bin_max2 = rho_plot_2->GetBinCenter(rho_plot_2->GetXaxis()->GetNbins()); 
    //find the histogram with largest range
    ifstream read1;
    if(bin_max1 > bin_max2){read1.open(filename1 + TString("_fixed.txt"));}
    else if(bin_max1 < bin_max2){read1.open(filename2 + TString("_fixed.txt"));}
    //read file corresponding to largest range and set xrange for the histograms
    while (read1 >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {if(EX>bin_max){bin_max=EX;}}
    read1.close();
    bin_max = ceil(bin_max+E_dif);
    //The new histogram for total level density
    TH1F* rho_tot = new TH1F("tot",TString("Graph of Level density vs. Excitation energy for ") + name , (bin_max)/BW, 0 , bin_max);
    gStyle->SetOptStat(kFALSE); 
    ofstream mm (name + TString("_tot.txt"));//print data to file
    mm<<"Bin no.     "<<"Bin center      "<<"value"<<endl; 
    for(int s=0; s<=max(rho_plot_2->GetXaxis()->GetNbins(),rho_plot_1->GetXaxis()->GetNbins()); ++s) //combine level densities of the 2 parities to get combined one
    {   
        double Bin_cont = rho_plot_1->GetBinContent( s ) + rho_plot_2->GetBinContent( s );
        rho_tot-> Fill(rho_plot_1->GetBinCenter(s), Bin_cont);
        mm<< s <<"     "<< rho_plot_1->GetBinCenter(s) <<"      "<<Bin_cont<<endl;
    }
    mm.close();
    rho_tot->GetXaxis()->SetTitle("Excitation energy (MeV)"); //plot environment
	rho_tot->GetYaxis()->SetTitle("Level density (MeV^{-1})");
    rho_tot->SetLineColor(kGreen);
    rho_tot->SetLineWidth(1);
    rho_tot->GetXaxis()->SetRangeUser(0,10);
    rho_tot->Draw("HIST,SAME"); 
    rho_plot_1->Draw("HIST,SAME"); 
    rho_plot_2->Draw("HIST,SAME");
    //selections on J_i
    int J_max1=0, J_max2=0;
    ifstream rr;
    rr.open(filename1 + TString("_fixed.txt"));
    while (rr >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)//get J range of either parities by reading files, and print them
    {if(J/2>J_max1){J_max1=J/2;} }
    rr.close();

    ifstream rrr;
    rrr.open(filename2 + TString("_fixed.txt"));
    while (rrr >> NZ >> N >> Z >> TJ >> J >> P >> NEQ >> N_eig >> CEQ >> C >> EXC >> EQ >> EX >> EEQ >> E)
    {if(J/2>J_max2){J_max2=J/2;} }
    rrr.close();
    cout<<"Jmax for (+) parity = "<<J_max1<<",  Jmax for (-) parity = "<<J_max2<<endl;
    int Jarr[max(J_max1, J_max2)], parity[max(J_max1, J_max2)], Jread, Pread;
    if(Jstat== 1)//if user wants J cuts, ask how many he wants
        {
            cout << "How many J values?"<<endl;
            cin >> Jnum;
            for(int g=0; g<Jnum; ++g)//read J and pi as per user request
            {
                cout<<"***** Enter J value no." + TString(Form("%d",g+1)) +" (below Jmax for non-zero outputs): *****"<<endl;
                cin >> Jread;
                Jarr[Jread]= Jread;
                cout <<"***** Enter parity for J value no." + TString(Form("%d",g+1)) +" ( (1 FOR +) or (-1 FOR -) ): *****"<<endl;
                cin >> Pread;
                parity[Jread] = Pread;
            }
        }
    TString rho_plot_p_str, rho_plot_n_str;
    int cc=5;
    TH1F *rho_plot_p_arr[max(J_max1, J_max2)];//arrays of histograms cut in J for either parity
    TH1F *rho_plot_n_arr[max(J_max1, J_max2)];
    for(int q=0; q<max(J_max1, J_max2)+1; ++q)//depending on user request, print only the ones user asked for
        {
            TString ftag,ptag;
            rho_plot_p_arr[q]   = plotter(filename1,BW,0,1,q);
            rho_plot_n_arr[q]   = plotter(filename2,BW,E_dif,1,q);
            if(q==Jarr[q])
            {
                if(parity[q] == 1)          
                {
                    rho_plot_p_arr[q] ->SetLineColor(kWhite + cc );
                    rho_plot_p_arr[q] ->Draw("HIST,SAME");
                    ftag = filename1;  
                    ptag = TString("+");
                    legend->AddEntry(filename1 + TString("1") + TString(Form("%d",q)) ,TString("J = ") + TString(Form("%d",Jarr[q])) + "^{" + ptag + "}","l");
                    ++cc;
                }
                if(parity[q] == -1)    
                {
                    rho_plot_n_arr[q] ->SetLineColor(kWhite + cc );
                    rho_plot_n_arr[q] ->Draw("HIST,SAME");
                    ftag = filename2;  
                    ptag = TString("-");
                    legend->AddEntry(filename2 + TString("1") + TString(Form("%d",q)),TString("J = ") + TString(Form("%d",Jarr[q])) + "^{" + ptag + "}","l");
                    ++cc;
                }
                
            }     
        }
    //the legend
    legend->AddEntry(rho_plot_1,"(+) parity states","l");
    legend->AddEntry(rho_plot_2,"(-) parity states","l");
    legend->AddEntry(rho_tot,"combined","l");
    legend->SetBorderSize(0);
    legend->Draw();
    c1-> SaveAs(TString("level_density_") + name + TString(".pdf"));
    int cutstat, cutnum=0, Jistat, J_inum=0, pstat;
    int J_i[10];
    for(int k=0; k<10; ++k){J_i[k]=0;}
    TString ptag;
    cout <<"******* Start of phase II: B(M1) plots *******"<<endl;
    if(name == "TE128")//user inputs
    {
        cout << "Which parity do you wish to analyze (input 1 for (+), -1 for (-))?"<<endl;
        cin >> pstat;  
        if(pstat==1){ptag = "pos";}//TStrings depending on user's parity choice
        else if(pstat==-1){ptag = "neg";}
    }
    cout << "Do you wish to produce plots with cuts in excitation energy (1 if YES, 0 if NO)?:"<<endl;
    cin >> cutstat;
    if(cutstat ==1)
    {
        cout << "How many cuts?:"<<endl;
        cin >> cutnum;
    }
    double UL[cutnum], LL[cutnum];//arrays storing upper and lower limits of cuts in energy
    if (cutstat== 1)//get the limits for each cut
    {
        for(int k=0; k<cutnum; ++k)
        {
            cout<<TString("*****For cut no.") + Form("%d",k+1) +TString(" : *****")<<endl;
            cout << "Input the UPPER LIMIT (in MeV s):"<<endl;
            cin >> UL[k];
            cout << "Input the LOWER LIMIT (in MeV s):"<<endl;
            cin >> LL[k];
        }
    }
    cout << "Do you wish to make a selection on J_i (1 if YES, 0 if NO)?:"<<endl;//similarly cuts in J
    cin >> Jistat;
    if(Jistat==1)
    {
        cout << "How many J_i values? :"<<endl;
        cin >> J_inum;
        
        for(int p=0; p<J_inum; ++p)
        {
            cout << "Input J_i :"<<endl;
            cin >> J_i[p];
        }
    }
    c1->Clear();
    auto nlegend = new TLegend(0.4,0.15,0.4,0.15);
    int colorcount=0, countnum=0, Jcnum=0;
    countnum=cutnum;
    Jcnum = J_inum;
    if(cutnum==0){countnum=1;}//if user asks no cuts in E_exc, set this value to 1 so as to go through the code uncut
    if(J_inum==0){Jcnum=1;}//if user asks no cuts in J, set this value to 1 so as to go through the code uncut
    if(name == TString("TE128"))
    {   //always do uncut plot as reference
        TString BM,BM1_ = "BM1_", sgn;
        if(pstat==1){sgn=TString("+");}
        else if(pstat==-1){sgn=TString("-");}
        TH1F *BM1 = Bplotter(name,BW,ptag,0,0,0,0,0);
        BM1 ->SetLineColor(kRed);
        BM1 ->SetMaximum(0.1);
        BM1 ->Draw("HIST");
        nlegend->AddEntry(BM1,"<B(M1,#pi^{" + sgn + "})> : uncut","l"); 
        for(int m=0; m<countnum; ++m)//cuts in J and E_exc printed or not as per user choice
        {
            for(int s=0; s<Jcnum; ++s)
            {
                BM = BM1_ + Form("%d",m) + Form("%d",s);
                TH1F *BM = Bplotter(name,BW,ptag,cutstat,LL[m],UL[m],Jistat,J_i[s]);
                BM ->SetLineColor(kWhite + colors[colorcount]);
                BM ->Draw("HIST,SAME");
                if(cutnum==0)   {if(Jistat==1) {nlegend->AddEntry(BM,TString("<B(M1,#pi^{" + sgn + "})> : All E_{EXC}, J_{i} = ") + TString(Form("%d",J_i[s])),"l");}} 
                else           {
                                    if(Jistat==1)   {nlegend->AddEntry(BM,TString("<B(M1,#pi^{" + sgn + "})> : ") + Form("%.1f",LL[m]) + TString("-") + Form("%.1f",UL[m]) + TString("MeV, J_{i} = ") + TString(Form("%d",J_i[s])),"l");}
                                    else            {nlegend->AddEntry(BM,TString("<B(M1,#pi^{" + sgn + "})> : ") + Form("%.1f",LL[m]) + TString("-") + Form("%.1f",UL[m]) + TString("MeV, all J_{i}") ,"l");}
                                }
                ++colorcount;
            }    
        }
        nlegend->SetBorderSize(0);
        nlegend->Draw();
        c1->Update();
        c1->WaitPrimitive();
        c1-> SaveAs("B(M1)_Te128.pdf");
        cout <<"******* Start of phase III: nuclear strength plots *******"<<endl;
        //creating fxl plots
        TH1F* fxlp = new TH1F("fxlp", TString("Graph of nuclear strength function vs. E_{#gamma} for B(M1) transitions of Te128"), 10/BW, 0, 10);//fxl for + parity 
        TH1F* fxln = new TH1F("fxln", TString("Graph of nuclear strength function vs. E_{#gamma} for B(M1) transitions of Te128"), 10/BW, 0, 10);//fxl for + parity
        TH1F* fxl = new TH1F("fxl", TString("Graph of nuclear strength function vs. E_{#gamma} for B(M1) transitions of Te128"), 10/BW, 0, 10);//fxl combined
        TH1F* ftp = new TH1F("testp", "testp", 5/BW, 0, 5); //dummy histogram to find bin index for Egamma, + parity
        TH1F* ftn = new TH1F("testn", "testn", 5/BW, 0, 5); //dummy histogram to find bin index for Egamma, - parity
        TH1F* fpcount = new TH1F("fpcount", name, (bin_max)/BW, 0 , bin_max); //dummy histogram to find bin index for E_exc, + parity
        TH1F* fncount = new TH1F("fncount", name, (bin_max)/BW, 0 , bin_max); //dummy histogram to find bin index for E_exc, - parity
        ifstream readpos;
        ifstream readpos1;
        ifstream readneg;
        ifstream readneg1;
        double Epmin=0, Exc=0, Ei, Egamma, Bif, Bfi;
        double Exc1=0, Ei1, Egamma1, Bif1, Bfi1, btot=0;
        double barrp[100][max(J_max1, J_max2)][100], barrn[100][max(J_max1, J_max2)][100];
        int bin=0, cbinp[120],cbinn[120],uv=0, Ji, Jf, temp=0,bcount=0;
        int bin1=0, Ji1, Jf1, countp[100][max(J_max1, J_max2)][100],countn[100][max(J_max1, J_max2)][100];
        for(int c=0; c<120; ++c){cbinp[c]=0;}
        readpos.open("te128pos_m1.txt");
        while (readpos >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
        {
            Exc = Ei - E_refp; 
            bin = fpcount->Fill(Exc);//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Exc
            bin1 = ftp->Fill(abs(Egamma));//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Egamma
            ++countp[bin][Ji/2][bin1];//store no. of B values of a given combo of (Exc,Ji,Egamma)
            if(Egamma>0){barrp[bin][Ji/2][bin1]=barrp[bin][Ji/2][bin1]+Bif;}//store total B of a given combo of (Exc,Ji,Egamma)
            if(Egamma<0){barrp[bin][Ji/2][bin1]=barrp[bin][Ji/2][bin1]+Bfi;}  
        }
        readpos.close();
        readpos.open("te128pos_m1.txt");
        while (readpos >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
        {
            Exc = Ei - E_refp;
            bin = fpcount->Fill(Exc);//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Exc
            bin1 = ftp->Fill(abs(Egamma));//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Egamma
            if(Egamma>0)         {uv = fxlp-> Fill(Egamma, 0.00000001154*barrp[bin][Ji/2][bin1]*rho_plot_p_arr[Ji/2]->GetBinContent(bin)/countp[bin][Ji/2][bin1]); ++cbinp[uv];  }
            else if(Egamma<0)    {uv = fxlp-> Fill(abs(Egamma), 0.00000001154*barrp[bin][Ji/2][bin1]*rho_plot_p_arr[Ji/2]->GetBinContent(bin)/countp[bin][Ji/2][bin1]); ++cbinp[uv]; }
        }
        readpos.close();

        ofstream pp ("te128pos_m1_f.txt");//print to data file, + parity
        pp<<"Bin no.    "<<"Bin center    "<<"value"<<"  No. of counts"<<endl;
        for(int e=0; e<=fxlp->GetXaxis()->GetNbins(); ++e)
        {
            fxlp->SetBinContent( e, fxlp->GetBinContent(e)/ cbinp[e] );
            if(cbinp[e]==0){fxlp->SetBinContent( e, 0 );}
            pp<<e<<"          "<<fxlp->GetBinCenter(e)<<"     "<<fxlp->GetBinContent(e)<<"   "<<cbinp[e]<<endl;
        }
        pp.close();
        
        readneg.open("te128neg_m1.txt");
        while (readneg >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
        {
            Exc = Ei - E_refp; 
            bin = fncount->Fill(Exc);//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Exc
            bin1 = ftn->Fill(abs(Egamma));//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Egamma
            ++countn[bin][Ji/2][bin1];//store no. of B values of a given combo of (Exc,Ji,Egamma)
            if(Egamma>0){barrn[bin][Ji/2][bin1]=barrn[bin][Ji/2][bin1]+Bif;}//store total B of a given combo of (Exc,Ji,Egamma)
            if(Egamma<0){barrn[bin][Ji/2][bin1]=barrn[bin][Ji/2][bin1]+Bfi;} 
        }
        readneg.close();
        readneg.open("te128neg_m1.txt"); //gnd state reference is E_refp
        for(int d=0; d<120; ++d){cbinn[d]=0;}
        while (readneg >> Ji >> Jf >> Ei >> Egamma >> Bif >> Bfi)
        {
            Exc = Ei - E_refp;
            bin = fncount->Fill(Exc);
            bin1 = ftn->Fill(abs(Egamma));
            uv = fxln-> Fill(abs(Egamma), 0.00000001154*barrn[bin][Ji/2][bin1]*rho_plot_n_arr[Ji/2]->GetBinContent(bin)/countn[bin][Ji/2][bin1]); ++cbinn[uv]; //calculation. the constant is the conversion factor.
        }
        readneg.close();
        ofstream nn ("te128neg_m1_f.txt");
        nn<<"Bin no.    "<<"Bin center    "<<"value"<<"  No. of counts"<<endl;
        for(int f=0; f<=fxln->GetXaxis()->GetNbins(); ++f)
        {
            fxln->SetBinContent( f, fxln->GetBinContent(f)/ cbinn[f] );
            if(cbinn[f]==0){fxln->SetBinContent( f, 0 );}
            nn<<f<<"          "<<fxln->GetBinCenter(f)<<"     "<<fxln->GetBinContent(f)<<"   "<<cbinn[f]<<endl; //print to data file, + parity
        }
        nn.close();
        ofstream tt ("te128_m1_f_tot.txt");
        tt<<"Bin no.    "<<"Bin center    "<<"value"<<endl;
        for(int g=0; g<fxl->GetXaxis()->GetNbins(); ++g)//extract values from either parities, compute the total fxl
        {
            fxl->SetBinContent( g, ( fxln->GetBinContent(g)*cbinn[g] + fxlp->GetBinContent(g)*cbinp[g] ) / ( cbinn[g] + cbinp[g] ) );
            if((cbinn[g] + cbinp[g]) ==0){fxl->SetBinContent( g, 0 );}
            tt<<g<<"          "<<fxln->GetBinCenter(g)<<"     "<<fxln->GetBinContent(g)<<endl; //print data
        }
        tt.close();
        auto flegend = new TLegend(0.4,0.2,0.4,0.2);
        flegend->AddEntry(fxlp,"f_{M1}, (+) parity states","l"); 
        flegend->AddEntry(fxln,"f_{M1}, (-) parity states","l"); 
        flegend->AddEntry(fxl,"f_{M1}, combined","l"); 
        flegend->SetBorderSize(0);
        fxlp->GetXaxis()->SetTitle("E_{#gamma} (MeV)"); 
        fxlp->GetYaxis()->SetTitle("f_{M1} (MeV^{-3})");
        fxlp->SetLineColor(kRed);
        fxlp->SetLineWidth(1);
        fxln->SetLineColor(kBlue);
        fxln->SetLineWidth(1);
        fxl->SetLineColor(kGreen);
        fxl->SetLineWidth(1);
        fxlp ->GetXaxis()->SetRangeUser(0,5);
        fxlp ->Draw("HIST");
        fxln ->GetXaxis()->SetRangeUser(0,5);
        fxln ->Draw("HIST, same");
        fxl ->GetXaxis()->SetRangeUser(0,5);
        fxl ->Draw("HIST, same");
        flegend->Draw();
        c1-> SaveAs("f(M1)_Te128.pdf");
    }
    else if(name == TString("SC44"))//same repeated for Sc44 (but algorithm slightly different because selection rules and file systems are different)
    {  
        TH1F *BE1 = Bplotter(name,BW,"",0,0,0,0,0);// the uncut plot in B(E1)
        BE1 ->SetLineColor(kRed);
        BE1->GetYaxis()->SetRangeUser(10e-9,10e-3);
        BE1 ->GetXaxis()->SetRangeUser(0,10);
        BE1 ->Draw("HIST");
        nlegend->AddEntry(BE1,"<B(E1)> : Uncut","l");
        TString BE, BE1_ = "BE1_";
        for(int m=0; m<countnum; ++m)//make histograms by cuts in J and E_exc depending on user choice
        {
            for(int s=0; s<Jcnum; ++s)
            {
                BE = BE1_ + Form("%d",m) + Form("%d",s);
                TH1F *BE = Bplotter(name,BW,"",cutstat,LL[m],UL[m],Jistat,J_i[s]);
                BE ->SetLineColor(kWhite + colors[colorcount]);
                BE ->Draw("HIST,SAME"); 
                if(cutnum==0)   {if(Jistat==1)  {nlegend->AddEntry(BE,TString("<B(E1)> : All E_{EXC}, J_{i} = ") + TString(Form("%d",J_i[s])),"l");}} 
                else            {
                                    if(Jistat==1)   {nlegend->AddEntry(BE,TString("<B(E1)> : ") + Form("%.1f",LL[m]) + TString("-") + Form("%.1f",UL[m]) + TString("MeV, J_{i} = ") + TString(Form("%d",J_i[s])),"l");}
                                    else            {nlegend->AddEntry(BE,TString("<B(E1)> : ") + Form("%.1f",LL[m]) + TString("-") + Form("%.1f",UL[m]) + TString("MeV, all J_{i}"),"l");}
                                }
                ++colorcount;
            }
        }
        nlegend->Draw();
        c1->Update();
        c1->WaitPrimitive();
        c1-> SaveAs("B(E1)_Sc44.pdf");//SAVE PDF
        cout <<"******* Start of phase III: nuclear strength plots *******"<<endl;
        TH1F* fxls = new TH1F("fxlps", TString("Graph of nuclear strength function vs. E_{#gamma} for B(E1) transitions of Sc44"), 12/BW, 0, 12);//HISTOGRAM STORING FXL. LIMITED TO 12 BECAUSE 10 IS WHAT WE NEED
        TH1F* fpcounts = new TH1F("fpcounts", name, (bin_max)/BW, 0 , bin_max);//DUMMY HISTOGRAM
        TH1F* fts = new TH1F("fts", "fts", 30/BW, 0, 30);
        ifstream readposs;
        double Epmins=0, Excs=0, Eis, Egammas, Bifs, Bfis;
        int bins=0, cbinps[120],uvs=0, Jis, Jfs, bin1s=0,T=0;
        int counts[120][max(J_max1, J_max2)][120];
        double barrs[120][max(J_max1, J_max2)][120];
        for(int c=0; c<120; ++c){cbinps[c]=0;}
        readposs.open("sc44e1.txt"); //OPEN SPECTRUM FILE
        while (readposs >> Jis >> Jfs >> Eis >> Egammas >> Bifs >> Bfis)
        {
            Excs = Eis - E_refp; 
            bins = fpcounts->Fill(Excs);//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Exc
            bin1s = fts->Fill(abs(Egammas));//DUMMY HISTOGRAMS TO EXTRACT BIN INDEX CORRESPONDING TO Egamma
            ++counts[bins][Jis/2][bin1s];//store no. of B values of a given combo of (Exc,Ji,Egamma)
            if( (Egammas>0) ){barrs[bins][Jis/2][bin1s]=barrs[bins][Jis/2][bin1s]+Bifs;}//store total B of a given combo of (Exc,Ji,Egamma)
            if( (Egammas<0) ){barrs[bins][Jis/2][bin1s]=barrs[bins][Jis/2][bin1s]+Bfis;}  
        }
        readposs.close();
        readposs.open("sc44e1.txt");
        while (readposs >> Jis >> Jfs >> Eis >> Egammas >> Bifs >> Bfis)
        {
            Excs = Eis - E_refp;
            bins = fpcounts->Fill(Excs);
            bin1s = fts->Fill(abs(Egammas));
            uvs = fxls-> Fill(abs(Egammas), 0.00000106*barrs[bins][Jis/2][bin1s]*rho_tot->GetBinContent(bins)/counts[bins][Jis/2][bin1s]); ++cbinps[uvs];   //calculation, same as in Te128
        }
        readposs.close();
        ofstream pps ("sc44e1_f.txt");
        pps<<"Bin no.    "<<"Bin center    "<<"value"<<"   counts"<<endl;//write data to file
        for(int es=0; es<fxls->GetXaxis()->GetNbins(); ++es)
        {
            fxls->SetBinContent( es, fxls->GetBinContent(es)/ cbinps[es] );
            if(cbinps[es]==0){fxls->SetBinContent( es, 0 );}
            pps<<es<<"          "<<fxls->GetBinCenter(es)<<"     "<<fxls->GetBinContent(es)<<"     "<<cbinps[es]<<endl;
        }
        pps.close();
        auto fslegend = new TLegend(0.4,0.2,0.4,0.2);
        fslegend->AddEntry(fxls,"f_{E1}","l"); 
        fxls->GetXaxis()->SetTitle("E_{#gamma} (MeV)"); 
        fxls->GetYaxis()->SetTitle("f_{E1} (MeV^{-3})");
        fxls->SetLineColor(kRed);
        fxls->SetLineWidth(1);
        fxls ->GetXaxis()->SetRangeUser(0,10);
        fxls ->Draw("HIST");
        fslegend->SetBorderSize(0);
        fslegend->Draw();
        c1-> SaveAs("f(E1)_Sc44.pdf");
    }
}
