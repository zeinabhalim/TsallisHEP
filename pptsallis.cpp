#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TGraph.h>
#include <TAxis.h>
#include <string>

//  Tsallis function as a return value

double Tsallis(double *x, double *par) 
{
    double mass = 0.13957018; // Constant mass 
    double xval = x[0];// transverse momentum variable
    double A = par[0];   // Normalization factor
    double q = par[1];   // non-equilibrium index
    double T = par[2];   // Temperature parameter
    double mu = par[3];   // chemical potential parameter
    double y = par[4];   // psaudo-rapidity parameter


    return par[0] * xval * sqrt(xval * xval + mass * mass) * 
           TMath::CosH(par[4]) * 
           TMath::Power(1 + ((par[1] - 1)*( 1/ par[2]) * 
           (sqrt(xval * xval + mass * mass) * TMath::CosH(par[4]) - par[3])), 
           (1 / (1-par[1])));
}
void pptsallis() {
    // Create a new canvas to hold the graphs together
    TCanvas *c1 = new TCanvas("canvas", "Multiple Graphs with Tsallis Fit Function", 1000, 800); 
     c1->Divide(2, 2); 

    // Define multiple sets of data files 
   std::vector<std::vector<const char*>> fileGroups = 
       { {"pp91.txt", "pp231.txt", "pp71.txt"}, 
        {"pp92.txt", "pp232.txt", "pp72.txt"}, 
        {"pp93.txt", "pp233.txt", "pp73.txt"}, 
        {"pp94.txt", "pp234.txt", "pp74.txt"}} ;
        
  // domain limits of fit function    
    double xmin = 0.15, xmax = 2.0;
    TF1 *tsallisFunp = new TF1("fitFunc", Tsallis, xmin, xmax, 5);

// setting Colors and marker styles
   int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan, kBlack, kViolet, kTeal, kPink, kGray, kSpring};
   int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33};
   
// Save parameters to a text file
          std::ofstream outfile("Tsallisfit_parametersppTev.txt");
              outfile << "fitting of the transverse momentum spectra of charged particles in proton-proton at Sqrt{s_NN} = 0.9,2.36,7 TeV, using Tsallis  distribution function."<<std::endl;
              outfile << std::endl;
              
    // Loop over each group of data sets  
    for (size_t groupIndex = 0; groupIndex < fileGroups.size(); groupIndex++) {
        c1->cd(groupIndex + 1); //to  Move around the canvas partitions
        TMultiGraph *mg = new TMultiGraph();
        TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.9);
        gPad->SetLogy(); // Apply log scale to Y-axis
        gPad->SetLogx(); // Apply log scale to x-axis
        // Loop over each file in each group in data set (3 files per group)
        for (size_t i = 0; i < fileGroups[groupIndex].size(); i++) {
            const char* filename = fileGroups[groupIndex][i];
            std::vector<double> x_vals, y_vals, x_errs, y_errs;

            // Read the data file
            std::ifstream file(filename);
            if (!file) {
               std::cerr << "Error: Cannot open file " << filename << std::endl;
                continue;
           }

            double x, y, xerr, yerr;
            while (file >> x >> y >> yerr >> xerr) {
                x_vals.push_back(x);
                y_vals.push_back(y);
                x_errs.push_back(xerr);
                y_errs.push_back(yerr);
            }
            file.close();
            
            // define error bars
            int dataIndex = groupIndex * 3 + i;
            TGraphErrors *gr = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], 0, &y_errs[0]);
            gr->SetMarkerStyle(markers[dataIndex]);
            gr->SetMarkerColor(colors[dataIndex]);
            gr->SetLineColor(colors[dataIndex]);

       
    //estimating the fitting parameters for each group of data set in each graph
    if (groupIndex == 0) {
        tsallisFunp->SetParLimits(0,0.1,3.0); // Initial Value for A (Normalization)
        tsallisFunp->SetParLimits(1,1.1, 1.17);//Initial Value for q (Non-Extensive index)
        tsallisFunp->SetParLimits(2, 0.1,0.35); // Initial Value for T (temperture)
        tsallisFunp->SetParLimits(3, -0.1,2.0); // Initial Value for mu (Chemical Pot.)
        tsallisFunp->FixParameter(4, 0.1); // Initial Value for eta (rapidity)
        tsallisFunp->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","rapidity");
        } else if (groupIndex == 1) {
      tsallisFunp->SetParLimits(0,0.1,3.0); 
        tsallisFunp->SetParLimits(1,1.1, 1.17);
        tsallisFunp->SetParLimits(2, 0.1,0.35); 
        tsallisFunp->SetParLimits(3, -0.1,2.0); 
        tsallisFunp->FixParameter(4, 0.3); 
        tsallisFunp->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","rapidity");
        } else if (groupIndex == 2) {
      tsallisFunp->SetParLimits(0,0.1,3.0); 
        tsallisFunp->SetParLimits(1,1.1, 1.165);
        tsallisFunp->SetParLimits(2, 0.1,0.3); 
        tsallisFunp->SetParLimits(3, -0.1,2.0); 
        tsallisFunp->FixParameter(4, 0.5); 
        tsallisFunp->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","rapidity");
        } 
        else {
       tsallisFunp->SetParLimits(0,0.1,3.0); 
        tsallisFunp->SetParLimits(1,1.0, 1.165);
        tsallisFunp->SetParLimits(2, 0.1,0.3); 
        tsallisFunp->SetParLimits(3, -0.1,2.0); 
        tsallisFunp->FixParameter(4, 0.7); 
        tsallisFunp->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","rapidity");
        
        }
  
          gr->Fit(tsallisFunp, "R");  

        const char *titles[4] = {
       "0.0< #eta<0.2",
        "0.2< #eta<0.4",
        "0.4< #eta<0.6",
        "0.6< #eta<0.8"
    };
   const char *energyLabels[3] = {
       "#sqrt{s_{NN}} = 0.9 TeV",
        "#sqrt{s_{NN}} = 2.36 TeV",
        "#sqrt{s_{NN}} = 7 TeV"
    };
       
       // define Chi-Square , Degrees of Freedom and p-value
        double chi2 =  tsallisFunp->GetChisquare();
        int ndf =  tsallisFunp->GetNDF();
        double chi2ndf =  chi2 / ndf ;
        double pValue = TMath::Prob(chi2, ndf);
        
        outfile << "Fitting results for p-p collision at : "<<energyLabels[i]<<", at pesudo-rapidity "<<titles[groupIndex]<<std::endl;
       outfile << "Fitting parameters:" << std::endl;
        outfile << "  A (Normalization): " <<  tsallisFunp->GetParameter(0) << " ± " <<  tsallisFunp->GetParError(0) << std::endl;
        outfile << "  q (non-equilibrium index): " <<  tsallisFunp->GetParameter(1)<< " ± " << tsallisFunp->GetParError(1) << std::endl;
       outfile <<"  T (Temperature): " <<  tsallisFunp->GetParameter(2) << " ± " <<  tsallisFunp->GetParError(2) << std::endl;
       outfile << "  mu (Chemical pot.): " <<  tsallisFunp->GetParameter(3) << " ± " <<  tsallisFunp->GetParError(2) << std::endl;
       
            outfile << "Estimating statistical accuracy"<<std::endl;
       outfile << "chi2 " << chi2<< std::endl;
        outfile << "ndf " << ndf<< std::endl;
       outfile << "chi2ndf " << chi2ndf<< std::endl;
        outfile << "p-Value: " << pValue << std::endl;
        outfile << "----------------------------------------------------------------------------------------------"<<std::endl;
              outfile << std::endl; 
              
         mg->SetTitle(titles[groupIndex]); // Add titles for each graph at the top
    
         leg->AddEntry(gr,  energyLabels[i] , "P");
         mg->Add(gr, "P"); 
         mg->GetXaxis()->SetLimits(0.1, 3.2);  // X-axis range
         mg->SetMinimum(0.05);  // Y-axis lower limit 
         mg->SetMaximum(1e2);   // Y-axis upper limit

        mg->GetXaxis()->SetTitle("P_{t}[GeV]");
        mg->GetYaxis()->SetTitle("(1/(N_{evt}))*d^2(N)/d\\eta /dP _{t }");  
             
        }
        

      tsallisFunp->SetLineColor(kRed);    // Set the color of fit function 
      tsallisFunp->SetLineWidth(2);
      leg->AddEntry(tsallisFunp, "Tsallis Fit", "l");
    
        mg->Draw("APL"); 

        leg->Draw();
    }
 
}
