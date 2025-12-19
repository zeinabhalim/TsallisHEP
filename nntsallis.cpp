#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// the Tsallis function as a return value

double TsallisFunction(double *x, double *par) 
{
    double mass = 0.13957018; // mass term of produced particle(consider the mass of pion particle)
    double xval = x[0]; //the transverse momentum 
    double A = par[0];   // Normalization factor
    double q = par[1];   // non-equilibrium tsallis index
    double T = par[2];   // average Temperature at freeze-out stage
    double mu = par[3];   // chemical potential
    double y = par[4];   // rapdity parameter(fixed to mid rapidity)


    return par[0] * xval * sqrt(xval * xval + mass * mass) * 
           TMath::CosH(par[4]) * 
           TMath::Power(1 + ((par[1] - 1)*( 1/ par[2]) * 
           (sqrt(xval * xval + mass * mass) * TMath::CosH(par[4]) - par[3])), 
           (1 / (1-par[1])));
}
void nntsallis() {
    
    TCanvas *canvas = new TCanvas("canvas", "xe-xe with Tsallis Fit Function", 1300, 1000);
      
    canvas->SetLogy(true);  
    canvas->SetLogx(true);  

 //locating the legend box       
    TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.9); 
   leg->SetBorderSize(1);  
    leg->SetTextSize(0.02); 
    leg->SetFillStyle(0.01);   
  
// loading the input data files
    std::vector<std::string> filenames = {
       "xxe2.txt" ,
       "xxe3.txt" ,
        "xxe4.txt" ,
        "xxe5.txt" ,
        "xxe6.txt",
        "xxe7.txt",
        "xxe8.txt" ,
        "xxe9.txt"  
         };
        

double xmin = 0.2, xmax = 3.03;
TF1 *tsallisFunc = new TF1("fitFunc", TsallisFunction, xmin, xmax, 5);
tsallisFunc->SetNpx(2000);
// Save parameters to a text file
          std::ofstream outfile("Tsallisfit_parametersxexe.txt");
              outfile << "fitting of the transverse momentum spectra of charged particles in Xe-Xe at Sqrt{s_NN} = 5.44 TeV, using Tsallis  distribution function."<<std::endl;
              outfile << std::endl;
       
    // Loop over each file and process the data
    for (size_t file_index = 0; file_index < filenames.size(); ++file_index) {

        std::vector<double> x_data, y_data, x_error, y_error;

        const std::string& filename = filenames[file_index];

        // Open the file
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;  
        }

        double x, y, ex, ey;
        while (input_file >> x >> y >> ex >> ey) {
            x_data.push_back(x);  // Add x-values to the vector
            y_data.push_back(y);  // Add y-values to the vector
            x_error.push_back(ex);  // Add x-error values to the vector
            y_error.push_back(ey);  // Add x-error values to the vector
        }


        input_file.close();

        // Create a TGraphErrors to store the data with error bars
        int n_points = x_data.size();
        TGraphErrors *graph = new TGraphErrors(n_points);

        // Set title and labels for the graph
        graph->SetTitle(" ;P_{t}[GeV];(1/(N_{evt}))*d^2(N)/d\\eta /dP _{t } ");
        tsallisFunc->SetLineColor(kBlack);
         
      
        const char* names[sizeof(filenames.size())]={"0-5%", "5-10%", "10-20%", "20-30%",
        "30-40%", "40-50%", "50-60%", "60-70%"};
        int colors[sizeof(filenames.size())] = {kRed, kBlue, kGreen + 2, kMagenta,
                            kOrange + 7, kCyan + 2, kViolet, kBlack};
            // define Chi-Square , Degrees of Freedom and p-value
        double chi2 =  tsallisFunc->GetChisquare();
        int ndf =  tsallisFunc->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 999;
// Print out the fitting parameters 

       // outfile << "Fitting results for centerality: "<<names[1]<<std::endl;
       outfile << "Fitting parameters:" << std::endl;
        outfile << "  A (Normalization): " << tsallisFunc->GetParameter(0) << " ± " << tsallisFunc->GetParError(0) << std::endl;
        outfile << "  q (non-equilibrium index): " << tsallisFunc->GetParameter(1)<< " ± " <<tsallisFunc->GetParError(1) << std::endl;
       outfile <<"  T (Temperature): " << tsallisFunc->GetParameter(2) << " ± " << tsallisFunc->GetParError(2) << std::endl;
       outfile << "  mu (Chemical pot.): " << tsallisFunc->GetParameter(3) << " ± " << tsallisFunc->GetParError(2) << std::endl;

            outfile << "Estimating statistical accuracy"<<std::endl;
       outfile << "chi2 " << chi2<< std::endl;
        outfile << "ndf " << ndf<< std::endl;
       outfile << "chi2ndf " << chi2ndf<< std::endl;
if (chi2ndf > 0.5 && chi2ndf < 8.0) {
         outfile << " VAILD FITTING   "<<std::endl;}
         else{ outfile << "                                          INVAILD FITTING                              "<<std::endl;}
        outfile << "----------------------------------------------------------------------------------------------"<<std::endl;
              outfile << std::endl;                  

         for (int i = 0; i <= filenames.size(); i++) {
         // Set color and marker style based on the file index
        if (file_index == i) {
            graph->SetMarkerColor(colors[i]);
            graph->SetMarkerStyle(20+i); //shapes 
            leg->AddEntry(graph, names[i], "p");
        } 
        }
         
        // Set the points and the errors (only y-error, x-error is 0)
        for (int i = 0; i <= n_points; ++i) {
            graph->SetPoint(i, x_data[i], y_data[i]);
            graph->SetPointError(i,0,y_error[i]);  // Only y-error, x-error is 0
        }

        // Draw the graph (with option "P" to plot points)
        if (file_index == 0) {
            graph->Draw("AP");  // "A" means Axis and "P" means Points (with markers)
        } else {
            graph->Draw("P");  // Draw subsequent graphs with only points (no axes)
            
        }

    graph->GetXaxis()->SetLimits(0.175, 3.2);  // X-axis range
    graph->SetMinimum(0.22);  // Y-axis lower limit (prevent log scale issues)
    graph->SetMaximum(1e5);   // Y-axis upper limit
   

 if (file_index == 0) {
        tsallisFunc->SetParLimits(0,0.01,4.0); // Initial guess for A 
        tsallisFunc->SetParLimits(1,1.11, 1.15);//Initial guess for q 
        tsallisFunc->SetParLimits(2, 0.3,0.55); // Initial guess for T 
        tsallisFunc->SetParLimits(3, 0.01,3.0); // Initial guess for mu 
        tsallisFunc->FixParameter(4, 0.0); // Initial guess for y 
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
        
        } else if (file_index == 1) {
       tsallisFunc->SetParLimits(0,0.1,4.0); 
        tsallisFunc->SetParLimits(1,1.11, 1.145);
        tsallisFunc->SetParLimits(2, 0.2,0.45); 
        tsallisFunc->SetParLimits(3, 0.01,3.0); 
        tsallisFunc->FixParameter(4, 0.0); 
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
        }
         else if(file_index == 2) {
        tsallisFunc->SetParLimits(0,0.1,4.0); 
        tsallisFunc->SetParLimits(1,1.1, 1.145);
        tsallisFunc->SetParLimits(2, 0.2,0.45); 
        tsallisFunc->SetParLimits(3, 0.01,3.0); 
        tsallisFunc->FixParameter(4, 0.0); 
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
           }
           else if(file_index == 3 ) {
        tsallisFunc->SetParLimits(0,0.1,4.0); 
        tsallisFunc->SetParLimits(1,1.1, 1.145);
        tsallisFunc->SetParLimits(2, 0.35,0.4); 
        tsallisFunc->SetParLimits(3, 0.1,3.0); 
        tsallisFunc->FixParameter(4, 0.0); 
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
           }else if( file_index == 4) {
        tsallisFunc->SetParLimits(0,0.1,4.0); 
        tsallisFunc->SetParLimits(1,1.1, 1.135);
        tsallisFunc->SetParLimits(2, 0.2,0.35); 
        tsallisFunc->SetParLimits(3, 0.1,3.0); 
        tsallisFunc->FixParameter(4, 0.0);
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
           }
            else {
        tsallisFunc->SetParLimits(0,0.1,4.0); 
        tsallisFunc->SetParLimits(1,1.1, 1.135);
        tsallisFunc->SetParLimits(2, 0.2,0.35); 
        tsallisFunc->SetParLimits(3, 0.1,3.0); 
        tsallisFunc->FixParameter(4, 0.0);
        tsallisFunc->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)","ripdity");
        }
        
         
    
       
     graph->Fit(tsallisFunc, "R");  // Fit the graph with the Tsallis function


        
       
       graph->SetTitle("Xe Xe --> Charged X ,  #sqrt{s_{NN}} =5.44 TeV "); 
       }
     
      
        
    leg->AddEntry(tsallisFunc, "Tsallis Fit", "l");
    leg->Draw();
}
