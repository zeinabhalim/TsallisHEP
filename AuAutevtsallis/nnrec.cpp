#include "TROOT.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TSystem.h"
#include "TApplication.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

// Define the Tsallis function as a return value
double TsallisFun(double *x, double *par) 
{
    double mass = 0.13957018; // Constant mass term
    double xval = x[0]; //transverse momentum term
    double A = par[0];   // Normalization factor
    double q = par[1];   // Tsallis non-equilibrium index
    double T = par[2];   // Temperature term
    double mu = par[3];   // Tsallis index
     double y = par[4];   // rapdity term


    return par[0] * xval * sqrt(xval * xval + mass * mass) * 
           TMath::CosH(par[4]) * 
           TMath::Power(1 + ((par[1] - 1)*( 1/ par[2]) * 
           (sqrt(xval * xval + mass * mass) * TMath::CosH(par[4]) - par[3])), 
           (1 / (1-par[1])));
}

void nnrec() {
    // Create a new canvas 
    TCanvas *canvas = new TCanvas("canvas", "Multiple Graphs with Tsallis Fit Function", 1300, 1000);
      
    canvas->SetLogy(true);  // Logarithmic Y-axis
    canvas->SetLogx(true);  // Logarithmic X-axis
        
    TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.85); // Left side
    leg->SetBorderSize(1);  
    leg->SetTextSize(0.02); 
  
    double xmin = 0.25, xmax = 1.98;
    TF1 *tsallisFun2 = new TF1("fitFunc", TsallisFun, xmin, xmax, 5);
  //  tsallisFun2->SetNpx(1000); 

    // List of input filenames
    std::vector<std::string> filenames = {
        "auau1.txt", //"auau2.txt", 
        "auau3.txt" , //"auau4.txt",
        "auau5.txt", "auau6.txt", "auau7.txt", "auau8.txt", "auau9.txt", "auau10.txt" //"auau11.txt"
    };

    // Arrays to store fitting results
    std::string centrality_labels[12];
    double norm_values[12], norm_errors[12];
    double q_values[12], q_errors[12];
    double T_values[12], T_errors[12];
    double mu_values[12], mu_errors[12];
    double chi2_values[12], chi2ndf_values[12];
    int ndf_values[12];

    // Save parameters to a text file
    std::ofstream outfile("Tsallisfit_parametersauau200.txt");
    outfile << "Fitting of the transverse momentum spectra of charged particles in Au-Au at √s_NN = 200 GeV, using Tsallis distribution function." << std::endl;
    outfile << std::endl;

    // Loop over each file and process the data
    for (size_t file_index = 0; file_index < filenames.size(); ++file_index) {
        // Vectors to hold data points and errors for the current dataset
        std::vector<double> x_data, y_data, y_stat_error;

        const std::string& filename = filenames[file_index];

        // Open the input file
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;  // Skip to the next file if the current one can't be opened
        }
        
        // Read data from the file (x, y, stat_error)
        double x, y, ey;
        std::cout << "Reading file: " << filename << std::endl;
        int point_count = 0;
        
        while (input_file >> x >> y >> ey) {
            x_data.push_back(x);
            y_data.push_back(y);
            y_stat_error.push_back(ey);
            point_count++;
        }
        
        // Close the input file after reading
        input_file.close();
        
        std::cout << "  Read " << point_count << " data points" << std::endl;

        // Check if data is available before creating the graph
        if (x_data.empty()) {
            std::cerr << "No data found in file: " << filename << std::endl;
            continue;
        }

        // Create a TGraphErrors to store the data with errors
        int n_points = x_data.size();
        TGraphErrors *graph = new TGraphErrors(n_points);

        // Set title and labels for the graph
        graph->SetTitle(" ;P_{t}[GeV];(1/(N_{evt}))*D^2(N)/d\\eta /dP_{t} ");
        tsallisFun2->SetLineColor(kBlack);
        
        const char* names[12] = {"0-5%", //"5-10%", 
        "10-15%", //"15-20%",
                               "20-30%", "30-40%", "40-50%", "50-60%","60-70%","70-80%","80-92%"};
        int colors[12] = {kRed, kBlue, kGreen + 2, kMagenta,
                        kOrange + 7, kCyan + 2, kViolet, kBlack, kRed + 2};
        
        // Set color and marker style based on the file index
        graph->SetMarkerColor(colors[file_index]);
        graph->SetMarkerStyle(20 + file_index);
        leg->AddEntry(graph, names[file_index], "p");
        
        // Set the points and the errors
        for (int i = 0; i < n_points; ++i) {
            graph->SetPoint(i, x_data[i], y_data[i]);
            graph->SetPointError(i, 0.0, y_stat_error[i]);  // x-error is 0, y-error is stat_error
        }
        
        // Draw the graph
        if (file_index == 0) {
            graph->Draw("AP");  // "A" means Axis and "P" means Points (with markers)
        } else {
            graph->Draw("P");  // Draw subsequent graphs with only points (no axes)
        }
        
        graph->GetXaxis()->SetLimits(0.22, 2.2);  // X-axis range
        graph->SetMinimum(1e-5);  // Y-axis lower limit (prevent log scale issues)
        graph->SetMaximum(1e4);   // Y-axis upper limit
        
        // Set parameter limits based on centrality
        if (file_index == 0) {
            // 0-5%
            tsallisFun2->SetParLimits(0, 0.01, 0.15);
            tsallisFun2->SetParLimits(1, 1.1, 1.114);
            tsallisFun2->SetParLimits(2, 0.15, 0.26);
            tsallisFun2->SetParLimits(3, 0.01, 2.0);
            tsallisFun2->FixParameter(4, 0.0);
        } else if (file_index == 1) {
            tsallisFun2->SetParLimits(0, 0.01, 0.25);
            tsallisFun2->SetParLimits(1, 1.0, 1.114);
            tsallisFun2->SetParLimits(2, 0.15, 0.26);
            tsallisFun2->SetParLimits(3, 0.1, 2.0);
            tsallisFun2->FixParameter(4, 0.0);
        } else if (file_index == 2) {
            tsallisFun2->SetParLimits(0, 0.01, 0.35);
            tsallisFun2->SetParLimits(1, 1.1, 1.114);
            tsallisFun2->SetParLimits(2, 0.15, 0.25);
            tsallisFun2->SetParLimits(3, 1.0, 2.0);
            tsallisFun2->FixParameter(4, 0.0);
        } else if (file_index == 3) {
            tsallisFun2->SetParLimits(0, 0.1, 0.35);
            tsallisFun2->SetParLimits(1, 1.1, 1.1145);
            tsallisFun2->SetParLimits(2, 0.15, 0.25);
            tsallisFun2->SetParLimits(3, 1.0, 2.0);
            tsallisFun2->FixParameter(4, 0.0);
        } else if (file_index == 4) {
            tsallisFun2->SetParLimits(0, 0.01, 0.3);
            tsallisFun2->SetParLimits(1, 1.1, 1.1145);
            tsallisFun2->SetParLimits(2, 0.15, 0.2);
            tsallisFun2->SetParLimits(3, 1.0, 2.0);
            tsallisFun2->FixParameter(4, 0.0);
        } else if (file_index == 5) {
            tsallisFun2->SetParLimits(0, 0.01, 0.3);
            tsallisFun2->SetParLimits(1, 1.1, 1.1145);
            tsallisFun2->SetParLimits(2, 0.15, 0.19);
            tsallisFun2->SetParLimits(3, 0.01, 1.5);
            tsallisFun2->FixParameter(4, 0.0);
            }
            else if (file_index == 6) {
            tsallisFun2->SetParLimits(0, 0.001, 0.2);
            tsallisFun2->SetParLimits(1, 1.1, 1.1155);
            tsallisFun2->SetParLimits(2, 0.15, 0.175);
            tsallisFun2->SetParLimits(3, 0.01, 1.5);
            tsallisFun2->FixParameter(4, 0.0);
        } else {
            tsallisFun2->SetParLimits(0, 0.001, 0.2);
            tsallisFun2->SetParLimits(1, 1.1, 1.1155);
            tsallisFun2->SetParLimits(2, 0.15, 0.17);
            tsallisFun2->SetParLimits(3, 0.01, 1.5);
            tsallisFun2->FixParameter(4, 0.0);
        }
        
        tsallisFun2->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)", "rapidity");
        
        // Perform the fit
        graph->Fit(tsallisFun2, "R");

        // Get Chi-Square and Degrees of Freedom
        double chi2 = tsallisFun2->GetChisquare();
        int ndf = tsallisFun2->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 999;

        // Store results in arrays
        centrality_labels[file_index] = names[file_index];
        norm_values[file_index] = tsallisFun2->GetParameter(0);
        norm_errors[file_index] = tsallisFun2->GetParError(0);
        q_values[file_index] = tsallisFun2->GetParameter(1);
        q_errors[file_index] = tsallisFun2->GetParError(1);
        T_values[file_index] = tsallisFun2->GetParameter(2);
        T_errors[file_index] = tsallisFun2->GetParError(2);
        mu_values[file_index] = tsallisFun2->GetParameter(3);
        mu_errors[file_index] = tsallisFun2->GetParError(3);
        chi2_values[file_index] = chi2;
        chi2ndf_values[file_index] = chi2ndf;
        ndf_values[file_index] = ndf;

        // Print out the fitting parameters 
        outfile << "==================================================================" << std::endl;
        outfile << "Centrality: " << names[file_index] << std::endl;
        outfile << "Data points: " << n_points << std::endl;
        outfile << "==================================================================" << std::endl;
        
        outfile << std::fixed << std::setprecision(6);
        outfile << "Fit Parameters:" << std::endl;
        outfile << "  A (Normalization): " << tsallisFun2->GetParameter(0) 
                << " ± " << tsallisFun2->GetParError(0) << std::endl;
        outfile << "  q (Non-extensive index): " << tsallisFun2->GetParameter(1)
                << " ± " << tsallisFun2->GetParError(1) << std::endl;
        outfile << "  T (Temperature): " << tsallisFun2->GetParameter(2)
                << " ± " << tsallisFun2->GetParError(2) << std::endl;
        outfile << "  μ (Chemical pot.): " << tsallisFun2->GetParameter(3)
                << " ± " << tsallisFun2->GetParError(3) << std::endl;
        outfile << "  Rapidity (fixed): " << tsallisFun2->GetParameter(4) << " (fixed)" << std::endl;
        
        outfile << std::endl << "Goodness of Fit:" << std::endl;
        outfile << "  χ²: " << chi2 << std::endl;
        outfile << "  NDF: " << ndf << std::endl;
        outfile << "  χ²/NDF: " << chi2ndf << std::endl;
        
       
        outfile << std::endl;
        
        graph->SetTitle("Au Au --> Charged X ,  #sqrt{s_{NN}} =200 GeV ");
    }
    
    // Define the legend
    leg->AddEntry(tsallisFun2, "Tsallis Fit", "l");
    
    // Draw legend
    leg->Draw("frame");
    outfile.close();
    
    // ========== CREATE SUMMARY TABLES ==========
    
    // 1. Create summary text table
    std::ofstream summaryfile("Tsallis_fit_summary_auau.txt");
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << "SUMMARY OF TSALLIS FIT PARAMETERS - Au-Au √s_NN = 200 GeV" << std::endl;
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << std::endl;
    summaryfile << std::left << std::setw(10) << "Centrality" 
                    << std::setw(7) << "  " 
                << std::setw(27) << "Norm" 
                << std::setw(20) << "q" 
                << std::setw(17) << "T (GeV)" 
                << std::setw(21) << "μ (GeV)"
                << std::setw(5) << "χ²/NDF" << std::endl;
    summaryfile << std::string(135, '-') << std::endl;
    
    for (size_t i = 0; i < filenames.size(); i++) {
        // Format each parameter with ± error
        std::stringstream norm_str, q_str, T_str, mu_str;
        norm_str << std::fixed << std::setprecision(6) << norm_values[i] << " ± " << norm_errors[i];
        q_str << std::fixed << std::setprecision(6) << q_values[i] << " ± " << q_errors[i];
        T_str << std::fixed << std::setprecision(4) << T_values[i] << " ± " << T_errors[i];
        mu_str << std::fixed << std::setprecision(4) << mu_values[i] << " ± " << mu_errors[i];
        
        summaryfile << std::left 
                    << std::setw(10) << centrality_labels[i]
                    << std::setw(25) << norm_str.str()
                    << std::setw(25) << q_str.str()
                    << std::setw(20) << T_str.str()
                    << std::setw(25) << mu_str.str()
                    << std::setw(30) << std::fixed << std::setprecision(4) << chi2ndf_values[i] 
                    << std::endl;
                
    }
    
    // 2. Create CSV file
    std::ofstream csvfile("Tsallis_fit_parameters_auau.csv");
    csvfile << "Centrality,Norm,Norm_error,q,q_error,T,T_error,mu,mu_error,"
            << "chi2,chi2ndf,NDF" << std::endl;
    
    for (size_t i = 0; i < filenames.size(); i++) {
        csvfile << std::fixed << std::setprecision(8);
        csvfile << centrality_labels[i] << ","
                << norm_values[i] << "," << norm_errors[i] << ","
                << q_values[i] << "," << q_errors[i] << ","
                << T_values[i] << "," << T_errors[i] << ","
                << mu_values[i] << "," << mu_errors[i] << ","
                << chi2_values[i] << "," << chi2ndf_values[i] << ","
                << ndf_values[i] << std::endl;
    }
    csvfile.close();

    // Save the canvas
    canvas->SaveAs("Tsallis_fit_AuAu_200GeV.png");
 
}
