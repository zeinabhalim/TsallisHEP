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

// Tsallis function as a return value
double Tsallis(double *x, double *par) 
{
    double mass = 0.13957018; // Constant mass 
    double xval = x[0]; // transverse momentum variable
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

// Function to create graph with combined errors
TGraphErrors* CreateGraphWithTotalErrors(const std::vector<double>& x_data,
                                       const std::vector<double>& y_data,
                                       const std::vector<double>& y_stat_error,
                                       const std::vector<double>& y_sys_error) {
    int n_points = x_data.size();
    TGraphErrors* graph = new TGraphErrors(n_points);
    
    for (int i = 0; i < n_points; ++i) {
        graph->SetPoint(i, x_data[i], y_data[i]);
        
        // Combine statistical and systematic errors in quadrature
        double stat_err = y_stat_error[i];
        double sys_err = y_sys_error[i];
        double total_err = sqrt(stat_err*stat_err + sys_err*sys_err);
        
        graph->SetPointError(i, 0.0, total_err);
    }
    
    return graph;
}

void ppgevtsallis() {
    // Create a new canvas to hold the graphs together
    TCanvas *c1 = new TCanvas("canvas", "Multiple Graphs with Tsallis Fit Function", 1000, 900); 
    c1->Divide(2, 2); 

    // Define multiple sets of data files 
    std::vector<std::vector<const char*>> fileGroups = 
       { { //"pp6gev1.txt",
       "pp7gev1.txt",
       //"pp8gev1.txt",
       "pp12gev1.txt", "pp17gev1.txt"}, 
        { //"pp6gev2.txt",
       "pp7gev2.txt",
       //"pp8gev2.txt",
       "pp12gev2.txt", "pp17gev2.txt"}, 
        { //"pp6gev3.txt",
       "pp7gev3.txt",
       //"pp8gev3.txt",
       "pp12gev3.txt", "pp17gev3.txt"}, 
        {//"pp6gev4.txt",
       "pp7gev4.txt",
       //"pp8gev4.txt",
       "pp12gev4.txt", "pp17gev4.txt"} };
        
    // domain limits of fit function    
    double xmin = 0.123, xmax = 1.4;
    TF1 *tsallisFunp = new TF1("fitFunc", Tsallis, xmin, xmax, 5);

    // setting Colors and marker styles
    int colors[] = {kBlue, kGreen+2, kBlack, kBlue, kGreen+2, kBlack, 
                   kBlue, kGreen+2, kBlack, kBlue, kGreen+2, kBlack};
    int markers[] = {20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30};
    
    // Arrays to store all fitting results
    const int total_entries = 4 * 3; // 4 eta bins × 3 energies
    std::string energy_labels[total_entries];
    std::string eta_labels[total_entries];
    double norm_values[total_entries], norm_errors[total_entries];
    double q_values[total_entries], q_errors[total_entries];
    double T_values[total_entries], T_errors[total_entries];
    double mu_values[total_entries], mu_errors[total_entries];
    double chi2_values[total_entries], chi2ndf_values[total_entries];
    int ndf_values[total_entries];
    
    // File streams for output
    std::ofstream outfile("Tsallisfit_parametersppgev.txt");
    std::ofstream summaryfile("Tsallis_fit_summary_pp.txt");
    std::ofstream csvfile("Tsallis_fit_parameters_pp.csv");
    
    outfile << "Fitting of the transverse momentum spectra of charged particles in proton-proton" << std::endl;
    outfile << "using Tsallis distribution function with combined errors (statistical + systematic)" << std::endl;
    outfile << "Data format: x y stat_error sys_error" << std::endl;
    outfile << std::endl;
    
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << "SUMMARY OF TSALLIS FIT PARAMETERS - PROTON-PROTON COLLISIONS" << std::endl;
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << std::endl;
    summaryfile << std::left << std::setw(15) << "Energy" 
                << std::setw(15) << "η Range" 
                << std::setw(20) << "Norm" 
                << std::setw(20) << "q" 
                << std::setw(20) << "T (GeV)" 
                << std::setw(20) << "μ (GeV)"
                << std::setw(15) << "χ²/NDF" << std::endl;
    summaryfile << std::string(130, '-') << std::endl;
    
    csvfile << "Energy,eta_range,Norm,Norm_error,q,q_error,T,T_error,mu,mu_error,chi2,chi2ndf,NDF" << std::endl;
    
    // Loop over each group of data sets  
    for (size_t groupIndex = 0; groupIndex < fileGroups.size(); groupIndex++) {
        c1->cd(groupIndex + 1); // Move to the canvas partition

        TMultiGraph *mg = new TMultiGraph();
        TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.9);
        gPad->SetLogy(); // Apply log scale to Y-axis
        gPad->SetLogx(); // Apply log scale to x-axis
        
        const char *titles[4] = {
       "0.0< #eta<0.2",
        "0.2< #eta<0.4",
        "0.4< #eta<0.6",
        "0.6< #eta<0.8"
    };
   const char *energyLabels[3] = {
        //"#sqrt{s_{NN}} = 6.3 GeV",
        "#sqrt{s_{NN}} = 7.7 GeV",
       // "#sqrt{s_{NN}} = 8.8 GeV",
        "#sqrt{s_{NN}} = 12.3 GeV",
        "#sqrt{s_{NN}} = 17.3 GeV"
    };
      
        
        // Loop over each file in each group in data set (3 files per group)
        for (size_t i = 0; i < fileGroups[groupIndex].size(); i++) {
            const char* filename = fileGroups[groupIndex][i];
            std::vector<double> x_vals, y_vals, y_stat_errs, y_sys_errs;

            // Read the data file (x, y, stat_error, sys_error)
            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Error: Cannot open file " << filename << std::endl;
                continue;
            }

            double x, y, yerr_stat, yerr_sys;
            while (file >> x >> y >> yerr_stat >> yerr_sys) {
                x_vals.push_back(x);
                y_vals.push_back(y ); 
                y_stat_errs.push_back(yerr_stat);
                y_sys_errs.push_back(yerr_sys);
            }
            file.close();
            
            // Create graph with combined errors
            TGraphErrors *gr = CreateGraphWithTotalErrors(x_vals, y_vals, y_stat_errs, y_sys_errs);
            
            int dataIndex = groupIndex * 3 + i;
            gr->SetMarkerStyle(markers[dataIndex]);
            gr->SetMarkerColor(colors[dataIndex]);
            gr->SetLineColor(colors[dataIndex]);

            // Set parameter limits based on group
            if (groupIndex == 0) {
                tsallisFunp->SetParLimits(0, 0.01, 3.0);
                tsallisFunp->SetParLimits(1, 0.95, 1.09);
                tsallisFunp->SetParLimits(2, 0.12, 0.35);
                tsallisFunp->SetParLimits(3, 0.4, 3.0);
                tsallisFunp->FixParameter(4, 0.0);
            } else if (groupIndex == 1) {
                tsallisFunp->SetParLimits(0, 0.01, 3.0);
                tsallisFunp->SetParLimits(1, 0.95, 1.09);
                tsallisFunp->SetParLimits(2, 0.12, 0.35);
                tsallisFunp->SetParLimits(3, 0.4, 3.0);
                tsallisFunp->FixParameter(4, 0.3);
            } else if (groupIndex == 2) {
                tsallisFunp->SetParLimits(0, 0.01, 3.0);
                tsallisFunp->SetParLimits(1, 0.95, 1.09);
                tsallisFunp->SetParLimits(2, 0.12, 0.35);
                tsallisFunp->SetParLimits(3, 0.4, 3.0);
                tsallisFunp->FixParameter(4, 0.5);
            } else {
                tsallisFunp->SetParLimits(0, 0.01, 3.0);
                tsallisFunp->SetParLimits(1, 0.95, 1.09);
                tsallisFunp->SetParLimits(2, 0.12, 0.35);
                tsallisFunp->SetParLimits(3, 0.4, 3.0);
                tsallisFunp->FixParameter(4, 0.7);
            }
            
            tsallisFunp->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)", "rapidity");
            
            // Perform the fit
            gr->Fit(tsallisFunp, "R");
            
            // Calculate χ² with total errors
            double chi2_total = 0.0;
            int ndf = tsallisFunp->GetNDF();
            
            for (int j = 0; j < gr->GetN(); ++j) {
                double x_val, y_val;
                gr->GetPoint(j, x_val, y_val);
                double model_val = tsallisFunp->Eval(x_val);
                double error = gr->GetErrorY(j);  // Combined stat + sys error
                chi2_total += ((y_val - model_val) * (y_val - model_val)) / (error * error);
            }
            
            double chi2ndf_total = chi2_total / ndf;
            
            // Store results
            int index = groupIndex * 3 + i;
            energy_labels[index] = energyLabels[i];
            eta_labels[index] = titles[groupIndex];
            norm_values[index] = tsallisFunp->GetParameter(0);
            norm_errors[index] = tsallisFunp->GetParError(0);
            q_values[index] = tsallisFunp->GetParameter(1);
            q_errors[index] = tsallisFunp->GetParError(1);
            T_values[index] = tsallisFunp->GetParameter(2);
            T_errors[index] = tsallisFunp->GetParError(2);
            mu_values[index] = tsallisFunp->GetParameter(3);
            mu_errors[index] = tsallisFunp->GetParError(3);
            chi2_values[index] = chi2_total;
            chi2ndf_values[index] = chi2ndf_total;
            ndf_values[index] = ndf;
            
            // Write to detailed output file
            outfile << "==================================================================" << std::endl;
            outfile << "Energy: " << energyLabels[i] << ", Pseudorapidity: " << titles[groupIndex] << std::endl;
            outfile << "Data points: " << gr->GetN() << std::endl;
            outfile << "==================================================================" << std::endl;
            
            outfile << std::fixed << std::setprecision(6);
            outfile << "Fit Parameters:" << std::endl;
            outfile << "  Norm:            " << tsallisFunp->GetParameter(0) 
                    << " ± " << tsallisFunp->GetParError(0) << std::endl;
            outfile << "  q (Non-extensive): " << tsallisFunp->GetParameter(1)
                    << " ± " << tsallisFunp->GetParError(1) << std::endl;
            outfile << "  T (Temperature):   " << tsallisFunp->GetParameter(2)
                    << " ± " << tsallisFunp->GetParError(2) << std::endl;
            outfile << "  μ (Chemical Pot.): " << tsallisFunp->GetParameter(3)
                    << " ± " << tsallisFunp->GetParError(3) << std::endl;
            outfile << "  Rapidity (fixed):  " << tsallisFunp->GetParameter(4) << std::endl;
            
            outfile << std::endl << "Goodness of Fit:" << std::endl;
            outfile << "  χ²/NDF:                " << chi2ndf_total << std::endl;
            outfile << "  NDF:                   " << ndf << std::endl;
            
      
            outfile << std::endl;
            
            // Add to multigraph
            mg->SetTitle(titles[groupIndex]);
            leg->AddEntry(gr, energyLabels[i], "P");
            mg->Add(gr, "P");
            mg->GetXaxis()->SetLimits(0.11, 1.5);
            mg->SetMinimum(0.0001);
            mg->SetMaximum(0.7e1);
            gr->SetMarkerSize(0.7);
            mg->GetXaxis()->SetTitle("P_{t}[GeV]");
        mg->GetYaxis()->SetTitle("(1/(N_{evt}))*d^2(N)/d\\eta /dP _{t }");  
        }
        
        tsallisFunp->SetLineColor(kRed);
        tsallisFunp->SetLineWidth(2);
        leg->AddEntry(tsallisFunp, "Tsallis Fit", "l");
        
        mg->Draw("APL");
        leg->Draw();
    }
    
    // Write summary table
    for (int i = 0; i < total_entries; i++) {
        // Format for summary text file
        std::stringstream norm_str, q_str, T_str, mu_str;
        norm_str << std::fixed << std::setprecision(4) << norm_values[i] << " ± " << norm_errors[i];
        q_str << std::fixed << std::setprecision(4) << q_values[i] << " ± " << q_errors[i];
        T_str << std::fixed << std::setprecision(4) << T_values[i] << " ± " << T_errors[i];
        mu_str << std::fixed << std::setprecision(4) << mu_values[i] << " ± " << mu_errors[i];
        
        summaryfile << std::left 
                    << std::setw(15) << energy_labels[i]
                    << std::setw(15) << eta_labels[i]
                    << std::setw(20) << norm_str.str()
                    << std::setw(20) << q_str.str()
                    << std::setw(20) << T_str.str()
                    << std::setw(20) << mu_str.str()
                    << std::setw(15) << std::fixed << std::setprecision(2) << chi2ndf_values[i] 
                    << std::endl;
        
        // Write to CSV file
        csvfile << std::fixed << std::setprecision(8);
        csvfile << energy_labels[i] << ","
                << eta_labels[i] << ","
                << norm_values[i] << "," << norm_errors[i] << ","
                << q_values[i] << "," << q_errors[i] << ","
                << T_values[i] << "," << T_errors[i] << ","
                << mu_values[i] << "," << mu_errors[i] << ","
                << chi2_values[i] << "," << chi2ndf_values[i] << ","
                << ndf_values[i] << std::endl;
    }
    
    // Close files
    outfile.close();
    summaryfile.close();
    csvfile.close();
    
    // Save canvas
    c1->SaveAs("Tsallis_fit_pp_multienergies.png");
    
}
