#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  

// the Tsallis function as a return value
double TsallisnFit(double *x, double *par) 
{
    double mass = 0.13957018; // mass term of produced particle(consider the mass of pion particle)
    double xval = x[0]; //the transverse momentum 
    double A = par[0];   // Normalization factor
    double q = par[1];   // non-equilibrium tsallis index
    double T = par[2];   // average Temperature at freeze-out stage
    double mu = par[3];   // chemical potential
    double y = par[4];   // rapdity parameter

    return par[0] * xval * sqrt(xval * xval + mass * mass) * 
           TMath::CosH(par[4]) * 
           TMath::Power(1 + ((par[1] - 1)*( 1/ par[2]) * 
           (sqrt(xval * xval + mass * mass) * TMath::CosH(par[4]) - par[3])), 
           (1 / (1-par[1])));
}

// Errors Function (statistical + systematic)
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

void auaugevfit () {
    
    // Create a canvas divided into 6 pads (3x2)
    TCanvas *c1 = new TCanvas("canvas", "Au-Au with Tsallis Fit Function", 1200, 800); 
    c1->Divide(3, 2); 

    // loading the input data files as 2D vector
    std::vector<std::vector<std::string>> fileGroups = {
        { "auau14gev1.txt", "auau19gev1.txt", "auau62gev1.txt" }, 
        { "auau14gev2.txt", "auau19gev2.txt", "auau62gev2.txt" }, 
        { "auau14gev3.txt", "auau19gev3.txt", "auau62gev3.txt" }, 
        { "auau14gev4.txt", "auau19gev4.txt", "auau62gev4.txt" } ,
        { "auau14gev5.txt", "auau19gev5.txt", "auau62gev5.txt" }, 
        { "auau14gev6.txt", "auau19gev6.txt", "auau62gev6.txt" }
    };
    
    // Energy and centrality labels
    const char* energy_names[3] = {"14.6 GeV", "19.6 GeV", "62.4 GeV"};
    const char* energy_labels_short[3] = {"14GeV", "19GeV", "62GeV"};
    
    const char* cent_ranges[6] = {
        "0-5%", "5-10%", "10-20%", 
        "20-30%", "30-40%", "40-50%"
    };
    
    // Colors and markers for different energies
    int colors[3] = {kRed, kBlue, kGreen + 2};
    int markers[3] = {20, 21, 22};
    
    // Arrays to store all fitting results
    const int total_files = 6 * 3; // 6 centrality bins × 3 energies = 18 files
    std::string cent_labels[total_files];
    std::string energy_labels[total_files];
    double norm_values[total_files], norm_errors[total_files];
    double q_values[total_files], q_errors[total_files];
    double T_values[total_files], T_errors[total_files];
    double mu_values[total_files], mu_errors[total_files];
    double chi2_total_values[total_files], chi2ndf_total_values[total_files];
    int ndf_values[total_files];
    
    double xmin = 0.55, xmax = 2.1;
    TF1 *TsallisFit = new TF1("fitFunc", TsallisnFit, xmin, xmax, 5);
    TsallisFit->SetNpx(1000);
    
    // File streams for output
    std::ofstream outfile("Tsallisfit_parameters_auau.txt");
    std::ofstream summaryfile("Tsallis_fit_summary_auau.txt");
    std::ofstream csvfile("Tsallis_fit_parameters_auau.csv");
    
    outfile << "Fitting of the transverse momentum spectra in Au-Au collisions" << std::endl;
    outfile << "Using Tsallis distribution function with combined errors (statistical + systematic)" << std::endl;
    outfile << std::endl;
    
    // Write summary header ONLY ONCE
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << "TSALLIS FIT PARAMETERS - Au-Au COLLISIONS" << std::endl;
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << "√s = 14.6, 19.6, 62.4 GeV" << std::endl;
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << std::endl;
    
    csvfile << "Energy,cent_range,Norm,Norm_error,q,q_error,T,T_error,mu,mu_error,chi2_total,chi2ndf_total,NDF" << std::endl;
    
    int file_counter = 0;
    
   
    // Loop over each centrality bin (6 groups = 6 pads)
    for (size_t cent_index = 0; cent_index < fileGroups.size(); ++cent_index) {
        c1->cd(cent_index + 1);
        
        // Create TMultiGraph for this pad
        TMultiGraph *mg = new TMultiGraph();
        TLegend *leg = new TLegend(0.5, 0.6, 0.85, 0.9);
        leg->SetBorderSize(2);
        leg->SetTextSize(0.05);
        
        // Set log scale for Y-axis only
        gPad->SetLogy();
        
        // Add centrality label INSIDE each graph
        TLatex *centLabel = new TLatex();
        centLabel->SetNDC();
        centLabel->SetTextAlign(22);
        centLabel->SetTextSize(0.08);
        centLabel->SetTextColor(kRed);
        centLabel->DrawLatex(0.5, 0.85, cent_ranges[cent_index]);
        
        // Loop over each energy in this centrality bin
        for (size_t energy_index = 0; energy_index < fileGroups[cent_index].size(); ++energy_index) {
            
            std::vector<double> x_data, y_data, y_stat_error, y_sys_error;
            const std::string& filename = fileGroups[cent_index][energy_index];

            // Open the file
            std::ifstream input_file(filename);
            if (!input_file.is_open()) {
                std::cerr << "Error opening file: " << filename << std::endl;
                continue;  
            }

            double x, y, ey_stat, ey_sys;
            std::cout << "Reading file: " << filename << std::endl;
            int point_count = 0;
            
            while (input_file >> x >> y >> ey_stat >> ey_sys) {
                x_data.push_back(x);
                y_data.push_back(y);
                y_stat_error.push_back(ey_stat);
                y_sys_error.push_back(ey_sys);
                point_count++;
            }
            input_file.close();
            
            std::cout << "  Read " << point_count << " data points" << std::endl;

            // Create graph with combined errors
            TGraphErrors *graph = CreateGraphWithTotalErrors(x_data, y_data, y_stat_error, y_sys_error);

            // Set style based on energy
            graph->SetMarkerColor(colors[energy_index]);
            graph->SetMarkerStyle(markers[energy_index]);
            graph->SetLineColor(colors[energy_index]);
            graph->SetMarkerSize(0.8);
            
            // Add to legend
            leg->AddEntry(graph, energy_names[energy_index], "p");
            
            // Set parameter limits based on centrality
            if (cent_index == 0) {
                TsallisFit->SetParLimits(0, 0.01, 2.0); 
                TsallisFit->SetParLimits(1, 1.0, 1.08);
                TsallisFit->SetParLimits(2, 0.1, 0.2); 
                TsallisFit->SetParLimits(3, 0.001, 1.5);  
                TsallisFit->FixParameter(4, 0.0);
            } else if (cent_index == 1) {
                TsallisFit->SetParLimits(0, 0.01, 2.0); 
                TsallisFit->SetParLimits(1, 1.0, 1.08);
                TsallisFit->SetParLimits(2, 0.1, 0.2); 
                TsallisFit->SetParLimits(3, 0.001, 1.5); 
                TsallisFit->FixParameter(4, 0.0);
            } else if (cent_index == 2) {
                TsallisFit->SetParLimits(0, 0.01, 2.0); 
                TsallisFit->SetParLimits(1, 1.0, 1.08);
                TsallisFit->SetParLimits(2, 0.1, 0.2); 
                TsallisFit->SetParLimits(3, 0.001, 1.5); 
                TsallisFit->FixParameter(4, 0.0);
            } else if (cent_index == 3) {
                TsallisFit->SetParLimits(0, 0.01, 2.0); 
                TsallisFit->SetParLimits(1, 1.0, 1.08);
                TsallisFit->SetParLimits(2, 0.1, 0.2); 
                TsallisFit->SetParLimits(3, 0.001, 1.5); 
                TsallisFit->FixParameter(4, 0.0);
            } else {
                TsallisFit->SetParLimits(0, 0.01, 2.0); 
                TsallisFit->SetParLimits(1, 1.0, 1.08);
                TsallisFit->SetParLimits(2, 0.1, 0.2); 
                TsallisFit->SetParLimits(3, 0.001, 1.5); 
                TsallisFit->FixParameter(4, 0.0);
            }
            
            TsallisFit->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)", "rapidity");
            
            // Perform the fit
            graph->Fit(TsallisFit, "R");
            
            // Calculate χ² with total errors
            double chi2_total = 0.0;
            int ndf = TsallisFit->GetNDF();
            
            for (int j = 0; j < graph->GetN(); ++j) {
                double x_val, y_val;
                graph->GetPoint(j, x_val, y_val);
                double model_val = TsallisFit->Eval(x_val);
                double error = graph->GetErrorY(j);
                chi2_total += ((y_val - model_val) * (y_val - model_val)) / (error * error);
            }
            
            double chi2ndf_total = chi2_total / ndf;
            
            // Store results
            int index = file_counter;
            cent_labels[index] = cent_ranges[cent_index];
            energy_labels[index] = energy_labels_short[energy_index];
            norm_values[index] = TsallisFit->GetParameter(0);
            norm_errors[index] = TsallisFit->GetParError(0);
            q_values[index] = TsallisFit->GetParameter(1);
            q_errors[index] = TsallisFit->GetParError(1);
            T_values[index] = TsallisFit->GetParameter(2);
            T_errors[index] = TsallisFit->GetParError(2);
            mu_values[index] = TsallisFit->GetParameter(3);
            mu_errors[index] = TsallisFit->GetParError(3);
            chi2_total_values[index] = chi2_total;
            chi2ndf_total_values[index] = chi2ndf_total;
            ndf_values[index] = ndf;
            
            // Write detailed output
            outfile << "==================================================================" << std::endl;
            outfile << "File: " << filename << std::endl;
            outfile << "Energy: " << energy_names[energy_index] << ", Centrality: " << cent_ranges[cent_index] << std::endl;
            outfile << "Data points: " << graph->GetN() << std::endl;
            outfile << "==================================================================" << std::endl;
            
            outfile << std::fixed << std::setprecision(6);
            outfile << "Fit Parameters:" << std::endl;
            outfile << "  Norm:            " << TsallisFit->GetParameter(0) 
                    << " ± " << TsallisFit->GetParError(0) << std::endl;
            outfile << "  q (Non-extensive): " << TsallisFit->GetParameter(1)
                    << " ± " << TsallisFit->GetParError(1) << std::endl;
            outfile << "  T (Temperature):   " << TsallisFit->GetParameter(2)
                    << " ± " << TsallisFit->GetParError(2) << std::endl;
            outfile << "  μ (Chemical Pot.): " << TsallisFit->GetParameter(3)
                    << " ± " << TsallisFit->GetParError(3) << std::endl;
            outfile << "  Rapidity (fixed):  " << TsallisFit->GetParameter(4) << " (fixed)" << std::endl;
            
            outfile << std::endl << "Goodness of Fit:" << std::endl;
            outfile << "  χ²: " << chi2_total << std::endl;
            outfile << "  χ²/NDF: " << chi2ndf_total << std::endl;
            outfile << "  NDF: " << ndf << std::endl;
            
            outfile << std::endl;
            
            // Add graph to multigraph
            mg->Add(graph, "P");
            
            file_counter++;
        }
        
        // Add fit function to legend (once per pad)
        TsallisFit->SetLineColor(kRed);
        TsallisFit->SetLineWidth(2);
        leg->AddEntry(TsallisFit, "Tsallis Fit", "l");
        
        // Set multigraph properties
        mg->SetTitle(cent_ranges[cent_index]);
        mg->GetXaxis()->SetLimits(0.45, 2.3);
        mg->SetMinimum(0.0001);
        mg->SetMaximum(1e2);
        mg->GetXaxis()->SetTitle("P_{t} [GeV]");
        mg->GetYaxis()->SetTitle("(1/(N_{evt}))*D^2(N)/d\\eta /dP_{t}");
        
        // Draw multigraph and legend
        mg->Draw("AP");
        leg->Draw();
    }
    
    // ========== WRITE SUMMARY TABLE (ONLY DATA, NO HEADERS) ==========
    // Write table headers ONCE
    summaryfile << std::left << std::setw(15) << "Energy" 
                << std::setw(15) << "Centrality" 
                << std::setw(25) << "Norm" 
                << std::setw(25) << "q" 
                << std::setw(20) << "T (GeV)" 
                << std::setw(25) << "μ (GeV)"
                << std::setw(15) << "χ²/NDF" << std::endl;

    
    // Write all data in one continuous table
    for (int i = 0; i < file_counter; i++) {
        std::stringstream norm_str, q_str, T_str, mu_str;
        norm_str << std::fixed << std::setprecision(4) << norm_values[i] << " ± " << norm_errors[i];
        q_str << std::fixed << std::setprecision(4) << q_values[i] << " ± " << q_errors[i];
        T_str << std::fixed << std::setprecision(4) << T_values[i] << " ± " << T_errors[i];
        mu_str << std::fixed << std::setprecision(4) << mu_values[i] << " ± " << mu_errors[i];
            summaryfile << "=========================================================================================================" << std::endl;
        summaryfile << std::left 
                    << std::setw(15) << energy_labels[i]
                    << std::setw(15) << cent_labels[i]
                    << std::setw(25) << norm_str.str()
                    << std::setw(25) << q_str.str()
                    << std::setw(20) << T_str.str()
                    << std::setw(25) << mu_str.str()
                    << std::setw(15) << std::fixed << std::setprecision(4) << chi2ndf_total_values[i] 
                    << std::endl;
    }
    
    // Write to CSV file
    for (int i = 0; i < file_counter; i++) {
        csvfile << std::fixed << std::setprecision(8);
        csvfile << energy_labels[i] << ","
                << "\"" << cent_labels[i] << "\"" << ","
                << norm_values[i] << "," << norm_errors[i] << ","
                << q_values[i] << "," << q_errors[i] << ","
                << T_values[i] << "," << T_errors[i] << ","
                << mu_values[i] << "," << mu_errors[i] << ","
                << chi2_total_values[i] << "," << chi2ndf_total_values[i] << ","
                << ndf_values[i] << std::endl;
    }
    
    // Close files
    outfile.close();
    summaryfile.close();
    csvfile.close();
    
    // Save canvas
    c1->SaveAs("Tsallis_fit_AuAu_multi.png");
    
}
