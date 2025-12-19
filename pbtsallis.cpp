#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
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

void pbtsallis() {
    
    TCanvas *canvas = new TCanvas("canvas", "Pb-Pb with Tsallis Fit Function", 1300, 1000);
      
    canvas->SetLogy(true);  
    canvas->SetLogx(true);  

    //locating the legend box       
    TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.9); 
    leg->SetBorderSize(1);  
    leg->SetTextSize(0.02); 
   
  
    // loading the input data files
    std::vector<std::string> filenames = {
        "pb1.txt", "pb2.txt", "pb3.txt", "pb4.txt",
        "pb5.txt", "pb6.txt", "pb7.txt", "pb8.txt"
    };
        
    double xmin = 0.175, xmax = 2.1;
    TF1 *TsallisFitpb = new TF1("fitFunc", TsallisnFit, xmin, xmax, 5);
    TsallisFitpb->SetNpx(1000);

    // Arrays to store all fitting results
    std::string centrality_labels[8];
    double norm_values[8], norm_errors[8];
    double q_values[8], q_errors[8];
    double T_values[8], T_errors[8];
    double mu_values[8], mu_errors[8];
    double chi2_total_values[8], chi2ndf_total_values[8];
    int ndf_values[8];
    
    // Save parameters to a text file
    std::ofstream outfile("Tsallisfit_parameterspbpb.txt");
    outfile << "Fitting of the transverse momentum spectra of charged particles in Pb-Pb at √s_NN = 2.76 TeV" << std::endl;
    outfile << "Using Tsallis distribution function with combined errors (statistical + systematic)" << std::endl;
    outfile << "Data format: x  y  stat_error  sys_error" << std::endl;
    outfile << std::endl;

    // Loop over each file and process the data
    for (size_t file_index = 0; file_index < filenames.size(); ++file_index) {
        std::vector<double> x_data, y_data, y_stat_error, y_sys_error;

        const std::string& filename = filenames[file_index];

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

        // Create graph 
        TGraphErrors *graph = CreateGraphWithTotalErrors(x_data, y_data, y_stat_error, y_sys_error);

        // Set title and labels
        graph->SetTitle(" ;P_{t}[GeV];(1/(N_{evt}))*d^2(N)/d\\eta /dP _{t } ");
        TsallisFitpb->SetLineColor(kBlack);
         
        const char* names[8]={"0-5%", "5-10%", "10-20%", "20-30%",
                             "30-40%", "40-50%", "50-60%", "60-70%"};
        int colors[8] = {kRed, kBlue, kGreen + 2, kMagenta,
                        kOrange + 7, kCyan + 2, kViolet, kBlack};
        
        // Set color and marker style
        graph->SetMarkerColor(colors[file_index]);
        graph->SetMarkerStyle(20 + file_index);
        graph->SetMarkerSize(0.9);
        leg->AddEntry(graph, names[file_index], "p");
         
        // Draw the graph
        if (file_index == 0) {
            graph->Draw("AP");
        } else {
            graph->Draw("P");
        }

        graph->GetXaxis()->SetLimits(0.155, 2.5);
        graph->SetMinimum(0.22);
        graph->SetMaximum(1e5);
    
        // Set parameter limits based on centrality
        if (file_index == 0) {
            // 0-5%
            TsallisFitpb->SetParLimits(0, 0.01, 11); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.145);
            TsallisFitpb->SetParLimits(2, 0.2, 0.46); 
            TsallisFitpb->SetParLimits(3,0.01, 3.5);  
            TsallisFitpb->FixParameter(4, 0.0); // Fix rapidity parameter
        } else if (file_index == 1) {
            // 5-10%
          TsallisFitpb->SetParLimits(0, 0.01, 10); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.145);
            TsallisFitpb->SetParLimits(2, 0.2, 0.46); 
            TsallisFitpb->SetParLimits(3,0.01, 3.5);  
            TsallisFitpb->FixParameter(4, 0.0);
        } else if (file_index == 2) {
            // 10-20%
             TsallisFitpb->SetParLimits(0, 0.01,10); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.145);
            TsallisFitpb->SetParLimits(2, 0.2, 0.45); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0);  
            TsallisFitpb->FixParameter(4, 0.0);
        } else if (file_index == 3) {
            // 20-30%
              TsallisFitpb->SetParLimits(0, 0.01, 10); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.146);
            TsallisFitpb->SetParLimits(2, 0.2, 0.45); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0);  
            TsallisFitpb->FixParameter(4, 0.0);
        } else if (file_index == 4) {
            // 30-40%
           TsallisFitpb->SetParLimits(0, 0.01, 10.0); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.146);
            TsallisFitpb->SetParLimits(2, 0.2, 0.45); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0); 
            TsallisFitpb->FixParameter(4, 0.0);
        } else if (file_index == 5) {
            // 40-50%
            TsallisFitpb->SetParLimits(0, 0.01, 10.0); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.146);
            TsallisFitpb->SetParLimits(2, 0.2, 0.45); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0); 
            TsallisFitpb->FixParameter(4, 0.0);
        } else if (file_index == 6) {
            // 50-60%
         TsallisFitpb->SetParLimits(0, 0.01, 3.0); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.146);
            TsallisFitpb->SetParLimits(2, 0.2, 0.4); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0); 
            TsallisFitpb->FixParameter(4, 0.0);
        } else {
            // 60-70%
            TsallisFitpb->SetParLimits(0, 0.01, 3.0); 
            TsallisFitpb->SetParLimits(1, 1.1, 1.146);
            TsallisFitpb->SetParLimits(2, 0.2, 0.4); 
            TsallisFitpb->SetParLimits(3,0.01, 3.0); 
            TsallisFitpb->FixParameter(4, 0.0);
        }
        
        // Set parameter names
        TsallisFitpb->SetParNames("Norm", "q (Non-Extensive)", "T (Temperature)", "mu (Chemical Pot.)", "rapidity");
        
        // Perform the fit
        graph->Fit(TsallisFitpb, "R");

        // Calculate χ² with total errors
        double chi2_total = 0.0;
        int ndf = TsallisFitpb->GetNDF();
        
        for (int i = 0; i < graph->GetN(); ++i) {
            double x_val, y_val;
            graph->GetPoint(i, x_val, y_val);
            double model_val = TsallisFitpb->Eval(x_val);
            double error = graph->GetErrorY(i);  // Combined stat + sys error
            chi2_total += ((y_val - model_val) * (y_val - model_val)) / (error * error);
        }
        
        double chi2ndf_total = chi2_total / ndf;

        // Store results in arrays
        centrality_labels[file_index] = names[file_index];
        norm_values[file_index] = TsallisFitpb->GetParameter(0);
        norm_errors[file_index] = TsallisFitpb->GetParError(0);
        q_values[file_index] = TsallisFitpb->GetParameter(1);
        q_errors[file_index] = TsallisFitpb->GetParError(1);
        T_values[file_index] = TsallisFitpb->GetParameter(2);
        T_errors[file_index] = TsallisFitpb->GetParError(2);
        mu_values[file_index] = TsallisFitpb->GetParameter(3);
        mu_errors[file_index] = TsallisFitpb->GetParError(3);
        chi2_total_values[file_index] = chi2_total;
        chi2ndf_total_values[file_index] = chi2ndf_total;
        ndf_values[file_index] = ndf;

        // Print results
        outfile << "==================================================================" << std::endl;
        outfile << "Centrality: " << names[file_index] << std::endl;
        outfile << "Data points: " << graph->GetN() << std::endl;
        outfile << "==================================================================" << std::endl;
        
        outfile << std::fixed << std::setprecision(6);
        outfile << "Fit Parameters:" << std::endl;
        outfile << "  Norm:            " << TsallisFitpb->GetParameter(0) 
                << " ± " << TsallisFitpb->GetParError(0) << std::endl;
        outfile << "  q (Non-extensive): " << TsallisFitpb->GetParameter(1)
                << " ± " << TsallisFitpb->GetParError(1) << std::endl;
        outfile << "  T (Temperature):   " << TsallisFitpb->GetParameter(2)
                << " ± " << TsallisFitpb->GetParError(2) << std::endl;
        outfile << "  μ (Chemical Pot.): " << TsallisFitpb->GetParameter(3)
                << " ± " << TsallisFitpb->GetParError(3) << std::endl;
        outfile << "  Rapidity (fixed):  " << TsallisFitpb->GetParameter(4) << " (fixed)" << std::endl;
        
        outfile << std::endl << "Goodness of Fit:" << std::endl;
        outfile << "  χ² (stat + sys errors):       " << chi2_total << std::endl;
        outfile << "  χ²/NDF (total errors):        " << chi2ndf_total << std::endl;
        outfile << "  NDF:                          " << ndf << std::endl;
        
        outfile << std::endl << "Fit Quality Assessment:" << std::endl;
        outfile << "  Using total errors (stat+sys): ";
        if (chi2ndf_total < 1.5) {
            outfile << "Excellent (χ²/NDF < 1.5)" << std::endl;
        } else if (chi2ndf_total < 3.0) {
            outfile << "Good (χ²/NDF < 3.0)" << std::endl;
        } else if (chi2ndf_total < 5.0) {
            outfile << "Acceptable (χ²/NDF < 5.0)" << std::endl;
        } else {
            outfile << "Poor (χ²/NDF ≥ 5.0)" << std::endl;
        }
        outfile << std::endl;
    }
     
    leg->AddEntry(TsallisFitpb, "Tsallis Fit", "l");
    leg->Draw();
    
    
    // Add title to the plot
    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextSize(0.035);
    title->DrawLatex(0.15, 0.96, "Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, |#eta| < 0.8");

    
    // Save the canvas
    canvas->SaveAs("Tsallis_fit_PbPb.png");
    canvas->Draw("APL");
    
    outfile.close();
    
    // ========== CREATE TABLES WITH ERROR VALUES ==========
    
    // 1. Create summary table with ± errors
    std::ofstream summaryfile("Tsallis_fit_summary.txt");
    summaryfile << "===============================================================================" << std::endl;
    summaryfile << "TSALLIS FIT PARAMETERS - Pb-Pb √s_NN = 2.76 TeV" << std::endl;
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
                    << std::setw(30) << std::fixed << std::setprecision(4) << chi2ndf_total_values[i] 
                    << std::endl;
    }
    summaryfile.close();
    
    // Create CSV file with separate columns for values and errors
    std::ofstream csvfile("Tsallis_fit_parameters.csv");
    csvfile << "Centrality,Norm,Norm_error,q,q_error,T,T_error,mu,mu_error,"
            << "chi2_total,chi2ndf_total,NDF" << std::endl;
    
    for (size_t i = 0; i < filenames.size(); i++) {
        csvfile << std::fixed << std::setprecision(8);
        csvfile << centrality_labels[i] << ","
                << norm_values[i] << "," << norm_errors[i] << ","
                << q_values[i] << "," << q_errors[i] << ","
                << T_values[i] << "," << T_errors[i] << ","
                << mu_values[i] << "," << mu_errors[i] << ","
                << chi2_total_values[i] << "," << chi2ndf_total_values[i] << ","
                << ndf_values[i] << std::endl;
    }
    csvfile.close();

}
