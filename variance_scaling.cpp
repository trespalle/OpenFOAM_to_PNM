#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <filesystem> 
#include <sstream> 
#include <ctime>   
#include "pore_networks.h"
#include "myfunctions.h" 
#include "Rand.h"  

namespace fs = std::filesystem;

// =======================================================================
// HELPER: Calculate Normalized Variance (CV^2 = Sigma^2 / Mu^2)
// =======================================================================
double calculate_normalized_variance(const std::vector<double>& data)
{
    if(data.empty()) return 0.0;
    
    double mean = 0.0, var = 0.0;
    // Compute mean and variance using the custom library (using 100% of data)
    myfun::med_var(data, mean, var, 99.5);

    // Avoid division by zero
    if(std::abs(mean) < 1e-20) return 0.0; 

    // Return squared coefficient of variation
    return var / (mean * mean);
}

// =======================================================================
// HELPER: Save Histogram with dynamic filename based on k
// =======================================================================
void save_tagged_histogram(const std::vector<double>& data, const std::string& base_name, double k_val, const std::string& out_dir)
{
    if(data.empty()) return;
    
    // Construct filename: e.g., "hist_flow_k_0.50.dat"
    std::stringstream ss;
    ss << base_name << "_k_" << std::fixed << std::setprecision(2) << k_val << ".dat";
    std::string filename = ss.str();

    // Work with absolute values for histograms
    std::vector<double> data_abs;
    data_abs.reserve(data.size());
    for(double v : data) data_abs.push_back(std::abs(v));

    double min_val, max_val, delta;
    // Build histogram with 40 bins for higher resolution
    auto h = myfun::buildHistogram(data_abs, 25, min_val, max_val, delta, 0.995); 
    myfun::guardaHistograma(h, min_val, delta, out_dir + filename);
}

// =======================================================================
// MAIN PROGRAM
// =======================================================================
int main(int argc, char** argv) 
{
    std::cout << "==========================================" << std::endl;
    std::cout << "      VARIANCE SCALING ANALYSIS (CV^2)    " << std::endl;
    std::cout << "==========================================" << std::endl;

    
    if (argc < 4) {
        std::cerr << "Error: Missing arguments.\n";
        std::cerr << "Use: " << argv[0] << " <output_directory> <links_file> <coords_file>\n";
        return 1;
    }

    std::string output_dir = argv[1];
    if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/'; 
    
    std::string link_file = argv[2];
    std::string coord_file = argv[3];

    if (!fs::exists(output_dir)) fs::create_directories(output_dir);

    // Define file schemas
    std::vector<LinkData> link_fmt = { LinkData::Weight, LinkData::HalfWidth, LinkData::FlowRate };
    std::vector<NodeData> node_fmt = { NodeData::Pressure };

    // -----------------------------------------------------------------------
    // PHASE 1: SETUP & BASELINE CALIBRATION
    // -----------------------------------------------------------------------
    std::cout << "[SETUP] Loading Topology and calibrating Baseline Mean..." << std::endl;
    
    // Load real network structure
    PoreNetwork PN_base(link_file, coord_file, link_fmt, node_fmt, 2);
    
    
    // Calculate effective conductances (inverse problem from exp data)
    PN_base.update_pressure_drops(); 
    PN_base.init_physics_from_experimental_data();

    // Back-calculate effective halfwidths to obtain the target physical mean
    PN_base.calculate_halfwidths_from_geometric_conductances(ModelType::Poiseuille, 0.001);

    // Compute baseline statistics (mean halfwidth)
    double mean_eff = 0.0, var_eff = 0.0;
    myfun::med_var(PN_base.effective_halfwidths, mean_eff, var_eff, 99.5);
    
    std::cout << "  Baseline mean halfwidth (mu): " << mean_eff << " m" << std::endl;
    
    // Calculate target total flow 
    // We compute the average experimental flow rate to ensure physical consistency 
    // across all synthetic simulations.
    double total_inlet_Q = 0.0;
    for (std::size_t id : PN_base.inlet) 
        for (std::size_t l_idx : PN_base.nodes[id].outlnks) total_inlet_Q += std::abs(PN_base.flows[l_idx]);
    
    double total_outlet_Q = 0.0;
    for (std::size_t id : PN_base.outlet) 
        for (std::size_t l_idx : PN_base.nodes[id].inlnks) total_outlet_Q += std::abs(PN_base.flows[l_idx]);

    double avg_Q_target = (total_inlet_Q + total_outlet_Q) / 2.0;
    std::cout << "  Target total flow rate: " << avg_Q_target << " m^3/s" << std::endl;


    // -----------------------------------------------------------------------
    // PHASE 2: THE LOOP (Variance scaling)
    // -----------------------------------------------------------------------
    // List of target shape parameters (k). Recall: NormVar = 1/k

    std::size_t nr_targets = 20;
    double half_width_lowest = 0.1;
    double half_width_increment = 0.1;
    std::vector<double> half_width_targets = { half_width_lowest };
    half_width_targets.reserve(nr_targets);
    for(std::size_t ii = 1; ii < nr_targets; ++ii)
        half_width_targets.push_back(half_width_targets[ii-1] + half_width_increment);

    std::vector<double> k_targets;
    k_targets.reserve(half_width_targets.size());
    for (auto val : half_width_targets)
        k_targets.push_back(1./val);


    for (auto val : k_targets)
        std::cout << val << std::endl;
    
    
    // Vectors to store the final curve points for Gnuplot
    std::vector<double> x_vals; // Normalized variance of halfwidths (input)
    std::vector<double> y_vals; // Normalized variance of pore flow rates (output)

    std::cout << "\n[LOOP] Starting simulation loop over k values..." << std::endl;

    // Use a robust Mersenne Twister RNG seeded with time
    RNG rng(123456789);

    for(double k_target : k_targets)
    {
        std::cout << "  -> Processing k_target = " << k_target << "..." << std::endl;

        // Create a fresh synthetic network structure (a copy of the real one)
        PoreNetwork PN_sim(link_file, coord_file, link_fmt, node_fmt, 2);

        // Calculate gamma parameters for halfwidth sampling
        // Keep mean = baseline, but vary variance via k_target
        // Relation: mu = k * theta  =>  theta = mu / k
        double theta = mean_eff / k_target;

        // Sample synthetic halfwidths from gamma distribution
        for(std::size_t i = 0; i < PN_sim.links.size(); ++i)
        {
            double w = rng.gamma(k_target, theta);
            // Physical clamp to avoid singularity (closed tubes)
            if(w < 1e-7) w = 1e-7; 
            PN_sim.effective_halfwidths[i] = w;
        }

        // Calculate input statistic (actual sampled normalized variance)
        // Note: we calculate what we actually sampled, rather than assuming 1/k_target
        double norm_var_in = calculate_normalized_variance(PN_sim.effective_halfwidths);
        x_vals.push_back(norm_var_in);


        // Solve direct physics 
        // Geometry -> Conductance
        PN_sim.set_geometric_conductances_from_geometry(ModelType::Poiseuille, 0.001);
        
        // Solve Pressure & Flow
        PN_sim.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);
        
        
        // Extract pore flows (throughput)
        std::vector<double> pore_masses;
        pore_masses.reserve(PN_sim.nodes.size());
        
        for(const auto& node : PN_sim.nodes) 
        {
            // Filter out nodes with negligible mass (disconnected or dead-ends)
            if(node.mass > 0.0) pore_masses.push_back(node.mass);
        }

        // Calculate output statistic (normalized variance of pore flows)
        double norm_var_out = calculate_normalized_variance(pore_masses);
        y_vals.push_back(norm_var_out);

        std::cout << "     [RESULT] In_NormVar (Geometry)=" << norm_var_in 
                  << " | Out_NormVar (Transport)=" << norm_var_out << std::endl;


        // Save histograms for this iteration (tagged with k value)
        save_tagged_histogram(PN_sim.flows, "hist_flow", k_target, output_dir);
        save_tagged_histogram(pore_masses, "hist_pore_flows", k_target, output_dir);
        save_tagged_histogram(PN_sim.pressure_drops, "hist_pressure_drop", k_target, output_dir);
        save_tagged_histogram(PN_sim.geometric_conductances, "hist_conductance", k_target, output_dir);
        save_tagged_histogram(PN_sim.effective_halfwidths, "hist_halfwidth", k_target, output_dir);
        
        // Extract current physical splitting fractions
        std::vector<double> current_weights;
        for(const auto& l : PN_sim.links) current_weights.push_back(l.weight);
        save_tagged_histogram(current_weights, "hist_splitting_fractions", k_target, output_dir);
    }

    // -----------------------------------------------------------------------
    // PHASE 3: SAVE SCALING RESULTS
    // -----------------------------------------------------------------------
    std::cout << "\n[OUTPUT] Saving variance scaling curve..." << std::endl;
    std::string curve_file = output_dir + "variance_scaling_curve.dat";
    std::ofstream fout(curve_file);
    
    if(fout.is_open())
    {
        fout << "# NormVar_Halfwidths(Input)\tNormVar_PoreMasses(Output)\n";
        for(size_t i = 0; i < x_vals.size(); ++i)
        {
            fout << x_vals[i] << "\t" << y_vals[i] << "\n";
        }
        fout.close();
        std::cout << "  Curve saved to: " << curve_file << std::endl;
    }
    else
    {
        std::cerr << "Error: Could not save scaling curve." << std::endl;
    }

    std::cout << "[DONE] Analysis complete." << std::endl;

    

    return 0;
}