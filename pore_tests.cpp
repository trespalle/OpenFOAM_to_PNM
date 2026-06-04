#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <random>
#include <cmath>
#include <iomanip>
#include <filesystem> 
#include "pore_networks.h"
#include "myfunctions.h"
#include "Rand.h"   

namespace fs = std::filesystem;

void randomize_network_weights(PoreNetwork& PN)
{
    std::cout << "  -> Randomizing splitting fractions (Null Model U(0,1))..." << std::endl;
    
    RNG rng(123456789); 


    for(std::size_t i = 0; i < PN.nodes.size(); ++i) 
    {
        auto& node = PN.nodes[i];
        std::size_t k = node.outlnks.size();
        
        if(k == 0) continue; // Outlet or dead-end

        std::vector<double> new_weights(k);

        // RANDOMIZATION LOGIC
        if(k == 1) 
        {
            new_weights[0] = 1.0;
        }
        else if(k == 2)
        {
            // Binary split: w and 1-w using your RNG
            double w = rng(); 
            new_weights[0] = w;
            new_weights[1] = 1.0 - w;
        }
        else 
        {
            // General case k > 2
            double sum = 0.0;
            for(std::size_t j=0; j<k; ++j) 
            {
                new_weights[j] = rng(); // Llamada al operator() de tu clase
                sum += new_weights[j];
            }
            // Normalize
            if(sum > 1e-9) for(std::size_t j=0; j<k; ++j) new_weights[j] /= sum;
            else for(std::size_t j=0; j<k; ++j) new_weights[j] = 1.0/k;   
        }

        // SYNCHRONIZATION (update node lists aND links list) 
        for(std::size_t j=0; j<k; ++j) 
        {
            double w = new_weights[j];
            std::size_t link_idx = node.outlnks[j];
            std::size_t neighbor = node.outnbrs[j];

            // Update node outgoing weight
            node.outwgs[j] = w;

            // Update global link list
            PN.links[link_idx].weight = w;

            // Update neighbor incoming weight
            bool found_link = false;
            for(std::size_t m=0; m < PN.nodes[neighbor].inlnks.size(); ++m) 
            {
                if(PN.nodes[neighbor].inlnks[m] == link_idx) 
                {
                    PN.nodes[neighbor].inwgs[m] = w;
                    found_link = true;
                    break;
                }
            }
            if(!found_link) std::cerr << "Warning: Sync issue at node " << i << " link " << link_idx << std::endl;
        }
    }
    
    // Refresh internal Eigen matrix for robustness
    PN.refresh_weights();
    std::cout << "Randomization complete and synchronized using Mersenne Twister." << std::endl;
}

// Helper to save Linear histograms (for Weights/Fractions 0-1)
void save_histogram(const std::vector<double>& data, const std::string& filename, const std::string& out_dir, bool use_abs = true) 
{
    if(data.empty()) return;
    
    // Create a local copy to process (abs) and sort (inside buildHistogram)
    std::vector<double> data_processed;
    data_processed.reserve(data.size());

    if (use_abs) 
    {
        // Apply absolute value to all elements
        for(double val : data) data_processed.push_back(std::abs(val));
    } 
    else 
    {
        // Just copy raw values
        data_processed = data;
    }

    double min_val, max_val, delta;
    
    // Assuming myfun::buildHistogram sorts the vector in-place and bins it
    auto h = myfun::buildHistogram(data_processed, 25, min_val, max_val, delta, 0.995);
    
    myfun::guardaHistograma(h, min_val, delta, out_dir + filename);
    std::cout << "  Saved " << filename << (use_abs ? " (Abs Values)" : " (Raw Values)") << std::endl;
}

void test_gamma_implementation(const std::string& out_dir)
{
    std::cout << "\n[TEST] Verifying Custom Gamma RNG (Marsaglia & Tsang)..." << std::endl;
    
    // Arbitrary Theoretical Parameters
    double k_shape = 0.5; 
    double theta_scale = 1.5;
    
    double theo_mean = k_shape * theta_scale;           // 2.5 * 1.5 = 3.75
    double theo_var  = k_shape * theta_scale * theta_scale; // 2.5 * 2.25 = 5.625

    std::cout << "  Parameters: Shape (k)=" << k_shape << ", Scale (theta)=" << theta_scale << std::endl;
    std::cout << "  Theoretical: Mean=" << theo_mean << ", Variance=" << theo_var << std::endl;

    // 2. Massive Sampling
    std::size_t N = 1000000; // 1 million points
    std::vector<double> samples;
    samples.reserve(N);

    RNG rng(12345); // Fixed seed for reproducibility

    for(std::size_t i = 0; i < N; ++i)
    {
        // Calling the manual gamma implementation
        samples.push_back(rng.gamma(k_shape, theta_scale)); 
    }

    // 3. Statistical Calculation
    double samp_mean = 0.0, samp_var = 0.0;
    myfun::med_var(samples, samp_mean, samp_var, 100.0);

    std::cout << "  Sample (N=1e6): Mean=" << samp_mean << ", Variance=" << samp_var << std::endl;
    
    // Relative Error
    double err_mean = std::abs(samp_mean - theo_mean) / theo_mean * 100.0;
    double err_var = std::abs(samp_var - theo_var) / theo_var * 100.0;
    std::cout << "  Error: Mean " << std::fixed << std::setprecision(4) << err_mean << "%, Var " << err_var << "%" << std::endl;

    // 4. Save Histogram for Gnuplot verification
    save_histogram(samples, "test_gamma_distribution.dat", out_dir, false); 
}

int main(int argc, char** argv) 
{
    std::cout << "==========================================" << std::endl;
    std::cout << "      PORE NETWORK COMPLETE ANALYSIS      " << std::endl;
    std::cout << "==========================================" << std::endl;

    if (argc < 4) {
        std::cerr << "Error: Arguments missing.\n";
        std::cerr << "Use: " << argv[0] << " <output_directory> <links_file> <coords_file>\n";
        return 1;
    }

    const double MU_DNS = 1e6;

    std::string output_dir = argv[1];
    if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/'; 
    
    std::string link_file = argv[2];
    std::string coord_file = argv[3];

    if (!fs::exists(output_dir)) fs::create_directories(output_dir);

    // Define file schemas
    std::vector<LinkData> link_fmt = { LinkData::Weight, LinkData::HalfWidth, LinkData::FlowRate };
    std::vector<NodeData> node_fmt = { NodeData::Pressure };

    // =======================================================================
    // PHASE 1: LOADING EXPERIMENTAL TOPOLOGY & "REAL/GEOMETRIC" DATA
    // =======================================================================
    std::cout << "\n[PHASE 1] Loading Experimental Topology & Data..." << std::endl;
    PoreNetwork PN(link_file, coord_file, link_fmt, node_fmt, 2);
    PN.viscosity = MU_DNS; 

    
    // Save REAL histograms (directly from file)
    // -----------------------------------------------------
    // Real half-widths (a_exp) 
    save_histogram(PN.effective_halfwidths, "hist_halfwidth_exp.dat", output_dir);

    // print some statistics:
    double mean_exp = 0.0, var_exp = 0.0;
    myfun::med_var(PN.effective_halfwidths, mean_exp, var_exp, 99.5); 
    std::cout << "  [STATS EXP] Halfwidth -> Mean: " << mean_exp << " | Variance: " << var_exp << std::endl;
    // -----------------------------
    
    // Real splitting fractions (weights)
    // We assume the weights loaded from file are the "real" ones.
    std::vector<double> real_weights;
    for(const auto& l : PN.links) real_weights.push_back(l.weight);
    save_histogram(real_weights, "hist_splitting_fractions_exp.dat", output_dir);

    // Q_exp and dP_exp (for reference)
    PN.update_pressure_drops(); // Ensure cache is updated
    save_histogram(PN.flows, "hist_flow_rate_exp.dat", output_dir);
    save_histogram(PN.pressure_drops, "hist_pressure_drop_exp.dat", output_dir);

    // Calculate GEOMETRIC conductances (g_geom) -> Requested in point (2)
    // ----------------------------------------------------------------------
    // We use Hele-Shaw model with H=1mm on the REAL apertures we just loaded.
    std::cout << "\n[PHASE 1b] Calculating geometric conductances (Hele-Shaw)..." << std::endl;
    PN.set_geometric_conductances_from_geometry(ModelType::Poiseuille, 0.001); // H = 1mm
    
    // Save histogram (to compare with g_eff later)
    save_histogram(PN.geometric_conductances, "hist_conductance_geom.dat", output_dir);


    // =======================================================================
    // PHASE 2: INVERSE PROBLEM (EFFECTIVE CALIBRATION)
    // =======================================================================
    std::cout << "\n[PHASE 2] Inverse problem: calculating effective conductances..." << std::endl;
    
    // IMPORTANT: This overwrites the geometric conductances calculated above
    // with the "effective" ones derived from Q_exp / dP_exp.
    // Ensure your function uses the Q and P currently loaded.
    PN.init_physics_from_experimental_data(); 

    // Save EFFECTIVE conductances (g_eff) 
    save_histogram(PN.geometric_conductances, "hist_conductance_eff.dat", output_dir);

    // Calculate EFFECTIVE half-widths (a_eff)
    // -------------------------------------------------------------------
    // We invert the Hele-Shaw model using the current effective conductances.
    std::cout << "          Back-calculating effective halfwidths..." << std::endl;
    PN.calculate_halfwidths_from_geometric_conductances(ModelType::Poiseuille, 0.001);

    // Save histogram a_eff
    save_histogram(PN.effective_halfwidths, "hist_halfwidth_eff.dat", output_dir);

    // print some statistics
    double mean_eff = 0.0, var_eff = 0.0;
    myfun::med_var(PN.effective_halfwidths, mean_eff, var_eff, 99.5);
    std::cout << "  [STATS EFF] Halfwidth -> Mean: " << mean_eff << " | Variance: " << var_eff << std::endl;
    // -----------------------------


    // =======================================================================
    // PHASE 3: FORWARD SIMULATION (DIRECT PROBLEM)
    // =======================================================================
    std::cout << "\n[PHASE 3] Forward simulation using EFFECTIVE conductances..." << std::endl;

    // Calculate target flow rate (average of experimental Inlet/Outlet)
    double total_inlet_Q = 0.0;
    for (std::size_t id : PN.inlet) 
        for (std::size_t l_idx : PN.nodes[id].outlnks) total_inlet_Q += std::abs(PN.flows[l_idx]);
    
    double total_outlet_Q = 0.0;
    for (std::size_t id : PN.outlet) 
        for (std::size_t l_idx : PN.nodes[id].inlnks) total_outlet_Q += std::abs(PN.flows[l_idx]);

    double avg_Q_target = (total_inlet_Q + total_outlet_Q) / 2.0;
    std::cout << "  Target Q: " << avg_Q_target << " m^3/s" << std::endl;

    // -----------------------------------------------------------------------
    // 3A. Solve the flow field using EFFECTIVE conductances (from Phase 2)
    // -----------------------------------------------------------------------
    PN.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);

    // Save Simulation Results (Effective)
    save_histogram(PN.flows, "hist_flow_rate_sim_eff.dat", output_dir);
    save_histogram(PN.pressure_drops, "hist_pressure_drop_sim_eff.dat", output_dir);

    // Save simulated splitting Fractions (Effective)
    std::vector<double> sim_weights_eff;
    for(const auto& l : PN.links) sim_weights_eff.push_back(l.weight);
    save_histogram(sim_weights_eff, "hist_splitting_fractions_sim_eff.dat", output_dir);


    // -----------------------------------------------------------------------
    // 3B. Solve the flow field using REAL/GEOMETRIC conductances
    // -----------------------------------------------------------------------
    std::cout << "\n[PHASE 3B] Forward simulation using REAL/GEOMETRIC conductances..." << std::endl;
    
    // Reload a fresh network to get the real (experimental) half-widths intact
    PoreNetwork PN_exp(link_file, coord_file, link_fmt, node_fmt, 2);
    PN_exp.viscosity = MU_DNS;

    // Set conductances based purely on geometry
    PN_exp.set_geometric_conductances_from_geometry(ModelType::Poiseuille, 0.001); // H = 1mm

    // Solve the flow field using the exact same target boundary condition
    PN_exp.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);

    // Save Simulation Results (Experimental/Geometric)
    save_histogram(PN_exp.flows, "hist_flow_rate_sim_exp.dat", output_dir);
    save_histogram(PN_exp.pressure_drops, "hist_pressure_drop_sim_exp.dat", output_dir);

    // Save simulated splitting Fractions (Experimental/Geometric)
    std::vector<double> sim_weights_exp;
    for(const auto& l : PN_exp.links) sim_weights_exp.push_back(l.weight);
    save_histogram(sim_weights_exp, "hist_splitting_fractions_sim_exp.dat", output_dir);


    // =============================================
    // PHASE 4: TOPOLOGICAL INVERSE PROBLEM (MaxEnt) 
    // =============================================
    std::cout << "\n[PHASE 4] Topological inverse problem (pure splitting fractions)..." << std::endl;
    
    // CRITICAL: We reload the network from scratch.
    // Why? Because 'solve_flow_field' in Phase 3 modified the link weights 
    // based on the physical flow physics. We want to use the "real" (file) 
    // topological weights for this prediction.
    PoreNetwork PN_topo(link_file, coord_file, link_fmt, node_fmt, 2);
    PN_topo.viscosity = MU_DNS;

    // We need to define boundary pressures for the topological method.
    // We will use the average of the experimental pressures (loaded at start in PN_topo 
    // from the file).
    double P_in_avg = 0.0;
    if(!PN_topo.inlet.empty()) 
    {
        for(auto id : PN_topo.inlet) P_in_avg += PN_topo.pressures[id];
        P_in_avg /= PN_topo.inlet.size();
    }

    double P_out_avg = 0.0;
    if(!PN_topo.outlet.empty()) 
    {
        for(auto id : PN_topo.outlet) P_out_avg += PN_topo.pressures[id];
        P_out_avg /= PN_topo.outlet.size();
    }

    std::cout << "  Using avg P_in=" << P_in_avg << " and P_out=" << P_out_avg << " for topology initialization." << std::endl;

    // Execute Topological Algorithm using the NEW method signature
    // 1. Calculates Q based on topological weights (Mass Balance)
    // 2. Calculates P based on Laplacian (Harmonic Field)
    // 3. Calibrates g based on Q/dP
    PN_topo.init_physics_from_splitting_fractions(avg_Q_target, P_in_avg, P_out_avg);

    // Save Topological Histograms
    // Note: The method already updates the pressure_drops cache internally, 
    // so the vectors are ready to use.
    save_histogram(PN_topo.geometric_conductances, "hist_conductance_topo.dat", output_dir);
    save_histogram(PN_topo.flows, "hist_flow_rate_topo.dat", output_dir);
    save_histogram(PN_topo.pressure_drops, "hist_pressure_drop_topo.dat", output_dir);
    
    // Optional: Calculate widths from these topological conductances if needed
    PN_topo.calculate_halfwidths_from_geometric_conductances(ModelType::Poiseuille, 0.001);
    save_histogram(PN_topo.effective_halfwidths, "hist_halfwidth_topo.dat", output_dir);


    // =======================================================================
    // PHASE 5: NULL MODEL (RANDOM SPLITTING FRACTIONS)
    // =======================================================================
    std::cout << "\n[PHASE 5] Null Model: Random splitting fractions (independent U(0,1))..." << std::endl;

    // Reload clean network structure
    PoreNetwork PN_rand(link_file, coord_file, link_fmt, node_fmt, 2);
    PN_rand.viscosity = MU_DNS;

    // Randomize weights using
    randomize_network_weights(PN_rand);

    // Save the generated random distribution
    std::vector<double> rand_weights;
    for(const auto& l : PN_rand.links) rand_weights.push_back(l.weight);
    save_histogram(rand_weights, "hist_splitting_fractions_rand.dat", output_dir);

    // Solve inverse problem (same boundary conditions as phase 4)
    PN_rand.init_physics_from_splitting_fractions(avg_Q_target, P_in_avg, P_out_avg);

    // Save null model histograms
    save_histogram(PN_rand.geometric_conductances, "hist_conductance_rand.dat", output_dir);
    save_histogram(PN_rand.flows, "hist_flow_rate_rand.dat", output_dir);
    save_histogram(PN_rand.pressure_drops, "hist_pressure_drop_rand.dat", output_dir);
    PN_rand.calculate_halfwidths_from_geometric_conductances(ModelType::Poiseuille, 0.001);
    save_histogram(PN_rand.effective_halfwidths, "hist_halfwidth_rand.dat", output_dir);

    // =======================================================================
    // PHASE 6: SYNTHETIC MODEL (SAMPLED GAMMA HALFWIDTHS)
    // =======================================================================
    /*std::cout << "\n[PHASE 6] Synthetic Model: Sampling halfwidths from Gamma..." << std::endl;

    // Calcular los parámetros Gamma usando la estadística REAL (Fase 2)
    if(mean_eff <= 1e-12) 
    {
        std::cerr << "Error: mean_eff is too close to zero. Cannot generate Gamma." << std::endl;
        return 1;
    }
    
    // Varianza normalizada real (CV al cuadrado)
    double real_normalized_var = var_eff / (mean_eff * mean_eff); 
    std::cout << "  Using Real Normalized Variance (1/k) = " << real_normalized_var << std::endl;
    
    // Parámetros de la Gamma
    double k_gamma = 1.0 / real_normalized_var;
    double theta_gamma = mean_eff / k_gamma;

    std::cout << "  Derived Gamma Params: k=" << k_gamma << ", theta=" << theta_gamma << std::endl;

    // Load topology
    PoreNetwork PN_synth(link_file, coord_file, link_fmt, node_fmt, 2);

    // Sample new halfwidths
    RNG rng_synth(123456789);
    
    for(std::size_t i = 0; i < PN_synth.links.size(); ++i)
    {
        // Sample from Gamma(k, theta)
        double w = rng_synth.gamma(k_gamma, theta_gamma);
        
        // Safety: width must be > 0. If 0 is sampled (rare but possible), clamp it.
        if(w < 1e-8) w = 1e-8;

        PN_synth.effective_halfwidths[i] = w;
    }
    save_histogram(PN_synth.effective_halfwidths, "hist_halfwidth_synthetic.dat", output_dir);

    // Calculate conductance from these synthetic widths
    std::cout << "  Calculating conductances from synthetic widths..." << std::endl;
    PN_synth.set_geometric_conductances_from_geometry(ModelType::Poiseuille, 0.001);
    save_histogram(PN_synth.geometric_conductances, "hist_conductance_synthetic.dat", output_dir);

    // Solve direct problem (forward simulation)
    // We reuse the same target flow rate (avg_Q_target)
    std::cout << "  Solving flow field for synthetic network..." << std::endl;
    PN_synth.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);

    // Save resulting histograms
    save_histogram(PN_synth.flows, "hist_flow_rate_synthetic.dat", output_dir);
    save_histogram(PN_synth.pressure_drops, "hist_pressure_drop_synthetic.dat", output_dir);

    // Extract resulting splitting fractions (generated by the physics of the random widths)
    std::vector<double> synth_weights;
    for(const auto& l : PN_synth.links) synth_weights.push_back(l.weight);
    save_histogram(synth_weights, "hist_splitting_fractions_synthetic.dat", output_dir);*/

    // THE END

    // =======================================================================
    // PHASE 6: SYNTHETIC MODEL (EMPIRICAL RESAMPLING / SHUFFLING)
    // =======================================================================
    std::cout << "\n[PHASE 6] Synthetic Model: Shuffling real halfwidths..." << std::endl;

    // Load topology (this automatically loads the real halfwidths from the file)
    PoreNetwork PN_synth(link_file, coord_file, link_fmt, node_fmt, 2);
    PN_synth.viscosity = MU_DNS;


    // Resampling directo: desordenamos aleatoriamente la lista de anchuras
    std::cout << "  Shuffling effective halfwidths across the network..." << std::endl;
    std::mt19937 g_shuffle(123456789); // Fixed seed for reproducibility
    std::shuffle(PN_synth.effective_halfwidths.begin(), PN_synth.effective_halfwidths.end(), g_shuffle);

    // Guardamos el histograma de las anchuras barajadas
    // (El histograma en sí se verá idéntico al experimental, pero su topología ha cambiado)
    save_histogram(PN_synth.effective_halfwidths, "hist_halfwidth_synthetic.dat", output_dir);

    // Calculate conductance from these shuffled widths
    std::cout << "  Calculating conductances from shuffled widths..." << std::endl;
    PN_synth.set_geometric_conductances_from_geometry(ModelType::Poiseuille, 0.001);
    save_histogram(PN_synth.geometric_conductances, "hist_conductance_synthetic.dat", output_dir);

    // Solve direct problem (forward simulation)
    // We reuse the same target flow rate (avg_Q_target) calculated in Phase 3
    std::cout << "  Solving flow field for synthetic network..." << std::endl;
    PN_synth.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);

    // Save resulting histograms
    save_histogram(PN_synth.flows, "hist_flow_rate_synthetic.dat", output_dir);
    save_histogram(PN_synth.pressure_drops, "hist_pressure_drop_synthetic.dat", output_dir);

    // Extract resulting splitting fractions (generated by the physics of the shuffled widths)
    std::vector<double> synth_weights;
    for(const auto& l : PN_synth.links) synth_weights.push_back(l.weight);
    save_histogram(synth_weights, "hist_splitting_fractions_synthetic.dat", output_dir);

    // THE END
    std::cout << "\n[SUCCESS] All requested histograms (Phase 1-6) saved to " << output_dir << std::endl;
    return 0;

   
}