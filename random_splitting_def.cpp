#include "myfunctions.h"
#include "Rand.h"
#include "graph_theory.h"

#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <functional>
#include <numeric>

struct Params
{
    std::string filename;
    int nreals = 1;
    int nbins = 50;
    std::string coords_file = "";
};

RNG Rand; // global RNG

Params read_params(int argc, char** argv);
void initialize_flows(Graph& G);
Graph create_big_tubes_graph(const Graph& G_original);
void create_output_directory(const std::string& dir_path);
void create_all_plots(const std::string& output_dir, const std::string& base_name_tubes, const std::string& base_name_pores, const std::string& analytical_base_name, double mu, double sigma2);
void save_inlet_outlet_coordinates(const Graph& G, const std::string& output_dir);
void save_node_layers(const Graph& G, const std::string& filename);


int main(int argc, char** argv){

    
    
    // Read parameters (number of histogram bins & number of realizations of the network:
    Params p = read_params(argc, argv);
    std::cout << "nbins = " << p.nbins << std::endl << "nreals = " << p.nreals << std::endl;

    // Create output directory:
    std::string base_folder = "output_def";

    const std::string output_dir = base_folder + "/";
    create_output_directory(output_dir);

    std::cout << "Running simulation with:" << std::endl;
    std::cout << "  Input File: " << p.filename << std::endl;
    std::cout << "  Output Dir: " << output_dir << std::endl;
    std::cout << "  Realizations: " << p.nreals << ", Bins: " << p.nbins << std::endl;

    // We read the graph from its list of links:
    Graph G(p.filename, p.coords_file);

    // IDENTIFY AGAIN INLETS AND OUTLETS (BASED ON GEOMETRY AND NOT TOPOLOGY)
    if (!p.coords_file.empty() && !G.nodes.empty()) 
    {
        double min_x = 1e9;
        double max_x = -1e9;
        
        // 1. Encontrar límites geométricos del dominio
        for (const auto& node : G.nodes) 
        {
            if (node.coordinates.empty()) continue;
            if (node.coordinates[0] < min_x) min_x = node.coordinates[0];
            if (node.coordinates[0] > max_x) max_x = node.coordinates[0];
        }

        double middle_point = (max_x + min_x) * 0.5; 

        std::vector<std::size_t> boundary_candidates;
        for (auto idx : G.inlet) boundary_candidates.push_back(idx);
        for (auto idx : G.outlet) boundary_candidates.push_back(idx);

        G.inlet.clear();
        G.outlet.clear();

        for (auto i : boundary_candidates) 
        {
            if (!G.nodes[i].coordinates.empty()) 
            {
                if (G.nodes[i].coordinates[0] <= middle_point) 
                {
                    G.inlet.push_back(i);
                    G.nodes[i].is_boundary = true; 
                }
                else 
                {
                    G.outlet.push_back(i);
                }
            }
        }
        std::cout << "-> Existing boundary nodes successfully re-sorted based on X coordinates." << std::endl;
    }

    // Check if the graph is a DAG
    bool is_dag = G.is_DAG();
    std::cout << "Graph is a DAG: " << (is_dag ? "true" : "false") << std::endl;

    

    // Save boundary coordinates if coordinates were loaded 
    if (!p.coords_file.empty() && !G.nodes.empty()) save_inlet_outlet_coordinates(G, output_dir);
    

    // We print some interesting information about the pore network
    std::size_t n_bulk = G.nodes.size() - G.inlet.size() - G.outlet.size();
    std::cout << "Number of inlets = " << G.inlet.size() << ", Number of outlets = " << G.outlet.size() << std::endl;
    std::cout << "Number of bulk nodes = " << n_bulk << std::endl;
    std::cout << "Among the bulk nodes there are: " << std::endl;
    

    std::size_t n_type_1, n_type_2;
    n_type_1 = n_type_2 = 0;
    for(auto& junction : G.nodes)
    {
        if((junction.outnbrs.size() == 2)&&(junction.innbrs.size() == 1)) 
        {
            ++ n_type_2;
            junction.state = 2;
        } 
        
        if((junction.outnbrs.size() == 1)&&(junction.innbrs.size() == 2))
        {
            ++ n_type_1;
            junction.state = 1;
        }
         
    }
    std::cout << n_type_1 << " type 1 junctions(merge) " << std::endl;
    std::cout << n_type_2 << " type 2 junctions(split) " << std::endl;
    std::cout << n_bulk - n_type_1 - n_type_2 << " others" << std::endl;

    // We compute the (k=1) degree distribution:
    std::vector<double> in_degrees, out_degrees;
    double min_in_degree, max_in_degree, min_out_degree, max_out_degree, delta;
    for(const auto& node : G.nodes)
    {
        in_degrees.push_back(static_cast<double>(node.innbrs.size()));
        out_degrees.push_back(static_cast<double>(node.outnbrs.size()));
    }
    
    auto degree_histogram = myfun::buildHistogram(in_degrees, p.nbins, min_in_degree, max_in_degree, delta, 1.0);
    myfun::guardaHistograma(degree_histogram, min_in_degree, delta, output_dir + "in_degree_dist.dat");
    std::cout << "Max in-degree = " << max_in_degree << std::endl;
    std::cout << "Min in-degree = " << min_in_degree << std::endl;

    degree_histogram = myfun::buildHistogram(out_degrees, p.nbins, min_out_degree, max_out_degree, delta, 1.0);
    myfun::guardaHistograma(degree_histogram, min_out_degree, delta, output_dir + "out_degree_dist.dat");
    std::cout << "Max out-degree = " << max_out_degree << std::endl;
    std::cout << "Min out-degree = " << min_out_degree << std::endl;

    G.compute_distances_from_inlet();
    save_node_layers(G, output_dir + "node_layers.dat");

    
    // We compute the stationary flows for p.nreals realizations of the disorder (weights):
    std::vector<double> mass_series, mass_1_series, mass_2_series;
    std::vector<double> evenness_series;
    std::vector<double> flow_series, big_in_flow_series, big_out_flow_series, small_in_flow_series, small_out_flow_series; 
    double minflow, maxflow;
    
    for(std::size_t e = 0; e < p.nreals; ++e)
    {
        
        initialize_flows(G); 
        G.compute_stationary_flows(Rand, false);
        
        //We collect flow-rate data at pores (junctions):
        for(auto& pore : G.nodes) if(!pore.is_boundary)
        {
            mass_series.push_back(pore.mass);
            if(pore.state == 1) mass_1_series.push_back(pore.mass);
            if(pore.state == 2) mass_2_series.push_back(pore.mass);
            if(pore.outnbrs.size() == 2) 
            {
                double ev = 2.0 * std::min(pore.outwgs[0], pore.outwgs[1]);
                evenness_series.push_back(ev);
            }
        }
        
        //We collect flow-rate data at tubes:
        for(auto& tube : G.links) 
        {
            if((!G.nodes[tube.in].is_boundary)&&(!G.nodes[tube.out].is_boundary))
            {
                double flow = G.nodes[tube.in].mass*tube.weight;
                flow_series.push_back(flow);
                if((G.nodes[tube.in].innbrs.size() == 2)) big_in_flow_series.push_back(flow);
                else small_in_flow_series.push_back(flow);
                if((G.nodes[tube.out].outnbrs.size() == 2)) big_out_flow_series.push_back(flow);
                else small_out_flow_series.push_back(flow);  
            }
                
        }         
    }

    // Calculate the mean and variance of the raw data
    double flow_mean, flow_variance;
    myfun::med_var(flow_series, flow_mean, flow_variance); 
    std::cout << "Calculated from data: mean = " << flow_mean << ", variance = " << flow_variance << std::endl;

    // Call Python script to generate the analytical curve
    std::cout << "Executing Python script to generate analytical fit..." << std::endl; 
    
    std::string  analytical_base_path = output_dir + "analytical";
    
    /* double mass_mean, mass_variance, mass_skewness;
    myfun::med_var(mass_series, mass_mean, mass_variance); 
    std::cout << "Calculated from data: mean = " << mass_mean << ", variance = " << mass_variance << std::endl;
    std::string command = "python calculate_fit.py " + std::to_string(mass_mean) + " " + std::to_string(mass_variance) + " " + analytical_base_path; 
     */

    std::string raw_data_path = output_dir + "raw_data_for_fit.tmp";
    std::ofstream raw_out(raw_data_path);
   
    for (double val : mass_series) raw_out << val << "\n";
    raw_out.close();
    std::cout << "Raw data saved to " << raw_data_path << std::endl;
    std::string command = "python calculate_fit.py " + raw_data_path + " " + analytical_base_path; 
    
    
    // Execute the command
    int result = system(command.c_str());
    if(result != 0) std::cerr << "Warning: Python script may have failed." << std::endl;

    
    std::cout << "\n--- Computing intra-layer statistics ---" << std::endl;

    // Group data by layer (topological distance to inlet)
    // layer_masses[d] = vector containing masses of all nodes at distance d
    // layer_evenness[d] = vector containing evenness of splitters at distance d
    std::map<int, std::vector<double>> layer_masses;
    std::map<int, std::vector<double>> layer_evenness;

    for(const auto& node : G.nodes)
    {
        int d = node.dist_to_inlet;
        if(d < 0) continue; // ignore unreachable nodes

        // Save mass (for all nodes in this layer)
        layer_masses[d].push_back(node.mass);

        // Save evenness (only if the node is a splitter)
        if(node.outnbrs.size() == 2) 
        {
            double ev = 2.0 * std::min(node.outwgs[0], node.outwgs[1]);
            layer_evenness[d].push_back(ev);
        }
    }

    // 2. Calculate "Correlation" (Homogeneity/CV)
    // NOTE: Calculating "correlation between nodes in the same layer" is tricky 
    // because there are no defined pairs. Instead, we measure HOMOGENEITY.
    // The standard measure is the Coefficient of Variation (CV = sigma / mu).
    // - If nodes behave identically (perfect correlation), CV = 0.
    // - If the layer is very heterogeneous (channeling), CV is high.
    
    std::vector<double> dist_vec, cv_mass_vec, cv_evenness_vec;
    std::vector<double> mean_mass_vec, mean_evenness_vec;

    for (auto const& [d, masses] : layer_masses)
    {
        // Mass Statistics
        double mu_m, var_m;
        if(masses.size() > 1) 
        {
            myfun::med_var(masses, mu_m, var_m); // Using your library
            // CV = std_dev / mean
            double cv_m = (mu_m > 1e-20) ? (std::sqrt(var_m) / mu_m) : 0.0;
            
            // Evenness Statistics (there might not be splitters in this layer)
            double mu_e = 0.0, var_e = 0.0, cv_e = 0.0;
            if(layer_evenness.count(d) && layer_evenness[d].size() > 1) 
            {
                myfun::med_var(layer_evenness[d], mu_e, var_e);
                cv_e = (mu_e > 1e-20) ? (std::sqrt(var_e) / mu_e) : 0.0;
            }

            dist_vec.push_back(static_cast<double>(d));
            cv_mass_vec.push_back(cv_m);
            cv_evenness_vec.push_back(cv_e);
            mean_mass_vec.push_back(mu_m);
            mean_evenness_vec.push_back(mu_e);
        }
    }

    // Save statistics to file
    std::string layer_stats_file = output_dir + "layer_statistics.dat";
    std::ofstream fout_layers(layer_stats_file);
    fout_layers << "# distance cv_mass cv_evenness mean_mass mean_evenness\n";
    for(std::size_t i = 0; i < dist_vec.size(); ++i) 
    {
        fout_layers << dist_vec[i] << "\t" 
                    << cv_mass_vec[i] << "\t" 
                    << cv_evenness_vec[i] << "\t"
                    << mean_mass_vec[i] << "\t"
                    << mean_evenness_vec[i] << "\n";
    }
    fout_layers.close();
    std::cout << "Layer statistics saved to " << layer_stats_file << std::endl;
   

    // Compute correspoinding histograms
    std::string tube_hist_filename = "tube_flow_dist.dat";
    auto histogram = myfun::buildHistogram(mass_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_flow_dist.dat");
    histogram = myfun::buildHistogram(mass_1_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_1_flow_dist.dat");
    histogram = myfun::buildHistogram(mass_2_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_2_flow_dist.dat");
    histogram = myfun::buildHistogram(flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + tube_hist_filename);
    histogram = myfun::buildHistogram(big_in_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "big_in_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(big_out_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "big_out_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(small_in_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "small_in_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(small_out_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "small_out_tube_flow_dist.dat");
    auto hist_even = myfun::buildHistogram(evenness_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(hist_even, minflow, delta, output_dir + "evenness_dist.dat");
    

   // MASS BALANCE
    double inlet_mass = 0.0, outlet_mass = 0.0;
    double reflected_mass = 0.0, trapped_mass = 0.0;

    for(auto& v : G.inlet) inlet_mass += G.nodes[v].mass;
    for(auto& v : G.outlet) outlet_mass += G.nodes[v].mass;

    // BACKFLOWS
    for(auto& v : G.inlet) 
    {
        for(std::size_t k = 0; k < G.nodes[v].innbrs.size(); ++k) 
        {
            std::size_t j = G.nodes[v].innbrs[k]; 
            double w_ji = G.nodes[v].inwgs[k];
            reflected_mass += G.nodes[j].mass * w_ji;
        }
        
        //false inlets
        if(G.nodes[v].outnbrs.empty()) 
        {
            reflected_mass += G.nodes[v].mass; 
        }
    }

    // trapped flow
    for(std::size_t i = 0; i < G.nodes.size(); ++i)
    {
        if(!G.nodes[i].is_boundary && G.nodes[i].outnbrs.empty())
        {
            trapped_mass += G.nodes[i].mass;
        }
    }

    std::cout << "\n--- TOPOLOGICAL MASS BALANCE ---" << std::endl;
    std::cout << "Mass Injected (inlets) = " << inlet_mass << std::endl;
    std::cout << "Transmitted (outlets)  = " << outlet_mass << std::endl;
    std::cout << "Reflected (backflow)   = " << reflected_mass << std::endl;
    std::cout << "Trapped (dead-ends)    = " << trapped_mass << std::endl;
    std::cout << "Total Check (Out+Ref+Trap) = " << outlet_mass + reflected_mass + trapped_mass << std::endl;
    std::cout << "--------------------------------\n" << std::endl;

    // Create the Gnuplot script for fitting
    myfun::med_var(flow_series, flow_mean, flow_variance, 99.5);
    create_all_plots(output_dir, tube_hist_filename, "pore_flow_dist.dat", "analytical", flow_mean, flow_variance);


    // We compute the fractions distribution
    std::vector<double> fractions;
    double fracmin, fracmax;
    for(auto& link : G.links ) fractions.push_back(link.weight);
    histogram = myfun::buildHistogram(fractions, p.nbins, fracmin, fracmax, delta);
    myfun::guardaHistograma(histogram, fracmin, delta, output_dir + "fractions_dist.dat");



    std::cout << "\nAll output files have been saved to the '" << output_dir << "' directory." << std::endl;

    return 0;

}



Params read_params(int argc, char** argv)
{
   Params p;

   if(argc < 2)
   {
    std::cerr << "Error: you must provide a filename as the first argument." << std::endl;
    std::cerr << "Usage: " << argv[0] << " <filename> [n_realizations] [n_bins] [coords_file]" << std::endl;
    std::cerr << "  <filename>       : Path to the links file." << std::endl;
    std::cerr << "  [n_realizations] : (Default 1) Number of realizations." << std::endl;
    std::cerr << "  [n_bins]         : (Default 50) Number of histogram bins." << std::endl;
    std::cerr << "  [coords_file]    : (Default \"\") Path to coordinates file. Use \"\" to skip." << std::endl;
    std::exit(1);
   }

   p.filename = argv[1];
   if(argc > 2) p.nreals = std::atoi(argv[2]);
   if(argc > 3) p.nbins = std::atoi(argv[3]);
   if(argc > 4) p.coords_file = argv[4];

   return p;
}

void initialize_flows( Graph& G)
{
    for(auto& v : G.nodes) v.mass = 0.0;
    double sum = 0.0;
    for(auto& i : G.inlet) 
    {
        double mass = Rand();
        G.nodes[i].mass = mass;
        sum += mass;
    }
    
    for(auto& i : G.inlet) G.nodes[i].mass /= sum;
}


void create_output_directory(const std::string& dir_path)
{
    try 
    {
        // Attempt to create the directory (and any parent directories needed)
        std::filesystem::create_directories(dir_path);
    } 
    catch (const std::exception& e) 
    {
        // If it fails, print an error message
        std::cerr << "Error creating directory '" << dir_path << "': " << e.what() << std::endl;
        // Exit the program, as the output directory is essential
        std::exit(1); 
    }
}

void create_all_plots(const std::string& output_dir, const std::string& base_name_tubes, const std::string& base_name_pores, const std::string& analytical_base_name, double mu, double sigma2) 
{
    // Derive the other histogram filenames automatically
    std::string big_in_tubes_hist_file = "big_in_" + base_name_tubes;
    std::string big_out_tubes_hist_file = "big_out_" + base_name_tubes;
    std::string small_in_tubes_hist_file = "small_in_" + base_name_tubes;
    std::string small_out_tubes_hist_file = "small_out_" + base_name_tubes;

    // SCRIPT 1: The main plot with P(q) and Alim fit 
    {
        std::string script_path = output_dir + "tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);
        
        gpout << "# --- Main Plot: Simulation vs. Full Analytical Fit and Gamma (Alim) Fit ---" << std::endl;
        gpout << "reset" << std::endl;
        
        // Alim Fit parameters (Standard Gamma from moments)
        gpout << "k_alim = (" << mu << "**2) / " << sigma2 << std::endl;
        gpout << "theta_alim = " << sigma2 << " / " << mu << std::endl;
        gpout << "GammaPDF(x) = (1.0/(gamma(k_alim)*theta_alim**k_alim)) * x**(k_alim-1) * exp(-x/theta_alim)" << std::endl;
        
        gpout << "set terminal wxt size 900,600 enhanced" << std::endl;
        gpout << "set title 'Full Distribution vs. Fits'" << std::endl;
        gpout << "set xlabel 'Flow Rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability Density'" << std::endl;
        
        gpout << "plot '" << base_name_tubes << "' u 2:3 with linespoints title 'All Tubes Data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_P_fit.dat' with lines lw 2 title 'P(q) Analytical Fit', \\" << std::endl;
        gpout << "     GammaPDF(x) with lines lw 2 dashtype 2 title 'Alim Fit (Gamma PDF)'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }

    // SCRIPT 2: Big tubes vs. f(q)
    {
        std::string script_path = output_dir + "big_tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);

        gpout << "# --- Big Tubes Plot: Simulation vs. f(q) component ---" << std::endl;
        gpout << "set xlabel 'Flow rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability density'" << std::endl;

        gpout << "plot '" << big_in_tubes_hist_file << "' u 2:3 with linespoints title 'Big-in tubes data', \\" << std::endl;
        gpout << "     '" << big_out_tubes_hist_file << "' u 2:3 with linespoints title 'Big-out tubes data', \\" << std::endl;
        gpout << "     '" << base_name_pores << "' u 2:3 with linespoints title 'Pore data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_g_fit.dat' with lines lw 2 title 'Model fit'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }

    // SCRIPT 3: Small tubes vs. g(q) 
    {
        std::string script_path = output_dir + "small_tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);

        gpout << "# Small Tubes Plot: Simulation vs. g(q) component" << std::endl;
        gpout << "set xlabel 'Flow rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability density'" << std::endl;

        gpout << "plot '" << small_in_tubes_hist_file << "' u 2:3 with linespoints title 'Small-in tubes data', \\" << std::endl;
        gpout << "     '" << small_out_tubes_hist_file << "' u 2:3 with linespoints title 'Small-out tubes data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_f_fit.dat' with lines lw 2 title 'Model prediction'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }
}

void save_inlet_outlet_coordinates(const Graph& G, const std::string& output_dir)
{
    // Save Inlet Coordinates
    std::string inlet_filename = output_dir + "inlet_coordinates.dat";
    std::ofstream fout_in(inlet_filename);
    if(!fout_in) std::cerr << "Error: Unable to open file: " << inlet_filename << std::endl;
    else 
    {
        for(std::size_t node_id : G.inlet) 
        {
            fout_in << node_id;
            for (double coord : G.nodes[node_id].coordinates) fout_in << " " << coord;
            fout_in << std::endl;
        }
        fout_in.close();
        std::cout << "Inlet coordinates saved to " << inlet_filename << std::endl;
    }

    // Save Outlet Coordinates
    std::string outlet_filename = output_dir + "outlet_coordinates.dat";
    std::ofstream fout_out(outlet_filename);
    if(!fout_out) std::cerr << "Error: Unable to open file: " << outlet_filename << std::endl;
    else 
    {
        for(std::size_t node_id : G.outlet) 
        {
            fout_out << node_id;
            for (double coord : G.nodes[node_id].coordinates) fout_out << " " << coord;
            fout_out << std::endl;
        }
        fout_out.close();
        std::cout << "Outlet coordinates saved to " << outlet_filename << std::endl;
    }
}

void save_node_layers(const Graph& G, const std::string& filename)
{
    std::ofstream fout(filename);
    if(!fout) 
    {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return;
    }

    // Header for Python (pandas)
    fout << "id x y distance\n";

    for(std::size_t i = 0; i < G.nodes.size(); ++i)
    {
        const auto& node = G.nodes[i];
        
        // Only save nodes that have valid coordinates
        // We also save nodes with distance -1 (unreachable) so Python can see them
        if (!node.coordinates.empty() && node.coordinates.size() >= 2)
        {
            fout << i << " " 
                 << node.coordinates[0] << " " 
                 << node.coordinates[1] << " " 
                 << node.dist_to_inlet << "\n";
        }
    }
    fout.close();
    std::cout << "Node layer data saved to " << filename << std::endl;
}

