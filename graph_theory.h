#ifndef GRAPH_THEORY_H
#define GRAPH_THEORY_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <map>           
#include <limits>        
#include "Rand.h"   
#include <Eigen/Sparse> 
#include <queue>
#include <Eigen/Dense>



struct Node
{
    double mass = 0.0;
    int state = 0;
    bool is_boundary = false;
    std::vector<double> coordinates;
    std::vector<std::size_t> outnbrs, innbrs, outlnks, inlnks; ///list of outnbrs and innbrs labels
    std::vector<double> outwgs, inwgs;
    int dist_to_inlet = -1;
    int dist_to_outlet = -1;
    
};

struct Link
{
    double weight;
    std::size_t in, out;
};

struct Graph
{
    Graph() = default; /// if the constructor I will define below is not used, then the default constructor is used 
    std::vector<Node> nodes;
    std::vector<Link> links;
    std::vector<std::size_t> inlet, outlet;

    // PBC requires knowledge of the spatial domain bounds
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

    Eigen::SparseMatrix<double, Eigen::RowMajor> weights_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> weights_matrix_powered;
    std::size_t current_k = 1;


private:
    std::size_t position(std::size_t guy, const std::vector<std::size_t>& people)
    {
        std::size_t pos = 0;
        for(pos = 0; pos < people.size(); ++pos) if(people[pos] == guy) return pos;
        std::cout << "The element was not found in the vector." << std::endl;
        std::exit(1);
    }

    // The following private method finds the min/max X and Y coordinates from all nodes. Necessary for PBC distance calculations.
     
    void update_domain_bounds()
    {
        if(nodes.empty()) return;

        bool first_node_found = false;
        // Find the first node with valid coordinates to initialize
        for(const auto& node : nodes) 
        {
            if(!node.coordinates.empty() && node.coordinates.size() >= 2) 
            {
                xmin = xmax = node.coordinates[0];
                ymin = ymax = node.coordinates[1];
                first_node_found = true;
                break;
            }
        }

        if(!first_node_found) 
        {
             std::cerr << "Warning: Cannot update domain bounds. No nodes with coordinates found." << std::endl;
             return;
        }

        // Iterate over the rest of the nodes
        for(const auto& node : nodes) 
        {
            if(!node.coordinates.empty() && node.coordinates.size() >= 2) 
            {
                xmin = std::min(xmin, node.coordinates[0]);
                xmax = std::max(xmax, node.coordinates[0]);
                ymin = std::min(ymin, node.coordinates[1]);
                ymax = std::max(ymax, node.coordinates[1]);
            }
        }
        // std::cout << "Domain bounds updated: X=[" << xmin << ", " << xmax << "], Y=[" << ymin << ", " << ymax << "]" << std::endl;
    }

    // The following private method calculates periodic distance (squared) between two nodes.
    // Assumes periodicity in X and hard walls in Y.
    
    double get_periodic_distance_sq(std::size_t node1_idx, std::size_t node2_idx) const
    {
        const auto& coords1 = nodes[node1_idx].coordinates;
        const auto& coords2 = nodes[node2_idx].coordinates;

        // If no coords, return a huge distance
        if(coords1.empty() || coords2.empty() || coords1.size() < 2 || coords2.size() < 2) 
            return std::numeric_limits<double>::max();

        double L_x = xmax - xmin;
        if(L_x <= 0) return std::numeric_limits<double>::max(); // Invalid bounds

        double dx = std::abs(coords1[0] - coords2[0]);
        double dx_pbc = std::min(dx, L_x - dx); // Periodic distance in X
        double dy = std::abs(coords1[1] - coords2[1]); // Normal distance in Y

        return (dx_pbc * dx_pbc) + (dy * dy);
    }

    // The following private method finds the N closest inlet nodes to a given outlet node.
    // Uses periodic distance.
    
    std::vector<std::size_t> find_closest_inlets(std::size_t outlet_node_idx, int n_closest) const
    {
        std::multimap<double, std::size_t> dist_map;
        for(std::size_t inlet_idx : inlet) 
        {
            double dist_sq = get_periodic_distance_sq(outlet_node_idx, inlet_idx);
            dist_map.insert({dist_sq, inlet_idx});
        }

        std::vector<std::size_t> closest;
        int count = 0;
        for(auto const& [dist, idx] : dist_map) 
        {
            closest.push_back(idx);
            count++;
            if (count == n_closest) break;
        }
        return closest;
    }

    bool is_DAG_util(std::size_t v, std::vector<int>& visited)
    {
        visited[v] = 1; // mark as "currently visiting" (on recursion stack)

        for(const auto& nbr_idx : nodes[v].outnbrs) 
        {
            if(visited[nbr_idx] == 0) // if not visited yet
            {
                if(is_DAG_util(nbr_idx, visited)) // recurse
                    return true; // cycle found downstream
            }
            else if (visited[nbr_idx] == 1) // if 'nbr' is already on our stack
            {
                // This is a "back edge"
                return true; // Cycle detected
            }
            // if visited[nbr_idx] == 2, it's fully processed, so we ignore it.
        }

        visited[v] = 2; // mark as "fully visited" (off recursion stack)
        return false; // no cycles found from this node
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix_power(std::size_t power) const
    {
        std::size_t N = nodes.size();
        // Ensure result and temp matrices are also RowMajor
        Eigen::SparseMatrix<double, Eigen::RowMajor> result(N, N);
        
        if(power == 0) 
        {
            result.setIdentity();
            return result;
        }

        result.setIdentity(); // Start with result = A^0 (Identity)
        
        // 'temp' will hold the powers of A: A^1, A^2, A^4, A^8, ...
        Eigen::SparseMatrix<double, Eigen::RowMajor> temp = weights_matrix; // Start with temp = A^1
        
        while(power > 0) 
        {
            // If the current bit of 'power' is 1 (i.e., power is odd)
            if(power % 2 == 1) result = (result * temp).pruned();
            
            // Square 'temp' for the next bit
            temp = (temp * temp).pruned();
            
            // Move to the next bit (integer division by 2)
            power /= 2;
        }
        return result;
    }

    // The following private method is a HELPER for compute_all_cumulative_degrees.
    // It runs BFS from a start node and returns a distance vector.
     
    Eigen::RowVectorXi run_bfs_distance_internal(std::size_t start_node, std::size_t N, bool use_innbrs)
    {
        Eigen::RowVectorXi distances(N);
        distances.setConstant(-1); // -1 = unvisited
        std::queue<std::size_t> q;

        q.push(start_node);
        distances(start_node) = 0;

        while(!q.empty())
        {
            std::size_t current_node_idx = q.front();
            q.pop();
            
            const auto& neighbors = use_innbrs ? nodes[current_node_idx].innbrs : nodes[current_node_idx].outnbrs;
            for(const auto& nbr_idx : neighbors)
            {
                if(distances(nbr_idx) == -1) // If unvisited
                {
                    distances(nbr_idx) = distances(current_node_idx) + 1;
                    q.push(nbr_idx);
                }
            }
        }
        return distances;
    }


     
    
public:
    void add_link(std::size_t i, std::size_t j, double weight=1.0)
    {
        if(i >= nodes.size()) nodes.resize(i+1);
        if(j >= nodes.size()) nodes.resize(j+1);

        nodes[i].outnbrs.push_back(j);
        nodes[i].outwgs.push_back(weight);
        nodes[i].outlnks.push_back(links.size());
        nodes[j].innbrs.push_back(i);
        nodes[j].inwgs.push_back(weight);
        nodes[j].inlnks.push_back(links.size());
        
        Link newlink;
        newlink.weight = weight;
        newlink.in = i;
        newlink.out = j;
        links.push_back(newlink);
    }

    void erase_link(std::size_t in, std::size_t out)
    {
        ///First, we delete "out" from the list of outnbrs of "in"
        std::size_t pos = position(out, nodes[in].outnbrs);
        nodes[in].outnbrs.erase(nodes[in].outnbrs.begin() + pos);
        nodes[in].outwgs.erase(nodes[in].outwgs.begin() + pos);
        std::size_t link_label = nodes[in].outlnks[pos];
        nodes[in].outlnks.erase(nodes[in].outlnks.begin() + pos);

        ///Then, we delete "in" from the list of innbrs of "out":
        pos = position(link_label, nodes[out].inlnks);
        nodes[out].innbrs.erase(nodes[out].innbrs.begin() + pos);
        nodes[out].inwgs.erase(nodes[out].inwgs.begin() + pos);
        nodes[out].inlnks.erase(nodes[out].inlnks.begin() + pos);

        ///Finally, we delete the link from the list of links, and we shift the labels of all the other links accordingly:
        links.erase(links.begin() + link_label);
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            for(std::size_t j = 0; j < nodes[i].outlnks.size(); ++j) 
                if(nodes[i].outlnks[j] > link_label) -- nodes[i].outlnks[j];
            for(std::size_t j = 0; j < nodes[i].inlnks.size(); ++j) 
                if(nodes[i].inlnks[j] > link_label) -- nodes[i].inlnks[j];
        }

    }

    void erase_link(std::size_t link_label)
    {
        ///First, we delete "out" from the list of outnbrs of "in"
        if(link_label >= links.size()) 
        {
            std::cout << "Error." << std::endl;
            std::exit(2);
        }
        std::size_t in, out, pos;
        in = links[link_label].in;
        out = links[link_label].out;
        pos = position(link_label, nodes[in].outlnks);
        nodes[in].outnbrs.erase(nodes[in].outnbrs.begin() + pos);
        nodes[in].outwgs.erase(nodes[in].outwgs.begin() + pos);
        nodes[in].outlnks.erase(nodes[in].outlnks.begin() + pos);

        ///Then, we delete "in" from the list of innbrs of "out":
        pos = position(link_label, nodes[out].inlnks);
        nodes[out].innbrs.erase(nodes[out].innbrs.begin() + pos);
        nodes[out].inwgs.erase(nodes[out].inwgs.begin() + pos);
        nodes[out].inlnks.erase(nodes[out].inlnks.begin() + pos);

        ///Finally, we delete the link from the list of links, and we shift the labels of all the other links accordingly:
        links.erase(links.begin() + link_label);
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            for(std::size_t j = 0; j < nodes[i].outlnks.size(); ++j) 
                if(nodes[i].outlnks[j] > link_label) -- nodes[i].outlnks[j];
            for(std::size_t j = 0; j < nodes[i].inlnks.size(); ++j) 
                if(nodes[i].inlnks[j] > link_label) -- nodes[i].inlnks[j];
        }

    }

    void erase_node(std::size_t label)
    {
        //First, we erase all the links in which it participates
        auto outlnks = nodes[label].outlnks;
        std::sort(outlnks.rbegin(), outlnks.rend()); //we order the vector in decreasing order, so that sequential elimination of links is not problematic
        for(std::size_t i = 0; i < outlnks.size(); ++i) erase_link(outlnks[i]);
        auto inlnks = nodes[label].inlnks;
        std::sort(inlnks.rbegin(), inlnks.rend());
        for(std::size_t i = 0; i < inlnks.size(); ++i) erase_link(inlnks[i]);

        // Then, we erase the node from the list of nodes
        nodes.erase(nodes.begin() + label);

        //Finally, we shift the labels of all the "bigger nodes" one unit down:
        for(std::size_t i = 0; i < nodes.size(); ++i)
        {
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j) 
                if(nodes[i].outnbrs[j] > label) -- nodes[i].outnbrs[j];

            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j) 
                if(nodes[i].innbrs[j] > label) -- nodes[i].innbrs[j];   
        }

        for(std::size_t i = 0; i < links.size(); ++i)
        {
            if(links[i].in > label) -- links[i].in;
            if(links[i].out > label) -- links[i].out;
            
        }
        
    }

    void check_and_renormalize_weights()
    {
        std::cout << "Starting K1L check and renormalization..." << std::endl;
        bool any_node_failed = false;
        std::size_t N = nodes.size();

        //Loop 1: check, renormalize outwgs, and update links ---
        for(std::size_t i = 0; i < N; ++i)
        {
            // We only care about nodes with outgoing links
            if(nodes[i].outwgs.empty()) continue; 

            // Calculate the sum of outgoing weights
            double sum = 0.0;
            for(double w : nodes[i].outwgs) sum += w;

            
            // Case 1: Sum is 0 (problematic)
            if (std::abs(sum) < 1e-10) 
            {
                if(!any_node_failed) 
                { 
                    // Print header only the first time
                    std::cerr << "--- K1L WARNINGS (Critical Errors) ---" << std::endl;
                    any_node_failed = true;
                }
                std::cerr << "  ERROR Node " << i << ": Has " << nodes[i].outwgs.size() << " outgoing links, but their weight sum is 0." << " Cannot renormalize! Flows will be incorrect." << std::endl;
                continue; // We cannot divide by zero, skip this node
            }
            
            // Case 2: Sum is not 1.0 (standard violation)
            if(std::abs(sum - 1.0) > 1e-9)
            {
                if(!any_node_failed) 
                { 
                    // Print header only the first time
                    std::cerr << "--- K1L WARNINGS (Corrections) ---" << std::endl;
                    any_node_failed = true;
                }
                std::cerr << "  Node " << i << ": Weight sum = " << sum << ". Forcing renormalization." << std::endl;
            }

            // Force renormalization (always, if sum != 0)
            // This corrects violations from Case 2 and also ensures
            // nodes that already summed to 1 are unaffected.
            for(std::size_t j = 0; j < nodes[i].outwgs.size(); ++j)
            {
                // 1. Renormalize the outwg of node i
                nodes[i].outwgs[j] /= sum;
                
                // 2. Update the weight in the global links list
                std::size_t link_idx = nodes[i].outlnks[j];
                links[link_idx].weight = nodes[i].outwgs[j];
            }
        }

        // Loop 2: update in-weights (copied from your logic in fix_disorder) 
        // This loop ensures that the 'inwgs' of neighbors reflect
        // the renormalized weights we just calculated.
        for (std::size_t i = 0; i < N; ++i) 
        {
            for (std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
            {
                auto outnbr  = nodes[i].outnbrs[j];
                auto outwg   = nodes[i].outwgs[j]; // The already renormalized weight
                auto& innbrs = nodes[outnbr].innbrs;
                auto& inwgs = nodes[outnbr].inwgs;
                
                // Find 'i' in the 'innbrs' list of 'outnbr'
                auto pos = std::find(innbrs.begin(), innbrs.end(), i);
                if (pos != innbrs.end()) 
                {
                    std::size_t k = std::size_t(pos - innbrs.begin());
                    if (k < inwgs.size()) inwgs[k] = outwg;           
                } 
                // No 'else' is needed; if not found, it's a graph
                // consistency error that 'add_link' should prevent.
            }
        }

        // Final Message 
        if (any_node_failed) 
        {
            std::cerr << "------------------------------------------" << std::endl;
            std::cout << "K1L check finished. Violations were found and corrected (or reported)." << std::endl;
        } 
        else std::cout << "K1L check passed. All nodes comply." << std::endl;
    
    }


    void update_boundaries()
    {
        // Reset the entire boundary state first
        inlet.clear();
        outlet.clear();
        for (auto& node : nodes) node.is_boundary = false;

        // Recalculate inlets and outlets in a single pass
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            // Check for inlet condition
            if (nodes[i].innbrs.empty()) 
            {
                inlet.push_back(i);
                nodes[i].is_boundary = true;
            }
            // Check for outlet condition
            if (nodes[i].outnbrs.empty()) 
            {
                outlet.push_back(i);
                nodes[i].is_boundary = true;
            }
        }
    }
    

    void read_coordinates(const std::string& filename)
    {
        std::ifstream fin(filename);
        if (!fin) throw std::runtime_error("Error: could not open coordinate file: " + filename);
    
        std::string line;
        std::size_t node_id;
        double coord_val;
        int line_num = 0;

        while (std::getline(fin, line)) 
        {
            line_num++;
            if(line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);

            // Read Node ID
            if (!(ss >> node_id)) continue;
            if (node_id >= nodes.size()) throw std::runtime_error("Error: Node ID " + std::to_string(node_id) + " from coordinate file is out of bounds (max is " + std::to_string(nodes.size() - 1) + ").");
        
            // Clear previous coordinates and read new ones
            nodes[node_id].coordinates.clear();
            while (ss >> coord_val) nodes[node_id].coordinates.push_back(coord_val);
            
            if (nodes[node_id].coordinates.empty()) throw std::runtime_error("Error: No coordinates found for node ID " + std::to_string(node_id) + " on line " + std::to_string(line_num) + " in " + filename);
            
        }
        fin.close();
        std::cout << "Successfully read coordinates for " << nodes.size() << " nodes from " << filename << std::endl;
        update_domain_bounds();
    }

    void check_mass_conservation()
    {
        std::cout << "Checking mass conservation for non-outlet nodes..." << std::endl;
        bool conservation_failed = false;

        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            // Solo nos importan los nodos que no son 'outlets'
            // (ya que los outlets por definición no tienen salida)
            if (nodes[i].outnbrs.empty()) continue; 

            double weight_sum = 0.0;
            for (double w : nodes[i].outwgs)
            {
                weight_sum += w;
            }

            // Comprobamos si la suma se desvía significativamente de 1.0
            if (std::abs(weight_sum - 1.0) > 1e-9)
            {
                std::cerr << "WARNING: Mass not conserved at node " << i << ". "
                        << "Sum of output weights is: " << weight_sum 
                        << " (Difference: " << (weight_sum - 1.0) << ")" << std::endl;
                conservation_failed = true;
                std::cout << "Outnbrs: " << std::endl;
                for( auto& v : nodes[i].outnbrs) std::cout << v << " ";
                std::cout << std::endl;

            }
        }

        if (!conservation_failed)
        {
            std::cout << "Mass conservation check passed for all nodes." << std::endl;
        }
    }

    void compute_distances_from_inlet()
    {
        std::cout << "Computing topological distances from inlets..." << std::endl;

        // Reset distances and prepare queue
        std::queue<std::size_t> q;
        
        for(auto& node : nodes) node.dist_to_inlet = -1; // reset to unvisited
        
        // Initialize BFS with all inlet nodes (distance 0)
        for(std::size_t inlet_idx : inlet) 
        {
            nodes[inlet_idx].dist_to_inlet = 0;
            q.push(inlet_idx);
        }

        // Run BFS
        while(!q.empty())
        {
            std::size_t u_idx = q.front();
            q.pop();

            int current_dist = nodes[u_idx].dist_to_inlet;

            // Explore neighbors (downstream)
            for (std::size_t v_idx : nodes[u_idx].outnbrs)
            {
                if (nodes[v_idx].dist_to_inlet == -1) // If unvisited
                {
                    nodes[v_idx].dist_to_inlet = current_dist + 1;
                    q.push(v_idx);
                }
            }
        }
        std::cout << "Distance computation complete." << std::endl;
    }

    void compute_distances_to_outlet()
    {
        std::cout << "Computing topological distances to outlet..." << std::endl;

        std::queue<std::size_t> q;
        
        // Reset and Initialize
        for(auto& node : nodes) node.dist_to_outlet = -1; 
        
        for(std::size_t outlet_idx : outlet) 
        {
            nodes[outlet_idx].dist_to_outlet = 0;
            q.push(outlet_idx);
        }

        // Run Backward BFS
        while(!q.empty())
        {
            std::size_t u = q.front();
            q.pop();
            int d = nodes[u].dist_to_outlet;

            // Traverse UPSTREAM (using in-neighbors)
            for(std::size_t v : nodes[u].innbrs)
            {
                if(nodes[v].dist_to_outlet == -1) // If unvisited
                {
                    nodes[v].dist_to_outlet = d + 1;
                    q.push(v);
                }
            }
        }
        std::cout << "Outlet distance computation complete." << std::endl;
    }



    void refresh_weights()
    {
        std::size_t N = nodes.size();
        if(N == 0) return;

        // Rebuild weights_matrix from updated links
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(links.size());
        for(const auto& link : links)
            triplets.push_back(Eigen::Triplet<double>(link.in, link.out, link.weight));
        
        weights_matrix.resize(N, N);
        
        weights_matrix.setZero(); // Clear old values
        weights_matrix.setFromTriplets(triplets.begin(), triplets.end());

        // Reset cache
        weights_matrix_powered = weights_matrix;
        current_k = 1;
        
        std::cout << "Graph weights refreshed from link list." << std::endl;
    }

    bool is_DAG()
    {
        std::size_t N = nodes.size();
        if(N == 0) return true; // an empty graph is a DAG

        // visited state: 0=unvisited, 1=visiting (on stack), 2=visited
        std::vector<int> visited(N, 0);

        for(std::size_t i = 0; i < N; ++i)
        {
            if(visited[i] == 0) // If we haven't processed this node yet
            {
                if (is_DAG_util(i, visited))
                {
                    std::cout << "Cycle detected starting from node " << i << std::endl;
                    return false; // Cycle found
                }
            }
        }
        return true; // no cycles found
    }

    // The following public method computes the full distance histograms for all nodes.
    // This is an expensive, one-time calculation.

    // max_k is the maximum distance to histogram.
    // out_dist_hist (Output) Vector [node_idx][dist] = count
    // in_dist_hist (Output) Vector [node_idx][dist] = count
     
    void compute_all_distance_histograms(std::size_t max_k, std::vector<std::vector<std::size_t>>& out_dist_hist, std::vector<std::vector<std::size_t>>& in_dist_hist)
    {
        std::size_t N = nodes.size();
        if(N == 0) return;

        std::cout << "Running one-time distance histogram calculation (N * BFS + O(N^2))..." << std::endl;

        // Resize output histograms
        // [node_idx][distance]
        out_dist_hist.assign(N, std::vector<std::size_t>(max_k + 1, 0));
        in_dist_hist.assign(N, std::vector<std::size_t>(max_k + 1, 0));

        // Run N * BFS and build the O(N^2) distance matrices locally
        Eigen::MatrixXi out_distance_matrix(N, N);
        Eigen::MatrixXi in_distance_matrix(N, N);

        for(std::size_t i = 0; i < N; ++i)
        {
            out_distance_matrix.row(i) = run_bfs_distance_internal(i, N, false);
            in_distance_matrix.row(i) = run_bfs_distance_internal(i, N, true);
        }

        // Scan the matrices ONCE (O(N^2)) to build histograms
        for(std::size_t i = 0; i < N; ++i)
        {
            for(std::size_t j = 0; j < N; ++j)
            {
                if(i == j) continue;

                // Out-degree: path from i to j
                int dist_out = out_distance_matrix(i, j);
                if(dist_out > 0 && dist_out <= max_k) out_dist_hist[i][dist_out]++;
                

                // In-degree: path from j to i
                int dist_in = in_distance_matrix(i, j); // This is dist(j -> i)
                if(dist_in > 0 && dist_in <= max_k) in_dist_hist[i][dist_in]++; // Histogram for node 'i'
            
            }
        }
        std::cout << "Distance histograms complete." << std::endl;
    }



    // The following method computes the stationary flows by solving the linear system directly.
    // Rand: a reference to the RNG object (needed for PBC logic).
    // PBC:  if true, solves the closed-loop stationary distribution.
    // If false (default), solves the open-loop flow propagation.
    void compute_stationary_flows(RNG& Rand, bool PBC = false)
    {
        std::size_t N = nodes.size();
        if (N == 0) return;

        // std::cout << "Building linear system for direct solver..." << std::endl;

        // The system is Ax = b
        // A is the matrix (I - T_internal - T_pbc)
        // x is the unknown mass vector (m)
        // b is the source vector (inlet masses, or 0 for PBC)
        
        Eigen::SparseMatrix<double> A(N, N);
        Eigen::VectorXd b(N);
        b.setZero(); // Initialize b vector to all zeros

        // We build the matrix A using a list of triplets for efficiency
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(N + links.size() + (PBC ? outlet.size() * 2 : 0));

        double total_inlet_mass = 0.0;

        // Build the (I - T_internal) part of the matrix 
        for(std::size_t i = 0; i < N; ++i) 
        {
            // Add the diagonal '1' (from the Identity matrix I)
            triplets.push_back(Eigen::Triplet<double>(i, i, 1.0));

            bool is_inlet_i = (std::find(inlet.begin(), inlet.end(), i) != inlet.end());
            if(is_inlet_i) 
            {
                // This is an INLET node.
                if(!PBC) b(i) = nodes[i].mass;
                total_inlet_mass += nodes[i].mass;
            }

            else 
            {
                // This is a BULK or OUTLET node.
                // Add the incoming flow contributions: -T
                // m_i = SUM(m_j * w_ji)  -->  m_i - SUM(m_j * w_ji) = 0
                for(std::size_t k = 0; k < nodes[i].innbrs.size(); ++k) 
                {
                    std::size_t j = nodes[i].innbrs[k]; // j is the source node
                    double w_ji = nodes[i].inwgs[k];    // w_ji is the weight from j to i
                    triplets.push_back(Eigen::Triplet<double>(i, j, -w_ji));
                }
            }
        }

        // Handle the PBC case (if applicable)
        if(PBC) 
        {
            // std::cout << "PBC enabled: Building feedback part of the matrix..." << std::endl;
            if(xmax - xmin <= 0) 
            {
                 std::cerr << "Fatal Error: Cannot run PBC without valid coordinates and domain bounds." << std::endl;
                 std::exit(1);
            }

            // Determine the fixed feedback rules ONCE
            struct PBCFeedbackRule { bool is_splitter = false; double fraction_to_first = 0.5; };
            std::map<size_t, PBCFeedbackRule> feedback_rules;
            for(std::size_t outlet_idx : outlet) 
            {
                PBCFeedbackRule rule;
                if(Rand() < 0.5) 
                { 
                    rule.is_splitter = true; 
                    rule.fraction_to_first = Rand(); 
                }
                feedback_rules[outlet_idx] = rule;
            }

            // Add the feedback links (-T_pbc) to the matrix
            // m_inlet = F(m_outlet)  -->  m_inlet - F(m_outlet) = 0
            for(size_t outlet_idx : outlet) 
            {
                const auto& rule = feedback_rules.at(outlet_idx);
                
                if(rule.is_splitter) 
                {
                    auto closest = find_closest_inlets(outlet_idx, 2);
                    if(closest.size() == 2) 
                    {
                        triplets.push_back(Eigen::Triplet<double>(closest[0], outlet_idx, -rule.fraction_to_first));
                        triplets.push_back(Eigen::Triplet<double>(closest[1], outlet_idx, -(1.0 - rule.fraction_to_first)));
                    } 
                    else if (closest.size() == 1) triplets.push_back(Eigen::Triplet<double>(closest[0], outlet_idx, -1.0));
            
                } 

                else 
                { 
                    // Not a splitter
                    auto closest = find_closest_inlets(outlet_idx, 1);
                    if (!closest.empty()) 
                        triplets.push_back(Eigen::Triplet<double>(closest[0], outlet_idx, -1.0));
                }
            }

            // Solve the stationary distribution (Homogeneous system) 
            // The system (I - T_int - T_pbc)m = 0 is singular (by design).
            // We must replace one equation with the conservation of mass: SUM(m_i) = 1.0 (or total_inlet_mass)
            
            std::vector<Eigen::Triplet<double>> final_triplets;

            if(inlet.empty()) 
            {
                std::cerr << "Fatal Error: PBC enabled but no inlet nodes found." << std::endl;
                std::exit(1);
            }
            std::size_t row_to_replace = inlet[0]; // Use the first inlet node
            // std::cout << "PBC: Replacing equation for node " << row_to_replace << " with conservation rule." << std::endl;

            // Copy all triplets *except* those for the row we're replacing
            for(const auto& t : triplets) if(t.row() != row_to_replace) final_triplets.push_back(t);
            
            // Add the NEW conservation equation: SUM(m_inlet) = total_inlet_mass
            for (std::size_t i = 0; i < N; ++i) 
            {
                bool is_inlet_i = (std::find(inlet.begin(), inlet.end(), i) != inlet.end());
                if (is_inlet_i) 
                    final_triplets.push_back(Eigen::Triplet<double>(row_to_replace, i, 1.0));
            }
    
            // Set the 'b' vector for this equation
            // b is already all zeros, so just set the one element
            b(row_to_replace) = total_inlet_mass;

            // Build the final matrix A
            A.setFromTriplets(final_triplets.begin(), final_triplets.end());

        } 
        
        else 
        {
            // Non-PBC case: just build the matrix from the first loop 
            A.setFromTriplets(triplets.begin(), triplets.end());
        }

        // Solve the system Ax = b 
        // std::cout << "Solving system with " << N << " nodes..." << std::endl;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        
        solver.compute(A);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Fatal Error: Eigen solver matrix decomposition failed." << std::endl;
            std::exit(1);
        }

        Eigen::VectorXd m = solver.solve(b);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Fatal Error: Eigen solver failed to solve the system." << std::endl;
            std::exit(1);
        }

        // std::cout << "System solved. Copying mass vector back to nodes." << std::endl;

        // Copy the solution vector 'm' back into the graph nodes
        for (std::size_t i = 0; i < N; ++i) nodes[i].mass = m(i);
    }

    void compute_weighted_power(std::size_t k)
    {
        std::size_t N = nodes.size();
        if(N == 0) return; // Empty graph

        // Case 0: Already computed
        if(k == current_k && weights_matrix_powered.size() > 0) 
        {
            std::cout << "Matrix A^" << k << " is already cached." << std::endl;
            return;
        }

        // Case 1: Compute k=0 (Identity)
        if(k == 0) 
        {
            std::cout << "Computing A^0 (Identity)..." << std::endl;
            weights_matrix_powered.resize(N, N);
            weights_matrix_powered.setIdentity();
            current_k = 0;
            return;
        }

        // Case 2: Compute k=1 (Base Adjacency Matrix)
        // This is the "bootstrapper" for all other calculations.
        if(k == 1) 
        {
            std::cout << "Computing A^1 (Base Adjacency Matrix)..." << std::endl;
            
            // Build A^1 (weights_matrix) if it's not already built
            if(weights_matrix.size() == 0) 
            {
                std::vector<Eigen::Triplet<double>> triplets;
                triplets.reserve(links.size());
                for(const auto& link : links) 
                {
                    triplets.push_back(Eigen::Triplet<double>(link.in, link.out, link.weight));
                }
                weights_matrix.resize(N, N);
                weights_matrix.setFromTriplets(triplets.begin(), triplets.end());
            }
            weights_matrix_powered = weights_matrix;
            current_k = 1;
            return;
        }
        
        // k>1:

        // Ensure A^1 (weights_matrix) is available to power up
        if(weights_matrix.size() == 0) 
        {
            std::cout << "Base matrix A^1 not found, computing it first..." << std::endl;
            compute_weighted_power(1); // This sets current_k=1 and weights_matrix_powered=A^1
        }
        
        // Case 3: k > current_k
        // We need to compute A^k. We have A^current_k.
        // We calculate A^(k - current_k) and multiply.
        if(k > current_k && current_k > 0)
        {
            std::size_t diff = k - current_k;
            std::cout << "Optimized path: Computing A^" << diff << " to multiply by cached A^" << current_k << "..." << std::endl;
            
            // Use the helper to get A^diff
            Eigen::SparseMatrix<double, Eigen::RowMajor> A_diff = matrix_power(diff);
            
            // A^k = A^current_k * A^diff
            weights_matrix_powered = (weights_matrix_powered * A_diff).pruned();
            current_k = k;
        }
        // Case 4: k < current_k OR we are starting from A^0 (Must recompute from scratch)
        // e.g., we have A^10 but want A^3. We can't divide.
        else 
        {
            std::cout << "Re-computing A^" << k << " from scratch..." << std::endl;
            
            // Use the helper to get A^k directly from A^1
            weights_matrix_powered = matrix_power(k);
            current_k = k;
        }
    }

    

    Graph(const std::string& file_name, const std::string& coords_filename = "") ///this is the constructor
    {
        //We read the list of links from a file:
        std::ifstream fin(file_name);
        if(!fin)
        {
            std::cerr << "Error: could not open file.";
            std::exit(5);
        }
        
        std::size_t i, j;
        double weight;

        std::string line;
        while(std::getline(fin, line)) 
        {
            if(line.empty() || line[0] == '#') continue;
            
            std::stringstream ss(line);
            // we try to read 3 values. If it doesnt work (because it is a header), next
            if(!(ss >> i >> j >> weight)) continue;

            add_link(i, j, weight);
        }
        
        fin.close();

        check_and_renormalize_weights(); // Enforce K1L

        // Read node coordinates (if provided)
        if (!coords_filename.empty()) 
        {
            try 
            {
                // This will also call update_domain_bounds() internally
                read_coordinates(coords_filename);
            } 
            catch (const std::runtime_error& e) 
            {
                 // Re-throw or handle more gracefully if needed
                 throw std::runtime_error("Error during coordinate reading: " + std::string(e.what()));
            }
        }

        update_boundaries();
        compute_weighted_power(1);
       
    }

    void print_links()
    {
        for(std::size_t i = 0; i < links.size(); ++i) 
            std::cout << links[i].in << " " << links[i].out << " " << links[i].weight << std::endl;
        std::cout << std::endl;
    }

    void save_links(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < links.size(); ++i) 
            fout << links[i].in << " " << links[i].out << " " << links[i].weight << std::endl;

    }

    void print_outnbrs()
    {
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            std::cout << i << " :";
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
                std::cout << " " << nodes[i].outnbrs[j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void save_outnbrs(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            fout << i << " ";
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
                fout << " " << nodes[i].outnbrs[j];
            fout << std::endl;
        }
    }

    void print_innbrs()
    {
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            std::cout << i << " :";
            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j)
                std::cout << " " << nodes[i].innbrs[j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void save_innbrs(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            fout << i << " ";
            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j)
                fout << " " << nodes[i].innbrs[j];
            fout << std::endl;
        }
    }

};

#endif