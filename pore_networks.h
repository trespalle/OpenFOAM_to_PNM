#ifndef PORE_NETWORK_H
#define PORE_NETWORK_H

#include "graph_theory.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <string>

// Schema for columns in the LINK file (after u, v)
enum class LinkData
{
    Weight,       // Graph weight
    HalfWidth,    // Geometric semi-aperture
    FlowRate,     // Experimental Flow Rate
    Ignored       // Skip column
};

// Schema for columns in the NODE file (after ID, coords...)
enum class NodeData
{
    Pressure,     // Node Pressure
    Ignored       // Skip column
};

enum class BoundaryCondition
{
    FixedFlowRate,    // user fixes Q (m3/s)
    FixedPressureDrop // user fixes Delta P (Pa)
};

enum class ModelType
{
    HeleShaw,
    Poiseuille
};

struct PoreNetwork : public Graph 
{
    // Physical properties of nodes
    // pressures[i] corresponds to node i
    std::vector<double> pressures; 

    // Physical properties of links
    // These vectors must have size equal to links.size()
    std::vector<double> geometric_conductances; 
    std::vector<double> flows;
    std::vector<double> effective_halfwidths; // Semi-aperture (a)
    std::vector<double> tube_lengths;
    std::vector<double> pressure_drops;

    // Global properties
    double viscosity = 0.001; // Pa*s
    double total_inlet_flow = 0.0;
    double equivalent_permeability = 0.0;
    double total_pressure_drop;

private:
    // Solves t - tanh(t) = y for t using Newton-Raphson
    // Robust against the singularity at t=0
    double solve_inverse_transcendental(double y_target)
    {
        // 1. Handle Sign (The function is odd: f(-t) = -f(t))
        // Work only with positive numbers to simplify logic
        double sign = (y_target < 0.0) ? -1.0 : 1.0;
        double y = std::abs(y_target);

        // 2. Small Value Approximation (Analytic)
        // Near 0, tanh(t) approx t - t^3/3
        // So: t - (t - t^3/3) = y  =>  t^3/3 = y  =>  t = cbrt(3*y)
        // This avoids division by zero in Newton-Raphson derivative
        if (y < 1e-4) 
        {
            return sign * std::cbrt(3.0 * y);
        }

        // 3. Initial Guess for Newton-Raphson
        // For large y, tanh(t) -> 1, so y approx t - 1 => t = y + 1
        double t = (y > 1.0) ? (y + 1.0) : std::cbrt(3.0 * y);

        // 4. Newton-Raphson Loop
        // f(t) = t - tanh(t) - y
        // f'(t) = 1 - sech^2(t) = tanh^2(t)
        const int max_iter = 20;
        const double tol = 1e-9;

        for (int k = 0; k < max_iter; ++k)
        {
            double th = std::tanh(t);
            double f_val = t - th - y;       // Function value
            double f_der = th * th;          // Derivative

            // Safety check for derivative
            if (f_der < 1e-10) break; 

            double dt = f_val / f_der;
            t -= dt;

            // Convergence check
            if (std::abs(dt) < tol) break;
        }

        return sign * t;
    }


public:
    // Helper: Updates the pressure_drops cache based on current node pressures
    // Should be called immediately after loading pressures or solving the system.
    void update_pressure_drops()
    {
        pressure_drops.resize(links.size());
        for(std::size_t i = 0; i < links.size(); ++i)
        {
            std::size_t u = links[i].in;
            std::size_t v = links[i].out;
            // Store norm of  drop consistent with link direction u->v
            pressure_drops[i] = pressures[u] - pressures[v];
        }
        std::cout << "Pressure drops (Delta P) cache updated." << std::endl;
    }

    // Parses extra node data (like Pressure) based on user schema
    // Assumes file format: ID  Coord_1 ... Coord_dim  [Column 0] [Column 1] ...
    void parse_custom_node_data(const std::string& filename, std::size_t dim, const std::vector<NodeData>& format)
    {
        std::ifstream fin(filename);
        if(!fin) { std::cerr << "Error: Cannot open node file: " << filename << std::endl; exit(1); }

        std::string line;
        std::size_t id;
        double val;
        double coord_skip;
        std::size_t line_count = 0;

        while(std::getline(fin, line))
        {
            line_count++;
            if(line.empty() || line[0] == '#') continue;
            
            std::stringstream ss(line);
            
            
            ss.imbue(std::locale::classic()); 
            
            // Try read ID. If fail (header), skip line.
            if(!(ss >> id)) continue;

            // Safety check
            if(id >= nodes.size()) 
            {
                // Si el ID es válido numéricamente pero está fuera de rango, avisamos
                // (Opcional: podríamos saltarlo silenciosamente si el archivo es de una red mayor)
                continue; 
            }

            // 2. Skip coordinates
            for(std::size_t k = 0; k < dim; ++k) 
            {
                if(!(ss >> coord_skip))
                {
                    std::cerr << "\n[PARSING ERROR] Failed reading coordinate." << std::endl;
                    std::cerr << "  File: " << filename << std::endl;
                    std::cerr << "  Line: " << line_count << " (Node ID " << id << ")" << std::endl;
                    std::cerr << "  Expected Coordinate index: " << k << std::endl;
                    std::cerr << "  Content: '" << line << "'" << std::endl;
                    exit(1);
                }
            }

            // Read Data Columns
            for(std::size_t k = 0; k < format.size(); ++k)
            {
                if(ss >> val) 
                {
                    switch(format[k]) {
                        case NodeData::Pressure: pressures[id] = val; break;
                        case NodeData::Ignored: break;
                    }
                }
                else
                {
                    // --- DEBUGGING BLOCK ---
                    std::cerr << "\n[PARSING ERROR] Missing node data column." << std::endl;
                    std::cerr << "  File: " << filename << std::endl;
                    std::cerr << "  Line: " << line_count << " (Node ID " << id << ")" << std::endl;
                    std::cerr << "  Expected Data Column: " << k << " (after coords)" << std::endl;
                    std::cerr << "  Content: '" << line << "'" << std::endl;
                    exit(1);
                    // -----------------------
                }
            }
        }
        fin.close();
        std::cout << "Custom node data parsed successfully." << std::endl;
        
        // Important: Update cache
        update_pressure_drops();
    }

    void parse_custom_link_data(const std::string& filename, const std::vector<LinkData>& format)
    {
        std::ifstream fin(filename);
        if(!fin) { std::cerr << "Error: Cannot open link file: " << filename << std::endl; exit(1); }

        std::string line;
        std::size_t link_idx = 0;
        std::size_t line_count = 0; 
        std::size_t u, v;
        double val;

        while(std::getline(fin, line))
        {
            line_count++;
            if(line.empty() || line[0] == '#') continue;
            
            std::stringstream ss(line);
            
            
            ss.imbue(std::locale::classic()); 
        

            if(!(ss >> u >> v)) continue; 

            if(link_idx >= links.size()) break; 

            // Check sync
            if(links[link_idx].in != u || links[link_idx].out != v) {
                 continue;
            }

            // Read Columns
            for(std::size_t k = 0; k < format.size(); ++k)
            {
                // Intentar leer
                if(!(ss >> val)) 
                {
                    // Mensaje de error detallado para que sepas qué pasa
                    std::cerr << "\n[PARSING ERROR] Fatal error reading link data." << std::endl;
                    std::cerr << "  File Line: " << line_count << std::endl;
                    std::cerr << "  Link Index: " << link_idx << std::endl;
                    std::cerr << "  Column Index: " << (k + 3) << std::endl;
                    std::cerr << "  Content: '" << line << "'" << std::endl;
                    exit(1);
                }

                switch(format[k]) {
                    //case LinkData::Weight: links[link_idx].weight = val; break;
                    case LinkData::HalfWidth: effective_halfwidths[link_idx] = val; break;
                    case LinkData::FlowRate: flows[link_idx] = val; break;
                    case LinkData::Ignored: break;
                }
            }
            link_idx++;
        }
        fin.close();
        std::cout << "Custom link data parsed. Processed " << link_idx << " links." << std::endl;
    }


    // Computes Euclidean length for all links once to avoid recalculation
    void compute_tube_lengths()
    {
        tube_lengths.resize(links.size());
        
        for(std::size_t i = 0; i < links.size(); ++i)
        {
            std::size_t u = links[i].in;
            std::size_t v = links[i].out;
            
            // Access coordinates safely
            if(nodes[u].coordinates.size() >= 2 && nodes[v].coordinates.size() >= 2)
            {
                double dx = nodes[u].coordinates[0] - nodes[v].coordinates[0];
                double dy = nodes[u].coordinates[1] - nodes[v].coordinates[1];
                double L = std::sqrt(dx*dx + dy*dy);
                
                // Physical protection against zero-length tubes (overlapping nodes)
                tube_lengths[i] = (L < 1e-13) ? 1e-13 : L;
            }
            else
            {
                // Fallback if coordinates are missing (should not happen based on Graph constructor)
                tube_lengths[i] = 1e-13; 
            }
        }
        std::cout << "Link lengths computed and cached." << std::endl;
    }
    
    // Master Constructor
    // link_schema: What are columns 3, 4... in link_file?
    // node_schema: What are columns (dim+2), (dim+3)... in coord_file?
    // dim: Number of coordinate dimensions (usually 2 or 3)
    PoreNetwork(const std::string& link_file, const std::string& coord_file, const std::vector<LinkData>& link_schema, const std::vector<NodeData>& node_schema, std::size_t dim = 2) 
        : Graph(link_file, coord_file) 
    {
        
        // --- IDENTIFY AGAIN INLETS AND OUTLETS (BASED ON GEOMETRY AND NOT TOPOLOGY) ---
        if(!nodes.empty() && !nodes[0].coordinates.empty()) 
        {
            double min_x = 1e9;
            double max_x = -1e9;
            
            for(const auto& n : nodes) {
                if(n.coordinates.empty()) continue;
                if(n.coordinates[0] < min_x) min_x = n.coordinates[0];
                if(n.coordinates[0] > max_x) max_x = n.coordinates[0];
            }
            
           
            double middle_point = (max_x + min_x) * 0.5;
            
            std::vector<std::size_t> boundary_candidates;
            for (auto idx : inlet) boundary_candidates.push_back(idx);
            for (auto idx : outlet) boundary_candidates.push_back(idx);
            
           
            inlet.clear(); 
            outlet.clear();
            
        
            for (auto i : boundary_candidates) 
            {
                if (!nodes[i].coordinates.empty()) 
                {
                    if (nodes[i].coordinates[0] <= middle_point) 
                    {
                        inlet.push_back(i);
                    }
                    else 
                    {
                        outlet.push_back(i);
                    }
                }
            }
            std::cout << "-> Existing boundary nodes successfully re-sorted based on X coordinates." << std::endl;
        }
        
        
        std::size_t N = nodes.size();
        std::size_t M = links.size();

        // Initialize vectors
        pressures.assign(N, 0.0);
        geometric_conductances.assign(M, 0.0); 
        flows.assign(M, 0.0);
        effective_halfwidths.assign(M, 0.0);
        pressure_drops.assign(M, 0.0);

        // Compute gometry (lengths)
        compute_tube_lengths();

        // Parse custom link data
        if(!link_schema.empty())
        {
            parse_custom_link_data(link_file, link_schema);
        }

        // Parse custom node data
        if(!node_schema.empty())
        {
            parse_custom_node_data(coord_file, dim, node_schema);
        }
    }

    // DIRECT PROBLEM METHODS (geometry to physics)

    // Calculates geometric conductances based on geometric width
    // Uses the 'effective_halfwidths' vector if filled, otherwise uses link weights as widths
    // H_meters is only used for HeleShaw
    void set_geometric_conductances_from_geometry(ModelType model, double H_meters = 0.001) 
    {
        double R_const = H_meters / std::sqrt(12.0);
        
        for(std::size_t i = 0; i < links.size(); ++i) 
        {
            
            double L = tube_lengths[i];
            double a = effective_halfwidths[i];

            if(a < 1e-12) geometric_conductances[i] = 0.0;
        
            else 
            {
                double S = 0.0;
                if(model == ModelType::HeleShaw)
                {
                    // Cubic law with tanh correction for rectangular cross-section
                    double ratio = R_const / a;
                    double tanh_term = std::tanh(a / R_const);
                    S = (H_meters * H_meters * H_meters * a * (1.0 - ratio * tanh_term)) / 6.0;
                }
                else if(model == ModelType::Poiseuille)
                {
                    // Planar Poiseuille law
                    S = 2.0 * H_meters * a * a * a / 3.0;
                }

                geometric_conductances[i] = S/L;
            }
        }
        std::cout << "Geometric conductances calculated from geometry." << std::endl;
    }

    // Inverse of the above: Calculates geometric semi-aperture (a) from permeability
    void calculate_halfwidths_from_geometric_conductances(ModelType model, double H_meters = 0.001)
    {
        double R_const = H_meters / std::sqrt(12.0);

        for(std::size_t i = 0; i < links.size(); ++i)
        {
            double S = geometric_conductances[i]*tube_lengths[i];
            if(S < 1e-25)
            {
                effective_halfwidths[i] = 0.0;
                continue;
            }

            if(model == ModelType::Poiseuille)
            {
                double term = 1.5 * S / H_meters;
                effective_halfwidths[i] = std::cbrt(term);
            }
            else if(model == ModelType::HeleShaw)
            {
                // INVERSE COMPLEX PLANAR (Newton-Raphson)
                // Y = (6 * S) / (H^3 * R) = t - tanh(t)
                // where t = a / R
                
                double numerator = 6.0 * S;
                double denominator = H_meters * H_meters * H_meters * R_const;
                double Y_target = numerator / denominator;

                // Call the numeric solver we built previously
                double t = solve_inverse_transcendental(Y_target);

                // Recover a
                effective_halfwidths[i] = t * R_const;
            }
        }
        std::cout << "Effective widths (semi-apertures) estimated from geometric conductances." << std::endl;
    }


    // SOLVER calculates pressures, flows and weights
    // 1. Solves linear pressure field
    // 2. Calculates flows
    // 3. Scales to match target total flow
    // 4. Updates topological weights in parent Graph
    void solve_flow_field(BoundaryCondition bc_type, double target_value = 1.0)
    {
        std::size_t N = nodes.size();
        Eigen::SparseMatrix<double> A(N, N);
        Eigen::VectorXd b(N);
        b.setZero();
        
        std::vector<Eigen::Triplet<double>> triplets;

        // First, we build the linear system for UNIT pressure drop (Pin=1, Pout=0) ---
        // We calculate the "conductance field". The absolute pressure values
        // don't matter yet, only the relative distribution.
        
        for(std::size_t i = 0; i < N; ++i) 
        {
            if(nodes[i].is_boundary) 
            {
                // Dirichlet condition
                triplets.push_back(Eigen::Triplet<double>(i, i, 1.0));
                
                // Arbitrary pressure drop P_in=1 P_out=0 (will be scaled later)
                bool is_inlet = (std::find(inlet.begin(), inlet.end(), i) != inlet.end());
                if(is_inlet) b(i) = 1.0;
                else b(i) = 0.0; 
            }
            else 
            {
                // Internal node mass balance
                double sum_conductance = 0.0;

                auto process_link = [&](std::size_t link_idx, std::size_t neighbor_idx)
                {
                    // geometric_conductance is g [m^3]
                    // Hydraulic Conductance = g / mu
                    double G_hyd = geometric_conductances[link_idx] / viscosity;
                    
                    triplets.push_back(Eigen::Triplet<double>(i, neighbor_idx, -G_hyd));
                    sum_conductance += G_hyd;
                };
                
                // Outgoing links
                for(std::size_t k = 0; k < nodes[i].outlnks.size(); ++k) 
                    process_link(nodes[i].outlnks[k], nodes[i].outnbrs[k]);

                // Incoming links
                for(std::size_t k = 0; k < nodes[i].inlnks.size(); ++k) 
                    process_link(nodes[i].inlnks[k], nodes[i].innbrs[k]);
                
                triplets.push_back(Eigen::Triplet<double>(i, i, sum_conductance));
            }
        }

        // Solve Pressure Field
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Error: Pressure solver decomposition failed." << std::endl; 
            exit(1);
        }
        Eigen::VectorXd x_base = solver.solve(b);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Error: Pressure solver failed." << std::endl; 
            exit(1);
        }

        // Calculate Base Flow (Q_base)
        // We compute the flow generated by DeltaP = 1.0 Pa
        double Q_base_inlet = 0.0;
        
        // Temporary storage for base flows to avoid recalculating
        std::vector<double> base_flows(links.size()); 

        for(std::size_t l = 0; l < links.size(); ++l) 
        {
            std::size_t u = links[l].in;
            std::size_t v = links[l].out;
            
            // Q = (g/mu) * dP
            double q = (geometric_conductances[l] / viscosity) * (x_base(u) - x_base(v));
            base_flows[l] = q;

            // Sum inlet flow
            bool is_inlet_u = (std::find(inlet.begin(), inlet.end(), u) != inlet.end());
            bool is_inlet_v = (std::find(inlet.begin(), inlet.end(), v) != inlet.end());
            
            if (is_inlet_u && !is_inlet_v) Q_base_inlet += q;        
            else if (is_inlet_v && !is_inlet_u) Q_base_inlet -= q;  
        }

        // Determine scaling factor 
        double scale_factor = 0.0;

        if (bc_type == BoundaryCondition::FixedFlowRate)
        {
            // We want Q_final = target_value
            // Q_final = Q_base * scale
            if(std::abs(Q_base_inlet) > 1e-20)
                scale_factor = target_value / Q_base_inlet;
            else
                std::cerr << "Warning: Network is impermeable (Q_base ~ 0)." << std::endl;
        }
        else if (bc_type == BoundaryCondition::FixedPressureDrop)
        {
            // We want DeltaP_final = target_value
            // Current DeltaP_base = 1.0 (Pin=1, Pout=0)
            // scale = target / 1.0
            scale_factor = target_value;
        }

        // Apply Scaling to Physics:
        // Scale Pressures
        // Note: Currently P_out is 0. If you wanted P_out != 0, you'd add an offset here.
        // But for DeltaP calculations, simple scaling works.
        for(std::size_t i = 0; i < N; ++i) pressures[i] = x_base(i) * scale_factor;

        // Scale Flows
        total_inlet_flow = 0.0;
        for(std::size_t l = 0; l < links.size(); ++l) 
        {
            flows[l] = base_flows[l] * scale_factor;
            
            std::size_t u = links[l].in;
            std::size_t v = links[l].out;
            bool is_inlet_u = (std::find(inlet.begin(), inlet.end(), u) != inlet.end());
            bool is_inlet_v = (std::find(inlet.begin(), inlet.end(), v) != inlet.end());
            
            if (is_inlet_u && !is_inlet_v) total_inlet_flow += flows[l];
            else if (is_inlet_v && !is_inlet_u) total_inlet_flow -= flows[l];
        }

        // 6. UPDATE WEIGHTS AND NODAL MASSES (OPTIMIZATION)
        // -------------------------------------------------
        // Reset masses first
        for(auto& node : nodes) node.mass = 0.0;

        for(std::size_t i = 0; i < N; ++i) 
        {
            double total_out_flow = 0.0;
            double total_in_flow = 0.0;

            // Sum outgoing flows
            for(std::size_t l_idx : nodes[i].outlnks) 
            {
                if(flows[l_idx] > 0) total_out_flow += flows[l_idx];
            }
            
            // Sum incoming flows
            for(std::size_t l_idx : nodes[i].inlnks)
            {
                if(flows[l_idx] > 0) total_in_flow += flows[l_idx];
            }

            // --- MASS CALCULATION ---
            bool is_inlet_i = (std::find(inlet.begin(), inlet.end(), i) != inlet.end());
            
            if (is_inlet_i) nodes[i].mass = total_out_flow; // Inlet real
            else            nodes[i].mass = total_in_flow;  // Bulk & Outlet
            
            // Normalize Mass (so sum(inlet_mass) = 1.0)
            if(total_inlet_flow > 1e-20) nodes[i].mass /= total_inlet_flow; 

            // Update Topological Weights (Splitting Fractions)
            for(std::size_t k = 0; k < nodes[i].outlnks.size(); ++k) 
            {
                std::size_t l_idx = nodes[i].outlnks[k];
                double w = 0.0;
                
                if(total_out_flow > 1e-20) 
                    w = std::max(0.0, flows[l_idx]) / total_out_flow;
                
                links[l_idx].weight = w;
                nodes[i].outwgs[k] = w;
            }
        }
        
        // Sync Eigen matrix for consistency
        refresh_weights();

        // Post-Processing 
        std::cout << "Flow field solved. Nodal masses and weights updated directly." << std::endl; 
        if(bc_type == BoundaryCondition::FixedFlowRate)
            std::cout << "  Condition: Fixed Flow = " << target_value << " m^3/s" << std::endl;
        else
            std::cout << "  Condition: Fixed DeltaP = " << target_value << " Pa" << std::endl;
            
        std::cout << "  Resulting Total Flow: " << total_inlet_flow << " m^3/s" << std::endl;
        // Important: Update caches
        update_pressure_drops();        // Update the DeltaP cache vector
    }
    

    void update_graph_weights_from_flows() 
    {
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            double total_out_flow = 0.0;
            for(std::size_t l_idx : nodes[i].outlnks) 
                if(flows[l_idx] > 0) total_out_flow += flows[l_idx];
            
            for(std::size_t k = 0; k < nodes[i].outlnks.size(); ++k) 
            {
                std::size_t l_idx = nodes[i].outlnks[k];
                double w = 0.0;
                
                if(total_out_flow > 1e-20) 
                    w = std::max(0.0, flows[l_idx]) / total_out_flow;
                
                links[l_idx].weight = w;
                nodes[i].outwgs[k] = w;
            }
        }
        
        // Important: update Eigen matrix for topological analysis
        refresh_weights(); 
        
        std::cout << "Graph weights updated from physical flows." << std::endl;
    }

    

    void compute_equivalent_permeability(double L_domain, double A_domain) 
    {
        // Get P_in from any inlet node
        double P_in = 0.0;
        if(!inlet.empty()) P_in = pressures[inlet[0]];

        // Assuming P_out scaled to 0 relative to P_in
        double delta_P = std::abs(P_in); 
        
        if(delta_P > 1e-20) 
        {
            equivalent_permeability = (total_inlet_flow * viscosity * L_domain) / (A_domain * delta_P);
            std::cout << "Equivalent Permeability (K_eq): " << equivalent_permeability << " m^2" << std::endl;
        } 
        else 
        {
            equivalent_permeability = 0.0;
        }
    }


    // INVERSE PROBLEM METHODS  (weights to physics)

    // Logic: We know Q (flow) and P (pressure) from files/experiments.
    // We want to find the conductance 'g' such that: Q = (g / mu) * DeltaP
    // Therefore: g_geom = (Q_exp * mu) / |P_u - P_v|
    void init_physics_from_experimental_data()
    {
        std::cout << "Calculating effective conductances from Q and Delta P..." << std::endl;

        // Ensure the pressure_drops cache is consistent with current node pressures
        update_pressure_drops();

        for(std::size_t i = 0; i < links.size(); ++i)
        {
            // We use absolute values because conductance is a positive scalar property
            // regardless of the flow direction.
            double dP = std::abs(pressure_drops[i]); 
            double q  = std::abs(flows[i]);          

            if(dP > 1e-20)
            {
                // Hydraulic Ohm's Law inverted:
                // g = (Q * mu) / dP
                geometric_conductances[i] = (q * viscosity) / dP;
            }
            else
            {
                // EDGE CASE: Zero Pressure Drop
                // If there is flow (Q > 0) but no pressure drop (dP ~ 0), 
                // it implies infinite conductance (super-conductive channel).
                // If there is no flow either, the tube is effectively closed (g = 0).
                if(q > 1e-20) 
                {
                    geometric_conductances[i] = 1e6; // Numeric representation of "Infinity"
                }
                else 
                {
                    geometric_conductances[i] = 0.0; // Stagnant / Clogged tube
                }
            }
        }
        
        std::cout << "Geometric conductances calibrated from experimental data." << std::endl;
        
        // OPTIONAL: Update the effective widths to match these new conductances.
        // This ensures geometric consistency (a <-> g) for visualization or further analysis.
        // We assume Planar Poiseuille (Slit) model as default for 2D.
        //calculate_widths_from_geometric_conductances(ModelType::Poiseuille);
    }

    // Calculates geometric conductances using splitting fractions (for Q) 
    // and harmonic/laplacian interpolation (for P).
    // This guarantees smooth pressure fields and avoids negative pressure drops.
    void init_physics_from_splitting_fractions(double Q_total, double P_in, double P_out) 
    {
        std::cout << "Initializing: Laplacian topological method..." << std::endl;
        
        // FIRST: CALCULATE FLOWS (Q) FROM WEIGHTS 
        // We propagate mass through the network using the original splitting fractions.

        if(inlet.empty()) 
        {
            std::cerr << "Error: No inlet nodes found for topological flow." << std::endl;
            return;
        }
        double mass_per_inlet = 1.0 / inlet.size();
        for(std::size_t id : inlet) nodes[id].mass = mass_per_inlet;
    
        
        // Use parent Graph mass solver:
        RNG temp_rng; 
        compute_stationary_flows(temp_rng, false); // Calculates 'nodes[i].mass'

        // Convert nodal mass to link flows
        // Q_link = Mass_node * Q_total * weight_link
        flows.assign(links.size(), 0.0);
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            double node_throughput = nodes[i].mass * Q_total;
            for(std::size_t k = 0; k < nodes[i].outlnks.size(); ++k) 
            {
                std::size_t link_idx = nodes[i].outlnks[k];
                double w = nodes[i].outwgs[k];
                flows[link_idx] = node_throughput * w;
            }
        }
        total_inlet_flow = Q_total;
        std::cout << "  -> Flows calculated from splitting fractions." << std::endl;


        // SECONDLY: CALCULATE PRESSURES (P) VIA LAPLACIAN 
        // We solve the linear system assuming homogeneous conductance (g=1).
        // This is equivalent to finding the "electric potential" in a resistor network
        // where all resistors are 1 Ohm. It creates a smooth, monotonic field.

        std::size_t N = nodes.size();
        Eigen::SparseMatrix<double> A(N, N);
        Eigen::VectorXd b(N);
        b.setZero();
        std::vector<Eigen::Triplet<double>> triplets;

        // Dummy conductance for the Laplacian solver
        double g_unity = 1.0; 

        for(std::size_t i = 0; i < N; ++i) 
        {
            if(nodes[i].is_boundary) 
            {
                // Dirichlet condition
                triplets.push_back(Eigen::Triplet<double>(i, i, 1.0));
                
                // Assign boundary pressure
                // We assume inlet nodes are at dist=0 (check 'dist_to_inlet' or 'inlet' vector)
                bool is_inlet = false;
                for(auto id : inlet) if(id == i) is_inlet = true;
                
                if(is_inlet) b(i) = P_in;
                else b(i) = P_out;
            }
            else 
            {
                // Internal nodes: Kirchhof law with g=1
                // sum(P_neighbor - P_i) = 0  =>  sum(P_neighbor) - degree * P_i = 0
                double sum_g = 0.0;

                auto add_connection = [&](std::size_t nbr_idx) 
                {
                    triplets.push_back(Eigen::Triplet<double>(i, nbr_idx, -g_unity));
                    sum_g += g_unity;
                };

                // Incoming neighbors
                for(std::size_t nbr : nodes[i].innbrs) add_connection(nbr);
                // Outgoing neighbors
                for(std::size_t nbr : nodes[i].outnbrs) add_connection(nbr);

                triplets.push_back(Eigen::Triplet<double>(i, i, sum_g));
            }
        }

        // Solve system
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Error: Laplacian decomposition failed." << std::endl; exit(1);
        }
        Eigen::VectorXd x_p = solver.solve(b);
        if(solver.info() != Eigen::Success) 
        {
            std::cerr << "Error: Laplacian solver failed." << std::endl; exit(1);
        }

        // Update pressures
        for(std::size_t i = 0; i < N; ++i) pressures[i] = x_p(i);
        
        // Update delta P cache
        update_pressure_drops();
        std::cout << "  -> Pressures interpolated using Laplacian (harmonic) field." << std::endl;


        // FINALLY: SOLVE INVERSE PROBLEM (FIND G)
        // g = Q_weights / dP_laplacian
        
        double min_g_cutoff = 1e-25; // Floor for conductance
        double max_g_clamp = 1e-3;   // Ceiling (safety)

        for(std::size_t l = 0; l < links.size(); ++l) 
        {
            double q = std::abs(flows[l]);
            double dP = std::abs(pressure_drops[l]); // Should be > 0 mostly, thanks to Laplacian

            if(dP > 1e-15) 
            {
                double g = (q * viscosity) / dP;
                
                // Safety clamps
                if(g > max_g_clamp) g = max_g_clamp;
                if(g < min_g_cutoff) g = 0.0;
                
                geometric_conductances[l] = g;
            }
            else 
            {
                // Zero pressure drop (very rare in Laplacian field unless disconnected)
                geometric_conductances[l] = (q > 1e-20) ? max_g_clamp : 0.0;
            }
        }
        
        std::cout << "Topological permeabilities calibrated." << std::endl;

        // Auto-calculate effective widths
        //calculate_halfwidths_from_geometric_conductances(ModelType::Poiseuille);
    }
};

#endif