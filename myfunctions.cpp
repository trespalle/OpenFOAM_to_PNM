#include "myfunctions.h"
#include <vector>
#include <stdexcept> 
#include <algorithm> 
#include <cstddef>   
#include <iostream>      
#include <fstream>       
#include <string>        
#include <numeric>       
#include <limits>      
#include <cmath>       

namespace myfun {

    std::vector<double> buildHistogram(std::vector<double>& input_data, std::size_t numBins, double& minValOut, double& maxValOut, double& delta, double percentile)
    {
        if(input_data.empty()) throw std::runtime_error("buildHistogram: Data vector is empty!");
        if(numBins == 0) throw std::runtime_error("buildHistogram: numBins cannot be 0!");

        // FILTRADO Y COPIA (Sanitize Data)
        // Ignoramos NaNs e Infs que rompen el cálculo del rango.
        std::vector<double> data;
        data.reserve(input_data.size());
        for(double v : input_data) {
            if(std::isfinite(v)) {
                data.push_back(v);
            }
        }

        if(data.empty()) {
            // Si todos eran NaN/Inf, devolvemos un histograma vacío o de ceros.
            minValOut = 0.0; maxValOut = 0.0; delta = 1.0;
            return std::vector<double>(numBins, 0.0);
        }

        // 2. ORDENAR
        std::sort(data.begin(), data.end());

        // 3. OBTENER RANGO (Percentiles)
        minValOut = data.front();
        //std::cout << "MIN = " << minValOut << std::endl;
        
        std::size_t p_index;
        if (percentile >= 1.0) p_index = data.size() - 1;
        else if (percentile <= 0.0) p_index = 0;
        else p_index = static_cast<std::size_t>(data.size() * percentile);

        // Safety clamp
        if (p_index >= data.size()) p_index = data.size() - 1;

        maxValOut = data[p_index];
        double range = maxValOut - minValOut;
        //std::cout << "max = " << maxValOut << std::endl;
        //std::cout << "range = " << range << std::endl;


        // GESTIÓN DE RANGO NULO O MUY PEQUEÑO (Singularidad)
        // Usamos una tolerancia basada en la magnitud de los datos
        double magnitude = std::max(std::abs(minValOut), std::abs(maxValOut));
        double tolerance = std::max(magnitude * 1e-14, 1e-100); // 1e-100 para evitar underflow absoluto

        // Si el rango es casi cero (ej: todos los datos son iguales), forzamos un bin ficticio.
        if(range <= tolerance) {
            delta = tolerance; // Evitar división por cero
            std::vector<double> histogram(1, 0.0);
            
            // Contamos cuántos datos caen en este "punto único"
            double count = 0.0;
            for(double v : data) {
                if(v <= maxValOut + tolerance && v >= minValOut - tolerance) count += 1.0;
            }
            
            // Densidad brutal (Dirac): Count / (N * dx)
            // Nota: Si dx es minúsculo, esto da un número gigante. Es matemáticamente correcto para una Delta.
            if(data.size() > 0) histogram[0] = count / (data.size() * delta);
            return histogram;

            //std::cout << "El rango es menor que la tolerancia. Al diablo " << std::endl;

        }

        // DETECCIÓN DE DATOS DISCRETOS ("Quantum")
        // Solo vale la pena buscar quantum si el rango es pequeño o los datos son enteros/similares.
        double quantum = std::numeric_limits<double>::max();
        // Heurística: revisamos una muestra o todo si es pequeño, para no perder tiempo O(N) si N es gigante.
        // Pero como ya ordenamos, es O(N) barato.
        for(std::size_t i = 1; i < data.size(); ++i) {
            double diff = data[i] - data[i-1];
            if(diff > tolerance && data[i] <= maxValOut) {
                if(diff < quantum) quantum = diff;
            }
        }

        //std::cout << "La diferencia minima entre dos valores es  = " << quantum << std::endl;


        // Decisión: ¿Es discreto o continuo?
        double standard_width = range / numBins;
        std::size_t effective_numbins = numBins;
        //std::cout << "La 'standard width' es = " << standard_width << std::endl;

        
        // Si el salto más pequeño entre datos es mayor que el ancho de bin estándar,
        // entonces el histograma "continuo" mentiría interpolando vacíos. Pasamos a modo discreto.
        if(quantum < std::numeric_limits<double>::max() && quantum > standard_width) {
            // Caso Discreto
            double bins_needed = std::round(range / quantum);
            // Limitamos para no crear millones de bins si hay un outlier raro
            if(bins_needed < 10000) { 
                effective_numbins = static_cast<std::size_t>(bins_needed) + 1;
                delta = quantum;
                //std::cout << "Deberiamos hacer el tratamiento cuantico y nos lo podemos permitir." << std::endl;
                //std::cout << "Usaremos este numero de bines: " << effective_numbins << std::endl;


            } else {
                delta = standard_width; // Fallback a continuo si pide demasiada memoria
                //std::cout << "Deberiamos hacer el tratamiento cuantico pero es muy caro." << std::endl;
            }
        } else {
            // Caso Continuo
            delta = standard_width;
            //std::cout << "Esto es continuo." << std::endl;

        }

        // LLENADO DEL HISTOGRAMA
        std::vector<double> histogram(effective_numbins, 0.0);
        double counted_elements = 0.0;

        for(double val : data) {
            // Importante: Si percentile < 1.0, ignoramos lo que esté fuera del rango superior
            if(val > maxValOut) break; 
            
            if(val >= minValOut) {
                // Cálculo del índice protegido contra errores de punto flotante
                double ratio = (val - minValOut) / delta;
                
                // floor es más seguro que cast directo para negativos (aunque aquí ordenamos)
                // Usamos un pequeño epsilon para corregir errores numéricos en el borde exacto
                std::size_t binIndex = static_cast<std::size_t>(std::floor(ratio + 1e-14));
                
                if(binIndex >= effective_numbins) binIndex = effective_numbins - 1;
                
                histogram[binIndex] += 1.0;
                counted_elements += 1.0;
            }
        }

        // NORMALIZACIÓN (Densidad de Probabilidad)
        // PDF: Sum(hist * delta) = 1.0
        // hist[i] = count / (N_total * delta)
        // Usamos data.size() como N_total para que la integral sobre TODO el dominio sea 1.
        // Si usáramos 'counted_elements', la integral sobre el rango recortado sería 1 (PDF condicional).
        // Generalmente en física queremos la PDF global.
        //std::cout << "delta = " << delta << std::endl;
        double norm_factor = data.size() * delta;
        if(norm_factor > 0) {
            for(auto& bin : histogram) bin /= norm_factor;
        }

        return histogram;
    }

    

    std::vector<double> twoDHistogram(
        const std::vector<double>& data1,
        const std::vector<double>& data2,
        const std::vector<double>& weights,
        std::size_t numBins1,
        std::size_t numBins2,
        double& min1, 
        double& max1,
        double& min2,
        double& max2 
    )
    {
        if(data1.empty()) throw std::runtime_error("twoDHistogram: El vector data1 está vacío!");
        if(data2.empty()) throw std::runtime_error("twoDHistogram: El vector data2 está vacío!");
        if(weights.empty()) throw std::runtime_error("twoDHistogram: El vector de pesos está vacío!");
        if(data1.size() != weights.size()) throw std::runtime_error("twoDHistogram: Los vectores data1 y weights no son compatibles!");

        if (numBins1 == 0 || numBins2 == 0) throw std::runtime_error("twoDHistogram: numBins no puede ser 0!");

        // Encontrar min y max para ambos conjuntos de datos
        min1 = *std::min_element(data1.begin(), data1.end());
        max1 = *std::max_element(data1.begin(), data1.end());
        min2 = *std::min_element(data2.begin(), data2.end());
        max2 = *std::max_element(data2.begin(), data2.end());

        std::vector<double> histogram(numBins1 * numBins2, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if((min1 == max1) || (min2 == max2))
        {
            histogram[0] = data1.size();
            return histogram;
        }

        // Llenar el histograma 2D
        for(std::size_t i = 0; i < data1.size(); ++i) 
        {
            double ratio1 = (data1[i] - min1)/(max1 - min1);
            std::size_t binIndex1 = static_cast<std::size_t>(ratio1 * numBins1);
            if(binIndex1 == numBins1) binIndex1--;

            double ratio2 = (data2[i] - min2)/(max2 - min2);
            std::size_t binIndex2 = static_cast<std::size_t>(ratio2 * numBins2);
            if(binIndex2 == numBins2) binIndex2--;

            histogram[binIndex1 * numBins2 + binIndex2] += weights[i];
        }

        // Normalización
        double delta1 = (max1 - min1)/numBins1;
        double delta2 = (max2 - min2)/numBins2;
        double suma_weights = 0.0;
        for(auto w : weights) suma_weights += w;
        for(auto &bin : histogram) bin /= (suma_weights * delta1 * delta2);

        return histogram;
    }

    void guardaHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double delta,
        const std::string& filename
    )
    {
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex\tbinCenter\tcount\n";

        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = minVal + delta * ((static_cast<double>(i) + 0.5));
            outFile << i << "\t" << binCenter << "\t" << histogram[i] << "\n";
        }

        outFile.close();
        std::cout << "Histograma guardado en " << filename << "\n";
    }

    void guardaWeightedHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    )
    {
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex\tbinCenter\tcount\n";

        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = minVal + (maxVal - minVal) * ((static_cast<double>(i) + 0.5)/numBins);
            outFile << i << "\t" << binCenter << "\t" << histogram[i] << "\n";
        }

        outFile.close();
        std::cout << "Histograma ponderado guardado en " << filename << "\n";
    }

    void med_var(
    const std::vector<double>& data,
    double& media,
    double& varianza,
    double percentile 
    )
    {
        // --- Initial sanity checks ---
        if (data.empty()) {
            throw std::invalid_argument("med_var: Data vector cannot be empty.");
        }
        if (percentile <= 0.0 || percentile > 100.0) {
            throw std::invalid_argument("med_var: Percentile must be between 0 (exclusive) and 100 (inclusive).");
        }

        // --- Case 1: Default behavior, use all data (most common and efficient) ---
        if (percentile == 100.0) {
            if (data.size() < 2) {
                throw std::invalid_argument("med_var: At least two elements are required to calculate variance.");
            }
            
            double media_acumulada = 0.0;
            double M2 = 0.0;
            std::size_t N = data.size();

            for(std::size_t i = 0; i < N; ++i) {
                double delta = data[i] - media_acumulada;
                media_acumulada += delta / static_cast<double>(i + 1);
                M2 += delta * (data[i] - media_acumulada);
            }

            media = media_acumulada;
            varianza = M2 / static_cast<double>(N - 1);
            return; // We are done
        }

        // --- Case 2: Filter data based on the specified percentile ---
        
        // 1. Find the value at the given percentile
        // We must make a copy because nth_element modifies the vector
        std::vector<double> data_copy = data;
        std::size_t percentile_index = static_cast<std::size_t>(data.size() * (percentile / 100.0));
        
        // Ensure the index is valid
        if (percentile_index >= data.size()) {
            percentile_index = data.size() - 1;
        }

        std::nth_element(data_copy.begin(), data_copy.begin() + percentile_index, data_copy.end());
        double percentile_value = data_copy[percentile_index];

        // 2. Apply Welford's algorithm only to elements below the percentile value
        double media_acumulada = 0.0;
        double M2 = 0.0;
        std::size_t count = 0; // Count of elements that meet the criteria

        for (const auto& value : data) {
            if (value < percentile_value) {
                count++;
                double delta = value - media_acumulada;
                media_acumulada += delta / static_cast<double>(count);
                M2 += delta * (value - media_acumulada);
            }
        }

        // 3. Finalize the calculation and check for errors on the filtered set
        if (count == 0) {
            throw std::runtime_error("med_var: No data found below the specified percentile. Mean and variance are undefined.");
        }
        if (count < 2) {
            throw std::runtime_error("med_var: Fewer than two data points found below the percentile. Variance is undefined.");
        }

        media = media_acumulada;
        varianza = M2 / static_cast<double>(count - 1);
    }

    void med_var_skew(
        const std::vector<double>& data,
        double& media,
        double& varianza,
        double& skewness,
        double percentile 
    )
    {
        // Data Filtering
        
        std::vector<double> filtered_data;

        if(percentile == 100.0) 
        {
            // If 100%, don't filter, just copy
            filtered_data = data;
        } 
        else 
        {
            // If there's a percentile, filter
            if(data.empty())
            {
                media = 0.0; 
                varianza = 0.0; 
                skewness = 0.0; 
                return;
            }

            
            if (percentile <= 0.0 || percentile > 100.0) throw std::invalid_argument("med_var_skew: Percentile must be (0, 100].");
            
            // Find the value at the percentile
            std::vector<double> data_copy = data;
            std::size_t percentile_index = static_cast<std::size_t>(data.size() * (percentile / 100.0));
            if(percentile_index >= data.size()) percentile_index = data.size() - 1;
            
            std::nth_element(data_copy.begin(), data_copy.begin() + percentile_index, data_copy.end());
            double percentile_value = data_copy[percentile_index];

            // Create the filtered data vector
            filtered_data.reserve(data.size());
            for(const auto& value : data) 
                if (value < percentile_value) filtered_data.push_back(value);
        
        }
        
        // Statistics Calculation (on filtered_data)

        std::size_t n = filtered_data.size();
        media = 0.0;
        varianza = 0.0;
        skewness = 0.0;

        if(n == 0) return; // No data, return 0

        // First pass: Calculate Mean
        for(const auto& value : filtered_data) media += value;
        media /= static_cast<double>(n);

        if(n < 2) return; // Variance and skewness are undefined, return 0

        // Second pass: Calculate Variance and 3rd Moment (m3) sums
        double m2_sum = 0.0; // Sum for variance
        double m3_sum = 0.0; // Sum for skewness
        for(const auto& value : filtered_data) 
        {
            double diff = value - media;
            m2_sum += diff*diff;
            m3_sum += diff*diff*diff;
        }
        
        // Unbiased sample variance
        varianza = m2_sum / (static_cast<double>(n) - 1.0);

        if(n < 3 || varianza == 0) return; // Skewness is undefined, return 0
        
        // Skewness
        double m3 = m3_sum / static_cast<double>(n); // 3rd central moment
        double s = std::sqrt(varianza); // Standard deviation
        
        // Apply bias correction (g1, the one SciPy uses)
        double n_dbl = static_cast<double>(n);
        double correction = std::sqrt(n_dbl * (n_dbl - 1.0)) / (n_dbl - 2.0);
        
        skewness = (m3 / (s*s*s)) * correction;
    }

    void guarda2DHistograma(
        const std::vector<double>& histogram,
        double min1,
        double max1,
        double min2,
        double max2,
        std::size_t numBins1,
        std::size_t numBins2,
        const std::string& filename
    )
    {
        //Chequiar el tamano
        if(histogram.size() != numBins1*numBins2)
        {
            std::cerr << "Error: histogram size does not match de number of bins";
            return;
        }
        
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex1\tbinCenter1\tbinIndex2\tbinCenter2\tcount\n";

        for(std::size_t i = 0; i < numBins1; ++i)
        {
            double binCenter1 = min1 + (max1 - min1) * ((static_cast<double>(i) + 0.5)/numBins1);
            for(std::size_t j = 0; j<numBins2; ++j)
            {
                double binCenter2 = min2 + (max2 - min2) * ((static_cast<double>(j) + 0.5)/numBins2);
                outFile << i << "\t" << binCenter1 << "\t" << j << "\t" << binCenter2 << "\t" << histogram[i*numBins2 + j] << "\n";
            }
            
        }

        outFile.close();
        std::cout << "Histograma 2D ponderado guardado en " << filename << "\n";
    }

    std::vector<double> weightedLogHistogram(
        const std::vector<double>& data,
        const std::vector<double>& weights,
        std::size_t numBins,
        double minValOut, 
        double maxValOut
    )
    {
        if(data.empty()) throw std::runtime_error("weightedHistogram: El vector de datos está vacío!");
        if(weights.empty()) throw std::runtime_error("weightedHistogram: El vector de pesos está vacío!");
        if(data.size() != weights.size()) throw std::runtime_error("weightedHistogram: Los vectores de datos y pesos deben tener el mismo tamaño!");

        if(numBins == 0) throw std::runtime_error("weightedHistogram: numBins no puede ser 0!");
        

        std::vector<double> histogram(numBins, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if(minValOut == maxValOut)
        {
            histogram[0] = data.size();
            return histogram;
        }

        if(minValOut<=0.0)
        {
            histogram[0] = data.size();
            return histogram;
        }

        double logf, f;

        logf = std::log(maxValOut/minValOut)/static_cast<double>(numBins);
        f = std::exp(logf);

        // Llenar el histograma logaritmico ponderado
        double suma_weights = 0.0;
        for(std::size_t i = 0; i < data.size(); ++i)
        {
            double ratio = std::log(data[i]/minValOut)/logf;
            if(ratio>=0.0)
            {
                std::size_t binIndex = static_cast<std::size_t>(ratio);
                if(binIndex >= numBins) binIndex = numBins-1;
                histogram[binIndex] += weights[i];
                suma_weights += weights[i];
            }  
        }

        // Normalización
        double delta_i, x_i;
        x_i = minValOut;
        for(std::size_t i=0; i<numBins; i++) 
        {
            delta_i = x_i*(f-1.0);
            histogram[i]/=(suma_weights*delta_i);
            x_i *= f;
        }
        return histogram;
    }

    void guardaWeightedLogHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    )
    {
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binCenter\tcount\n";

        double x_i = minVal;
        double logf, sqf;
        
        logf = std::log(maxVal/minVal)/static_cast<double>(numBins);
        sqf = std::exp(0.5*logf);
    


        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = x_i*sqf;
            outFile << binCenter << "\t" << histogram[i] << "\n";
            x_i = binCenter*sqf;
        }

        outFile.close();
        std::cout << "Histograma logaritmico ponderado guardado en " << filename << "\n";
    }

    std::vector<double> twoDLogHistogram(
        const std::vector<double>& data1,
        const std::vector<double>& data2,
        const std::vector<double>& weights,
        std::vector<double>& true_means,
        std::vector<double>& errors,
        std::size_t numBins1,
        std::size_t numBins2,
        double min1, 
        double max1,
        double min2,
        double max2 
    )
    {
        if(data1.empty()) throw std::runtime_error("twoDHistogram: El vector data1 está vacío!");
        if(data2.empty()) throw std::runtime_error("twoDHistogram: El vector data2 está vacío!");
        if(weights.empty()) throw std::runtime_error("twoDHistogram: El vector de pesos está vacío!");
        if(data1.size() != weights.size()) throw std::runtime_error("twoDHistogram: Los vectores data1 y weights no son compatibles!");
        if (data2.size() != data1.size())
            throw std::runtime_error("twoDHistogram: data2 and data1 must have the same size!");
        if (data2.size() != weights.size())
            throw std::runtime_error("twoDHistogram: data2 and weights must have the same size!");
        if (numBins1 == 0 || numBins2 == 0) throw std::runtime_error("twoDHistogram: numBins no puede ser 0!");
        std::cout << "numbins1 = " << numBins1 << ", " << "numbins2 = " << numBins2 << "\n";
        

        std::cout << "Las velocidades maxima y minima son " << max1 << " y " << min1 << "\n"; 
        std::cout << "Los shears maximo y minimo son " << max2 << " y " << min2 << "\n";
        std::vector<double> histogram(numBins1 * numBins2, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if((min1 == max1) || (min2 == max2))
        {
            histogram[0] = data1.size();
            return histogram;
        }

        double logf1, logf2;
        double eps = 1.0E-30;
        logf1 = std::log(max1/min1)/static_cast<double>(numBins1);
        logf2 = std::log(max2/min2)/static_cast<double>(numBins2);
        std::cout<<"logf1 = "<<logf1<<"logf2 ="<<logf2<<"\n";
       
        // Llenar el histograma 2D
        for(std::size_t i = 0; i < data1.size(); ++i) 
        {
            double ratio1 = std::log(std::fabs(data1[i]/(min1+eps)))/logf1;
            if(ratio1>=0.0) 
            {
                std::size_t binIndex1 = static_cast<std::size_t>(ratio1);
                if(binIndex1 >= numBins1) binIndex1  = numBins1-1;

                double ratio2 = std::log(std::fabs(data2[i]/(min2+eps)))/logf2;
                if(ratio2>=0.0) 
                {
                    std::size_t binIndex2 = static_cast<std::size_t>(ratio2);
                    if(binIndex2 >= numBins2) binIndex2 = numBins2-1;

                    histogram[binIndex1 * numBins2 + binIndex2] += weights[i];
                    true_means[binIndex1 * numBins2 + binIndex2] += weights[i]*data2[i];
                    errors[binIndex1 * numBins2 + binIndex2] += weights[i]*data2[i]*data2[i];
                }
                 
            }
            
        }

        
        return histogram;
    }

    void guarda2DLogHistograma(
        const std::vector<double>& histogram,
        double min1,
        double max1,
        double min2,
        double max2,
        std::size_t numBins1,
        std::size_t numBins2,
        const std::string& filename
    )
    {
        //Chequiar el tamano
        if(histogram.size() != numBins1*numBins2)
        {
            std::cerr << "Error: histogram size does not match de number of bins";
            return;
        }
        
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex1\tbinCenter1\tbinIndex2\tbinCenter2\tcount\n";

        double v_i, s_i;
        double logf1, sqf1, logf2, sqf2;
        logf1 = std::log(max1/min1)/static_cast<double>(numBins1);
        sqf1 = std::exp(0.5*logf1);
        logf2 = std::log(max2/min2)/static_cast<double>(numBins2);
        sqf2 = std::exp(0.5*logf2);

        v_i = min1;
        for(std::size_t i = 0; i < numBins1; ++i)
        {
            s_i = min2;
            double binCenter1 = v_i*sqf1;
            for(std::size_t j = 0; j<numBins2; ++j)
            {
                double binCenter2 = s_i*sqf2;
                outFile << i << "\t" << binCenter1 << "\t" << j << "\t" << binCenter2 << "\t" << histogram[i*numBins2 + j] << "\n";
                s_i = binCenter2*sqf2;
            }
            v_i = binCenter1*sqf1;  
        }

        outFile.close();
        std::cout << "Histograma 2D logaritmico ponderado guardado en " << filename << "\n";
    }





} // namespace myfun
