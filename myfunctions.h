// myfunctions.h
#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <vector>
#include <algorithm> // Para min_element y max_element
#include <stdexcept> // Para runtime_error
#include <cmath>
#include <cstddef> // Para size_t
#include <string>
#include <fstream> // Para ofstream
#include <iostream> // Para cerr
#include <iomanip> //para setprecision

namespace myfun {

    // Función para construir un histograma
    std::vector<double> buildHistogram(
        std::vector<double>& data,
        std::size_t numBins,
        double& minValOut, 
        double& maxValOut,
        double& delta,
        double percentile = 0.99
    );

    // Función para construir un histograma ponderado
    std::vector<double> weightedHistogram(
        const std::vector<double>& data,
        const std::vector<double>& weights,
        std::size_t numBins,
        double& minValOut, 
        double& maxValOut
    );

    // Función para construir un histograma 2D
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
    );

    // Función para guardar un histograma
    void guardaHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double delta,
        const std::string& filename
    );

    // Función para guardar un histograma ponderado
    void guardaWeightedHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    );

    // Función para calcular media y varianza utilizando el algoritmo de Welford
    void med_var(
        const std::vector<double>& data,
        double& media,
        double& varianza,
        double percentile = 100.0
    );

    void med_var_skew(
        const std::vector<double>& data,
        double& media,
        double& varianza,
        double& skewness,
        double percentile = 99.6 
    );

    //Funcion para guardar un histograma 2D ponderado
    void guarda2DHistograma(
        const std::vector<double>& histogram,
        double min1,
        double max1,
        double min2,
        double max2,
        std::size_t numBins1,
        std::size_t numBins2,
        const std::string& filename
    );

    std::vector<double> weightedLogHistogram(
        const std::vector<double>& data,
        const std::vector<double>& weights,
        std::size_t numBins,
        double minValOut, 
        double maxValOut
    );

    void guardaWeightedLogHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    );

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
    );

    void guarda2DLogHistograma(
        const std::vector<double>& histogram,
        double min1,
        double max1,
        double min2,
        double max2,
        std::size_t numBins1,
        std::size_t numBins2,
        const std::string& filename
    );
    


} // namespace myfun

#endif // MYFUNCTIONS_H
