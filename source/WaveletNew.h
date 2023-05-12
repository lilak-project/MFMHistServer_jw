#ifndef WAVELETNEW_H
#define WAVELETNEW_H

#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <fstream>
#include <inttypes.h>
#include <complex>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <typeinfo>
#include <vector>

#define PI_1_4 1.3313353638

class WaveletNew {
public:
    WaveletNew(const std::vector<double>, int);
    WaveletNew(const std::vector<double>, int, bool);
    WaveletNew(const std::vector<double>, int, int, int);
    WaveletNew(const std::vector<double>, std::vector<double>);
    WaveletNew(const std::vector<double>, std::vector<double>, bool);
// Getters and Setters
public:
    std::vector<double> GetScale() {return scale;};
    std::vector<std::vector<double> > GetCWT() {return cwt;};
    std::vector<double> GetPa() {return pa;};

    void CalcCWT();
    void CalcCWTFast();

// Functions
private:
    std::vector<double> linspace(double a, double b, int n);

    // Wavelet
    void MakeNumWavelets();
    double Ricker(double f, double width, double length);
    std::pair<double, double> Morlet(int n, double width, double length);
    std::vector<double> fftw_convolve(const std::vector<double>, const std::vector<std::pair<double, double> >);
    std::vector<std::vector<double> > fftw_convolve_fast(const std::vector<double>, const std::vector<std::vector<std::pair<double, double> > >);

// Variables
private:
    std::vector<double> time;
    std::vector<double> volt;
    bool useMorlet;

    size_t size;

    // Wavelet
    int optimalSpace[18] = {1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048};
    int startNumber;
    int stopNumber;
    int numWavelets;
    std::vector<double> scale;
    std::vector<std::vector<double> > cwt;
    std::vector<double> pa;
};

#endif
