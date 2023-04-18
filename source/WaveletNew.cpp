#include "WaveletNew.h"

WaveletNew::WaveletNew(const std::vector<double> voltArray, int numWaveletsInput) : startNumber(32), stopNumber(2048) {
    volt = voltArray;
    size = volt.size();
    numWavelets = numWaveletsInput;

    useMorlet = false;

    MakeNumWavelets();
}

WaveletNew::WaveletNew(const std::vector<double> voltArray, int numWaveletsInput, bool Morlet) : startNumber(32), stopNumber(2048) {
    volt = voltArray;
    size = volt.size();
    numWavelets = numWaveletsInput;

    useMorlet = true;

    MakeNumWavelets();
}

WaveletNew::WaveletNew(const std::vector<double> voltArray, int numWaveletsInput, int start, int stop) {
    volt = voltArray;
    size = volt.size();
    numWavelets = numWaveletsInput;
    startNumber = start;
    stopNumber = stop;

    useMorlet = false;

    MakeNumWavelets();
}

WaveletNew::WaveletNew(const std::vector<double> voltArray, std::vector<double> range) {
    volt = voltArray;
    size = volt.size();

    scale = range;
    numWavelets = scale.size();

    useMorlet = false;
}

WaveletNew::WaveletNew(const std::vector<double> voltArray, std::vector<double> range, bool Morlet) {
    volt = voltArray;
    size = volt.size();

    scale = range;
    numWavelets = scale.size();

    useMorlet = true;
}


void WaveletNew::MakeNumWavelets() {
    scale = linspace(startNumber, stopNumber, numWavelets);
}

double WaveletNew::Ricker(double f, double width, double length) {
  double x = f - length/2. + 1./2.;
  double func1 = 2./(sqrt(3.*width)*PI_1_4);
  double func2 = (1. - (x*x)/(width*width));
  double func3 = exp(-(x*x)/(2.*width*width));
  double total = func1 * func2 * func3;
  return total;
}

std::pair<double, double> WaveletNew::Morlet(int n, double width, double length) {
  double stepSize = 4.*M_PI/(length - 1.);
  double step = n*stepSize - 2.*M_PI;
  const std::complex<double> i(0., 1.);
  std::complex<double> func1(1./PI_1_4, 0.);
  std::complex<double> func2 = exp(i*width*step); // width > 5 doesn't need correction
  std::complex<double> func3 = exp(-0.5*step*step);
  std::complex<double> sum = func1*func2*func3;
  return std::make_pair<double, double>(real(sum), imag(sum));
}

std::vector<double> WaveletNew::linspace(double a, double b, int n) {
    std::vector<double> array;
    double step = static_cast<double>((b - a)/(n - 1));

    while(a <= b + 1) {
        array.push_back(a);
        a += step;
    }
    return array;
}

void WaveletNew::CalcCWT() {
    std::vector<std::vector<std::pair<double, double> > > waveletVector;

    for(size_t j = 0; j < numWavelets; j++) {
        size_t wavelet_size = std::min(10*scale[j], (double)volt.size());
        std::vector<std::pair<double, double> > wavelet;
        for(size_t i = 0; i < wavelet_size; i++) {
            if(!useMorlet) {
                std::pair<double, double> rickerWavelet = {Ricker(i, scale[j], wavelet_size), 0.};
                // rickerWavelet[0] = Ricker(i, scale[j], wavelet_size);
                // rickerWavelet[1] = 0.;
                wavelet.push_back(rickerWavelet);
            }
            else {
                std::pair<double, double> morletPair = Morlet(i, scale[j], wavelet_size);
                // fftw_complex morlet;
                // morlet[0] = morletPair.first;
                // morlet[1] = morletPair.second;
                wavelet.push_back(morletPair);
            }
        }
        waveletVector.push_back(wavelet);
        std::vector<double> cwt = fftw_convolve(volt, wavelet);

        double sum = 0.;
        for(uint i = 0; i < cwt.size(); i++) {
            sum += cwt[i]*cwt[i];
        }
        pa.push_back(sum/cwt.size());
    }
}

void WaveletNew::CalcCWTFast() {
    std::vector<std::vector<std::pair<double, double> > > waveletVector;

    for(size_t j = 0; j < numWavelets; j++) {
        size_t wavelet_size = volt.size();
        std::vector<std::pair<double, double> > wavelet;
        for(size_t i = 0; i < wavelet_size; i++) {
            if(!useMorlet) {
                std::pair<double, double> rickerWavelet = {Ricker(i, scale[j], wavelet_size), 0.};
                wavelet.push_back(rickerWavelet);
            }
            else {
                std::pair<double, double> morletWavelet = Morlet(i, scale[j], wavelet_size);
                wavelet.push_back(morletWavelet);
            }
        }
        waveletVector.push_back(wavelet);
    }

    std::vector<std::vector<double> > cwt = fftw_convolve_fast(volt, waveletVector);

    for(uint i = 0; i < cwt.size(); i++) {
        double sum = 0;
        for(uint j = 0; j < cwt[i].size(); j++) {
            sum += cwt[i][j]*cwt[i][j];
        }
        pa.push_back(sum/cwt[i].size());
    }
}

std::vector<double> WaveletNew::fftw_convolve(const std::vector<double> a, const std::vector<std::pair<double, double> > b) {
    size_t n_a = a.size();
    size_t n_b = b.size();
    int padded_length = n_a + n_b - 1;

    int optimal_length = 0;
    if(padded_length == 2047) {
        optimal_length = 2048;
    }
    else {
        for(uint i = 0; i < 17; i++) {
            if(padded_length > optimalSpace[i] && padded_length < optimalSpace[i + 1]) {
                optimal_length = optimalSpace[i + 1];
                break;
            }
        }
    }

    fftw_complex xData[optimal_length], yData[optimal_length];
    fftw_complex xWavelet[optimal_length], yWavelet[optimal_length];
    for(uint i = 0; i < optimal_length; i++) {
        xData[i][0] = 0.; xData[i][1] = 0.;
        xWavelet[i][0] = 0.; xWavelet[i][1] = 0.;
    }
    for(uint i = 0; i < n_a; i++) {
        xData[i][0] = a[i];
    }
    fftw_plan plan = fftw_plan_dft_1d(optimal_length, xData, yData, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    for(uint i = 0; i < n_b; i++) {
        xWavelet[i][0] = b[i].first;
        xWavelet[i][1] = b[i].second;
    }

    fftw_execute_dft(plan, xWavelet, yWavelet);
    fftw_destroy_plan(plan);

    fftw_complex xCombined[optimal_length], yCombined[optimal_length];
    for(uint i = 0; i < optimal_length; i++) {
        xCombined[i][0] = yData[i][0]*yWavelet[i][0] - yData[i][1]*yWavelet[i][1];
        xCombined[i][1] = yData[i][0]*yWavelet[i][1] + yData[i][1]*yWavelet[i][0];
    }
    fftw_plan plan_c = fftw_plan_dft_1d(optimal_length, xCombined, yCombined, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_c);
    fftw_destroy_plan(plan_c);

    fftw_complex yPadded[padded_length];
    for(uint i = 0; i < padded_length; i++) {
        yPadded[i][0] = yCombined[i][0]/optimal_length;
        yPadded[i][1] = yCombined[i][1]/optimal_length;
    }

    int startInd = (padded_length - n_a)/2.;

    std::vector<double> result;
    for(uint i = startInd; i < n_a + startInd; i++) {
        result.push_back(yPadded[i][0]);
    }

    return result;
}

// If all of 10*widths are over len(data), use this for speed up
std::vector<std::vector<double> > WaveletNew::fftw_convolve_fast(const std::vector<double> a, const std::vector<std::vector<std::pair<double, double> > > b) {
    size_t n_a = a.size();
    int padded_length = n_a + n_a - 1;

    int optimal_length = 2048;

    fftw_complex xData[optimal_length], yData[optimal_length];
    fftw_complex xWavelet[numWavelets][optimal_length], yWavelet[numWavelets][optimal_length];

    for(uint i = 0; i < optimal_length; i++) {
        xData[i][0] = 0.; xData[i][1] = 0.;
        for(uint j = 0; j < numWavelets; j++) {
            xWavelet[j][i][0] = 0.; xWavelet[j][i][1] = 0.;
        }
    }

    for(uint i = 0; i < n_a; i++) {
        xData[i][0] = a[i];
        for(uint j = 0; j < numWavelets; j++) {
            xWavelet[j][i][0] = b[j][i].first;
            xWavelet[j][i][1] = b[j][i].second;
        }
    }

    fftw_plan plan = fftw_plan_dft_1d(optimal_length, xData, yData, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    for(uint j = 0; j < numWavelets; j++) {
        fftw_execute_dft(plan, xWavelet[j], yWavelet[j]);
    }
    fftw_destroy_plan(plan);

    fftw_complex xCombined[numWavelets][optimal_length], yCombined[numWavelets][optimal_length];
    for(uint i = 0; i < optimal_length; i++) {
        for(uint j = 0; j < numWavelets; j++) {
            xCombined[j][i][0] = yData[i][0]*yWavelet[j][i][0] - yData[i][1]*yWavelet[j][i][1];
            xCombined[j][i][1] = yData[i][0]*yWavelet[j][i][1] + yData[i][1]*yWavelet[j][i][0];
        }
    }

    plan = fftw_plan_dft_1d(optimal_length, xCombined[0], yCombined[0], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    for(uint j = 1; j < numWavelets; j++) {
        fftw_execute_dft(plan, xCombined[j], yCombined[j]);
    }
    fftw_destroy_plan(plan);
    fftw_cleanup();

    fftw_complex yPadded[numWavelets][padded_length];
    for(uint i = 0; i < padded_length; i++) {
        for(uint j = 0; j < numWavelets; j++) {
            yPadded[j][i][0] = yCombined[j][i][0]/optimal_length;
            yPadded[j][i][1] = yCombined[j][i][1]/optimal_length;
        }
    }

    int startInd = (padded_length - n_a)/2.;

    std::vector<std::vector<double> > result;
    for(uint j = 0; j < numWavelets; j++) {
        std::vector<double> resultInd;
        for(uint i = startInd; i < n_a + startInd; i++) {
            resultInd.push_back(yPadded[j][i][0]);
        }
        result.push_back(resultInd);
    }

    return result;
}
