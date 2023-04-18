std::vector<Double_t> widths;widths.push_back(8); widths.push_back(16); widths.push_back(32);
     WaveletNew* wavelet = new WaveletNew(waveletY, widths, true);
     wavelet->CalcCWTFast();
     std::vector<Double_t> scale = wavelet->GetScale();
     std::vector<Double_t> pa = wavelet->GetPa();
     delete wavelet;
