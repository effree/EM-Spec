#include <napi.h>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

const double PI = 3.14159265358979323846;

// Bit reversal for FFT
size_t reverseBits(size_t num, size_t bits) {
    size_t result = 0;
    for (size_t i = 0; i < bits; i++) {
        if (num & (1 << i)) {
            result |= 1 << (bits - 1 - i);
        }
    }
    return result;
}

// Cooley-Tukey FFT (in-place, power of 2)
void fft(std::vector<std::complex<double>>& x) {
    const size_t N = x.size();
    if (N <= 1) return;
    
    // Check if power of 2
    if ((N & (N - 1)) != 0) {
        throw std::runtime_error("FFT size must be power of 2");
    }
    
    // Bit reversal
    size_t bits = 0;
    size_t temp = N;
    while (temp > 1) {
        bits++;
        temp >>= 1;
    }
    
    for (size_t i = 0; i < N; i++) {
        size_t j = reverseBits(i, bits);
        if (j > i) {
            std::swap(x[i], x[j]);
        }
    }
    
    // Cooley-Tukey decimation-in-time
    for (size_t len = 2; len <= N; len *= 2) {
        double angle = -2.0 * PI / len;
        std::complex<double> wlen(std::cos(angle), std::sin(angle));
        
        for (size_t i = 0; i < N; i += len) {
            std::complex<double> w(1.0, 0.0);
            
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<double> u = x[i + j];
                std::complex<double> v = x[i + j + len / 2] * w;
                
                x[i + j] = u + v;
                x[i + j + len / 2] = u - v;
                
                w *= wlen;
            }
        }
    }
}

// Apply Hann-like window (cosine squared)
void applyWindow(std::vector<double>& data) {
    size_t N = data.size();
    double norm = 0.0;
    
    std::vector<double> window(N);
    for (size_t i = 0; i < N; i++) {
        double x = (2.0 * i / (N - 1)) - 1.0;
        double cosVal = std::cos(x * PI / 2.0);
        window[i] = cosVal * cosVal;
        norm += window[i];
    }
    
    // Apply window with normalization
    double normFactor = N / (norm * std::sqrt(2.0));
    for (size_t i = 0; i < N; i++) {
        data[i] *= window[i] * normFactor;
    }
}

// Apply delta window (for reassignment)
void applyDeltaWindow(std::vector<double>& data) {
    size_t N = data.size();
    double norm = 0.0;
    
    std::vector<double> window(N);
    for (size_t i = 0; i < N; i++) {
        window[i] = std::cos((2.0 * i / (N - 1) - 1.0) * PI / 2.0);
        window[i] = window[i] * window[i];
        norm += window[i];
    }
    
    // Apply delta window
    for (size_t i = 0; i < N; i++) {
        double x1 = (2.0 * (i - 0.5) / (N - 1)) - 1.0;
        double x2 = (2.0 * (i + 0.5) / (N - 1)) - 1.0;
        double cos1 = std::cos(x1 * PI / 2.0);
        double cos2 = std::cos(x2 * PI / 2.0);
        double w1 = cos1 * cos1;
        double w2 = cos2 * cos2;
        
        data[i] *= (w1 - w2);
    }
    
    // Normalize
    double normFactor = N / (norm * std::sqrt(2.0));
    for (size_t i = 0; i < N; i++) {
        data[i] *= normFactor;
    }
}

// Pink noise weighting
void applyPinkNoiseWeighting(std::vector<double>& magnitudes, double sampleRate, size_t fftSize) {
    double freqPerBin = sampleRate / fftSize;
    
    for (size_t i = 0; i < magnitudes.size(); i++) {
        double freq = i * freqPerBin;
        if (freq >= 20.0) {
            double octaves = std::log2(freq / 1000.0);
            double dbAdjustment = octaves * 6.0;
            magnitudes[i] *= std::pow(10.0, dbAdjustment / 20.0);
        } else {
            magnitudes[i] = 0.0;
        }
    }
}

// Low-end boost
void applyLowEndBoost(std::vector<double>& magnitudes, double sampleRate, size_t fftSize, double boostAmount) {
    if (boostAmount <= 1.0) return;
    
    double freqPerBin = sampleRate / fftSize;
    size_t bin20 = (size_t)std::ceil(20.0 / freqPerBin);
    size_t bin60 = (size_t)std::ceil(60.0 / freqPerBin);
    size_t bin250 = (size_t)std::ceil(250.0 / freqPerBin);
    size_t bin500 = (size_t)std::ceil(500.0 / freqPerBin);
    
    double boostFactor = boostAmount - 1.0;
    
    // Sub-bass: 20-60 Hz
    for (size_t i = bin20; i < bin60 && i < magnitudes.size(); i++) {
        double t = (double)(i - bin20) / (bin60 - bin20);
        magnitudes[i] *= 1.0 + boostFactor * 1.3 * (1.0 - t * 0.3);
    }
    
    // Bass: 60-250 Hz
    double bassBoost = 1.0 + boostFactor * 1.5;
    for (size_t i = bin60; i < bin250 && i < magnitudes.size(); i++) {
        magnitudes[i] *= bassBoost;
    }
    
    // Low mids: 250-500 Hz
    for (size_t i = bin250; i < bin500 && i < magnitudes.size(); i++) {
        double t = (double)(i - bin250) / (bin500 - bin250);
        magnitudes[i] *= 1.0 + boostFactor * 0.8 * (1.0 - t);
    }
}

// Noise gate
void applyNoiseGate(std::vector<double>& magnitudes, double thresholdDB) {
    if (thresholdDB <= -120.0) return;
    
    double threshold = std::pow(10.0, thresholdDB / 20.0);
    for (size_t i = 0; i < magnitudes.size(); i++) {
        if (magnitudes[i] < threshold) {
            magnitudes[i] = 0.0;
        }
    }
}

// Main spectrum processing function
Napi::Value ProcessSpectrum(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    
    try {
        // Parse arguments
        if (info.Length() < 2) {
            Napi::TypeError::New(env, "Expected 2 arguments").ThrowAsJavaScriptException();
            return env.Null();
        }
        
        Napi::Float32Array audioData = info[0].As<Napi::Float32Array>();
        Napi::Object options = info[1].As<Napi::Object>();
        
        double gain = options.Get("gain").As<Napi::Number>().DoubleValue();
        double sampleRate = options.Get("sampleRate").As<Napi::Number>().DoubleValue();
        size_t fftSize = options.Get("fftSize").As<Napi::Number>().Uint32Value();
        double lowEndBoost = options.Get("lowEndBoost").As<Napi::Number>().DoubleValue();
        double noiseGate = options.Get("noiseGate").As<Napi::Number>().DoubleValue();
        bool naturalWeighting = options.Get("naturalWeighting").As<Napi::Boolean>().Value();
        bool deltaWindow = options.Get("deltaWindow").As<Napi::Boolean>().Value();
        
        // Copy and apply gain
        std::vector<double> samples(fftSize, 0.0);
        size_t copyLen = std::min(fftSize, (size_t)audioData.ElementLength());
        for (size_t i = 0; i < copyLen; i++) {
            samples[i] = audioData[i] * gain;
        }
        
        // Apply window
        if (deltaWindow) {
            applyDeltaWindow(samples);
        } else {
            applyWindow(samples);
        }
        
        // Prepare FFT input
        std::vector<std::complex<double>> fftData(fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            fftData[i] = std::complex<double>(samples[i], 0.0);
        }
        
        // Compute FFT
        fft(fftData);
        
        // Calculate magnitudes
        std::vector<double> magnitudes(fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            magnitudes[i] = std::abs(fftData[i]) / (fftSize / 2.0);
        }
        
        // Apply natural weighting
        if (naturalWeighting) {
            applyPinkNoiseWeighting(magnitudes, sampleRate, fftSize);
        }
        
        // Apply noise gate
        applyNoiseGate(magnitudes, noiseGate);
        
        // Apply low-end boost
        applyLowEndBoost(magnitudes, sampleRate, fftSize, lowEndBoost);
        
        // Create result object with both magnitudes and complex data
        Napi::Object result = Napi::Object::New(env);
        
        // Magnitudes array
        Napi::Float32Array magArray = Napi::Float32Array::New(env, fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            magArray[i] = (float)magnitudes[i];
        }
        result.Set("magnitudes", magArray);
        
        // Real and imaginary parts (for reassignment)
        Napi::Float32Array realArray = Napi::Float32Array::New(env, fftSize);
        Napi::Float32Array imagArray = Napi::Float32Array::New(env, fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            realArray[i] = (float)fftData[i].real();
            imagArray[i] = (float)fftData[i].imag();
        }
        result.Set("real", realArray);
        result.Set("imag", imagArray);
        
        return result;
        
    } catch (const std::exception& e) {
        Napi::Error::New(env, e.what()).ThrowAsJavaScriptException();
        return env.Null();
    }
}

// Initialize module
Napi::Object Init(Napi::Env env, Napi::Object exports) {
    exports.Set("processSpectrum", Napi::Function::New(env, ProcessSpectrum));
    return exports;
}

NODE_API_MODULE(spectrogram_native, Init)