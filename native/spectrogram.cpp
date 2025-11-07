#include <napi.h>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

const double PI = 3.14159265358979323846;

// ===== FFT Implementation (unchanged but with inline optimization) =====

inline size_t reverseBits(size_t num, size_t bits) {
    size_t result = 0;
    for (size_t i = 0; i < bits; i++) {
        if (num & (1 << i)) {
            result |= 1 << (bits - 1 - i);
        }
    }
    return result;
}

void fft(std::vector<std::complex<double>>& x) {
    const size_t N = x.size();
    if (N <= 1) return;
    
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

// ===== Optimized Window Functions =====

// Pre-compute and cache window coefficients
class WindowCache {
private:
    std::vector<double> cosSquaredWindow;
    std::vector<double> deltaWindow;
    size_t cachedSize = 0;
    
public:
    void ensureSize(size_t N) {
        if (cachedSize == N) return;
        
        cosSquaredWindow.resize(N);
        deltaWindow.resize(N);
        
        double norm = 0.0;
        for (size_t i = 0; i < N; i++) {
            double x = (2.0 * i / (N - 1)) - 1.0;
            double cosVal = std::cos(x * PI / 2.0);
            cosSquaredWindow[i] = cosVal * cosVal;
            norm += cosSquaredWindow[i];
        }
        
        double normFactor = N / (norm * std::sqrt(2.0));
        for (size_t i = 0; i < N; i++) {
            cosSquaredWindow[i] *= normFactor;
        }
        
        // Compute delta window
        for (size_t i = 0; i < N; i++) {
            double x1 = (2.0 * (i - 0.5) / (N - 1)) - 1.0;
            double x2 = (2.0 * (i + 0.5) / (N - 1)) - 1.0;
            double cos1 = std::cos(x1 * PI / 2.0);
            double cos2 = std::cos(x2 * PI / 2.0);
            deltaWindow[i] = (cos1 * cos1 - cos2 * cos2) * normFactor;
        }
        
        cachedSize = N;
    }
    
    const std::vector<double>& getCosSquared() const { return cosSquaredWindow; }
    const std::vector<double>& getDelta() const { return deltaWindow; }
};

static WindowCache windowCache;

inline void applyWindowFast(std::vector<double>& data, const std::vector<double>& window) {
    const size_t N = data.size();
    for (size_t i = 0; i < N; i++) {
        data[i] *= window[i];
    }
}

// ===== Optimized Processing Functions =====

inline void applyPinkNoiseWeighting(std::vector<double>& magnitudes, double freqPerBin) {
    // Vectorized loop with early exit
    for (size_t i = 0; i < magnitudes.size(); i++) {
        double freq = i * freqPerBin;
        if (freq < 20.0) {
            magnitudes[i] = 0.0;
        } else {
            double octaves = std::log2(freq / 1000.0);
            magnitudes[i] *= std::pow(10.0, octaves * 0.3); // 6dB/20 = 0.3
        }
    }
}

inline void applyLowEndBoost(std::vector<double>& magnitudes, double freqPerBin, double boostAmount) {
    if (boostAmount <= 1.0) return;
    
    size_t bin20 = (size_t)std::ceil(20.0 / freqPerBin);
    size_t bin60 = (size_t)std::ceil(60.0 / freqPerBin);
    size_t bin250 = (size_t)std::ceil(250.0 / freqPerBin);
    size_t bin500 = (size_t)std::ceil(500.0 / freqPerBin);
    
    double boostFactor = boostAmount - 1.0;
    double bassBoost = 1.0 + boostFactor * 1.5;
    
    // Optimized loops
    const size_t maxIdx = magnitudes.size();
    
    for (size_t i = bin20; i < std::min(bin60, maxIdx); i++) {
        double t = (double)(i - bin20) / (bin60 - bin20);
        magnitudes[i] *= 1.0 + boostFactor * 1.3 * (1.0 - t * 0.3);
    }
    
    for (size_t i = bin60; i < std::min(bin250, maxIdx); i++) {
        magnitudes[i] *= bassBoost;
    }
    
    for (size_t i = bin250; i < std::min(bin500, maxIdx); i++) {
        double t = (double)(i - bin250) / (bin500 - bin250);
        magnitudes[i] *= 1.0 + boostFactor * 0.8 * (1.0 - t);
    }
}

inline void applyNoiseGate(std::vector<double>& magnitudes, double thresholdDB) {
    if (thresholdDB <= -120.0) return;
    
    double threshold = std::pow(10.0, thresholdDB / 20.0);
    for (size_t i = 0; i < magnitudes.size(); i++) {
        if (magnitudes[i] < threshold) {
            magnitudes[i] = 0.0;
        }
    }
}

// ===== Optimized Reassignment (Single Pass) =====

inline double logScale(double freq, double scale) {
    const double MIN_FREQ = 20.0;
    const double MAX_FREQ = 20000.0;
    
    freq = std::max(MIN_FREQ, std::min(freq, MAX_FREQ));
    
    if (scale < 0.5) {
        double linearAmount = 1.0 - (scale * 2.0);
        double logAmount = scale * 2.0;
        double linearX = (freq - MIN_FREQ) / (MAX_FREQ - MIN_FREQ);
        double logX = std::log2(freq / MIN_FREQ) / std::log2(MAX_FREQ / MIN_FREQ);
        return linearX * linearAmount + logX * logAmount;
    }
    return std::log2(freq / MIN_FREQ) / std::log2(MAX_FREQ / MIN_FREQ);
}

// ===== OPTIMIZED: Single-pass reassigned spectrum =====

Napi::Value ProcessReassignedSpectrumOptimized(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    
    try {
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
        size_t numBands = options.Get("numBands").As<Napi::Number>().Uint32Value();
        double frequencyScale = options.Get("frequencyScale").As<Napi::Number>().DoubleValue();
        
        // Ensure windows are cached
        windowCache.ensureSize(fftSize);
        
        // Process main window
        std::vector<double> mainSamples(fftSize, 0.0);
        size_t copyLen = std::min(fftSize, (size_t)audioData.ElementLength());
        
        // Combined copy and gain application
        for (size_t i = 0; i < copyLen; i++) {
            double sample = audioData[i];
            mainSamples[i] = std::isfinite(sample) ? sample * gain : 0.0;
        }
        
        applyWindowFast(mainSamples, windowCache.getCosSquared());
        
        std::vector<std::complex<double>> mainFFT(fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            mainFFT[i] = std::complex<double>(mainSamples[i], 0.0);
        }
        fft(mainFFT);
        
        // Process delta window (reuse mainSamples buffer)
        for (size_t i = 0; i < copyLen; i++) {
            double sample = audioData[i];
            mainSamples[i] = std::isfinite(sample) ? sample * gain : 0.0;
        }
        for (size_t i = copyLen; i < fftSize; i++) {
            mainSamples[i] = 0.0;
        }
        
        applyWindowFast(mainSamples, windowCache.getDelta());
        
        std::vector<std::complex<double>> auxFFT(fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            auxFFT[i] = std::complex<double>(mainSamples[i], 0.0);
        }
        fft(auxFFT);
        
        // Calculate magnitudes and reassigned frequencies in single pass
        std::vector<double> spectrum(numBands, 0.0);
        double freqPerBin = sampleRate / fftSize;
        double correction = sampleRate / (2.0 * PI);
        
        // Pre-calculate processing parameters
        double pinkFactor = naturalWeighting ? 1.0 : 0.0;
        double gateThreshold = std::pow(10.0, noiseGate / 20.0);
        
        for (size_t i = 0; i < fftSize; i++) {
            // Calculate magnitude
            double mag = std::abs(mainFFT[i]) / (fftSize / 2.0);
            
            // Apply noise gate early (skip expensive calculations)
            if (mag < gateThreshold) continue;
            
            // Apply natural weighting
            if (naturalWeighting) {
                double freq = i * freqPerBin;
                if (freq < 20.0) continue;
                double octaves = std::log2(freq / 1000.0);
                mag *= std::pow(10.0, octaves * 0.3);
            }
            
            // Calculate reassigned frequency
            double mainRe = mainFFT[i].real();
            double mainIm = mainFFT[i].imag();
            double auxRe = auxFFT[i].real();
            double auxIm = auxFFT[i].imag();
            
            double denom = mainRe * mainRe + mainIm * mainIm;
            double reassignedFreq = i * freqPerBin;
            
            if (denom > 1e-10) {
                double corr = -(mainRe * auxIm - auxRe * mainIm) / denom;
                reassignedFreq -= corr * correction;
            }
            
            // Bin into spectrum
            double x = logScale(reassignedFreq, frequencyScale);
            size_t binIdx = (size_t)(x * numBands);
            
            if (binIdx < numBands) {
                spectrum[binIdx] = std::max(spectrum[binIdx], mag);
            }
        }
        
        // Apply low-end boost (on final spectrum, much faster)
        applyLowEndBoost(spectrum, sampleRate / numBands, lowEndBoost);
        
        // Return result
        Napi::Float32Array result = Napi::Float32Array::New(env, numBands);
        for (size_t i = 0; i < numBands; i++) {
            result[i] = (float)spectrum[i];
        }
        
        return result;
        
    } catch (const std::exception& e) {
        Napi::Error::New(env, e.what()).ThrowAsJavaScriptException();
        return env.Null();
    }
}

// ===== OPTIMIZED: Regular spectrum with caching =====

Napi::Value ProcessSpectrumOptimized(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    
    try {
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
        
        windowCache.ensureSize(fftSize);
        
        // Combined copy, gain, and window
        std::vector<double> samples(fftSize, 0.0);
        size_t copyLen = std::min(fftSize, (size_t)audioData.ElementLength());
        
        const auto& window = windowCache.getCosSquared();
        for (size_t i = 0; i < copyLen; i++) {
            double sample = audioData[i];
            samples[i] = (std::isfinite(sample) ? sample * gain : 0.0) * window[i];
        }
        
        std::vector<std::complex<double>> fftData(fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            fftData[i] = std::complex<double>(samples[i], 0.0);
        }
        
        fft(fftData);
        
        // Calculate and process magnitudes in single pass
        std::vector<double> magnitudes(fftSize);
        double freqPerBin = sampleRate / fftSize;
        double gateThreshold = std::pow(10.0, noiseGate / 20.0);
        double normFactor = 2.0 / fftSize;
        
        for (size_t i = 0; i < fftSize; i++) {
            double mag = std::abs(fftData[i]) * normFactor;
            
            if (mag < gateThreshold) {
                magnitudes[i] = 0.0;
                continue;
            }
            
            if (naturalWeighting) {
                double freq = i * freqPerBin;
                if (freq < 20.0) {
                    magnitudes[i] = 0.0;
                    continue;
                }
                double octaves = std::log2(freq / 1000.0);
                mag *= std::pow(10.0, octaves * 0.3);
            }
            
            magnitudes[i] = mag;
        }
        
        applyLowEndBoost(magnitudes, freqPerBin, lowEndBoost);
        
        // Create result
        Napi::Object result = Napi::Object::New(env);
        Napi::Float32Array magArray = Napi::Float32Array::New(env, fftSize);
        for (size_t i = 0; i < fftSize; i++) {
            magArray[i] = (float)magnitudes[i];
        }
        result.Set("magnitudes", magArray);
        
        // Only return real/imag if needed (save memory)
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

// ===== OPTIMIZED: Smoothing =====

Napi::Value ApplySmoothing(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    
    try {
        Napi::Float32Array currentSpectrum = info[0].As<Napi::Float32Array>();
        Napi::Float32Array previousSpectrum = info[1].As<Napi::Float32Array>();
        double smoothFactor = info[2].As<Napi::Number>().DoubleValue();
        
        size_t len = currentSpectrum.ElementLength();
        Napi::Float32Array result = Napi::Float32Array::New(env, len);
        
        double invSmoothFactor = 1.0 - smoothFactor;
        
        // Vectorized smoothing
        for (size_t i = 0; i < len; i++) {
            float curr = currentSpectrum[i];
            float prev = previousSpectrum[i];
            float smoothed = curr * invSmoothFactor + prev * smoothFactor;
            result[i] = std::max(curr, smoothed * 0.95f);
        }
        
        return result;
        
    } catch (const std::exception& e) {
        Napi::Error::New(env, e.what()).ThrowAsJavaScriptException();
        return env.Null();
    }
}

// ===== Module Initialization =====

Napi::Object Init(Napi::Env env, Napi::Object exports) {
    exports.Set("processSpectrum", Napi::Function::New(env, ProcessSpectrumOptimized));
    exports.Set("processReassignedSpectrum", Napi::Function::New(env, ProcessReassignedSpectrumOptimized));
    exports.Set("applySmoothing", Napi::Function::New(env, ApplySmoothing));
    return exports;
}

NODE_API_MODULE(spectrogram_native, Init)