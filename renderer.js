const { ipcRenderer } = require('electron');
const { getLoopbackAudioMediaStream } = require('electron-audio-loopback');

// ===== Utility Functions =====
const map = (x, min, max, targetMin, targetMax) => 
  (x - min) / (max - min) * (targetMax - targetMin) + targetMin;

const clamp = (x, min, max) => Math.min(Math.max(x, min), max);

const idxWrapOver = (x, length) => (x % length + length) % length;

const hertzToFFTBin = (x, y = 'round', bufferSize = 4096, sampleRate = 44100) => {
  const bin = x * bufferSize / sampleRate;
  const func = ['floor','ceil','trunc'].includes(y) ? y : 'round';
  return Math[func](bin);
};

const fftBinToHertz = (x, bufferSize = 4096, sampleRate = 44100) => 
  x * sampleRate / bufferSize;

// Consolidated ascale with proper dbRange
const ascale = (x, dbRange = 90) => {
  const db = 20 * Math.log10(Math.max(x, 1e-10));
  return clamp(map(db, -dbRange, 0, 0, 1), 0, 1);
};

// ===== Audio Processing Functions =====

// Improved low-end boost - optimized with pre-calculated frequency bins
function applyLowEndBoost(complexFFT, sampleRate, fftSize, boostAmount) {
  const freqPerBin = sampleRate / fftSize;
  
  // Pre-calculate bin ranges (faster than calculating frequency per iteration)
  const bin20 = Math.ceil(20 / freqPerBin);
  const bin60 = Math.ceil(60 / freqPerBin);
  const bin250 = Math.ceil(250 / freqPerBin);
  const bin500 = Math.ceil(500 / freqPerBin);
  
  const boostFactor = boostAmount - 1.0;
  
  // Sub-bass: 20-60 Hz
  for (let i = bin20; i < bin60; i++) {
    const t = (i - bin20) / (bin60 - bin20);
    complexFFT[i].magnitude *= 1.0 + boostFactor * 1.3 * (1 - t * 0.3);
  }
  
  // Bass: 60-250 Hz
  const bassBoost = 1.0 + boostFactor * 1.5;
  for (let i = bin60; i < bin250; i++) {
    complexFFT[i].magnitude *= bassBoost;
  }
  
  // Low mids: 250-500 Hz
  for (let i = bin250; i < bin500; i++) {
    const t = (i - bin250) / (bin500 - bin250);
    complexFFT[i].magnitude *= 1.0 + boostFactor * 0.8 * (1 - t);
  }
  
  return complexFFT;
}

// Noise gate - optimized
function applyNoiseGate(complexFFT, thresholdDB = -80) {
  const thresholdLinear = 10 ** (thresholdDB / 20);
  
  for (let i = 0; i < complexFFT.length; i++) {
    if (complexFFT[i].magnitude < thresholdLinear) {
      complexFFT[i].magnitude = 0;
    }
  }
  
  return complexFFT;
}

// Pink noise weighting - cached log2 calculations
const LOG2_CACHE = new Map();
function applyPinkNoiseWeighting(freq) {
  if (freq < 20) return 0;
  
  if (!LOG2_CACHE.has(freq)) {
    // Limit cache size
    if (LOG2_CACHE.size > 10000) {
      const firstKey = LOG2_CACHE.keys().next().value;
      LOG2_CACHE.delete(firstKey);
    }
    
    const octaves = Math.log2(freq / 1000);
    const dbAdjustment = octaves * 6;
    LOG2_CACHE.set(freq, 10 ** (dbAdjustment / 20));
  }
  
  return LOG2_CACHE.get(freq);
}

// Frequency scaling
const fscale = (x) => Math.log2(x);

function logScale(x, scale = 1.0) {
  if (scale < 0.5) {
    const linearAmount = 1 - (scale * 2);
    const logAmount = scale * 2;
    return map(x, 20, 20000, 0, 1) * linearAmount + 
           map(fscale(x), fscale(20), fscale(20000), 0, 1) * logAmount;
  }
  return map(fscale(x), fscale(20), fscale(20000), 0, 1);
}

// ===== FFT Functions (optimized) =====

function calcComplexFFT(input) {
  const fft = new Float32Array(input);
  const fft2 = new Float32Array(input);
  transform(fft, fft2);
  
  const norm = input.length / 2;
  const result = new Array(input.length);
  
  for (let i = 0; i < input.length; i++) {
    const re = fft[i] / norm;
    const im = fft2[i] / norm;
    result[i] = {
      re,
      im,
      magnitude: Math.hypot(re, im),
      phase: Math.atan2(im, re)
    };
  }
  
  return result;
}

function transform(real, imag) {
  const n = real.length;
  if (n != imag.length) throw "Mismatched lengths";
  if (n <= 0) return;
  
  if ((2 ** Math.trunc(Math.log2(n))) === n) {
    transformRadix2(real, imag);
  } else {
    transformBluestein(real, imag);
  }
}

// Pre-allocate trig tables for common FFT sizes
const TRIG_TABLES = new Map();

function getTrigTables(n) {
  if (!TRIG_TABLES.has(n)) {
    const halfN = n / 2;
    const cosTable = new Float32Array(halfN);
    const sinTable = new Float32Array(halfN);
    const angle = 2 * Math.PI / n;
    
    for (let i = 0; i < halfN; i++) {
      cosTable[i] = Math.cos(angle * i);
      sinTable[i] = Math.sin(angle * i);
    }
    
    TRIG_TABLES.set(n, { cosTable, sinTable });
  }
  
  return TRIG_TABLES.get(n);
}

function transformRadix2(real, imag) {
  const n = real.length;
  if (n != imag.length) throw "Mismatched lengths";
  if (n <= 1) return;
  
  const logN = Math.log2(n);
  if ((2 ** Math.trunc(logN)) !== n) throw "Length is not a power of 2";
  
  const { cosTable, sinTable } = getTrigTables(n);
  
  // Bit reversal
  for (let i = 0; i < n; i++) {
    let j = 0;
    for (let k = 0, x = i; k < logN; k++, x >>>= 1) {
      j = (j << 1) | (x & 1);
    }
    
    if (j > i) {
      [real[i], real[j]] = [real[j], real[i]];
      [imag[i], imag[j]] = [imag[j], imag[i]];
    }
  }
  
  // FFT
  for (let size = 2; size <= n; size *= 2) {
    const halfsize = size / 2;
    const tablestep = n / size;
    
    for (let i = 0; i < n; i += size) {
      for (let j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
        const l = j + halfsize;
        const cos = cosTable[k];
        const sin = sinTable[k];
        const tpre = real[l] * cos + imag[l] * sin;
        const tpim = -real[l] * sin + imag[l] * cos;
        
        real[l] = real[j] - tpre;
        imag[l] = imag[j] - tpim;
        real[j] += tpre;
        imag[j] += tpim;
      }
    }
  }
}

function transformBluestein(real, imag) {
  const n = real.length;
  if (n != imag.length) throw "Mismatched lengths";
  
  const m = 2 ** Math.trunc(Math.log2(n * 2) + 1);
  
  const cosTable = new Float32Array(n);
  const sinTable = new Float32Array(n);
  const piOverN = Math.PI / n;
  
  for (let i = 0; i < n; i++) {
    const j = (i * i) % (n * 2);
    const angle = piOverN * j;
    cosTable[i] = Math.cos(angle);
    sinTable[i] = Math.sin(angle);
  }
  
  const areal = new Float32Array(m);
  const aimag = new Float32Array(m);
  
  for (let i = 0; i < n; i++) {
    const cos = cosTable[i];
    const sin = sinTable[i];
    areal[i] = real[i] * cos + imag[i] * sin;
    aimag[i] = -real[i] * sin + imag[i] * cos;
  }
  
  const breal = new Float32Array(m);
  const bimag = new Float32Array(m);
  breal[0] = cosTable[0];
  bimag[0] = sinTable[0];
  
  for (let i = 1; i < n; i++) {
    breal[i] = breal[m - i] = cosTable[i];
    bimag[i] = bimag[m - i] = sinTable[i];
  }
  
  const creal = new Float32Array(m);
  const cimag = new Float32Array(m);
  convolveComplex(areal, aimag, breal, bimag, creal, cimag);
  
  for (let i = 0; i < n; i++) {
    const cos = cosTable[i];
    const sin = sinTable[i];
    real[i] = creal[i] * cos + cimag[i] * sin;
    imag[i] = -creal[i] * sin + cimag[i] * cos;
  }
}

function convolveComplex(xreal, ximag, yreal, yimag, outreal, outimag) {
  const n = xreal.length;
  
  // Create copies
  const xr = new Float32Array(xreal);
  const xi = new Float32Array(ximag);
  const yr = new Float32Array(yreal);
  const yi = new Float32Array(yimag);
  
  transform(xr, xi);
  transform(yr, yi);
  
  for (let i = 0; i < n; i++) {
    const temp = xr[i] * yr[i] - xi[i] * yi[i];
    xi[i] = xi[i] * yr[i] + xr[i] * yi[i];
    xr[i] = temp;
  }
  
  inverseTransform(xr, xi);
  
  for (let i = 0; i < n; i++) {
    outreal[i] = xr[i] / n;
    outimag[i] = xi[i] / n;
  }
}

function inverseTransform(real, imag) {
  transform(imag, real);
}

// Pre-calculated window function
const WINDOW_CACHE = new Map();

function getWindow(size) {
  if (!WINDOW_CACHE.has(size)) {
    const window = new Float32Array(size);
    for (let i = 0; i < size; i++) {
      const x = map(i, 0, size - 1, -1, 1);
      const cos = Math.cos(x * Math.PI / 2);
      window[i] = cos * cos;
    }
    WINDOW_CACHE.set(size, window);
  }
  return WINDOW_CACHE.get(size);
}

function generateSpectrum(data, deltaWindow, gain, naturalWeighting, sampleRate, fftSize, lowEndBoost, noiseGateDB) {
  const window = getWindow(data.length);
  const dataArray = new Float32Array(data.length);
  let norm = 0;
  
  if (deltaWindow) {
    // Delta window mode
    for (let i = 0; i < data.length; i++) {
      const x1 = map(i - 0.5, 0, data.length - 1, -1, 1);
      const x2 = map(i + 0.5, 0, data.length - 1, -1, 1);
      const cos1 = Math.cos(x1 * Math.PI / 2);
      const cos2 = Math.cos(x2 * Math.PI / 2);
      const w1 = cos1 * cos1;
      const w2 = cos2 * cos2;
      const sample = isFinite(data[i]) ? data[i] * gain : 0;
      dataArray[i] = sample * (w1 - w2);
      norm += window[i];
    }
  } else {
    // Standard window mode
    for (let i = 0; i < data.length; i++) {
      const sample = isFinite(data[i]) ? data[i] * gain : 0;
      dataArray[i] = sample * window[i];
      norm += window[i];
    }
  }
  
  // Normalize and calculate FFT
  const normFactor = norm > 0 ? (data.length / (norm * Math.SQRT2)) : 0;
  for (let i = 0; i < dataArray.length; i++) {
    dataArray[i] *= normFactor;
  }
  
  let complexFFT = calcComplexFFT(dataArray);
  
  // Apply natural weighting
  if (naturalWeighting) {
    const freqPerBin = sampleRate / fftSize;
    for (let i = 0; i < complexFFT.length; i++) {
      const freq = i * freqPerBin;
      complexFFT[i].magnitude *= applyPinkNoiseWeighting(freq);
    }
  }
  
  // Apply noise gate
  if (noiseGateDB > -120) {
    complexFFT = applyNoiseGate(complexFFT, noiseGateDB);
  }
  
  // Apply low-end boost
  if (lowEndBoost > 1.0) {
    complexFFT = applyLowEndBoost(complexFFT, sampleRate, fftSize, lowEndBoost);
  }

  return complexFFT;
}

// ===== Spectrum Generation Functions =====

// Cache spectrum bar data
const SPECTRUM_BAR_CACHE = new Map();

function generateSpectrumBarData(size, length, sampleRate, frequencyScale) {
  const key = `${size}_${length}_${sampleRate}_${frequencyScale}`;
  
  if (SPECTRUM_BAR_CACHE.has(key)) {
    return SPECTRUM_BAR_CACHE.get(key);
  }
  
  const min = 20;
  const max = 20000;
  const spectrogramBars = [];
  
  for (let pixel = 0; pixel < size; pixel++) {
    const pixelStart = pixel / size;
    const pixelEnd = (pixel + 1) / size;
    
    let freqStart, freqEnd;
    
    if (frequencyScale < 0.5) {
      const linearAmount = 1 - (frequencyScale * 2);
      const logAmount = frequencyScale * 2;
      
      freqStart = (min + pixelStart * (max - min)) * linearAmount + 
                  (min * Math.pow(max / min, pixelStart)) * logAmount;
      freqEnd = (min + pixelEnd * (max - min)) * linearAmount + 
                (min * Math.pow(max / min, pixelEnd)) * logAmount;
    } else {
      freqStart = min * Math.pow(max / min, pixelStart);
      freqEnd = min * Math.pow(max / min, pixelEnd);
    }
    
    const binStart = freqStart * length / sampleRate;
    const binEnd = freqEnd * length / sampleRate;
    
    spectrogramBars.push({
      lo: Math.floor(binStart),
      hi: Math.ceil(binEnd),
      loFrac: binStart - Math.floor(binStart),
      hiFrac: Math.ceil(binEnd) - binEnd,
      start: pixel,           // <-- Changed: use integer pixel positions
      end: pixel + 1          // <-- Changed: exactly 1 pixel height
    });
  }
  
  SPECTRUM_BAR_CACHE.set(key, spectrogramBars);
  return spectrogramBars;
}

function calcReassignedFreqs(data, auxData, bufferSize, sampleRate) {
  const freqs = new Float32Array(data.length);
  const correction = sampleRate / (2 * Math.PI);
  
  for (let i = 0; i < freqs.length; i++) {
    const denom = data[i].re * data[i].re + data[i].im * data[i].im;
    const corr = denom > 0 ? 
      -(data[i].re * auxData[i].im - auxData[i].re * data[i].im) / denom : 0;
    freqs[i] = fftBinToHertz(i, bufferSize, sampleRate) - corr * correction;
  }
  
  return freqs;
}

function calcReassignedSpectrum(fft, freqs, bands, frequencyScale) {
  const buckets = new Float32Array(bands);
  const bandsFactor = bands;
  
  for (let i = 0; i < fft.length; i++) {
    if (!isFinite(fft[i])) continue;
    
    const x = logScale(freqs[i], frequencyScale);
    const posX = Math.round(x * bandsFactor);
    
    if (posX >= 0 && posX < bands) {
      buckets[posX] = Math.max(buckets[posX], fft[i]);
    }
  }
  
  return buckets;
}

function renderOscilloscope(data, oscilloscopeIdx) {
  const dataset = new Float32Array(data.length);
  const len = data.length;
  
  for (let i = 0; i < len; i++) {
    dataset[i] = data[idxWrapOver(i + oscilloscopeIdx - len, len)];
  }
  
  return dataset;
}

// ===== Note/Frequency Functions =====

const NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'];
const C0_FREQ = 440 * 2 ** -4.75;

function frequencyToNote(freq) {
  if (freq < 16.35) return null;
  
  const halfSteps = 12 * Math.log2(freq / C0_FREQ);
  const octave = Math.floor(halfSteps / 12);
  const noteIndex = Math.round(halfSteps % 12);
  const cents = Math.round((halfSteps % 1) * 100);
  
  return {
    note: `${NOTE_NAMES[noteIndex % 12]}${octave}`,
    freq: freq.toFixed(2),
    cents: cents !== 0 ? cents : null
  };
}

function yPositionToFrequency(yPos, canvasHeight, frequencyScale) {
  const normalizedY = 1 - (yPos / canvasHeight);
  const minFreq = 20;
  const maxFreq = 20000;
  
  if (frequencyScale < 0.5) {
    const linearAmount = 1 - (frequencyScale * 2);
    const logAmount = frequencyScale * 2;
    const linearFreq = minFreq + normalizedY * (maxFreq - minFreq);
    const logFreq = minFreq * (maxFreq / minFreq) ** normalizedY;
    return linearFreq * linearAmount + logFreq * logAmount;
  }
  
  return minFreq * (maxFreq / minFreq) ** normalizedY;
}

// ===== AGC Functions =====

function calculateRMS(data) {
  let sum = 0;
  for (let i = 0; i < data.length; i++) {
    sum += data[i] * data[i];
  }
  return Math.sqrt(sum / data.length);
}

function calculatePeakLevel(data) {
  let peak = 0;
  for (let i = 0; i < data.length; i++) {
    const abs = Math.abs(data[i]);
    if (abs > peak) peak = abs;
  }
  return peak;
}

// ===== Colormap System =====

const COLORMAPS = {
  aging: [[0.00, 0, 30, 70], [0.50, 130, 120, 90], [1.00, 255, 235, 160]],
  arctic: [[0.00, 0, 0, 0], [0.00, 0, 0, 50], [0.50, 100, 150, 200], [1.00, 200, 240, 255]],
  ayeye: [[0.00, 5, 0, 30], [0.33, 80, 0, 150], [0.66, 150, 50, 255], [1.00, 200, 200, 255]],
  bubblegum: [[0.00, 0, 0, 0], [0.00, 30, 0, 30], [0.50, 200, 100, 180], [1.00, 255, 180, 240]],
  chroma: [[0.00, 30, 0, 60], [0.25, 150, 0, 200], [0.50, 255, 100, 150], [0.75, 255, 200, 0], [1.00, 150, 255, 100]],
  cold: [[0.00, 0, 0, 30], [0.25, 50, 80, 150], [0.50, 100, 180, 230], [0.75, 150, 220, 255], [1.00, 200, 240, 255]],
  craft: [[0.00, 0, 0, 0], [0.50, 100, 0, 150], [1.00, 220, 150, 255]],
  dusk: [[0.00, 10, 10, 30], [0.50, 150, 80, 100], [1.00, 255, 180, 150]],
  ember: [[0.00, 20, 0, 0], [0.50, 200, 60, 0], [1.00, 255, 180, 100]],
  fall: [[0.00, 30, 0, 0], [0.33, 180, 100, 0], [0.66, 255, 180, 50], [1.00, 255, 220, 180]],
  flamangos: [[0.00, 30, 0, 20], [0.50, 255, 120, 180], [1.00, 255, 200, 220]],
  frogger: [[0.00, 0, 20, 10], [0.33, 0, 100, 80], [0.66, 100, 200, 100], [1.00, 200, 255, 150]],
  gemini: [[0.00, 20, 0, 40], [0.33, 120, 50, 200], [0.66, 180, 150, 255], [1.00, 230, 220, 255]],
  gothic: [[0.00, 20, 0, 30], [0.50, 120, 50, 150], [1.00, 200, 150, 230]],
  granny: [[0.00, 0, 0, 0], [0.50, 100, 180, 50], [1.00, 200, 255, 150]],
  greyscale: [[0.00, 0, 0, 0], [1.00, 200, 200, 200]],
  inferno: [[0.00, 0, 0, 0], [0.25, 80, 0, 80], [0.50, 200, 50, 50], [0.75, 255, 150, 0], [1.00, 255, 255, 200]],
  jurassic: [[0.00, 0, 0, 0], [0.50, 180, 60, 0], [1.00, 255, 200, 100]],
  light: [[0.00, 30, 20, 40], [0.50, 180, 150, 200], [1.00, 240, 220, 255]],
  magma: [[0.00, 0, 0, 10], [0.25, 80, 20, 80], [0.50, 200, 70, 100], [0.75, 255, 150, 120], [1.00, 255, 240, 200]],
  neon: [[0.00, 0, 0, 0], [0.33, 255, 0, 150], [0.66, 0, 255, 200], [1.00, 150, 255, 255]],
  original: [[0.00, 30, 30, 120], [0.25, 50, 150, 255], [0.50, 100, 255, 150], [0.75, 255, 230, 50], [1.00, 200, 50, 0]],
  pepper: [[0.00, 20, 0, 0], [0.50, 150, 100, 50], [1.00, 255, 220, 150]],
  plasma: [[0.00, 10, 0, 80], [0.25, 150, 0, 150], [0.50, 255, 80, 100], [0.75, 255, 180, 50], [1.00, 255, 255, 100]],
  sepia: [[0.00, 20, 15, 10], [0.50, 150, 120, 90], [1.00, 240, 220, 200]],
  sharp: [[0.00, 0, 0, 50], [0.33, 100, 0, 200], [0.66, 200, 100, 255], [1.00, 255, 230, 255]],
  slimer: [[0.00, 10, 20, 10], [0.50, 150, 200, 150], [1.00, 240, 255, 240]],
  smooth: [[0.00, 0, 0, 20], [0.33, 100, 50, 100], [0.66, 255, 120, 100], [1.00, 255, 230, 180]],
  soft: [[0.00, 0, 0, 0], [0.25, 60, 50, 100], [0.50, 150, 120, 140], [0.75, 230, 180, 160], [1.00, 255, 255, 240]],
  space: [[0.00, 0, 0, 0], [0.50, 50, 40, 80], [1.00, 255, 200, 100]],
  sunrise: [[0.00, 0, 20, 40], [0.33, 50, 100, 180], [0.66, 255, 150, 120], [1.00, 255, 200, 200]],
  torch: [[0.00, 20, 0, 0], [0.50, 255, 100, 0], [1.00, 255, 240, 180]],
  toxic: [[0.00, 0, 0, 0], [0.50, 150, 255, 0], [1.00, 230, 255, 150]],
  tropical: [[0.00, 0, 30, 30], [0.33, 255, 100, 150], [0.66, 255, 200, 100], [1.00, 255, 255, 200]],
  twentyfive: [[0.00, 70, 15, 90], [0.25, 60, 90, 140], [0.50, 40, 160, 140], [0.75, 140, 210, 100], [1.00, 255, 230, 60]]
};

function getColormapColor(value, colormapName) {
  value = clamp(value, 0, 1);
  
  const colormap = COLORMAPS[colormapName] || COLORMAPS.greyscale || COLORMAPS.inferno;
  
  if (!colormap || !colormap.length) {
    console.error(`Invalid colormap: ${colormapName}, falling back to greyscale`);
    return 'rgb(128,128,128)';
  }
  
  // Find interpolation range
  let i = 0;
  while (i < colormap.length - 1 && value > colormap[i + 1][0]) {
    i++;
  }
  
  if (i === colormap.length - 1) {
    const [, r, g, b] = colormap[i];
    return `rgb(${r},${g},${b})`;
  }
  
  // Interpolate
  const [pos1, r1, g1, b1] = colormap[i];
  const [pos2, r2, g2, b2] = colormap[i + 1];
  const t = (value - pos1) / (pos2 - pos1);
  
  const r = (r1 + (r2 - r1) * t) | 0;
  const g = (g1 + (g2 - g1) * t) | 0;
  const b = (b1 + (b2 - b1) * t) | 0;
  
  return `rgb(${r},${g},${b})`;
}

// ===== Spectrogram Rendering =====

function printSpectrogram(auxCtx, auxCanvas, auxCtx2, auxCanvas2, data, linearBinData, colormap, dbRange, shouldShift) {
  const height = auxCanvas.height;
  const width = auxCanvas.width;
  
  if (shouldShift && width > 0 && height > 0) {
    // Copy from canvas1 to canvas2, shifted left by 1 pixel
    auxCtx2.drawImage(auxCanvas, 1, 0, width - 1, height, 0, 0, width - 1, height);
    
    // Clear rightmost column on canvas2
    auxCtx2.fillStyle = '#000000';
    auxCtx2.fillRect(width - 1, 0, 1, height);
    
    // Swap: copy canvas2 back to canvas1
    auxCtx.drawImage(auxCanvas2, 0, 0);
  }
  
  const xPos = width - 1;
  
  if (linearBinData) {
    for (let pixelY = 0; pixelY < height; pixelY++) {
      let mag = 0;
      
      const binIndex = Math.floor((pixelY / height) * linearBinData.length);
      if (binIndex < linearBinData.length) {
        const bin = linearBinData[binIndex];
        
        for (let idx = bin.lo; idx <= bin.hi; idx++) {
          const binIdx = idxWrapOver(idx, data.length);
          mag = Math.max(mag, data[binIdx]);
        }
      }
      
      const amp = ascale(mag, dbRange);
      const color = getColormapColor(isFinite(amp) ? amp : 0, colormap);
      
      auxCtx.fillStyle = color;
      auxCtx.fillRect(xPos, height - pixelY - 1, 1, 1);
    }
  } else {
    for (let pixelY = 0; pixelY < height; pixelY++) {
      const binIndex = Math.floor((pixelY / height) * data.length);
      const mag = data[binIndex] || 0;
      
      const amp = ascale(mag, dbRange);
      const color = getColormapColor(isFinite(amp) ? amp : 0, colormap);
      
      auxCtx.fillStyle = color;
      auxCtx.fillRect(xPos, height - pixelY - 1, 1, 1);
    }
  }
}
// ===== Main AudioSpectrogram Class =====

class AudioSpectrogram {
  constructor() {
    this.canvas = document.getElementById('canvas');
    this.ctx = this.canvas.getContext('2d', { 
      alpha: false,
      desynchronized: this.isMac
    });
    this.ctx.imageSmoothingEnabled = false;

    this.canvas.addEventListener('webglcontextlost', (e) => {
      e.preventDefault();
      console.warn('Canvas context lost');
      this.isPaused = true;
    });

    this.canvas.addEventListener('webglcontextrestored', () => {
      console.log('Canvas context restored');
      this.resizeCanvas();
      this.isPaused = false;
    });

    this.statusEl = document.getElementById('status');

    this.audioCtx = new (window.AudioContext || window.webkitAudioContext)();
    this.analyser = this.audioCtx.createAnalyser();
    // Detect platform for optimizations
    this.isMac = navigator.platform.toLowerCase().includes('mac');

    const defaultSettings = {
      fftSize: 4096,
      hopSize: 64,
      // hopSizeRatio: 0.125,
      scrollSpeed: this.isMac ? (window.devicePixelRatio > 1.5 ? 3.0 : 2.0) : 1.0,  // Higher on Retina
      useReassignment: true,
      colormap: 'inferno',
      dbRange: 50,
      gain: 10.0,
      naturalWeighting: true,
      frequencyScale: 1.0,
      lowEndBoost: 10.0,
      smoothing: 0.35,
      noiseGate: -60,
      alwaysOnTop: true,
      enableAGC: true,
      agcStrength: 1.0,
      brightness: 0.44,
      sampleRate: this.audioCtx.sampleRate
    };

    this.settings = { ...defaultSettings };

    this.isPaused = false;
    this.previewColormap = null;

    this.state = {
      isRunning: false,
      spectrum: [],
      prevSpectrum: [],
      oscilloscopeBuffer: new Float32Array(this.settings.fftSize),
      oscilloscopeIdx: 0,
      sampleCounter: 0,
      staticSpectrogramIdx: 0,
      isReassigned: false,
      targetRMS: 0.1,
      currentGain: 1.0,
      gainSmoothingFactor: 0.95,
      frameBuffer: [],           // Stores computed frames
      frameSkipCounter: 0.0,     // Accumulates scroll speed
      skipFrameCounter: 0,
      lastRenderTime: null
    };

    this.auxCanvas = new OffscreenCanvas(this.canvas.width, this.canvas.height);
    this.auxCtx = this.auxCanvas.getContext('2d', { 
      alpha: false,
      willReadFrequently: true,  // Optimizes for frequent readback
      desynchronized: this.isMac
    });
    this.auxCtx.imageSmoothingEnabled = false;

    // ADD THESE LINES HERE - Second buffer for double-buffering
    this.auxCanvas2 = new OffscreenCanvas(this.canvas.width, this.canvas.height);
    this.auxCtx2 = this.auxCanvas2.getContext('2d', { 
      alpha: false,
      willReadFrequently: true,
      desynchronized: this.isMac
    });
    this.auxCtx2.imageSmoothingEnabled = false;

    this.setupUI();
    this.resizeCanvas();
    this.animate();
    
    this.initSettings(defaultSettings);

    // Listen for setting changes from settings window
    ipcRenderer.on('setting-changed', (event, { key, value }) => {
      this.settings[key] = value;
      
      // Handle specific setting changes
      if (key === 'fftSize') {
        this.updateFFTSize();
      } else if (key === 'scrollSpeed') {
        // Adjust hopSize based on scroll speed for fast scrolling
        this.settings.hopSize = Math.max(16, Math.floor(64 / Math.max(1, value)));
      } else if (key === 'frequencyScale') {
        SPECTRUM_BAR_CACHE.clear();
      } else if (key === 'brightness') {
        this.updateBrightness();
      } else if (key === 'enableAGC' && !value) {
        this.state.currentGain = 1.0;
      }
    });

    // Handle colormap preview
    ipcRenderer.on('preview-colormap', (event, colormap) => {
      this.previewColormap = colormap;
    });

    // Handle settings reset
    ipcRenderer.on('settings-reset', (event, defaultSettings) => {
      this.settings = { ...defaultSettings };
      this.updateFFTSize();
      this.updateBrightness();
      SPECTRUM_BAR_CACHE.clear();
    });
  }

  updateFFTSize() {
    // Update analyser FFT size (must be power of 2, between 32 and 32768)
    const analyserSize = Math.max(32, Math.min(32768, this.settings.fftSize * 2));
    this.analyser.fftSize = analyserSize;
    
    // Resize oscilloscope buffer
    const oldBuffer = this.state.oscilloscopeBuffer;
    this.state.oscilloscopeBuffer = new Float32Array(this.settings.fftSize);
    
    // Copy over old data if possible
    const copyLength = Math.min(oldBuffer.length, this.state.oscilloscopeBuffer.length);
    for (let i = 0; i < copyLength; i++) {
      this.state.oscilloscopeBuffer[i] = oldBuffer[i];
    }
    
    // Reset index if it's out of bounds
    if (this.state.oscilloscopeIdx >= this.state.oscilloscopeBuffer.length) {
      this.state.oscilloscopeIdx = 0;
    }
    
    // Clear spectrum cache
    this.state.spectrum = [];
    this.state.prevSpectrum = [];
    SPECTRUM_BAR_CACHE.clear();
    
    console.log(`FFT size updated to ${this.settings.fftSize}, analyser size: ${analyserSize}`);
  }

  async initSettings(defaultSettings) {
    await this.loadSettings(defaultSettings);
    this.updateUIFromSettings();
    this.updateFFTSize();
    setTimeout(() => this.start(), 500);
  }

  updateUIFromSettings() {
    const updateElement = (id, value, textId, formatter) => {
      const el = document.getElementById(id);
      const textEl = document.getElementById(textId);
      if (el) el.value = value;
      if (textEl) textEl.textContent = formatter ? formatter(value) : value;
    };

    // Update dropdowns (FFT size)
    const fftSelect = document.getElementById('fftSizeSelect');
    if (fftSelect) {
      fftSelect.value = this.settings.fftSize;
      document.getElementById('fftSizeValue').textContent = this.settings.fftSize;
    }

    // Update colormap
    updateElement('colormapSelect', this.settings.colormap, 'colormapValue', 
      v => v.charAt(0).toUpperCase() + v.slice(1));
    
    // Update sliders
    updateElement('dbRange', this.settings.dbRange, 'dbRangeValue');
    updateElement('gain', this.settings.gain, 'gainValue', v => v.toFixed(1));
    updateElement('freqScale', this.settings.frequencyScale, 'freqScaleValue', v => v.toFixed(1));
    updateElement('lowEnd', this.settings.lowEndBoost, 'lowEndValue', v => v.toFixed(1) + 'x');
    updateElement('smoothing', this.settings.smoothing, 'smoothingValue', v => v.toFixed(2));
    updateElement('noiseGate', this.settings.noiseGate, 'noiseGateValue', v => v + ' dB');
    updateElement('agcStrength', this.settings.agcStrength, 'agcStrengthValue', v => v.toFixed(2));
    updateElement('scrollSpeed', this.settings.scrollSpeed, 'scrollSpeedValue', v => v.toFixed(1) + 'x');
    updateElement('brightness', this.settings.brightness, 'brightnessValue', v => Math.round(v * 100) + '%');

    // Update checkboxes
    const checkboxes = [
      ['reassignedCheck', 'useReassignment'],
      ['naturalCheck', 'naturalWeighting'],
      ['alwaysOnTopCheck', 'alwaysOnTop'],
      ['agcCheck', 'enableAGC']
    ];
    
    checkboxes.forEach(([id, setting]) => {
      const el = document.getElementById(id);
      if (el) el.checked = this.settings[setting];
    });

    // Update colormap dropdown selected display
    const dropdownSelected = document.getElementById('colormapSelected');
    if (dropdownSelected) {
      const selectedName = dropdownSelected.querySelector('.colormap-name');
      const selectedGradient = dropdownSelected.querySelector('.colormap-gradient');
      if (selectedName) {
        selectedName.textContent = this.settings.colormap.charAt(0).toUpperCase() + this.settings.colormap.slice(1);
      }
      if (selectedGradient) {
        selectedGradient.style.background = this.generateColormapGradient(this.settings.colormap);
      }
    }

    // Update which option has the selected class
    const dropdownOptions = document.getElementById('colormapOptions');
    if (dropdownOptions) {
      const options = dropdownOptions.querySelectorAll('.dropdown-option');
      options.forEach(option => {
        if (option.dataset.colormap === this.settings.colormap) {
          option.classList.add('selected');
        } else {
          option.classList.remove('selected');
        }
      });
    }

    this.updateBrightness();
  }

  async loadSettings(defaultSettings) {
    try {
      const savedSettings = await ipcRenderer.invoke('load-settings');
      if (savedSettings) {
        this.settings = { ...defaultSettings, ...savedSettings };
        console.log('Settings loaded successfully');
      } else {
        console.log('No saved settings found, using defaults');
      }
    } catch (err) {
      console.error('Error loading settings:', err);
      this.settings = { ...defaultSettings };
    }
  }

  async saveSettings() {
    try {
      const settingsToSave = {
        fftSize: this.settings.fftSize,
        hopSize: this.settings.hopSize,
        scrollSpeed: this.settings.scrollSpeed,
        useReassignment: this.settings.useReassignment,
        colormap: this.settings.colormap,
        dbRange: this.settings.dbRange,
        gain: this.settings.gain,
        naturalWeighting: this.settings.naturalWeighting,
        frequencyScale: this.settings.frequencyScale,
        lowEndBoost: this.settings.lowEndBoost,
        smoothing: this.settings.smoothing,
        noiseGate: this.settings.noiseGate,
        alwaysOnTop: this.settings.alwaysOnTop,
        enableAGC: this.settings.enableAGC,
        agcStrength: this.settings.agcStrength,
        brightness: this.settings.brightness
      };
      
      const result = await ipcRenderer.invoke('save-settings', settingsToSave);
      if (!result.success) {
        console.error('Failed to save settings:', result.error);
      }
    } catch (err) {
      console.error('Error saving settings:', err);
    }
  }

  setupUI() {
    const exitBtn = document.getElementById('exitBtn');
    const settingsBtn = document.getElementById('settingsBtn');
    const pauseBtn = document.getElementById('pauseBtn');
    
    exitBtn?.addEventListener('click', () => window.close());
    settingsBtn?.addEventListener('click', () => {
      ipcRenderer.send('open-settings');
    });
    pauseBtn?.addEventListener('click', () => this.togglePause());
    
    window.addEventListener('resize', () => this.resizeCanvas());
    
    this.freqHover = document.getElementById('freqHover');
    this.freqLine = document.getElementById('freqLine');
    this.freqLabel = document.getElementById('freqLabel');
    
    const canvas = this.canvas;
    let isOverButton = false;
    let shiftPressed = false;
    
    const buttons = [exitBtn, settingsBtn, pauseBtn];
    
    buttons.forEach(btn => {
      btn?.addEventListener('mouseenter', () => {
        isOverButton = true;
        this.freqHover?.classList.remove('active');
      });
      
      btn?.addEventListener('mouseleave', () => {
        isOverButton = false;
      });
    });
    
    // Track Shift key
    window.addEventListener('keydown', (e) => {
      if (e.key === 'Shift') {
        shiftPressed = true;
        // Disable drag when shift is held
        canvas.style.webkitAppRegion = 'no-drag';
      }
    });
    
    window.addEventListener('keyup', (e) => {
      if (e.key === 'Shift') {
        shiftPressed = false;
        // Re-enable drag when shift is released
        canvas.style.webkitAppRegion = 'drag';
        this.freqHover?.classList.remove('active');
      }
    });
    
    canvas.addEventListener('mousemove', (e) => {
      if (shiftPressed && !isOverButton) {
        this.freqHover?.classList.add('active');
        this.updateFrequencyDisplay(e);
      } else {
        this.freqHover?.classList.remove('active');
      }
    });
    
    canvas.addEventListener('mouseleave', () => {
      this.freqHover?.classList.remove('active');
    });

    this.setupResizeHandles();
  }

  setupResizeHandles() {
    const resizeHandles = document.querySelectorAll('.resize-handle');
    
    resizeHandles.forEach(handle => {
      const direction = handle.dataset.direction;
      
      handle.addEventListener('mousedown', (e) => {
        e.preventDefault();
        e.stopPropagation();
        
        // Get initial window bounds via IPC
        ipcRenderer.invoke('get-window-bounds').then(startBounds => {
          const startX = e.screenX;
          const startY = e.screenY;
          
          const onMouseMove = (moveEvent) => {
            const deltaX = moveEvent.screenX - startX;
            const deltaY = moveEvent.screenY - startY;
            
            let newBounds = { ...startBounds };
            
            // Calculate new bounds based on direction
            if (direction.includes('n')) {
              newBounds.y = startBounds.y + deltaY;
              newBounds.height = startBounds.height - deltaY;
            }
            if (direction.includes('s')) {
              newBounds.height = startBounds.height + deltaY;
            }
            if (direction.includes('e')) {
              newBounds.width = startBounds.width + deltaX;
            }
            if (direction.includes('w')) {
              newBounds.x = startBounds.x + deltaX;
              newBounds.width = startBounds.width - deltaX;
            }
            
            // Enforce minimum size
            if (newBounds.width < 150) {
              if (direction.includes('w')) {
                newBounds.x = startBounds.x + startBounds.width - 150;
              }
              newBounds.width = 150;
            }
            if (newBounds.height < 150) {
              if (direction.includes('n')) {
                newBounds.y = startBounds.y + startBounds.height - 150;
              }
              newBounds.height = 150;
            }
            
            // Send resize command via IPC
            ipcRenderer.send('resize-window-bounds', newBounds);
          };
          
          const onMouseUp = () => {
            document.removeEventListener('mousemove', onMouseMove);
            document.removeEventListener('mouseup', onMouseUp);
          };
          
          document.addEventListener('mousemove', onMouseMove);
          document.addEventListener('mouseup', onMouseUp);
        });
      });
    });
  }

  togglePause() {
    this.isPaused = !this.isPaused;
    const pauseBtn = document.getElementById('pauseBtn');
    if (pauseBtn) {
      if (this.isPaused) {
        pauseBtn.classList.add('paused');
      } else {
        pauseBtn.classList.remove('paused');
      }
    }
  }

  updateFrequencyDisplay(e) {
    if (!this.freqHover || !this.freqLine || !this.freqLabel) return;
    
    const rect = this.canvas.getBoundingClientRect();
    const y = e.clientY;
    
    this.freqLine.style.top = y + 'px';
    
    const freq = yPositionToFrequency(y, rect.height, this.settings.frequencyScale);
    const noteInfo = frequencyToNote(freq);
    
    if (noteInfo) {
      let labelText = `<span class="note">${noteInfo.note}</span> ${noteInfo.freq} Hz`;
      if (noteInfo.cents) {
        const centsSign = noteInfo.cents > 0 ? '+' : '';
        labelText += `<span class="cents">(${centsSign}${noteInfo.cents}¢)</span>`;
      }
      this.freqLabel.innerHTML = labelText;
      
      const labelHeight = 30;
      const padding = 10;
      let labelY = y;
      
      if (y < labelHeight / 2 + padding) {
        labelY = labelHeight / 2 + padding;
      } else if (y > rect.height - labelHeight / 2 - padding) {
        labelY = rect.height - labelHeight / 2 - padding;
      }
      
      this.freqLabel.style.top = labelY + 'px';
    }
  }

  updateBrightness() {
    const blurLayer = document.getElementById('blur-layer');
    if (blurLayer) {
      const brightnessHex = Math.round(this.settings.brightness * 255).toString(16).padStart(2, '0');
      blurLayer.style.backgroundImage = `linear-gradient(to right, #00000017 0 3rem, #ffffff${brightnessHex} 15rem 3rem)`;
    }
  }

  resizeCanvas() {
    const dpr = window.devicePixelRatio || 1;
    
    // Main canvas at device pixel ratio
    this.canvas.width = window.innerWidth * dpr;
    this.canvas.height = window.innerHeight * dpr;
    this.canvas.style.width = window.innerWidth + 'px';
    this.canvas.style.height = window.innerHeight + 'px';
    
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    
    // Save old content
    const oldAuxCanvas = this.auxCanvas;
    const oldWidth = oldAuxCanvas.width;
    const oldHeight = oldAuxCanvas.height;
    
    const tempCanvas = new OffscreenCanvas(oldWidth, oldHeight);
    const tempCtx = tempCanvas.getContext('2d', { alpha: false });
    if (oldWidth > 0 && oldHeight > 0) {
      tempCtx.drawImage(oldAuxCanvas, 0, 0);
    }
    
    // Resize BOTH aux canvases
    this.auxCanvas.width = this.canvas.width;
    this.auxCanvas.height = this.canvas.height;
    this.auxCtx.imageSmoothingEnabled = false;
    
    this.auxCanvas2.width = this.canvas.width;  // <-- Added
    this.auxCanvas2.height = this.canvas.height;  // <-- Added
    this.auxCtx2.imageSmoothingEnabled = false;  // <-- Added
    
    // Restore old content
    if (oldWidth > 0 && oldHeight > 0) {
      this.auxCtx.drawImage(tempCanvas, 
        0, 0, oldWidth, oldHeight,
        0, 0, this.canvas.width, this.canvas.height
      );
    }
    
    SPECTRUM_BAR_CACHE.clear();
  }

  async start() {
    try {
      this.statusEl.textContent = 'Starting audio capture...';
      
      const stream = await getLoopbackAudioMediaStream();
      
      if (!stream) {
        this.statusEl.textContent = 'Error: Could not get audio stream';
        return;
      }

      const source = this.audioCtx.createMediaStreamSource(stream);
      source.connect(this.analyser);

      this.state.isRunning = true;
      this.statusEl.textContent = 'Capturing audio ▮';
      
      this.captureLoop();
    } catch (err) {
      this.statusEl.textContent = `Error: ${err.message}`;
      console.error('Start capture error:', err);
    }
  }

  captureLoop() {
    if (!this.state.isRunning) return;

    const timeDomainData = new Float32Array(this.analyser.fftSize);
    this.analyser.getFloatTimeDomainData(timeDomainData);

    // Copy all samples into oscilloscope buffer
    for (let i = 0; i < timeDomainData.length; i++) {
      this.state.oscilloscopeBuffer[this.state.oscilloscopeIdx] = timeDomainData[i];
      this.state.oscilloscopeIdx = idxWrapOver(this.state.oscilloscopeIdx + 1, this.settings.fftSize);
      this.state.sampleCounter++;

      // Analyze at hop intervals
      if (this.state.sampleCounter >= this.settings.hopSize) {
        this.state.sampleCounter = 0;
        
        // Calculate time since last render
        const now = performance.now();
        const timeSinceLastRender = this.state.lastRenderTime ? (now - this.state.lastRenderTime) : 0;
        
        // Calculate desired time between renders based on scroll speed
        // At scroll speed 1.0, render approximately every 16ms (60fps)
        // Adjust based on FFT size - larger FFTs need less frequent renders
        const baseRenderInterval = 16; // ms
        const fftSizeFactor = this.settings.fftSize / 4096; // normalize to 4096
        const targetRenderInterval = (baseRenderInterval * fftSizeFactor) / this.settings.scrollSpeed;
        
        // Only render if enough time has passed
        if (!this.state.lastRenderTime || timeSinceLastRender >= targetRenderInterval) {
          this.analyzeFrame();
          this.state.lastRenderTime = now;
        }
      }
    }

    if (this.isMac) {
      setTimeout(() => this.captureLoop(), 8);
    } else {
      requestAnimationFrame(() => this.captureLoop());
    }
  }

  analyzeFrame() {
    if (this.isPaused) return;
    
    this.state.isReassigned = this.settings.useReassignment;
    const pcmData = renderOscilloscope(this.state.oscilloscopeBuffer, this.state.oscilloscopeIdx);

    let effectiveGain = this.settings.gain;
    
    if (this.settings.enableAGC) {
      const rms = calculateRMS(pcmData);
      const peak = calculatePeakLevel(pcmData);
      const currentLevel = rms * 2 + peak * 0.5;
      
      if (currentLevel > 0.001) {
        const desiredGain = this.state.targetRMS / currentLevel;
        const agcGain = 1.0 + (desiredGain - 1.0) * this.settings.agcStrength;
        
        this.state.currentGain = this.state.currentGain * this.state.gainSmoothingFactor + 
                                  agcGain * (1 - this.state.gainSmoothingFactor);
        this.state.currentGain = clamp(this.state.currentGain, 0.1, 10);
        effectiveGain = this.settings.gain * this.state.currentGain;
      }
    } else {
      this.state.currentGain = 1.0;
    }
    
    const spectrumData = generateSpectrum(
      pcmData, false, effectiveGain, this.settings.naturalWeighting, 
      this.audioCtx.sampleRate, this.settings.fftSize, 
      this.settings.lowEndBoost, this.settings.noiseGate
    );

    let spectrum;

    if (this.state.isReassigned) {
      const auxSpectrum = generateSpectrum(
        pcmData, true, effectiveGain, this.settings.naturalWeighting, 
        this.audioCtx.sampleRate, this.settings.fftSize, 
        this.settings.lowEndBoost, this.settings.noiseGate
      );
      
      const magnitudes = new Float32Array(spectrumData.length);
      for (let i = 0; i < spectrumData.length; i++) {
        magnitudes[i] = spectrumData[i].magnitude;
      }
      
      spectrum = calcReassignedSpectrum(
        magnitudes,
        calcReassignedFreqs(spectrumData, auxSpectrum, this.state.oscilloscopeBuffer.length, this.audioCtx.sampleRate),
        this.canvas.width,
        this.settings.frequencyScale
      );
    } else {
      spectrum = new Float32Array(spectrumData.length);
      for (let i = 0; i < spectrumData.length; i++) {
        spectrum[i] = spectrumData[i].magnitude;
      }
    }

    // Apply smoothing
    if (this.settings.smoothing > 0 && this.state.prevSpectrum.length === spectrum.length) {
      const smoothFactor = this.settings.smoothing;
      const invSmoothFactor = 1 - smoothFactor;
      
      for (let i = 0; i < spectrum.length; i++) {
        const smoothed = spectrum[i] * invSmoothFactor + 
                        this.state.prevSpectrum[i] * smoothFactor;
        spectrum[i] = Math.max(spectrum[i], smoothed * 0.95);
      }
    }
    this.state.prevSpectrum = new Float32Array(spectrum);

    // Render immediately (timing-controlled in captureLoop)
    this.renderSpectrumFrame({
      spectrum: spectrum,
      isReassigned: this.state.isReassigned
    });
  }

  renderSpectrumFrame(frame) {
    const activeColormap = this.previewColormap || this.settings.colormap;

    if (frame.isReassigned) {
      printSpectrogram(
        this.auxCtx, this.auxCanvas, this.auxCtx2, this.auxCanvas2,  // <-- Added both canvases
        frame.spectrum, undefined,
        activeColormap, this.settings.dbRange, true
      );
    } else {
      const linearBinData = generateSpectrumBarData(
        this.canvas.height, this.settings.fftSize, 
        this.audioCtx.sampleRate, this.settings.frequencyScale
      );
      printSpectrogram(
        this.auxCtx, this.auxCanvas, this.auxCtx2, this.auxCanvas2,  // <-- Added both canvases
        frame.spectrum, linearBinData,
        activeColormap, this.settings.dbRange, true
      );
    }

    this.state.staticSpectrogramIdx = idxWrapOver(
      this.state.staticSpectrogramIdx + 1, 
      this.auxCanvas.width
    );
  }

  animate() {
    if (this.canvas.width > 0 && this.canvas.height > 0) {
      this.ctx.fillStyle = '#000';
      this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

      if (this.auxCanvas.width > 0 && this.auxCanvas.height > 0) {
        // 1:1 copy, no scaling
        this.ctx.drawImage(this.auxCanvas, 0, 0);  // <-- Removed scaling parameters
      }
    }

    if (this.isMac) {
      setTimeout(() => this.animate(), 16);
    } else {
      requestAnimationFrame(() => this.animate());
    }
  }
}

document.addEventListener('DOMContentLoaded', () => {
  new AudioSpectrogram();
});