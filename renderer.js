const { getLoopbackAudioMediaStream } = require('electron-audio-loopback');

function map(x, min, max, targetMin, targetMax) {
  return (x - min) / (max - min) * (targetMax - targetMin) + targetMin;
}

function clamp(x, min, max) {
  return Math.min(Math.max(x, min), max);
}

function idxWrapOver(x, length) {
  return (x % length + length) % length;
}

function hertzToFFTBin(x, y = 'round', bufferSize = 4096, sampleRate = 44100) {
  const bin = x * bufferSize / sampleRate;
  let func = y;
  if (!['floor','ceil','trunc'].includes(func))
    func = 'round';
  return Math[func](bin);
}

function fftBinToHertz(x, bufferSize = 4096, sampleRate = 44100) {
  return x * sampleRate / bufferSize;
}

// Consolidated ascale function with proper dbRange parameter
function ascale(x, dbRange = 90) {
  const db = 20 * Math.log10(Math.max(x, 1e-10)); // Prevent log(0)
  return clamp(map(db, -dbRange, 0, 0, 1), 0, 1);
}

// Improved low-end boost with better frequency targeting
function applyLowEndBoost(complexFFT, sampleRate, fftSize, boostAmount) {
  // More precise frequency bands
  for (let i = 0; i < complexFFT.length; i++) {
    const freq = fftBinToHertz(i, fftSize, sampleRate);
    
    if (freq < 20) continue; // Skip sub-bass noise
    
    let boostFactor = 1.0;
    
    // Sub-bass: 20-60 Hz - Moderate boost
    if (freq >= 20 && freq < 60) {
      const t = (freq - 20) / 40;
      boostFactor = 1.0 + (boostAmount - 1.0) * 1.3 * (1 - t * 0.3);
    }
    // Bass: 60-250 Hz - Strong boost
    else if (freq >= 60 && freq < 250) {
      boostFactor = 1.0 + (boostAmount - 1.0) * 1.5;
    }
    // Low mids: 250-500 Hz - Gradual reduction
    else if (freq >= 250 && freq < 500) {
      const t = (freq - 250) / 250;
      boostFactor = 1.0 + (boostAmount - 1.0) * 0.8 * (1 - t);
    }
    
    complexFFT[i].magnitude *= boostFactor;
  }
  
  return complexFFT;
}

// Noise gate to reduce low-level clutter
function applyNoiseGate(complexFFT, thresholdDB = -80) {
  const thresholdLinear = Math.pow(10, thresholdDB / 20);
  
  for (let i = 0; i < complexFFT.length; i++) {
    if (complexFFT[i].magnitude < thresholdLinear) {
      complexFFT[i].magnitude = 0;
    }
  }
  
  return complexFFT;
}

// Pink noise weighting with better low-frequency handling
function applyPinkNoiseWeighting(freq, sampleRate) {
  if (freq < 20) return 0; // Eliminate sub-20Hz
  
  // Pink noise: -3dB per octave from reference
  const referenceFreq = 1000; // 1kHz reference
  const octaves = Math.log2(freq / referenceFreq);
  const dbAdjustment = octaves * 3;
  
  return Math.pow(10, dbAdjustment / 20);
}

function fscale(x) {
  return Math.log2(x);
}

function logScale(x, scale = 1.0) {
  if (scale < 0.5) {
    const linearAmount = 1 - (scale * 2);
    const logAmount = scale * 2;
    return map(x, 20, 20000, 0, 1) * linearAmount + map(fscale(x), fscale(20), fscale(20000), 0, 1) * logAmount;
  } else {
    return map(fscale(x), fscale(20), fscale(20000), 0, 1);
  }
}

function calcComplexFFT(input) {
  let fft = input.map(x => x);
  let fft2 = input.map(x => x);
  transform(fft, fft2);
  return input.map((_, i, arr) => {
    return {
      re: fft[i]/(arr.length/2),
      im: fft2[i]/(arr.length/2),
      magnitude: Math.hypot(fft[i], fft2[i])/(arr.length/2),
      phase: Math.atan2(fft2[i], fft[i])
    };
  });
}

function transform(real, imag) {
  const n = real.length;
  if (n != imag.length)
    throw "Mismatched lengths";
  if (n <= 0)
    return;
  else if ((2 ** Math.trunc(Math.log2(n))) === n)
    transformRadix2(real, imag);
  else
    transformBluestein(real, imag);
}

function transformRadix2(real, imag) {
  const n = real.length;
  if (n != imag.length)
    throw "Mismatched lengths";
  if (n <= 1)
    return;
  const logN = Math.log2(n);
  if ((2 ** Math.trunc(logN)) !== n)
    throw "Length is not a power of 2";
  
  let cosTable = new Array(n / 2);
  let sinTable = new Array(n / 2);
  for (let i = 0; i < n / 2; i++) {
    cosTable[i] = Math.cos(2 * Math.PI * i / n);
    sinTable[i] = Math.sin(2 * Math.PI * i / n);
  }
  
  for (let i = 0; i < n; i++) {
    let j = reverseBits(i, logN);
    if (j > i) {
      let temp = real[i];
      real[i] = real[j];
      real[j] = temp;
      temp = imag[i];
      imag[i] = imag[j];
      imag[j] = temp;
    }
  }
  
  for (let size = 2; size <= n; size *= 2) {
    let halfsize = size / 2;
    let tablestep = n / size;
    for (let i = 0; i < n; i += size) {
      for (let j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
        const l = j + halfsize;
        const tpre = real[l] * cosTable[k] + imag[l] * sinTable[k];
        const tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
        real[l] = real[j] - tpre;
        imag[l] = imag[j] - tpim;
        real[j] += tpre;
        imag[j] += tpim;
      }
    }
  }

  function reverseBits(x, bits) {
    let y = 0;
    for (let i = 0; i < bits; i++) {
      y = (y << 1) | (x & 1);
      x >>>= 1;
    }
    return y;
  }
}

function transformBluestein(real, imag) {
  const n = real.length;
  if (n != imag.length)
    throw "Mismatched lengths";
  const m = 2 ** Math.trunc(Math.log2(n*2)+1);
  
  let cosTable = new Array(n);
  let sinTable = new Array(n);
  for (let i = 0; i < n; i++) {
    let j = i * i % (n * 2);
    cosTable[i] = Math.cos(Math.PI * j / n);
    sinTable[i] = Math.sin(Math.PI * j / n);
  }
  
  let areal = newArrayOfZeros(m);
  let aimag = newArrayOfZeros(m);
  for (let i = 0; i < n; i++) {
    areal[i] = real[i] * cosTable[i] + imag[i] * sinTable[i];
    aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
  }
  
  let breal = newArrayOfZeros(m);
  let bimag = newArrayOfZeros(m);
  breal[0] = cosTable[0];
  bimag[0] = sinTable[0];
  for (let i = 1; i < n; i++) {
    breal[i] = breal[m - i] = cosTable[i];
    bimag[i] = bimag[m - i] = sinTable[i];
  }
  
  let creal = new Array(m);
  let cimag = new Array(m);
  convolveComplex(areal, aimag, breal, bimag, creal, cimag);
  
  for (let i = 0; i < n; i++) {
    real[i] = creal[i] * cosTable[i] + cimag[i] * sinTable[i];
    imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
  }
}

function convolveComplex(xreal, ximag, yreal, yimag, outreal, outimag) {
  const n = xreal.length;
  xreal = xreal.slice();
  ximag = ximag.slice();
  yreal = yreal.slice();
  yimag = yimag.slice();
  
  transform(xreal, ximag);
  transform(yreal, yimag);
  
  for (let i = 0; i < n; i++) {
    const temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
    ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
    xreal[i] = temp;
  }
  
  inverseTransform(xreal, ximag);
  
  for (let i = 0; i < n; i++) {
    outreal[i] = xreal[i] / n;
    outimag[i] = ximag[i] / n;
  }
}

function inverseTransform(real, imag) {
  transform(imag, real);
}

function newArrayOfZeros(n) {
  let result = new Array(n).fill(0);
  return result;
}

function applyWindow(x) {
  return Math.cos(x*Math.PI/2) ** 2;
}

function generateSpectrumBarData(size = 1920, length = 4096, sampleRate = 48000, frequencyScale = 1.0) {
  const min = 20,
    max = 20000,
    isFlipped = min > max,
    minIdx = hertzToFFTBin(min, isFlipped ? 'ceil' : 'floor', length, sampleRate),
    maxIdx = hertzToFFTBin(max, isFlipped ? 'floor' : 'ceil', length, sampleRate);

  const spectrogramBars = [];
  for (let i = Math.min(minIdx, maxIdx); i <= Math.max(minIdx, maxIdx); i++) {
    const lowerBound = logScale(fftBinToHertz(i-0.5, length, sampleRate), frequencyScale),
      higherBound = logScale(fftBinToHertz(i+0.5, length, sampleRate), frequencyScale),
      lowerVisible = clamp(Math.round(lowerBound * size), 0, size),
      higherVisible = clamp(Math.round(higherBound * size), 0, size);
    
    if (lowerVisible !== higherVisible) {
      spectrogramBars.push({
        lo: i,
        hi: i,
        start: lowerVisible,
        end: higherVisible
      });
    }
    else if (spectrogramBars.length > 0) {
      const lastBin = spectrogramBars[spectrogramBars.length-1];
      lastBin.lo = Math.min(lastBin.lo, i);
      lastBin.hi = Math.max(lastBin.hi, i);
    }
  }
  return spectrogramBars;
}

function calcReassignedFreqs(data, auxData, bufferSize, sampleRate) {
  if (data.length !== auxData.length)
    throw 'Mismatched lengths';
  const freqs = new Array(data.length);
  for (let i = 0; i < freqs.length; i++) {
    const denom = data[i].re ** 2 + data[i].im ** 2,
      correction = denom > 0 ? -(data[i].re * auxData[i].im - auxData[i].re * data[i].im) / denom : 0;
    freqs[i] = fftBinToHertz(i, bufferSize, sampleRate) - correction * ((sampleRate/2)/Math.PI);
  }
  return freqs;
}

function calcReassignedSpectrum(fft, freqs, bands = 1920, frequencyScale = 1.0) {
  const buckets = new Array(bands).fill(0);
  for (let i = 0; i < fft.length; i++) {
    const x = map(logScale(freqs[i], frequencyScale), 0, 1, 0, buckets.length),
      posX = Math.round(x);
    if (posX >= 0 && posX < buckets.length)
      buckets[posX] = isFinite(buckets[posX]) ? Math.max(buckets[posX], isFinite(fft[i]) ? fft[i] : 0) : fft[i];
  }
  return buckets;
}

function renderOscilloscope(data, oscilloscopeIdx) {
  const dataset = [];
  for (let i = 0; i < data.length; i++) {
    dataset[i] = data[idxWrapOver(i+oscilloscopeIdx-data.length, data.length)];
  }
  return dataset;
}

// Frequency to note conversion
function frequencyToNote(freq) {
  if (freq < 16.35) return null; // Below C0
  
  const noteNames = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'];
  const a4 = 440;
  const c0 = a4 * Math.pow(2, -4.75); // C0 frequency
  
  const halfSteps = 12 * Math.log2(freq / c0);
  const octave = Math.floor(halfSteps / 12);
  const noteIndex = Math.round(halfSteps % 12);
  
  const noteName = noteNames[noteIndex % 12];
  const cents = Math.round((halfSteps % 1) * 100);
  
  return {
    note: `${noteName}${octave}`,
    freq: freq.toFixed(2),
    cents: cents !== 0 ? cents : null
  };
}

// Convert Y position to frequency based on current settings
function yPositionToFrequency(yPos, canvasHeight, frequencyScale) {
  // Invert Y (bottom = low freq, top = high freq)
  const normalizedY = 1 - (yPos / canvasHeight);
  
  // Convert from normalized position back to frequency
  // Reverses the logScale function
  const minFreq = 20;
  const maxFreq = 20000;
  
  if (frequencyScale < 0.5) {
    // Handle linear/log blend
    const linearAmount = 1 - (frequencyScale * 2);
    const logAmount = frequencyScale * 2;
    
    // Approximate inverse (this is simplified but works well enough)
    const linearFreq = minFreq + normalizedY * (maxFreq - minFreq);
    const logFreq = minFreq * Math.pow(maxFreq / minFreq, normalizedY);
    
    return linearFreq * linearAmount + logFreq * logAmount;
  } else {
    // Pure logarithmic
    return minFreq * Math.pow(maxFreq / minFreq, normalizedY);
  }
}

// Automatic Gain Control
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
    peak = Math.max(peak, Math.abs(data[i]));
  }
  return peak;
}

// Reordered processing pipeline for better results
function generateSpectrum(data, deltaWindow, gain = 1.0, naturalWeighting = true, sampleRate = 44100, fftSize = 4096, lowEndBoost = 1.0, noiseGateDB = -80) {
  const dataArray = [];
  let norm = 0;
  
  // Apply window function
  for (let i = 0; i < data.length; i++) {
    const x = map(i, 0, data.length-1, -1, 1),
      x1 = map(i-0.5, 0, data.length-1, -1, 1),
      x2 = map(i+0.5, 0, data.length-1, -1, 1),
      w = applyWindow(x),
      w1 = applyWindow(x1),
      w2 = applyWindow(x2),
      sample = isFinite(data[i]) ? data[i] * gain : 0,
      amp = sample*(deltaWindow ? w1-w2 : w);
    dataArray[i] = amp;
    norm += w;
  }
  
  // Calculate FFT
  let complexFFT = calcComplexFFT(dataArray.map(x => x/norm*data.length/Math.SQRT2));
  
  // Apply natural weighting FIRST (before boost)
  if (naturalWeighting) {
    for (let i = 0; i < complexFFT.length; i++) {
      const freq = fftBinToHertz(i, fftSize, sampleRate);
      const weight = applyPinkNoiseWeighting(freq, sampleRate);
      complexFFT[i].magnitude *= weight;
    }
  }
  
  // Apply noise gate to reduce clutter
  if (noiseGateDB > -120) {
    complexFFT = applyNoiseGate(complexFFT, noiseGateDB);
  }
  
  // Apply low-end boost AFTER weighting and gating
  if (lowEndBoost > 1.0) {
    complexFFT = applyLowEndBoost(complexFFT, sampleRate, fftSize, lowEndBoost);
  }
  
  return complexFFT;
}

// Colormap functions
function getColormapColor(value, colormap) {
  value = clamp(value, 0, 1);
  
  switch(colormap) {
    case 'flame':
      if (value < 0.25) {
        const t = value / 0.25;
        return `rgb(${Math.round(100 * t)}, 0, ${Math.round(150 * t)})`;
      } else if (value < 0.5) {
        const t = (value - 0.25) / 0.25;
        return `rgb(${Math.round(100 + 155 * t)}, 0, ${Math.round(150 - 100 * t)})`;
      } else if (value < 0.75) {
        const t = (value - 0.5) / 0.25;
        return `rgb(255, ${Math.round(100 * t)}, 0)`;
      } else {
        const t = (value - 0.75) / 0.25;
        return `rgb(255, ${Math.round(100 + 155 * t)}, 0)`;
      }
    
    case 'cool':
      if (value < 0.2) {
        const t = value / 0.2;
        return `rgb(0, 0, ${Math.round(255 * t)})`;
      } else if (value < 0.4) {
        const t = (value - 0.2) / 0.2;
        return `rgb(0, ${Math.round(255 * t)}, 255)`;
      } else if (value < 0.6) {
        const t = (value - 0.4) / 0.2;
        return `rgb(0, 255, ${Math.round(255 - 255 * t)})`;
      } else if (value < 0.8) {
        const t = (value - 0.6) / 0.2;
        return `rgb(${Math.round(255 * t)}, 255, 0)`;
      } else {
        const t = (value - 0.8) / 0.2;
        return `rgb(255, ${Math.round(255 - 100 * t)}, 0)`;
      }
    
    case 'twilight':
      if (value < 0.33) {
        const t = value / 0.33;
        return `rgb(${Math.round(100 + 50 * t)}, 0, ${Math.round(200 - 100 * t)})`;
      } else if (value < 0.66) {
        const t = (value - 0.33) / 0.33;
        return `rgb(${Math.round(150 - 50 * t)}, ${Math.round(100 * t)}, ${Math.round(100 + 155 * t)})`;
      } else {
        const t = (value - 0.66) / 0.34;
        return `rgb(${Math.round(100 + 155 * t)}, ${Math.round(100 - 50 * t)}, 255)`;
      }
    
    case 'hot':
      if (value < 0.33) {
        const t = value / 0.33;
        return `rgb(${Math.round(255 * t)}, 0, 0)`;
      } else if (value < 0.66) {
        const t = (value - 0.33) / 0.33;
        return `rgb(255, ${Math.round(255 * t)}, 0)`;
      } else {
        const t = (value - 0.66) / 0.34;
        return `rgb(255, 255, ${Math.round(255 * t)})`;
      }
    
    default:
      const gray = Math.round(255 * value);
      return `rgb(${gray}, ${gray}, ${gray})`;
  }
}

function printSpectrogram(auxCtx, auxCanvas, data, linearBinData, colormap, dbRange, staticSpectrogramIdx, isStatic) {
  const length = auxCanvas.height,
    isLinearBins = linearBinData !== undefined,
    actualLength = isLinearBins ? linearBinData.length : data.length;
  
  for (let i = 0; i < actualLength; i++) {
    let mag = 0;
    if (isLinearBins) {
      for (let idx = linearBinData[i].lo; idx <= linearBinData[i].hi; idx++) {
        const binIdx = idxWrapOver(idx, data.length);
        mag = Math.max(mag, data[binIdx]);
      }
    }
    
    const start = isLinearBins ? isNaN(linearBinData[i].start) ? 0 : clamp(linearBinData[i].start, 0, auxCanvas.height) : Math.trunc(i/data.length*length),
      end = isLinearBins ? isNaN(linearBinData[i].end) ? 0 : clamp(linearBinData[i].end, 0, auxCanvas.height) : Math.trunc((i+1)/data.length*length),
      delta = end-start,
      amp = ascale(isLinearBins ? mag : data[i], dbRange);
    
    auxCtx.fillStyle = getColormapColor(isFinite(amp) ? amp : 0, colormap);
    auxCtx.fillRect(isStatic ? staticSpectrogramIdx+1 : auxCanvas.width, auxCanvas.height-start, -1, -delta);
  }
  
  if (auxCanvas.width > 0 && auxCanvas.height > 0)
    auxCtx.drawImage(auxCanvas, -1, 0);
}

class AudioSpectrogram {
  constructor() {
    this.canvas = document.getElementById('canvas');
    this.ctx = this.canvas.getContext('2d');
    this.ctx.imageSmoothingEnabled = false;
    this.statusEl = document.getElementById('status');
    this.settingsPanel = document.getElementById('settingsPanel');
    this.settingsContent = document.getElementById('settingsContent');

    this.audioCtx = new (window.AudioContext || window.webkitAudioContext)();
    this.analyser = this.audioCtx.createAnalyser();
    this.analyser.fftSize = 8192;

    // Default settings with noise gate
    const defaultSettings = {
      fftSize: 4096,
      hopSize: 8192,
      useReassignment: true,
      colormap: 'flame',
      dbRange: 120,
      gain: 4.6,
      naturalWeighting: true,
      frequencyScale: 1.0,
      lowEndBoost: 5.0,
      smoothing: 0.5,
      noiseGate: -80,
      alwaysOnTop: true,
      enableAGC: true,
      agcStrength: 0.7,
      sampleRate: this.audioCtx.sampleRate
    };

    this.settings = { ...defaultSettings };

    this.state = {
      isRunning: false,
      spectrum: [],
      prevSpectrum: [],
      oscilloscopeBuffer: new Float32Array(this.settings.fftSize),
      oscilloscopeIdx: 0,
      sampleCounter: 0,
      staticSpectrogramIdx: 0,
      isReassigned: false,
      targetRMS: 0.1,           // Target RMS level
      currentGain: 1.0,         // Current AGC gain multiplier
      gainSmoothingFactor: 0.95 // Smoothing for gain changes
    };

    this.auxCanvas = new OffscreenCanvas(this.canvas.width, this.canvas.height);
    this.auxCtx = this.auxCanvas.getContext('2d');
    this.auxCtx.imageSmoothingEnabled = false;

    this.setupUI();
    this.resizeCanvas();
    this.animate();
    
    this.initSettings(defaultSettings);
  }

  async initSettings(defaultSettings) {
    await this.loadSettings(defaultSettings);
    this.updateUIFromSettings();
    setTimeout(() => this.start(), 500);
  }

  updateUIFromSettings() {
    const updateElement = (id, value, textId, formatter) => {
      const el = document.getElementById(id);
      const textEl = document.getElementById(textId);
      if (el) el.value = value;
      if (textEl) textEl.textContent = formatter ? formatter(value) : value;
    };

    updateElement('colormapSelect', this.settings.colormap, 'colormapValue', 
      v => v.charAt(0).toUpperCase() + v.slice(1));
    updateElement('dbRange', this.settings.dbRange, 'dbRangeValue');
    updateElement('gain', this.settings.gain, 'gainValue', v => v.toFixed(1));
    updateElement('freqScale', this.settings.frequencyScale, 'freqScaleValue', v => v.toFixed(1));
    updateElement('lowEnd', this.settings.lowEndBoost, 'lowEndValue', v => v.toFixed(1) + 'x');
    updateElement('smoothing', this.settings.smoothing, 'smoothingValue', v => v.toFixed(2));
    updateElement('noiseGate', this.settings.noiseGate, 'noiseGateValue', v => v + ' dB');

    const reassignedCheck = document.getElementById('reassignedCheck');
    const naturalCheck = document.getElementById('naturalCheck');
    const alwaysOnTopCheck = document.getElementById('alwaysOnTopCheck');
    const agcCheck = document.getElementById('agcCheck');
    if (reassignedCheck) reassignedCheck.checked = this.settings.useReassignment;
    if (naturalCheck) naturalCheck.checked = this.settings.naturalWeighting;
    if (alwaysOnTopCheck) alwaysOnTopCheck.checked = this.settings.alwaysOnTop;
    if (agcCheck) agcCheck.checked = this.settings.enableAGC;

    updateElement('agcStrength', this.settings.agcStrength, 'agcStrengthValue', v => v.toFixed(2));
  }

  async loadSettings(defaultSettings) {
    const { ipcRenderer } = require('electron');
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
    const { ipcRenderer } = require('electron');
    try {
      const settingsToSave = {
        fftSize: this.settings.fftSize,
        hopSize: this.settings.hopSize,
        useReassignment: this.settings.useReassignment,
        colormap: this.settings.colormap,
        dbRange: this.settings.dbRange,
        gain: this.settings.gain,
        naturalWeighting: this.settings.naturalWeighting,
        frequencyScale: this.settings.frequencyScale,
        lowEndBoost: this.settings.lowEndBoost,
        smoothing: this.settings.smoothing,
        noiseGate: this.settings.noiseGate
      };
      
      const result = await ipcRenderer.invoke('save-settings', settingsToSave);
      if (result.success) {
        console.log('Settings saved successfully');
      } else {
        console.error('Failed to save settings:', result.error);
      }
    } catch (err) {
      console.error('Error saving settings:', err);
    }
  }

  setupUI() {
    const exitBtn = document.getElementById('exitBtn');
    const settingsBtn = document.getElementById('settingsBtn');
    
    if (exitBtn) {
      exitBtn.addEventListener('click', () => {
        window.close();
      });
    }

    if (settingsBtn) {
      settingsBtn.addEventListener('click', () => {
        this.toggleSettings();
      });
    }

    if (this.settingsPanel) {
      this.settingsPanel.addEventListener('click', (e) => {
        if (e.target === this.settingsPanel) {
          this.toggleSettings();
        }
      });
    }

    window.addEventListener('resize', () => this.resizeCanvas());
    this.createSettingsUI();

    // Setup frequency hover display
    this.freqHover = document.getElementById('freqHover');
    this.freqLine = document.getElementById('freqLine');
    this.freqLabel = document.getElementById('freqLabel');
    
    const { ipcRenderer } = require('electron');
    
    // Start mouse tracking from main process
    ipcRenderer.send('start-mouse-tracking');
    
    ipcRenderer.on('mouse-position', (event, mousePos) => {
      const container = document.getElementById('container');
      const rect = container.getBoundingClientRect();
      
      const clientX = mousePos.x;
      const clientY = mousePos.y;
      
      // Check if mouse is inside container
      const isInside = clientX >= 0 && clientX <= rect.width &&
                       clientY >= 0 && clientY <= rect.height;
      
      // Check if settings panel is open
      const settingsOpen = this.settingsPanel && this.settingsPanel.classList.contains('open');
      
      // Check if mouse is over button areas
      let isOverButton = false;
      const exitButton = document.getElementById('exitBtn');
      const settingsButton = document.getElementById('settingsBtn');
      
      if (exitButton) {
        const btnRect = exitButton.getBoundingClientRect();
        isOverButton = isOverButton || (
          clientX >= btnRect.left && clientX <= btnRect.right &&
          clientY >= btnRect.top && clientY <= btnRect.bottom
        );
      }
      if (settingsButton) {
        const btnRect = settingsButton.getBoundingClientRect();
        isOverButton = isOverButton || (
          clientX >= btnRect.left && clientX <= btnRect.right &&
          clientY >= btnRect.top && clientY <= btnRect.bottom
        );
      }
      
      if (isInside && !settingsOpen && !isOverButton) {
        if (this.freqHover) this.freqHover.classList.add('active');
        this.updateFrequencyDisplay({ clientX, clientY });
      } else {
        if (this.freqHover) this.freqHover.classList.remove('active');
      }
    });
  }

  toggleSettings() {
    if (this.settingsPanel) {
      this.settingsPanel.classList.toggle('open');
    }
  }

  updateFrequencyDisplay(e) {
    if (!this.freqHover || !this.freqLine || !this.freqLabel) return;
    
    const rect = this.canvas.getBoundingClientRect();
    const y = e.clientY;
    
    // Update line position
    this.freqLine.style.top = y + 'px';
    
    // Calculate frequency from Y position using actual canvas height
    // Use the visual height (rect.height), not the internal canvas height
    const freq = yPositionToFrequency(y, rect.height, this.settings.frequencyScale);
    const noteInfo = frequencyToNote(freq);
    
    if (noteInfo) {
      // Format label text
      let labelText = `<span class="note">${noteInfo.note}</span> ${noteInfo.freq} Hz`;
      if (noteInfo.cents) {
        const centsSign = noteInfo.cents > 0 ? '+' : '';
        labelText += `<span class="cents">(${centsSign}${noteInfo.cents}¢)</span>`;
      }
      this.freqLabel.innerHTML = labelText;
      
      // Position label with padding to keep it visible
      // Estimate label height as ~30px
      const labelHeight = 30;
      const padding = 10;
      
      let labelY = y;
      
      // If too close to top, move it down
      if (y < labelHeight / 2 + padding) {
        labelY = labelHeight / 2 + padding;
      }
      // If too close to bottom, move it up
      else if (y > rect.height - labelHeight / 2 - padding) {
        labelY = rect.height - labelHeight / 2 - padding;
      }
      
      this.freqLabel.style.top = labelY + 'px';
    }
  }

  createSettingsUI() {
    this.settingsContent.innerHTML = `
      <div class="setting-group">
        <button id="resetSettingsBtn" class="reset-btn">Reset to Default</button>
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Colormap</label>
          <span class="setting-value" id="colormapValue">Flame</span>
        </div>
        <select id="colormapSelect">
          <option value="flame">Flame</option>
          <option value="cool">Cool</option>
          <option value="twilight">Twilight</option>
          <option value="hot">Hot</option>
          <option value="grayscale">Grayscale</option>
        </select>
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>dB Range</label>
          <span class="setting-value" id="dbRangeValue">120</span>
        </div>
        <input type="range" id="dbRange" min="30" max="120" value="120">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Gain</label>
          <span class="setting-value" id="gainValue">4.6</span>
        </div>
        <input type="range" id="gain" min="0.1" max="10" step="0.1" value="4.6">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Freq Scale</label>
          <span class="setting-value" id="freqScaleValue">1.0</span>
        </div>
        <input type="range" id="freqScale" min="0" max="1" step="0.1" value="1">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Low End Boost</label>
          <span class="setting-value" id="lowEndValue">5.0x</span>
        </div>
        <input type="range" id="lowEnd" min="1" max="10" step="0.1" value="5.0">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Noise Gate</label>
          <span class="setting-value" id="noiseGateValue">-80 dB</span>
        </div>
        <input type="range" id="noiseGate" min="-120" max="-40" step="5" value="-80">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Smoothing</label>
          <span class="setting-value" id="smoothingValue">0.50</span>
        </div>
        <input type="range" id="smoothing" min="0" max="0.5" step="0.05" value="0.5">
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <input type="checkbox" id="reassignedCheck" checked>
          <label for="reassignedCheck">Enhanced</label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <input type="checkbox" id="naturalCheck" checked>
          <label for="naturalCheck">Natural Audio</label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <input type="checkbox" id="alwaysOnTopCheck" checked>
          <label for="alwaysOnTopCheck">Always On Top</label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <input type="checkbox" id="agcCheck" checked>
          <label for="agcCheck">Auto Gain Control</label>
        </div>
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>AGC Strength</label>
          <span class="setting-value" id="agcStrengthValue">0.70</span>
        </div>
        <input type="range" id="agcStrength" min="0" max="1" step="0.05" value="0.7">
      </div>
    `;

    // Attach event listeners with auto-save
    document.getElementById('colormapSelect').addEventListener('change', (e) => {
      this.settings.colormap = e.target.value;
      document.getElementById('colormapValue').textContent = e.target.value.charAt(0).toUpperCase() + e.target.value.slice(1);
      this.saveSettings();
    });

    document.getElementById('dbRange').addEventListener('input', (e) => {
      this.settings.dbRange = parseFloat(e.target.value);
      document.getElementById('dbRangeValue').textContent = e.target.value;
      this.saveSettings();
    });

    document.getElementById('gain').addEventListener('input', (e) => {
      this.settings.gain = parseFloat(e.target.value);
      document.getElementById('gainValue').textContent = parseFloat(e.target.value).toFixed(1);
      this.saveSettings();
    });

    document.getElementById('freqScale').addEventListener('input', (e) => {
      this.settings.frequencyScale = parseFloat(e.target.value);
      document.getElementById('freqScaleValue').textContent = parseFloat(e.target.value).toFixed(1);
      this.saveSettings();
    });

    document.getElementById('lowEnd').addEventListener('input', (e) => {
      this.settings.lowEndBoost = parseFloat(e.target.value);
      document.getElementById('lowEndValue').textContent = parseFloat(e.target.value).toFixed(1) + 'x';
      this.saveSettings();
    });

    document.getElementById('noiseGate').addEventListener('input', (e) => {
      this.settings.noiseGate = parseFloat(e.target.value);
      document.getElementById('noiseGateValue').textContent = e.target.value + ' dB';
      this.saveSettings();
    });

    document.getElementById('smoothing').addEventListener('input', (e) => {
      this.settings.smoothing = parseFloat(e.target.value);
      document.getElementById('smoothingValue').textContent = parseFloat(e.target.value).toFixed(2);
      this.saveSettings();
    });

    document.getElementById('reassignedCheck').addEventListener('change', (e) => {
      this.settings.useReassignment = e.target.checked;
      this.saveSettings();
    });

    document.getElementById('naturalCheck').addEventListener('change', (e) => {
      this.settings.naturalWeighting = e.target.checked;
      this.saveSettings();
    });

    document.getElementById('alwaysOnTopCheck').addEventListener('change', (e) => {
      this.settings.alwaysOnTop = e.target.checked;
      const { ipcRenderer } = require('electron');
      ipcRenderer.send('set-always-on-top', e.target.checked);
      this.saveSettings();
    });

    document.getElementById('agcCheck').addEventListener('change', (e) => {
      this.settings.enableAGC = e.target.checked;
      if (!e.target.checked) {
        // Reset gain when disabling AGC
        this.state.currentGain = 1.0;
      }
      this.saveSettings();
    });

    document.getElementById('agcStrength').addEventListener('input', (e) => {
      this.settings.agcStrength = parseFloat(e.target.value);
      document.getElementById('agcStrengthValue').textContent = parseFloat(e.target.value).toFixed(2);
      this.saveSettings();
    });

    document.getElementById('resetSettingsBtn').addEventListener('click', () => {
      // Reset to default values
      const defaultSettings = {
        fftSize: 4096,
        hopSize: 8192,
        useReassignment: true,
        colormap: 'flame',
        dbRange: 120,
        gain: 4.6,
        naturalWeighting: true,
        frequencyScale: 1.0,
        lowEndBoost: 5.0,
        smoothing: 0.5,
        noiseGate: -80,
        alwaysOnTop: true,
        enableAGC: true,
        agcStrength: 0.7,
        sampleRate: this.audioCtx.sampleRate
      };
      
      this.settings = { ...defaultSettings };
      this.updateUIFromSettings();
      this.saveSettings();
      
      // Update always on top
      const { ipcRenderer } = require('electron');
      ipcRenderer.send('set-always-on-top', this.settings.alwaysOnTop);
    });
  }

  resizeCanvas() {
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = window.innerWidth * dpr;
    this.canvas.height = window.innerHeight * dpr;
    this.canvas.style.width = window.innerWidth + 'px';
    this.canvas.style.height = window.innerHeight + 'px';
    
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    
    this.auxCanvas.width = this.canvas.width;
    this.auxCanvas.height = this.canvas.height;
    this.auxCtx.imageSmoothingEnabled = false;
    
    this.state.staticSpectrogramIdx = 0;
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

    for (let i = 0; i < timeDomainData.length; i++) {
      this.state.oscilloscopeBuffer[this.state.oscilloscopeIdx] = timeDomainData[i];
      this.state.oscilloscopeIdx = idxWrapOver(this.state.oscilloscopeIdx + 1, this.settings.fftSize);
      this.state.sampleCounter++;

      if (this.state.sampleCounter >= this.settings.hopSize) {
        this.analyzeFrame();
        this.state.sampleCounter = 0;
      }
    }

    requestAnimationFrame(() => this.captureLoop());
  }

  analyzeFrame() {
    this.state.isReassigned = this.settings.useReassignment;
    const pcmData = renderOscilloscope(this.state.oscilloscopeBuffer, this.state.oscilloscopeIdx);

    let effectiveGain = this.settings.gain;
    
    if (this.settings.enableAGC) {
      // Calculate current audio level
      const rms = calculateRMS(pcmData);
      const peak = calculatePeakLevel(pcmData);
      
      // Use RMS for gentle normalization, but prevent clipping with peak
      const currentLevel = rms * 2 + peak * 0.5; // Weighted combination
      
      if (currentLevel > 0.001) { // Only adjust if there's actual signal
        // Calculate desired gain
        const desiredGain = this.state.targetRMS / currentLevel;
        
        // Apply AGC strength (0 = no AGC, 1 = full AGC)
        const agcGain = 1.0 + (desiredGain - 1.0) * this.settings.agcStrength;
        
        // Smooth the gain changes to avoid sudden jumps
        this.state.currentGain = this.state.currentGain * this.state.gainSmoothingFactor + 
                                  agcGain * (1 - this.state.gainSmoothingFactor);
        
        // Clamp gain to reasonable range (0.1x to 10x)
        this.state.currentGain = clamp(this.state.currentGain, 0.1, 10);
        
        // Combine manual gain with AGC gain
        effectiveGain = this.settings.gain * this.state.currentGain;
      }
    } else {
      // Reset to manual gain when AGC is off
      this.state.currentGain = 1.0;
    }
    
    // Pass noiseGate parameter
    const spectrumData = generateSpectrum(
      pcmData, 
      false, 
      effectiveGain,
      this.settings.naturalWeighting, 
      this.audioCtx.sampleRate, 
      this.settings.fftSize, 
      this.settings.lowEndBoost,
      this.settings.noiseGate
    );

    if (this.state.isReassigned) {
      const auxSpectrum = generateSpectrum(
        pcmData, 
        true, 
        effectiveGain,
        this.settings.naturalWeighting, 
        this.audioCtx.sampleRate, 
        this.settings.fftSize, 
        this.settings.lowEndBoost,
        this.settings.noiseGate
      );
      this.state.spectrum = calcReassignedSpectrum(
        spectrumData.map(x => x.magnitude),
        calcReassignedFreqs(spectrumData, auxSpectrum, this.state.oscilloscopeBuffer.length, this.audioCtx.sampleRate),
        this.canvas.width,
        this.settings.frequencyScale
      );
    } else {
      this.state.spectrum = spectrumData.map(x => x.magnitude);
    }

    // Apply smoothing if enabled
    if (this.settings.smoothing > 0 && this.state.prevSpectrum.length === this.state.spectrum.length) {
      for (let i = 0; i < this.state.spectrum.length; i++) {
        this.state.spectrum[i] = this.state.spectrum[i] * (1 - this.settings.smoothing) + this.state.prevSpectrum[i] * this.settings.smoothing;
      }
    }
    this.state.prevSpectrum = [...this.state.spectrum];

    if (this.state.isReassigned) {
      printSpectrogram(
        this.auxCtx,
        this.auxCanvas,
        this.state.spectrum,
        undefined,
        this.settings.colormap,
        this.settings.dbRange,
        this.state.staticSpectrogramIdx,
        false
      );
    } else {
      const linearBinData = generateSpectrumBarData(this.canvas.height, this.settings.fftSize, this.audioCtx.sampleRate, this.settings.frequencyScale);
      printSpectrogram(
        this.auxCtx,
        this.auxCanvas,
        this.state.spectrum,
        linearBinData,
        this.settings.colormap,
        this.settings.dbRange,
        this.state.staticSpectrogramIdx,
        false
      );
    }

    this.state.staticSpectrogramIdx = idxWrapOver(this.state.staticSpectrogramIdx + 1, this.auxCanvas.width);
  }

  animate() {
    this.ctx.fillStyle = '#000';
    this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

    if (this.auxCanvas.width > 0 && this.auxCanvas.height > 0) {
      this.ctx.drawImage(this.auxCanvas, 0, 0, this.canvas.width, this.canvas.height);
    }

    requestAnimationFrame(() => this.animate());
  }
}

document.addEventListener('DOMContentLoaded', () => {
  new AudioSpectrogram();
});