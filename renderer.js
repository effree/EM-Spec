const { ipcRenderer } = require('electron');
const { getLoopbackAudioMediaStream } = require('electron-audio-loopback');
const nativeBridge = require('./native-bridge');
const { COLORMAPS } = require('./colormaps');

window.onerror = function(message, source, lineno, colno, error) {
  console.error('Window Error:', message, source, lineno, colno, error);
  require('electron-log').error('Window Error:', message, source, lineno, colno, error);
};

const map = (x, min, max, targetMin, targetMax) => 
  (x - min) / (max - min) * (targetMax - targetMin) + targetMin;

const clamp = (x, min, max) => Math.min(Math.max(x, min), max);

const idxWrapOver = (x, length) => (x % length + length) % length;

// Consolidated ascale with proper dbRange
const ascale = (x, dbRange = 90) => {
  const db = 20 * Math.log10(Math.max(x, 1e-10));
  return clamp(map(db, -dbRange, 0, 0, 1), 0, 1);
};

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
      const exponent = 1.5 + (frequencyScale * 0.5); // 1.5 to 2.0
      freqStart = min * Math.pow(max / min, Math.pow(pixelStart, exponent));
      freqEnd = min * Math.pow(max / min, Math.pow(pixelEnd, exponent));
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

function renderOscilloscope(data, oscilloscopeIdx) {
  const len = data.length;
  const dataset = new Float32Array(len);
  
  // More efficient copy
  const startIdx = oscilloscopeIdx >= len ? 0 : oscilloscopeIdx;
  const firstPartLen = len - startIdx;
  
  dataset.set(data.subarray(startIdx, len), 0);
  if (startIdx > 0) {
    dataset.set(data.subarray(0, startIdx), firstPartLen);
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

function getColormapColorRGB(value, colormapName) {
  value = clamp(value, 0, 1);
  
  const colormap = COLORMAPS[colormapName] || COLORMAPS.inferno;
  
  if (!colormap || !colormap.length) {
    return { r: 128, g: 128, b: 128 };
  }
  
  // Find interpolation range
  let i = 0;
  while (i < colormap.length - 1 && value > colormap[i + 1][0]) {
    i++;
  }
  
  if (i === colormap.length - 1) {
    const [, r, g, b] = colormap[i];
    return { r, g, b };
  }
  
  // Interpolate
  const [pos1, r1, g1, b1] = colormap[i];
  const [pos2, r2, g2, b2] = colormap[i + 1];
  const t = (value - pos1) / (pos2 - pos1);
  
  return {
    r: (r1 + (r2 - r1) * t) | 0,
    g: (g1 + (g2 - g1) * t) | 0,
    b: (b1 + (b2 - b1) * t) | 0
  };
}

// ===== Spectrogram Rendering =====

function printSpectrogram(auxCtx, auxCanvas, data, linearBinData, colormap, dbRange, shouldShift, spectrogramState) {
  const height = auxCanvas.height;
  const width = auxCanvas.width;
  
  if (width === 0 || height === 0) return;
  
  // Draw new column at current position
  const xPos = spectrogramState.xPos;
  
  const columnData = auxCtx.createImageData(1, height);
  const columnPixels = columnData.data;
  
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
      const color = getColormapColorRGB(isFinite(amp) ? amp : 0, colormap);
      
      const pixelIndex = (height - pixelY - 1) * 4;
      columnPixels[pixelIndex] = color.r;
      columnPixels[pixelIndex + 1] = color.g;
      columnPixels[pixelIndex + 2] = color.b;
      columnPixels[pixelIndex + 3] = 255;
    }
  } else {
    for (let pixelY = 0; pixelY < height; pixelY++) {
      const binIndex = Math.floor((pixelY / height) * data.length);
      const mag = data[binIndex] || 0;
      
      const amp = ascale(mag, dbRange);
      const color = getColormapColorRGB(isFinite(amp) ? amp : 0, colormap);
      
      const pixelIndex = (height - pixelY - 1) * 4;
      columnPixels[pixelIndex] = color.r;
      columnPixels[pixelIndex + 1] = color.g;
      columnPixels[pixelIndex + 2] = color.b;
      columnPixels[pixelIndex + 3] = 255;
    }
  }
  
  // Draw new column
  auxCtx.putImageData(columnData, xPos, 0);

  // Draw a black line just 1 pixels wide (enough to clear old data)
  const nextXPos = (xPos + 1) % width;
  auxCtx.fillStyle = '#000000';
  auxCtx.fillRect(nextXPos, 0, 1, height);  // Only 1 pixel

  // Advance position
  spectrogramState.xPos = nextXPos;
}

// ===== Main AudioSpectrogram Class =====

class AudioSpectrogram {

  constructor() {
    this.canvas = document.getElementById('canvas');
    this.ctx = this.canvas.getContext('2d', { 
      alpha: false,
      desynchronized: true
    });
    this.ctx.imageSmoothingEnabled = false;
    this.ctx.imageSmoothingQuality = 'low';

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

    this.MIN_SPECTROGRAM_WIDTH = 2048;
    this.MIN_SPECTROGRAM_HEIGHT = 1024;
    this.MAX_SPECTROGRAM_WIDTH = 4096;
    this.MAX_SPECTROGRAM_HEIGHT = 2160;

    const defaultSettings = {
      fftSize: 4096,
      hopSize: 64,
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
      frameSkipCounter: 0.0,
      lastRenderTime: null,
      lastAnalysisTime: null
    };

    this.spectrogramState = { 
      xPos: 0,
      lastCleanup: 0,
      maxAge: 10000  // Keep max 10 seconds of offscreen data (at 1 col/frame, ~600 frames/sec)
    };

    // Call updateSpectrogramResolution to create canvases:
    this.updateSpectrogramResolution();

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
        // Adjust hopSize based on scroll speed for fast scrolling - disabled for now
        // this.settings.hopSize = Math.max(16, Math.floor(64 / Math.max(1, value)));
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

  updateSpectrogramResolution() {
    const dpr = window.devicePixelRatio || 1;
    
    // Try forcing DPR to 1 for sharper pixels
    const effectiveDPR = 1; // Instead of using actual DPR
    
    const targetWidth = window.innerWidth * effectiveDPR;
    const targetHeight = window.innerHeight * effectiveDPR;
    
    const needsResize = !this.auxCanvas || 
                        Math.abs(this.auxCanvas.width - targetWidth) > 0 ||
                        Math.abs(this.auxCanvas.height - targetHeight) > 0;
    
    if (needsResize) {
      let oldCanvas = null;
      let oldXPos = 0;
      let oldWidth = 0;
      let oldHeight = 0;
      
      if (this.auxCanvas && this.auxCanvas.width > 0) {
        oldCanvas = new OffscreenCanvas(this.auxCanvas.width, this.auxCanvas.height);
        const tempCtx = oldCanvas.getContext('2d', { alpha: false });
        tempCtx.drawImage(this.auxCanvas, 0, 0);
        oldXPos = this.spectrogramState.xPos;
        oldWidth = this.auxCanvas.width;
        oldHeight = this.auxCanvas.height;
      }
      
      // Create new canvases
      this.auxCanvas = new OffscreenCanvas(targetWidth, targetHeight);
      this.auxCtx = this.auxCanvas.getContext('2d', { 
        alpha: false,
        willReadFrequently: false,
        desynchronized: true
      });
      this.auxCtx.imageSmoothingEnabled = false;
      this.auxCtx.imageSmoothingQuality = 'low';
      
      // Copy with proper wrapping
      if (oldCanvas && oldWidth > 0 && oldHeight > 0) {
        const copyHeight = Math.min(oldHeight, targetHeight);
        const SEPARATOR_WIDTH = 1;
        
        // Skip just the separator area (3px around xPos)
        const safeStart = (oldXPos + SEPARATOR_WIDTH + 1) % oldWidth;
        const safeCopyWidth = oldWidth - SEPARATOR_WIDTH - 1;
        
        // First section: from safe start to end of old buffer
        const firstSectionWidth = Math.min(safeCopyWidth, oldWidth - safeStart);
        if (firstSectionWidth > 0) {
          this.auxCtx.drawImage(oldCanvas,
            safeStart, 0, firstSectionWidth, copyHeight,
            0, 0, firstSectionWidth, copyHeight
          );
        }
        
        // Second section: wrap around if needed
        const remainingWidth = safeCopyWidth - firstSectionWidth;
        if (remainingWidth > 0) {
          this.auxCtx.drawImage(oldCanvas,
            0, 0, remainingWidth, copyHeight,
            firstSectionWidth, 0, remainingWidth, copyHeight
          );
        }
        
        // Set xPos to continue from where we left off
        this.spectrogramState.xPos = safeCopyWidth % targetWidth;
      } else {
        this.spectrogramState.xPos = 0;
      }
      
      SPECTRUM_BAR_CACHE.clear();
    }
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
    
    // console.log(`FFT size updated to ${this.settings.fftSize}, analyser size: ${analyserSize}`);
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
    
    // Update spectrogram resolution to match new size
    this.updateSpectrogramResolution();
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
        this.state.sampleCounter = 0;
        
        // Always analyze at hop intervals
        const now = performance.now();
        const timeSinceLastAnalysis = this.state.lastAnalysisTime ? (now - this.state.lastAnalysisTime) : 1000;
        
        // Adaptive: only analyze if enough time has passed (prevent overload)
        // At 2048 FFT: min 2ms between analyses
        // At 4096 FFT: min 3ms between analyses  
        // At 8192 FFT: min 5ms between analyses
        const minAnalysisInterval = (this.settings.fftSize / 2048) * 2;
        
        if (timeSinceLastAnalysis >= minAnalysisInterval) {
          this.analyzeFrame();
          this.state.lastAnalysisTime = now;
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
    
    const profile = {};
    const totalStart = performance.now();
    
    // Oscilloscope
    let t = performance.now();
    const pcmData = renderOscilloscope(this.state.oscilloscopeBuffer, this.state.oscilloscopeIdx);
    profile.oscilloscope = performance.now() - t;

    let effectiveGain = this.settings.gain;
    
    // AGC
    t = performance.now();
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
    profile.agc = performance.now() - t;
    
    let spectrum;
    
    // Spectrum processing
    t = performance.now();
    if (this.settings.useReassignment) {
    const options = {
        gain: effectiveGain,
        sampleRate: this.audioCtx.sampleRate,
        fftSize: this.settings.fftSize,
        lowEndBoost: this.settings.lowEndBoost,
        noiseGate: this.settings.noiseGate,
        naturalWeighting: this.settings.naturalWeighting,
        numBands: this.auxCanvas.width,
        frequencyScale: this.settings.frequencyScale
    };
    
    spectrum = nativeBridge.processReassignedSpectrum(pcmData, options);
    } else {
    const options = {
        gain: effectiveGain,
        sampleRate: this.audioCtx.sampleRate,
        fftSize: this.settings.fftSize,
        lowEndBoost: this.settings.lowEndBoost,
        noiseGate: this.settings.noiseGate,
        naturalWeighting: this.settings.naturalWeighting,
        deltaWindow: false
    };
    
    const result = nativeBridge.processSpectrum(pcmData, options);
    spectrum = result.magnitudes;
    }
    profile.processing = performance.now() - t;

    // Smoothing
    t = performance.now();
    if (this.settings.smoothing > 0 && this.state.prevSpectrum.length === spectrum.length) {
    spectrum = nativeBridge.applySmoothing(spectrum, this.state.prevSpectrum, this.settings.smoothing);
    }
    this.state.prevSpectrum = new Float32Array(spectrum);
    profile.smoothing = performance.now() - t;

    // SCROLL SPEED LOGIC: Render multiple columns if scroll speed is high
    t = performance.now();
    
    // Initialize accumulator
    if (!this.renderAccumulator) this.renderAccumulator = 0;
    
    // Add scroll speed to accumulator
    this.renderAccumulator += this.settings.scrollSpeed;
    
    // Should we render this frame?
    let didRender = false;
    if (this.renderAccumulator >= 1.0) {
        this.renderAccumulator = this.renderAccumulator % 1.0;
      // YES - render exactly ONE column
      this.renderSpectrumFrame({
        spectrum: spectrum,
        isReassigned: this.settings.useReassignment
      });
      
      // Subtract 1.0 from accumulator (keep fractional part)
      this.renderAccumulator -= 1.0;
      didRender = true;
    }
    
    profile.render = performance.now() - t;
    profile.didRender = didRender;

    profile.total = performance.now() - totalStart;

    // if (!this.speedDebugCounter) this.speedDebugCounter = 0;
    // this.speedDebugCounter++;
    
    // if (this.speedDebugCounter % 100 === 0) {
    //   console.log(`Speed: ${this.settings.scrollSpeed.toFixed(1)}x | Rendered: ${didRender ? 'YES' : 'NO'} | Accum: ${this.renderAccumulator.toFixed(2)}`);
    // }

  }

  renderSpectrumFrame(frame) {
    const activeColormap = this.previewColormap || this.settings.colormap;

    if (frame.isReassigned) {
      printSpectrogram(
        this.auxCtx, this.auxCanvas,
        frame.spectrum, undefined,
        activeColormap, this.settings.dbRange, true,
        this.spectrogramState
      );
    } else {
      const linearBinData = generateSpectrumBarData(
        this.auxCanvas.height,  // Use current buffer height
        this.settings.fftSize, 
        this.audioCtx.sampleRate, 
        this.settings.frequencyScale
      );
      printSpectrogram(
        this.auxCtx, this.auxCanvas,
        frame.spectrum, linearBinData,
        activeColormap, this.settings.dbRange, true,
        this.spectrogramState
      );
    }

    this.state.staticSpectrogramIdx = idxWrapOver(
      this.state.staticSpectrogramIdx + 1, 
      this.auxCanvas.width
    );
  }

  animate() {
    if (this.canvas.width > 0 && this.canvas.height > 0) {
      this.ctx.imageSmoothingEnabled = false; // Force it every frame

      this.ctx.fillStyle = '#000';
      this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

      if (this.auxCanvas.width > 0 && this.auxCanvas.height > 0) {
        const xPos = this.spectrogramState.xPos;
        const srcWidth = this.auxCanvas.width;
        const srcHeight = this.auxCanvas.height;
        
        const displayWidth = this.canvas.width;
        const displayHeight = this.canvas.height;
        
        this.ctx.imageSmoothingEnabled = false;
        
        if (displayWidth >= srcWidth) {
          // Display wider - tile buffer from RIGHT edge backwards
          const SEPARATOR_WIDTH = 1;
          const usableWidth = srcWidth - SEPARATOR_WIDTH;
          
          // Calculate how many tiles we need
          const tilesNeeded = Math.ceil(displayWidth / usableWidth);
          
          // Start drawing from the RIGHT edge, going LEFT
          for (let r = 0; r < tilesNeeded; r++) {
            // Calculate destination X (from right edge, going left)
            const drawX = displayWidth - (r + 1) * usableWidth;
            
            // Skip the separator area (3px after write position)
            const srcStart = (xPos + SEPARATOR_WIDTH + 1) % srcWidth;
            const copyWidth = Math.min(usableWidth, displayWidth - Math.max(0, drawX));
            
            if (copyWidth <= 0) break;
            
            const actualDrawX = Math.max(0, drawX);
            const skipFromCopy = actualDrawX - drawX; // If drawX is negative
            
            if (srcStart + copyWidth - skipFromCopy <= srcWidth) {
              // No wrap
              this.ctx.drawImage(this.auxCanvas, 
                srcStart + skipFromCopy, 0, copyWidth - skipFromCopy, srcHeight,
                actualDrawX, 0, copyWidth - skipFromCopy, srcHeight
              );
            } else {
              // Wrap around
              const firstPart = srcWidth - (srcStart + skipFromCopy);
              if (firstPart > 0) {
                this.ctx.drawImage(this.auxCanvas,
                  srcStart + skipFromCopy, 0, firstPart, srcHeight,
                  actualDrawX, 0, firstPart, srcHeight
                );
              }
              
              const secondPart = (copyWidth - skipFromCopy) - firstPart;
              if (secondPart > 0) {
                this.ctx.drawImage(this.auxCanvas,
                  0, 0, secondPart, srcHeight,
                  actualDrawX + firstPart, 0, secondPart, srcHeight
                );
              }
            }
          }
        } else {
          // Display narrower - show recent portion
          const visibleWidth = displayWidth;
          let startX = xPos - visibleWidth;
          
          if (startX < 0) {
            const endWidth = Math.abs(startX);
            const endStart = srcWidth + startX;
            
            this.ctx.drawImage(this.auxCanvas,
              endStart, 0, endWidth, srcHeight,
              0, 0, endWidth, srcHeight
            );
            
            const beginWidth = visibleWidth - endWidth;
            if (beginWidth > 0) {
              this.ctx.drawImage(this.auxCanvas,
                0, 0, beginWidth, srcHeight,
                endWidth, 0, beginWidth, srcHeight
              );
            }
          } else {
            this.ctx.drawImage(this.auxCanvas,
              startX, 0, visibleWidth, srcHeight,
              0, 0, visibleWidth, srcHeight
            );
          }
        }
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