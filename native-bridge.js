// native-bridge.js
let nativeAddon = null;
let useNative = false;

try {
  nativeAddon = require('./build/Release/spectrogram_native.node');
  useNative = true;
  console.log('✓ Native C++ addon loaded successfully');
} catch (err) {
  console.warn('⚠ Native addon not available, falling back to JavaScript:', err.message);
  useNative = false;
}

module.exports = {
  isNativeAvailable: () => useNative,
  
  processSpectrum: (audioData, options) => {
    if (!useNative || !nativeAddon) {
      throw new Error('Native addon not available');
    }
    return nativeAddon.processSpectrum(audioData, options);
  },
  
  processReassignedSpectrum: (audioData, options) => {
    if (!useNative || !nativeAddon) {
      throw new Error('Native addon not available');
    }
    return nativeAddon.processReassignedSpectrum(audioData, options);
  },
  
  // ADD THIS:
  applySmoothing: (currentSpectrum, previousSpectrum, smoothFactor) => {
    if (!useNative || !nativeAddon) {
      throw new Error('Native addon not available');
    }
    return nativeAddon.applySmoothing(currentSpectrum, previousSpectrum, smoothFactor);
  }
};