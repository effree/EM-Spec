// native-bridge.js
const nativeAddon = require('./build/Release/spectrogram_native.node');

module.exports = {
  isNativeAvailable: () => true,
  
  processSpectrum: (audioData, options) => {
    return nativeAddon.processSpectrum(audioData, options);
  },
  
  processReassignedSpectrum: (audioData, options) => {
    return nativeAddon.processReassignedSpectrum(audioData, options);
  },
  
  applySmoothing: (currentSpectrum, previousSpectrum, smoothFactor) => {
    return nativeAddon.applySmoothing(currentSpectrum, previousSpectrum, smoothFactor);
  }
};