const path = require('path');

// Get the correct path whether in dev or production
let nativeAddon;
try {
  // In development, __dirname is the project root
  // In production with asar disabled, __dirname is inside resources/app/
  const addonPath = path.join(__dirname, 'build', 'Release', 'spectrogram_native.node');
  nativeAddon = require(addonPath);
} catch (err) {
  console.error('Failed to load native addon:', err);
  throw err;
}

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