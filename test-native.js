const native = require('./build/Release/spectrogram_native.node');

// Test data
const testData = new Float32Array(4096);
for (let i = 0; i < testData.length; i++) {
    testData[i] = Math.sin(2 * Math.PI * 440 * i / 44100); // 440 Hz sine wave
}

const options = {
    gain: 1.0,
    sampleRate: 44100,
    fftSize: 4096,
    lowEndBoost: 1.0,
    noiseGate: -80,
    naturalWeighting: false,
    deltaWindow: false
};

console.log('Testing native addon...');
const result = native.processSpectrum(testData, options);
console.log('Success! Magnitudes length:', result.magnitudes.length);
console.log('First 10 magnitudes:', Array.from(result.magnitudes.slice(0, 10)));