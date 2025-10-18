# EM-Spec 0.1.1 beta

A real-time audio spectrogram visualizer built with Electron that captures desktop audio.

## Features
- Real-time spectrogram visualization
- Frequency reassignment for enhanced clarity
- Low-end boost and noise gate controls
- Multiple colormap options
- Hover to see frequency and musical note information
- Customizable settings that persist between sessions
- Ableton Live companion M4L device to hide when info box is not visible

## Requirements
- Node.js (v16 or higher recommended)
- npm or yarn

## Installation

1. Clone the repository:
```bash
git clone https://github.com/effree/em-spec.git
cd em-spec
```

2. Install dependencies:
```bash
npm install
```

3. Run the application:
```bash
npm start
```

## Usage

- The application will capture your system audio automatically
- Hover over the spectrogram to see frequency and note information
- Click the settings gear icon to adjust visualization parameters
- Drag the window from anywhere to reposition

## Settings

- **Colormap**: Choose from Flame, Cool, Twilight, Hot, or Grayscale
- **dB Range**: Adjust dynamic range (30-120 dB)
- **Gain**: Overall amplitude boost
- **Freq Scale**: Adjust between linear and logarithmic frequency scaling
- **Low End Boost**: Enhance bass frequencies (1-10x)
- **Noise Gate**: Filter out low-level noise (-120 to -40 dB)
- **Smoothing**: Temporal smoothing of the display
- **Enhanced**: Toggle frequency reassignment method
- **Natural Audio**: Apply pink noise weighting

## Development

Built with:
- Electron
- electron-audio-loopback for system audio capture
- FFT implementation with frequency reassignment
