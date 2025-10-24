# EM-Spec 0.2.5

A real-time audio spectrogram visualizer built with Electron that captures desktop audio without additional audio routing software.  Using the reassignment method this spectrogram can enhance low-end clarity.  Logorithmic frequency scaling to emulate human hearing and linear frequency scaling to use standard values.

An Ableton Live companion Max for Live device is available.  This device can be placed in your main chain within your Abelton Live Projects to hide and show the device when the "information box" is and isn't visible.  This allows you to place the spectrogram over the often unused information box space in the lower left of the appliaction, to appear as a built in device within Ableton Live.

## Features
- Real-time spectrogram visualization
- Adjustable FFT and hop sizes
- Frequency reassignment for enhanced low-end clarity
- Low-end boost and noise gate controls
- Freeze button to stop writing new data to the screen
- 35 colormap options
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

- **FFT Size**: Select the FFT size manually *Default: 4096
- **Hop Size**: Select the hop size manually *Default: 8192
- **Colormap**: Choose from 35 options *Default: Inferno
- **dB Range**: Adjust dynamic range (30-120 dB)
- **Gain**: Overall amplitude boost
- **Freq Scale**: Adjust between linear and logarithmic frequency scaling
- **Low End Boost**: Enhance bass frequencies (1-10x)
- **Noise Gate**: Filter out low-level noise (-120 to -40 dB)
- **Smoothing**: Temporal smoothing of the display
- **Enhanced**: Toggle frequency reassignment method
- **Natural Audio**: Apply pink noise weighting

## Screenshots

- Spectrogram / Colormap / Hover on Corners
![Spectrogram](/screenshots/spectrogram.png?raw=true "Spectrogram")
![With colormap option](/screenshots/spectrogram-2.png?raw=true "With colormap option")
![With colormap option](/screenshots/spectrogram-3.png?raw=true "With colormap option")
![On corners hover](/screenshots/spectrogram-hover.png?raw=true "On corners hover")
- Settings
![Settings](/screenshots/spectrogram-settings.png?raw=true "Settings")
![Settings colormaps](/screenshots/spectrogram-settings-colormap.png?raw=true "Settings displaying colormaps")
- Hover Notes / Frequency
![Notes and Frequency](/screenshots/spectrogram-notes-frequency.png?raw=true "Hover shows notes and frequency")

## Development

Built with:
- Electron
- electron-audio-loopback for system audio capture