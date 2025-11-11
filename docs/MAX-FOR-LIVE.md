# EM-Spec Max for Live Integration

EM-Spec includes a companion Max for Live (M4L) device for seamless integration with Ableton Live.

## What It Does

The EM-Spec M4L device creates a smart connection between Ableton Live and the EM-Spec application:

- **Automatic Window Management** - EM-Spec automatically hides and shows based on Ableton Live's Info View visibility
- **Seamless Workflow** - Keep your screen organized without manually managing windows
- **Stay Focused** - The spectrogram appears when you need it and gets out of the way when you don't

## How It Works

When the EM-Spec M4L device is loaded in your Ableton Live set:

1. **Info View Shown** → EM-Spec window **restores** (becomes visible)
2. **Info View Hidden** → EM-Spec window **minimizes** (gets out of the way)

This synchronization happens automatically in real-time, allowing you to toggle both the Info View and EM-Spec with Ableton's built-in Info View shortcut.

## Installation

### 1. Install EM-Spec Application

First, download and install the EM-Spec desktop application from the [Releases](https://github.com/effree/EM-Spec/releases) page.

### 2. Get the M4L Device

[Download .amxd file](https://maxforlive.com/library/device.php?id=13898)

### 3. Install in Ableton Live

1. Download the `EM-Spec.amxd` file
2. Place it in one of these locations:
   - **User Library**: `~/Music/Ableton/User Library/Presets/Audio Effects/Max Audio Effect/`
   - **Any Ableton Project**: Just drag it into your project and save
3. Restart Ableton Live (if it was running)

### 4. Add to Your Set

1. Open Ableton Live
2. Launch the EM-Spec application
3. Drag the EM-Spec M4L device onto any track in your Live set
4. Toggle Ableton's Info View (default shortcut: `?` or `Shift + ?`) to test the integration
5. Edit the username section of the M4L device and toggle PC or Mac for your system

## Usage Tips

- **Always Running**: Keep the EM-Spec application running in the background while using Ableton
- **One Device**: You only need one instance of the M4L device in your entire Live set
- **Any Track**: The device can be placed on any track - it doesn't process audio
- **Save with Set**: Save the device in your default template so it's always available

## How the Integration Works Technically

The M4L device monitors Ableton Live's Info View state and writes updates to a state file that EM-Spec watches:

- **State File Location**: `%APPDATA%/EM-Spec/live_state.json` (Windows) or `~/Library/Application Support/EM-Spec/live_state.json` (macOS)
- **States**: `"minimized"` or `"restored"`
- **Update Frequency**: Real-time monitoring with minimal CPU usage

## Troubleshooting

### EM-Spec doesn't respond to Info View changes

1. **Check EM-Spec is running**: Make sure the application is launched
2. **Verify M4L device is active**: The device should show "Monitoring On" or a similar status
3. **Check file permissions**: Ensure EM-Spec has write access to its user data folder
4. **Restart both applications**: Close both Ableton and EM-Spec, then restart

### State file not found

The state file is automatically created by EM-Spec on first launch. If you're seeing errors:

1. Launch EM-Spec at least once to create the file
2. Check your user data directory exists
3. Verify no antivirus software is blocking file creation

### EM-Spec minimizes but won't restore

- Try manually toggling the Info View a few times
- Check if EM-Spec is set to "Always On Top" in settings
- Restart the M4L device by removing and re-adding it to your track

## Requirements

- **EM-Spec**: Version 0.3.4 or later
- **Ableton Live**: Version 10 or later with Max for Live
- **Operating System**: Windows 10+ or macOS 10.13+

## Standalone Use

EM-Spec works perfectly fine without the M4L device - you can use it as a standalone spectrogram visualizer for any DAW or audio application. The M4L integration is optional and only adds the automatic window management feature.

---

**Questions or Issues?** Open an issue on the [GitHub Issues](https://github.com/effree/EM-Spec/issues) page.