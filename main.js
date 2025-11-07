const { app, BrowserWindow, ipcMain, Menu } = require('electron');
const { initMain } = require('electron-audio-loopback');
const { autoUpdater } = require('electron-updater');
autoUpdater.logger = require('electron-log');
autoUpdater.logger.transports.file.level = 'info';
const windowStateKeeper = require('electron-window-state');
const path = require('path');
const fs = require('fs');

let mainWindow;
let settingsWindow = null;

autoUpdater.on('error', (err) => {
  console.error('Update error:', err);
  // Send to settings window instead of main window
  settingsWindow?.webContents.send('update-error', err.message);
});

// Path to store settings
const userDataPath = app.getPath('userData');
const settingsPath = path.join(userDataPath, 'settings.json');

// Ableton Live State
const liveStateFile = path.join(app.getPath('userData'), 'live_state.json');

// Initialize the state file if it doesn't exist
function ensureStateFile() {
  if (!fs.existsSync(liveStateFile)) {
    fs.writeFileSync(liveStateFile, JSON.stringify({ state: 'restored' }), 'utf8');
  }
}

// Watch the file for changes
function watchStateFile() {
  fs.watchFile(liveStateFile, { interval: 300 }, () => {
    try {
      const data = JSON.parse(fs.readFileSync(liveStateFile, 'utf8'));
      if (!data || !data.state) return;

      if (data.state === 'minimized') {
        if (!mainWindow.isMinimized()) {
          console.log('Minimizing window...');
          mainWindow.minimize();
        }
      } else if (data.state === 'restored') {
        if (mainWindow.isMinimized()) {
          console.log('Restoring window...');
          mainWindow.restore();
        }
      }
    } catch (err) {
      console.error('Failed to read window state file:', err);
    }
  });
}

// IPC handlers for settings
ipcMain.handle('load-settings', () => {
  try {
    if (fs.existsSync(settingsPath)) {
      const data = fs.readFileSync(settingsPath, 'utf8');
      return JSON.parse(data);
    }
  } catch (err) {
    console.error('Failed to load settings:', err);
  }
  return null; // Return null if file doesn't exist or is corrupt
});

ipcMain.handle('save-settings', (event, settings) => {
  try {
    fs.writeFileSync(settingsPath, JSON.stringify(settings, null, 2), 'utf8');
    return { success: true };
  } catch (err) {
    console.error('Failed to save settings:', err);
    return { success: false, error: err.message };
  }
});

function createSettingsWindow() {
  if (settingsWindow) {
    settingsWindow.focus();
    return;
  }

  let settingsWindowState = windowStateKeeper({
    file: 'settings-window-state.json',
    defaultWidth: 500,
    defaultHeight: 500,
    fullScreen: false
  });

  settingsWindow = new BrowserWindow({
    x: settingsWindowState.x,
    y: settingsWindowState.y,
    width: 500,
    height: 500,
    title: 'Settings',
    parent: mainWindow,
    alwaysOnTop: true,
    autoHideMenuBar: true,
    frame: false,
    fullscreen: false,
    minimizable: false,
    resizable: false,
    icon: path.join(__dirname, 'assets', process.platform === 'win32' ? 'icon.ico' : process.platform === 'darwin' ? 'icon.icns' : 'icon.png'),
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false
    }
  });

  settingsWindowState.manage(settingsWindow);

  settingsWindow.loadFile('settings.html');
  // settingsWindow.webContents.openDevTools({mode:'detach'});

  settingsWindow.on('closed', () => {
    settingsWindow = null;
  });
}

app.on('ready', () => {
  console.log('Initializing electron-audio-loopback in main process...');
  initMain();
  
  let mainWindowState = windowStateKeeper({
    file: 'window-state.json',
    defaultWidth: 1400,
    defaultHeight: 900,
    fullScreen: false
  });

  mainWindow = new BrowserWindow({
    width: mainWindowState.width,
    height: mainWindowState.height,
    x: mainWindowState.x,
    y: mainWindowState.y,
    minWidth: 150,
    minHeight: 150,
    frame: false,
    transparent: true,
    resizable: false,
    icon: path.join(__dirname, 'assets', process.platform === 'win32' ? 'icon.ico' : process.platform === 'darwin' ? 'icon.icns' : 'icon.png'),
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false
    }
  });

  mainWindow.setAlwaysOnTop(true, "floating");
  mainWindow.setVisibleOnAllWorkspaces(true);
  mainWindow.setFullScreenable(false);

  mainWindowState.manage(mainWindow);

  mainWindow.loadFile('index.html');
  // mainWindow.webContents.openDevTools({mode:'detach'});

  mainWindow.webContents.on('did-finish-load', () => {
    ensureStateFile();
    watchStateFile();
    
    // Check for updates after a short delay
    setTimeout(() => {
      console.log('Checking for updates...');
      autoUpdater.checkForUpdatesAndNotify();
    }, 2000);
  });

  const WM_INITMENU = 0x0116;
  const menu = Menu.buildFromTemplate([{label: "🛈 Hold shift and move mouse over the spectrogram to view note and frequency information"}]);
  mainWindow.hookWindowMessage(WM_INITMENU, () => {
    mainWindow.setEnabled(false);
    mainWindow.setEnabled(true);

    menu.popup();
  });

});

// Open settings window
ipcMain.on('open-settings', () => {
  createSettingsWindow();
});

// Get current settings
ipcMain.handle('get-settings', async () => {
  try {
    if (fs.existsSync(settingsPath)) {
      const data = fs.readFileSync(settingsPath, 'utf8');
      return JSON.parse(data);
    }
  } catch (err) {
    console.error('Failed to load settings:', err);
  }
  return null;
});

// Update a single setting
ipcMain.on('update-setting', (event, { key, value }) => {
  try {
    let settings = {};
    if (fs.existsSync(settingsPath)) {
      const data = fs.readFileSync(settingsPath, 'utf8');
      settings = JSON.parse(data);
    }
    
    settings[key] = value;
    fs.writeFileSync(settingsPath, JSON.stringify(settings, null, 2), 'utf8');
    
    // Notify main window of the change
    if (mainWindow && !mainWindow.isDestroyed()) {
      mainWindow.webContents.send('setting-changed', { key, value });
    }
    
    // Handle special cases
    if (key === 'alwaysOnTop' && mainWindow) {
      mainWindow.setAlwaysOnTop(value, "floating");
    }
  } catch (err) {
    console.error('Failed to update setting:', err);
  }
});

// Preview colormap
ipcMain.on('preview-colormap', (event, colormap) => {
  if (mainWindow && !mainWindow.isDestroyed()) {
    mainWindow.webContents.send('preview-colormap', colormap);
  }
});

// Reset settings
ipcMain.on('reset-settings', () => {
  const isMac = process.platform === 'darwin';
  const defaultSettings = {
    fftSize: 4096,
    hopSize: 64,
    scrollSpeed: isMac ? 2.9 : 1.7,
    useReassignment: true,
    colormap: 'inferno',
    dbRange: 47,
    gain: 8.1,
    naturalWeighting: true,
    frequencyScale: 1.0,
    lowEndBoost: 10.0,
    smoothing: 0.20,
    noiseGate: -90,
    alwaysOnTop: true,
    enableAGC: true,
    agcStrength: 1.0,
    brightness: 0.70
  };
  
  try {
    fs.writeFileSync(settingsPath, JSON.stringify(defaultSettings, null, 2), 'utf8');
    
    // Notify both windows
    if (mainWindow && !mainWindow.isDestroyed()) {
      mainWindow.webContents.send('settings-reset', defaultSettings);
    }
    if (settingsWindow && !settingsWindow.isDestroyed()) {
      settingsWindow.webContents.send('settings-updated', defaultSettings);
    }
    
    if (mainWindow) {
      mainWindow.setAlwaysOnTop(defaultSettings.alwaysOnTop, "floating");
    }
  } catch (err) {
    console.error('Failed to reset settings:', err);
  }
});

// Listen for manual update checks from renderer
ipcMain.handle('check-for-updates', async () => {
  try {
    const result = await autoUpdater.checkForUpdates();
    return result; // will include version info
  } catch (err) {
    return { error: err.message };
  }
});

// Send update events to settings window instead of main window
autoUpdater.on('update-available', (info) => {
  settingsWindow?.webContents.send('update-available', info);
});

autoUpdater.on('update-not-available', () => {
  settingsWindow?.webContents.send('update-not-available');
});

autoUpdater.on('download-progress', (progress) => {
  settingsWindow?.webContents.send('update-download-progress', progress);
});

autoUpdater.on('update-downloaded', () => {
  settingsWindow?.webContents.send('update-downloaded');
});

autoUpdater.on('checking-for-update', () => {
  console.log('Checking for update...');
  settingsWindow?.webContents.send('checking-for-update');
});

// Install update when user accepts
ipcMain.handle('install-update', async () => {
  autoUpdater.quitAndInstall();
});

ipcMain.on('set-always-on-top', (event, shouldBeOnTop) => {
  if (mainWindow) {
    mainWindow.setAlwaysOnTop(shouldBeOnTop, "floating");
  }
});

// Get window bounds
ipcMain.handle('get-window-bounds', () => {
  if (mainWindow && !mainWindow.isDestroyed()) {
    return mainWindow.getBounds();
  }
  return { x: 0, y: 0, width: 800, height: 600 };
});

// Handle window resize with position (for dragging edges)
ipcMain.on('resize-window-bounds', (event, bounds) => {
  if (mainWindow && !mainWindow.isDestroyed()) {
    mainWindow.setBounds(bounds);
  }
});

app.on('window-all-closed', () => {
 
  fs.unwatchFile(liveStateFile);
  
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('before-quit', () => {
  console.log('Cleaning up before quit...');
  
  // Stop file watching
  fs.unwatchFile(liveStateFile);
  
  // Close windows gracefully
  if (settingsWindow && !settingsWindow.isDestroyed()) {
    settingsWindow.destroy();
  }
  if (mainWindow && !mainWindow.isDestroyed()) {
    mainWindow.destroy();
  }
});