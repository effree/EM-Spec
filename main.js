const { app, BrowserWindow, ipcMain } = require('electron');
const { initMain } = require('electron-audio-loopback');
const windowStateKeeper = require('electron-window-state');
const path = require('path');
const fs = require('fs');

let mainWindow;

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

app.on('ready', () => {
  console.log('Initializing electron-audio-loopback in main process...');
  initMain();
  
  let mainWindowState = windowStateKeeper({
    defaultWidth: 1400,
    defaultHeight: 900
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
    icon: path.join(__dirname, 'assets', process.platform === 'win32' ? 'icon.ico' : process.platform === 'darwin' ? 'icon.icns' : 'icon.png'),
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false,
      enableRemoteModule: true
    }
  });
  mainWindowState.manage(mainWindow);

  mainWindow.setAlwaysOnTop(true, "floating");
  mainWindow.setVisibleOnAllWorkspaces(true);
  mainWindow.setFullScreenable(false);

  mainWindow.loadFile('index.html');
  mainWindow.webContents.openDevTools();

  ensureStateFile();
  watchStateFile();
});

// Mouse position tracking for renderer
let mouseTrackingInterval;

ipcMain.on('start-mouse-tracking', () => {
  if (mouseTrackingInterval) return;
  
  const { screen } = require('electron');
  mouseTrackingInterval = setInterval(() => {
    if (mainWindow && !mainWindow.isDestroyed()) {
      const mousePos = screen.getCursorScreenPoint();
      const windowBounds = mainWindow.getBounds();
      
      mainWindow.webContents.send('mouse-position', {
        x: mousePos.x - windowBounds.x,
        y: mousePos.y - windowBounds.y
      });
    }
  }, 16);
});

ipcMain.on('stop-mouse-tracking', () => {
  if (mouseTrackingInterval) {
    clearInterval(mouseTrackingInterval);
    mouseTrackingInterval = null;
  }
});

ipcMain.on('set-always-on-top', (event, shouldBeOnTop) => {
  if (mainWindow) {
    mainWindow.setAlwaysOnTop(shouldBeOnTop, "floating");
  }
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    fs.unwatchFile(liveStateFile);
    app.quit();
  }
});