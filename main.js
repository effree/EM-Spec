const { app, BrowserWindow, ipcMain } = require('electron');
const { initMain } = require('electron-audio-loopback');
const { autoUpdater } = require('electron-updater');
autoUpdater.logger = require('electron-log');
autoUpdater.logger.transports.file.level = 'info';
const windowStateKeeper = require('electron-window-state');
const path = require('path');
const fs = require('fs');

let mainWindow;

autoUpdater.on('error', (err) => {
  console.error('Update error:', err);
  mainWindow?.webContents.send('update-error', err.message);
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
  // mainWindow.webContents.openDevTools();

  mainWindow.webContents.on('did-finish-load', () => {
    ensureStateFile();
    watchStateFile();
    
    // Check for updates after a short delay
    setTimeout(() => {
      console.log('Checking for updates...');
      autoUpdater.checkForUpdatesAndNotify();
    }, 2000);
  });

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

// Send update events to renderer
autoUpdater.on('update-available', (info) => {
  mainWindow?.webContents.send('update-available', info);
});

autoUpdater.on('update-not-available', () => {
  mainWindow?.webContents.send('update-not-available');
});

autoUpdater.on('download-progress', (progress) => {
  mainWindow?.webContents.send('update-download-progress', progress);
});

autoUpdater.on('update-downloaded', () => {
  mainWindow?.webContents.send('update-downloaded');
});

autoUpdater.on('checking-for-update', () => {
  console.log('Checking for update...');
  mainWindow?.webContents.send('checking-for-update');
});

// Install update when user accepts
ipcMain.handle('install-update', async () => {
  autoUpdater.quitAndInstall();
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