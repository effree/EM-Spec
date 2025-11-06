const { ipcRenderer } = require('electron');
const { COLORMAPS } = require('./colormaps');

// ===== Auto Update Functions =====
const updateStatus = document.getElementById('updateStatus');
const updateButton = document.getElementById('updateButton');

updateButton.style.display = 'none';

ipcRenderer.on('update-available', (event, info) => {
  updateStatus.innerText = `Update Available: ${info.version}`;
  updateButton.style.display = 'inline';
});

ipcRenderer.on('update-not-available', () => {
  updateStatus.innerText = 'Already up to date';
  updateButton.style.display = 'none';
});

ipcRenderer.on('update-download-progress', (event, progress) => {
  updateStatus.innerText = `Downloading: ${Math.round(progress.percent)}%`;
});

ipcRenderer.on('update-downloaded', () => {
  updateStatus.innerText = 'Update downloaded! Click to restart.';
  updateButton.style.display = 'inline';
  updateButton.innerText = 'Install Update';
});

updateButton.addEventListener('click', async () => {
  ipcRenderer.invoke('install-update');
});

ipcRenderer.on('update-error', (event, error) => {
  updateStatus.innerText = `Update check failed: ${error}`;
  updateButton.style.display = 'none';
  console.error('Update error:', error);
});

// ===== Colormap System =====

function generateColormapGradient(colormapName) {
  const colormap = COLORMAPS[colormapName];
  if (!colormap) return 'linear-gradient(to right, #000, #fff)';
  
  const stops = colormap.map(([pos, r, g, b]) => {
    return `rgb(${r},${g},${b}) ${pos * 100}%`;
  }).join(', ');
  
  return `linear-gradient(to right, ${stops})`;
}

// ===== Settings Management =====
let currentSettings = null;

async function loadSettings() {
  currentSettings = await ipcRenderer.invoke('get-settings');
  if (currentSettings) {
    updateUIFromSettings();
  }
}

function updateUIFromSettings() {
  if (!currentSettings) return;

  const updateElement = (id, value, textId, formatter) => {
    const el = document.getElementById(id);
    const textEl = document.getElementById(textId);
    if (el) el.value = value;
    if (textEl) textEl.textContent = formatter ? formatter(value) : value;
  };

  // Update FFT size
  const fftSelect = document.getElementById('fftSizeSelect');
  if (fftSelect) {
    fftSelect.value = currentSettings.fftSize;
    document.getElementById('fftSizeValue').textContent = currentSettings.fftSize;
  }

  // Update colormap
  updateElement('colormapSelect', currentSettings.colormap, 'colormapValue', 
    v => v.charAt(0).toUpperCase() + v.slice(1));
  
  // Update sliders
  updateElement('dbRange', currentSettings.dbRange, 'dbRangeValue');
  updateElement('gain', currentSettings.gain, 'gainValue', v => v.toFixed(1));
  updateElement('freqScale', currentSettings.frequencyScale, 'freqScaleValue', v => v.toFixed(1));
  updateElement('lowEnd', currentSettings.lowEndBoost, 'lowEndValue', v => v.toFixed(1) + 'x');
  updateElement('smoothing', currentSettings.smoothing, 'smoothingValue', v => v.toFixed(2));
  updateElement('noiseGate', currentSettings.noiseGate, 'noiseGateValue', v => v + ' dB');
  updateElement('agcStrength', currentSettings.agcStrength, 'agcStrengthValue', v => v.toFixed(2));
  updateElement('scrollSpeed', currentSettings.scrollSpeed, 'scrollSpeedValue', v => v.toFixed(1) + 'x');
  updateElement('brightness', currentSettings.brightness, 'brightnessValue', v => Math.round(v * 100) + '%');

  // Update checkboxes
  const checkboxes = [
    ['reassignedCheck', 'useReassignment'],
    ['naturalCheck', 'naturalWeighting'],
    ['alwaysOnTopCheck', 'alwaysOnTop'],
    ['agcCheck', 'enableAGC']
  ];
  
  checkboxes.forEach(([id, setting]) => {
    const el = document.getElementById(id);
    if (el) el.checked = currentSettings[setting];
  });

  // Update colormap dropdown
  const dropdownSelected = document.getElementById('colormapSelected');
  if (dropdownSelected) {
    const selectedName = dropdownSelected.querySelector('.colormap-name');
    const selectedGradient = dropdownSelected.querySelector('.colormap-gradient');
    if (selectedName) {
      selectedName.textContent = currentSettings.colormap.charAt(0).toUpperCase() + currentSettings.colormap.slice(1);
    }
    if (selectedGradient) {
      selectedGradient.style.background = generateColormapGradient(currentSettings.colormap);
    }
  }

  const dropdownOptions = document.getElementById('colormapOptions');
  if (dropdownOptions) {
    const options = dropdownOptions.querySelectorAll('.dropdown-option');
    options.forEach(option => {
      if (option.dataset.colormap === currentSettings.colormap) {
        option.classList.add('selected');
      } else {
        option.classList.remove('selected');
      }
    });
  }
}

let saveTimeout = null;

function updateSetting(key, value) {
  if (!currentSettings) return;
  currentSettings[key] = value;
  
  // Send to main window immediately for live updates (NO DEBOUNCE)
  ipcRenderer.send('update-setting', { key, value });
  
  // Debounce only the file save
  if (saveTimeout) {
    clearTimeout(saveTimeout);
  }
  saveTimeout = setTimeout(() => {
    // This comment was here but no actual save - that's the bug!
    saveTimeout = null;
  }, 500);
}
// ===== UI Creation =====
function createSettingsUI() {
  const settingsContent = document.getElementById('settingsContent');
  
  // Generate colormap options
  const colormapOptions = Object.keys(COLORMAPS).sort().map(name => {
    const gradient = generateColormapGradient(name);
    const isSelected = currentSettings?.colormap === name ? 'selected' : '';
    return `
      <div class="dropdown-option ${isSelected}" data-colormap="${name}">
        <div class="colormap-gradient" style="background: ${gradient}"></div>
        <span class="colormap-name">${name.charAt(0).toUpperCase() + name.slice(1)}</span>
      </div>
    `;
  }).join('');

  const currentGradient = generateColormapGradient(currentSettings?.colormap || 'inferno');
  const currentName = (currentSettings?.colormap || 'inferno').charAt(0).toUpperCase() + (currentSettings?.colormap || 'inferno').slice(1);

  // Generate FFT size options
  const fftSizes = [1024, 2048, 4096, 8192];
  const fftSizeOptions = fftSizes.map(size => 
    `<option value="${size}" ${currentSettings?.fftSize === size ? 'selected' : ''}>${size}</option>`
  ).join('');

  settingsContent.innerHTML = `
    <div class="fft-colormap-group">
      <div class="setting-group">
        <div class="setting-label">
          <label>FFT Size</label>
          <span class="setting-value" id="fftSizeValue">${currentSettings?.fftSize || 4096}</span>
        </div>
        <select id="fftSizeSelect" class="setting-dropdown">
          ${fftSizeOptions}
        </select>
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Colormap</label>
          <span class="setting-value" id="colormapValue">${currentName}</span>
        </div>
        <div class="custom-dropdown">
          <div class="dropdown-selected" id="colormapSelected">
            <div style="display: flex; align-items: center; gap: 10px;">
              <div class="colormap-gradient" style="background: ${currentGradient}"></div>
              <span class="colormap-name">${currentName}</span>
            </div>
            <div class="dropdown-arrow"></div>
          </div>
          <div class="dropdown-options" id="colormapOptions">
            ${colormapOptions}
          </div>
        </div>
      </div>
    </div>

    <div class="brightness">
      <div class="setting-group">
        <div class="setting-label">
          <label>Brightness</label>
          <span class="setting-value" id="brightnessValue">44%</span>
        </div>
        <input type="range" id="brightness" min="0.44" max="0.99" step="0.01" value="0.44">
      </div>
    </div>

    <div class="slider-group">
      <div class="setting-group">
        <div class="setting-label">
          <label>dB Range</label>
          <span class="setting-value" id="dbRangeValue">120</span>
        </div>
        <input type="range" id="dbRange" min="30" max="120" value="120">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Gain</label>
          <span class="setting-value" id="gainValue">4.6</span>
        </div>
        <input type="range" id="gain" min="0.1" max="10" step="0.1" value="4.6">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Freq Scale</label>
          <span class="setting-value" id="freqScaleValue">1.0</span>
        </div>
        <input type="range" id="freqScale" min="0" max="1" step="0.1" value="1">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Low End Boost</label>
          <span class="setting-value" id="lowEndValue">5.0x</span>
        </div>
        <input type="range" id="lowEnd" min="1" max="10" step="0.1" value="5.0">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Noise Gate</label>
          <span class="setting-value" id="noiseGateValue">-80 dB</span>
        </div>
        <input type="range" id="noiseGate" min="-120" max="-40" step="5" value="-80">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>AGC Strength</label>
          <span class="setting-value" id="agcStrengthValue">0.70</span>
        </div>
        <input type="range" id="agcStrength" min="0" max="1" step="0.05" value="0.7">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Smoothing</label>
          <span class="setting-value" id="smoothingValue">0.50</span>
        </div>
        <input type="range" id="smoothing" min="0" max="0.5" step="0.05" value="0.5">
      </div>

      <div class="setting-group">
        <div class="setting-label">
          <label>Scroll Speed</label>
          <span class="setting-value" id="scrollSpeedValue">1.0x</span>
        </div>
        <input type="range" id="scrollSpeed" min="0.1" max="2.9" step="0.1" value="1.0">
      </div>
    </div>

    <div class="button-group">
      <div class="setting-group">
        <div class="checkbox-group">
          <label for="reassignedCheck">
            <input type="checkbox" id="reassignedCheck" checked>
            <span>Enhanced</span>
          </label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <label for="naturalCheck">
            <input type="checkbox" id="naturalCheck" checked>
            <span>Natural</span>
          </label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <label for="alwaysOnTopCheck">  
            <input type="checkbox" id="alwaysOnTopCheck" checked>
            <span>On Top</span>
          </label>
        </div>
      </div>

      <div class="setting-group">
        <div class="checkbox-group">
          <label for="agcCheck">
            <input type="checkbox" id="agcCheck" checked>
            <span>Auto Gain</span>
          </label>
        </div>
      </div>
    </div>
  `;

  attachSettingsListeners();
}

function setupColormapDropdown() {
  const dropdownSelected = document.getElementById('colormapSelected');
  const dropdownOptions = document.getElementById('colormapOptions');
  
  if (!dropdownSelected || !dropdownOptions) return;

  dropdownSelected.addEventListener('click', (e) => {
    e.stopPropagation();
    dropdownSelected.classList.toggle('open');
    dropdownOptions.classList.toggle('open');
  });

  document.addEventListener('click', (e) => {
    if (!dropdownSelected.contains(e.target) && !dropdownOptions.contains(e.target)) {
      dropdownSelected.classList.remove('open');
      dropdownOptions.classList.remove('open');
    }
  });

  const options = dropdownOptions.querySelectorAll('.dropdown-option');
  options.forEach(option => {
    const colormapValue = option.dataset.colormap;
    
    option.addEventListener('mouseenter', () => {
      ipcRenderer.send('preview-colormap', colormapValue);
    });
    
    option.addEventListener('click', (e) => {
      e.stopPropagation();
      updateSetting('colormap', colormapValue);
      currentSettings.colormap = colormapValue;
      
      const selectedName = dropdownSelected.querySelector('.colormap-name');
      const selectedGradient = dropdownSelected.querySelector('.colormap-gradient');
      if (selectedName) {
        selectedName.textContent = colormapValue.charAt(0).toUpperCase() + colormapValue.slice(1);
      }
      if (selectedGradient) {
        selectedGradient.style.background = generateColormapGradient(colormapValue);
      }
      
      document.getElementById('colormapValue').textContent = 
        colormapValue.charAt(0).toUpperCase() + colormapValue.slice(1);
      
      options.forEach(opt => opt.classList.remove('selected'));
      option.classList.add('selected');
      
      dropdownSelected.classList.remove('open');
      dropdownOptions.classList.remove('open');
    });
  });

  dropdownOptions.addEventListener('mouseleave', () => {
    ipcRenderer.send('preview-colormap', null);
  });
}

function attachSettingsListeners() {
  const addListener = (id, event, handler) => {
    const el = document.getElementById(id);
    if (el) el.addEventListener(event, handler);
  };

  setupColormapDropdown();

  addListener('fftSizeSelect', 'change', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('fftSize', value);
    document.getElementById('fftSizeValue').textContent = value;
  });

  addListener('brightness', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('brightness', value);
    document.getElementById('brightnessValue').textContent = Math.round(value * 100) + '%';
  });

  addListener('scrollSpeed', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('scrollSpeed', value);
    document.getElementById('scrollSpeedValue').textContent = value.toFixed(1) + 'x';
  });

  addListener('dbRange', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('dbRange', value);
    document.getElementById('dbRangeValue').textContent = value;
  });

  addListener('gain', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('gain', value);
    document.getElementById('gainValue').textContent = value.toFixed(1);
  });

  addListener('freqScale', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('frequencyScale', value);
    document.getElementById('freqScaleValue').textContent = value.toFixed(1);
  });

  addListener('lowEnd', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('lowEndBoost', value);
    document.getElementById('lowEndValue').textContent = value.toFixed(1) + 'x';
  });

  addListener('noiseGate', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('noiseGate', value);
    document.getElementById('noiseGateValue').textContent = value + ' dB';
  });

  addListener('smoothing', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('smoothing', value);
    document.getElementById('smoothingValue').textContent = value.toFixed(2);
  });

  addListener('reassignedCheck', 'change', (e) => {
    updateSetting('useReassignment', e.target.checked);
  });

  addListener('naturalCheck', 'change', (e) => {
    updateSetting('naturalWeighting', e.target.checked);
  });

  addListener('alwaysOnTopCheck', 'change', (e) => {
    updateSetting('alwaysOnTop', e.target.checked);
  });

  addListener('agcCheck', 'change', (e) => {
    updateSetting('enableAGC', e.target.checked);
  });

  addListener('agcStrength', 'input', (e) => {
    const value = parseFloat(e.target.value);
    updateSetting('agcStrength', value);
    document.getElementById('agcStrengthValue').textContent = value.toFixed(2);
  });

  addListener('resetSettingsBtn', 'click', () => {
    ipcRenderer.send('reset-settings');
  });
}

// ===== Initialize =====
document.addEventListener('DOMContentLoaded', async () => {
  console.log('DOM Loaded');
  await loadSettings();
  createSettingsUI();
  console.log('Create setting ui called from DOM Loaded');
  updateUIFromSettings();
  console.log('Update settings');
  
  // Listen for settings updates from main window
  ipcRenderer.on('settings-updated', (event, settings) => {
    currentSettings = settings;
    updateUIFromSettings();
  });
});