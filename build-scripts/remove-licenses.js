const fs = require('fs');
const path = require('path');

exports.default = async function(context) {
  const appOutDir = context.appOutDir;
  const platform = context.electronPlatformName;
  
  console.log('🧹 Starting aggressive cleanup...');
  console.log('Platform:', platform);
  console.log('Output dir:', appOutDir);

  let totalSaved = 0;

  // Remove unused locales (can save 50-80MB)
  const localesPath = path.join(appOutDir, 'locales');
  if (fs.existsSync(localesPath)) {
    const locales = fs.readdirSync(localesPath);
    const keepLocales = ['en-US.pak'];
    
    locales.forEach(locale => {
      if (!keepLocales.includes(locale)) {
        const localePath = path.join(localesPath, locale);
        const stats = fs.statSync(localePath);
        fs.unlinkSync(localePath);
        totalSaved += stats.size;
        console.log(`  Removed locale: ${locale}`);
      }
    });
  }

  function getFileSize(filePath) {
    try {
      return fs.statSync(filePath).size;
    } catch {
      return 0;
    }
  }

  function removeMatchingFiles(dir, depth = 0) {
    if (!fs.existsSync(dir) || depth > 8) return;
    
    try {
      const files = fs.readdirSync(dir);
      
      files.forEach(file => {
        const filePath = path.join(dir, file);
        
        try {
          const stat = fs.statSync(filePath);
          
          if (stat.isDirectory()) {
            // Remove entire unnecessary directories
            const dirsToRemove = [
              '.github', 'docs', 'doc', 'test', 'tests', 'testing',
              'examples', 'example', '__tests__', 'coverage', '.nyc_output',
              'benchmark', 'demo', 'samples', '.vscode', '.idea'
            ];
            
            if (dirsToRemove.some(d => file.toLowerCase() === d)) {
              const size = getDirectorySize(filePath);
              fs.rmSync(filePath, { recursive: true, force: true });
              totalSaved += size;
              console.log(`  Removed dir: ${file} (${formatBytes(size)})`);
              return;
            }
            
            removeMatchingFiles(filePath, depth + 1);
          } else {
            const fileName = file.toLowerCase();
            const ext = path.extname(fileName);
            
            // Files to remove
            const shouldRemove = 
              fileName.includes('license') || 
              fileName.includes('licence') ||
              fileName.includes('changelog') || 
              fileName.includes('readme') ||
              fileName.includes('authors') ||
              fileName.includes('contributors') ||
              fileName === 'licenses.chromium.html' ||
              ext === '.md' ||
              ext === '.markdown' ||
              ext === '.txt' ||
              ext === '.map' ||
              ext === '.ts' && !fileName.includes('.d.ts') ||
              ext === '.yml' ||
              ext === '.yaml' ||
              fileName === '.npmignore' ||
              fileName === '.gitignore' ||
              fileName === '.eslintrc' ||
              fileName === '.prettierrc' ||
              fileName === 'makefile' ||
              fileName === 'gruntfile.js' ||
              fileName === 'gulpfile.js';
            
            if (shouldRemove) {
              const size = stat.size;
              fs.unlinkSync(filePath);
              totalSaved += size;
              
              if (size > 500000) { // Log files > 500KB
                console.log(`  Removed large file: ${file} (${formatBytes(size)})`);
              }
            }
          }
        } catch (err) {
          // Silently skip files we can't access
        }
      });
    } catch (err) {
      console.warn(`  Error reading ${dir}`);
    }
  }

  function getDirectorySize(dirPath) {
    let size = 0;
    try {
      const files = fs.readdirSync(dirPath);
      files.forEach(file => {
        const filePath = path.join(dirPath, file);
        const stats = fs.statSync(filePath);
        if (stats.isDirectory()) {
          size += getDirectorySize(filePath);
        } else {
          size += stats.size;
        }
      });
    } catch {}
    return size;
  }

  function formatBytes(bytes) {
    return (bytes / 1024 / 1024).toFixed(2) + ' MB';
  }

  // Clean resources directory
  const resourcesPath = path.join(appOutDir, 'resources');
  if (fs.existsSync(resourcesPath)) {
    console.log('🔍 Cleaning resources...');
    removeMatchingFiles(resourcesPath);
  }

  // Clean main directory
  console.log('🔍 Cleaning main directory...');
  removeMatchingFiles(appOutDir);

  console.log(`✅ Cleanup complete! Saved: ${formatBytes(totalSaved)}`);
};