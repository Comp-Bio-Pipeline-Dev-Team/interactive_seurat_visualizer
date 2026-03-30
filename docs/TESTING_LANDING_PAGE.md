# Landing Page Testing Guide

## What's Been Implemented

### ✅ Completed
1. **Landing Page UI** (`ui_landing_page.R`)
   - Clean, centered upload interface
   - File upload for .rds files only
   - Placeholder sections for object preview, metadata types, and Ensembl conversion

2. **Conditional UI Rendering** (`app.R`)
   - App shows landing page when no object is loaded
   - App shows main interface when object is loaded and user clicks "Proceed"
   - "Load New Data" button in navbar to return to landing page

3. **Upload Logic**
   - Removed .h5ad support
   - Only accepts .rds Seurat objects
   - Progress indicators during upload

4. **Server Logic** (`server_landing_page.R`)
   - Object preview with stats (cells, features, assays, reductions)
   - Metadata type inference table (displays but editing not fully functional yet)
   - Ensembl conversion UI (displays but not functional yet)
   - Proceed to App button

### ⚠️ Known Limitations (To Be Fixed)
1. **Metadata type editing** - Table displays but changes aren't applied yet
2. **Ensembl conversion** - UI shows but conversion doesn't work on landing page yet
3. **Species selection** - Not yet connected to Ensembl conversion

## Testing Steps

### 1. Start the Container
```bash
sudo docker build -t seurat-shiny-app .
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

### 2. Access the App
Open browser to: **http://localhost:3838**

### 3. Test Landing Page
**Expected:** You should see:
- ✅ Clean landing page with upload area
- ✅ "Upload Your Seurat Object" heading
- ✅ File input accepting .rds files

**Test:**
- Try uploading a .rds Seurat object (e.g., `test_seurat.rds`)

### 4. Test Object Preview
**Expected after upload:**
- ✅ Object preview box showing:
  - Number of cells
  - Number of features
  - Number of assays
  - Number of reductions
  - List of assay names
  - List of reduction names

### 5. Test Metadata Type Table
**Expected:**
- ✅ Table showing all metadata columns
- ✅ "Inferred Type" column (numeric or categorical)
- ✅ "User Type" column (editable but changes not applied yet)

### 6. Test Ensembl Conversion Section
**Expected:**
- ✅ Section appears with species selection
- ✅ "Convert Ensembl IDs to Symbols" button
- ⚠️ Button doesn't work yet (to be implemented)

### 7. Test Proceed Button
**Expected:**
- ✅ "Proceed to App" button appears
- ✅ Clicking it transitions to main app interface
- ✅ Main app shows all tabs (Visualization, DE, Enrichment, Heatmap)

### 8. Test "Load New Data" Button
**Expected:**
- ✅ "Load New Data" button appears in navbar
- ✅ Clicking it returns to landing page
- ✅ Previous object is cleared

## Common Issues

### Issue: Landing page doesn't appear
**Cause:** App might be starting with an error
**Solution:** Check Docker logs:
```bash
sudo docker logs seurat-app
```

### Issue: Upload fails
**Cause:** File might not be a valid Seurat object or wrong format
**Solution:** 
- Ensure file is .rds format
- Ensure file contains a Seurat object
- Check Docker logs for error message

### Issue: Proceed button doesn't work
**Cause:** Server logic might have an error
**Solution:** Check Docker logs for R errors

## Next Steps After Testing

Based on test results, we need to:
1. ✅ Verify landing page appears
2. ✅ Verify upload works
3. ✅ Verify object preview displays correctly
4. ⚠️ Implement metadata type application logic
5. ⚠️ Move Ensembl conversion to landing page (currently still in main app)
6. ⚠️ Connect species selection to conversion
7. ✅ Verify proceed button transitions to main app
8. ✅ Verify "Load New Data" returns to landing page

## Docker Commands Reference

**View logs:**
```bash
sudo docker logs -f seurat-app
```

**Stop container:**
```bash
sudo docker stop seurat-app
```

**Remove container:**
```bash
sudo docker rm seurat-app
```

**Rebuild and restart:**
```bash
sudo docker stop seurat-app && sudo docker rm seurat-app
sudo docker build -t seurat-shiny-app .
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```
