// Emmanuel Makinde
// For counting and calculating fraction of cells that are HA/mRuby+

macro "LNGFR_FractionPositive_by_DAPI_ROIs" {
	
	// I believe these are the channels
    haCh = 1;
    dapiCh = 3;
    mrubyCh = 4;

    // DAPI segmentation params
    rbRadius = 50;
    thrLow = 60;
    minSize = 50;
    maxSize = 1e12;
    minCirc = 0.20;
    maxCirc = 1.00;

    // positivity rule parameters which seem to work (I think?)
    ratioThreshHA = 1.30;
    ratioThreshMR = 1.30;

    inputFolder = getDirectory("");
    if (inputFolder=="") exit("No folder chosen.");

    outputFolder = inputFolder + "LNGFR_Output" + File.separator;
    File.makeDirectory(outputFolder);

    outCsv = outputFolder + "fraction_positive_by_DAPI.csv";
    File.saveString("file,dapi_cells,ha_pos,mruby_pos,frac_ha,frac_mruby\n", outCsv);

    setBatchMode(true);
    processDir(inputFolder, outCsv, haCh, dapiCh, mrubyCh,
               rbRadius, thrLow, minSize, maxSize, minCirc, maxCirc,
               ratioThreshHA, ratioThreshMR);
    setBatchMode(false);

    print("DONE. Saved in: " + outputFolder);
}

function processDir(dir, outCsv, haCh, dapiCh, mrubyCh,
                    rbRadius, thrLow, minSize, maxSize, minCirc, maxCirc,
                    ratioHA, ratioMR) {
    items = getFileList(dir);
    
    for (i=0; i<lengthOf(items); i++) {
        name = items[i];

        if (endsWith(name, "/")) { processDir(dir+name, outCsv, haCh, dapiCh, mrubyCh,
                                              rbRadius, thrLow, minSize, maxSize, minCirc, maxCirc,
                                              ratioHA, ratioMR); continue; }

        if (!(endsWith(name,".tif")||endsWith(name,".tiff")||endsWith(name,".TIF")||endsWith(name,".TIFF"))) continue;

        path = dir + name;
        print("Processing: " + path);

        open(path);
        orig = getTitle();
        getDimensions(w,h,channels,slices,frames);

        // --- Build DAPI channel image as CHAN ---
        selectWindow(orig);
        run("Duplicate...", "title=CHAN duplicate channels=" + dapiCh + " slices=1-" + slices + " frames=1-" + frames);

        // Z max + segment + put nuclei ROIs into ROI Manager
        dapiCells = buildDapiRois(rbRadius, thrLow, minSize, maxSize, minCirc, maxCirc);

        // If no cells, write zeros and move on / avoids division by zero
        if (dapiCells == 0) {
            File.append(path + ",0,0,0,0,0\n", outCsv);
            cleanupAll();
            continue;
        }

        // --- Compute global background levels for HA and mRuby ---
        haBg = channelBackground(orig, haCh, slices, frames);
        mrBg = channelBackground(orig, mrubyCh, slices, frames);

        // --- Count positive cells by ROI intensity ratio ---
        haPos = countPositiveCells(orig, haCh, slices, frames, haBg, ratioHA);
        mrPos = countPositiveCells(orig, mrubyCh, slices, frames, mrBg, ratioMR);

        fracHA = haPos / dapiCells;
        fracMR = mrPos / dapiCells;

        File.append(path + "," + dapiCells + "," + haPos + "," + mrPos + "," +
                    d2s(fracHA,4) + "," + d2s(fracMR,4) + "\n", outCsv);

        cleanupAll();
    }
}


// Uses active window "CHAN" (DAPI) -> Zmax -> threshold -> analyze particles -> ROIs
function buildDapiRois(rbRadius, thrLow, minSize, maxSize, minCirc, maxCirc) {

    // Z max projection
    getDimensions(w,h,channels,slices,frames);
    if (slices > 1) { run("Z Project...", "projection=[Max Intensity]"); rename("DAPI_ZMAX"); }
    else { rename("DAPI_ZMAX"); }

    run("Duplicate...", "title=WORK");
    selectWindow("WORK");
	
	// Cleaning up the signal and binarizing!
    run("Subtract Background...", "rolling=" + rbRadius);
    run("Smooth");
    setAutoThreshold("Otsu dark");
    setThreshold(thrLow, 255);
    run("Convert to Mask");
    run("Fill Holes");
    run("Watershed");

    roiManager("Reset");
    run("Clear Results");

    opts = "size=" + minSize + "-" + maxSize +
           " circularity=" + minCirc + "-" + maxCirc +
           " show=Nothing add";
    run("Analyze Particles...", opts);

    n = roiManager("Count");

    close(); // WORK
    // keep DAPI_ZMAX open until cleanupAll
    return n;
}


// crude background function to normalize signal later on with mean intensity
function channelBackground(origTitle, chIdx, slices, frames) {
    selectWindow(origTitle);
    run("Duplicate...", "title=BGCHAN duplicate channels=" + chIdx + " slices=1-" + slices + " frames=1-" + frames);
    getDimensions(w,h,c,s,f);
    if (s > 1) { run("Z Project...", "projection=[Max Intensity]"); rename("BGZ"); } else { rename("BGZ"); }

    // simple background estimate: mean of the whole image (works if most pixels are background)
    getStatistics(area, mean, min, max);
    // close temps
    if (isOpen("BGCHAN")) { selectWindow("BGCHAN"); close(); }
    if (isOpen("BGZ")) { selectWindow("BGZ"); close(); }
    return mean;
}


// For each DAPI ROI, measure mean intensity in channel; positive if mean/bg >= ratioThresh
function countPositiveCells(origTitle, chIdx, slices, frames, bg, ratioThresh) {

    selectWindow(origTitle);
    run("Duplicate...", "title=MEASCHAN duplicate channels=" + chIdx + " slices=1-" + slices + " frames=1-" + frames);
    getDimensions(w,h,c,s,f);
    if (s > 1) { run("Z Project...", "projection=[Max Intensity]"); rename("MEASZ"); } else { rename("MEASZ"); }

    run("Set Measurements...", "mean redirect=None decimal=3");
    run("Clear Results");

    n = roiManager("Count");
    pos = 0;

    for (r=0; r<n; r++) {
        roiManager("Select", r);
        run("Measure");
        row = nResults - 1;
        m = getResult("Mean", row);
        if (bg > 0 && (m / bg) >= ratioThresh) pos++;
    }

    // close temps
    if (isOpen("MEASCHAN")) { selectWindow("MEASCHAN"); close(); }
    if (isOpen("MEASZ")) { selectWindow("MEASZ"); close(); }

    return pos;
}

function cleanupAll() {
    roiManager("Reset");
    run("Clear Results");
    run("Close All");
}

function isOpen(title) {
    titles = getList("window.titles");
    for (j=0; j<lengthOf(titles); j++) if (titles[j] == title) return 1;
    return 0;
}
