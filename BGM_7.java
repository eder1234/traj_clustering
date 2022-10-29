import ij.plugin.filter.PlugInFilter;
import java.awt.Color;
import java.util.*;
import java.io.*;
import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.process.*;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.Analyzer;
import ij.measure.*;

//computer assisted sperm analysis plugin, based on mTrack2r

public class BGM_7 implements PlugInFilter,Measurements  {

    int		nParticles;
    float[][]	ssx;
    float[][]	ssy;
    String directory,filename;
    //minimum sperm size
    float	minSize = 0;

    //maximum sperm size
    float	maxSize = 40;


    ImagePlus	imp;

    //cutOFF parameters inserted for boar sperm analysis
    float cutOffMinVAPforMotile = 0;
    float cutOffMinMotileSTR = 0;
    //minimum length of sperm track in frames
    float 	minTrackLength = 30;
    //show or do not show all paths of tracked sperm
    boolean	bShowPaths = true;
    //a left over from the original program, I used their length calculation area to do my initial sperm velocity parameter analysis for each track (this should not be a choice it is just a hold over from the original)
    boolean	bShowPathLengths = true;
    //maximum velocity that is allowable between frames (this is the search radius for the next sperm location in a track... it will only look w/in this distance)
    float 	maxVelocity = 40;
    //minimum first point to end point velocity to be called motile
    float	minVSL = 10;
    //minimum curvilinear velocity to be termed motile - this is the point to point (all all them) velocity
    float	minVCL = 25;
    float	maxVAPSlow = 30;
    float	maxVAPMedium = 50;
    float	maxVAPRapid = 100;
    //min velocity on the average path to be termed motile - this is the path that has been averaged with a roaming avg
    float	minVAP = 15;

    //low vcl and VAP values used in finding floaters
    //high wob / lin values for floater finding

    //this specifies the video frame rate / second
    int frameRate = 59;

    //this is the ratio of microns per pixels - size standard used in final velocity calculations
    float	microPerPixel = 500;
    //crap from old program
    float 	maxColumns=75;
    //print the xy co-ordinates for all tracks?
    int printXY = 0;
    //print the motion character for all motile sperm?
    int printSperm = 0;
    //print median values for characteristics?
    int printMedian = 0;

    public class particle {
        float	x;
        float	y;
        int	z;
        int	trackNr;
        boolean inTrack=false;
        boolean flag=false;

        public void copy(particle source) {
            this.x=source.x;
            this.y=source.y;
            this.z=source.z;
            this.inTrack=source.inTrack;
            this.flag=source.flag;
        }

        public float distance (particle p) {
            return (float) Math.sqrt(sqr(this.x-p.x) + sqr(this.y-p.y));
        }
    }

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (IJ.versionLessThan("1.17y"))
            return DONE;
        else
            return DOES_8G+NO_CHANGES;
    }

    public void run(ImageProcessor ip) {
        //the stuff below is the box that pops up to ask for pertinant values - why doesn't it remember the values entered????

        GenericDialog gd = new GenericDialog("Sperm Tracker");
        gd.addNumericField("a, Minimum sperm area (pixels^2):", minSize,0);
        gd.addNumericField("b, Maximum sperm area (pixels^2):", maxSize,40);
        gd.addNumericField("c, Minimum track length (frames):", minTrackLength,100);
        gd.addNumericField("d, Maximum sperm velocity between frames (pixels):", maxVelocity,8);
        gd.addNumericField("e, Minimum VSL for motile (um/s):", minVSL,3);
        gd.addNumericField("f, Minimum VAP for motile (um/s):", minVAP,20);
        gd.addNumericField("g, Minimum VCL for motile (um/s):", minVCL,25);
        gd.addNumericField("h, Maximum VAP for SLOW (um/s):", maxVAPSlow,25);
        gd.addNumericField("i, Maximum VAP for MEDIUM (um/s):", maxVAPMedium,25);
        gd.addNumericField("l, Frame Rate (frames per second):", frameRate,97);
        gd.addNumericField("m, Microns per 1000 pixels:", microPerPixel, 1075);
        gd.addNumericField("n, Print xy co-ordinates for all tracked sperm?(1 Yes, 0 No)", printXY, 0);
        gd.addNumericField("o, Print motion characteristics for all motile sperm?(1 Yes, 0 No)", printSperm,0);
        gd.addNumericField("p, Print mean and median values for motion characteristics?(1 Yes, 0 No)", printMedian,0);
        gd.addMessage("                   ---------    ADVANCED PARAMETERS    ---------");
        gd.addNumericField("q, Minimum VAP for PM", cutOffMinVAPforMotile, 15);
        gd.addNumericField("r, Minimum STR for PM (%)", cutOffMinMotileSTR, 75);
        
        gd.showDialog();
        if (gd.wasCanceled())
            return;

        minSize = (float)gd.getNextNumber();
        maxSize = (float)gd.getNextNumber();
        minTrackLength = (float)gd.getNextNumber();
        maxVelocity = (float)gd.getNextNumber();
        minVSL = (float)gd.getNextNumber();
        minVAP = (float)gd.getNextNumber();
        minVCL = (float)gd.getNextNumber();
        maxVAPSlow = (float)gd.getNextNumber();
        maxVAPMedium = (float)gd.getNextNumber();
        frameRate = (int)gd.getNextNumber();
        microPerPixel = (float)gd.getNextNumber();
        printXY = (int)gd.getNextNumber();
        printSperm = (int)gd.getNextNumber();
        printMedian = (int)gd.getNextNumber();
        cutOffMinVAPforMotile = (float)gd.getNextNumber();
        cutOffMinMotileSTR = (float)gd.getNextNumber();


        //below I am trying to convert integer values required for float to decimal for percent calculations

        microPerPixel = microPerPixel / 1000;
        track(imp, minSize, maxSize, maxVelocity);
    }


    public double myAngle(double dX, double dY){
        double aRad = 0;
        double aFinal = 0;
        if(dY > 0)
        {
            //gives us 0->90, the angle is fine
            if(dX > 0)
            {
                aRad = Math.atan((dY) / (dX));
                aFinal = ((aRad * 180) / Math.PI);
            }
            else if(dX < 0)
            //gives us 0->-90, we need to add 180 (90->180)
            {
                aRad = Math.atan((dY) / (dX));
                aFinal =180 + ((aRad * 180) / Math.PI);
            }
            else if(dX == 0)
            {
                aFinal = 90;
            }
        }
        else if(dY < 0)
        {
            //gives us 0->-90, add 360 and we've got it (270->360)
            if(dX > 0)
            {
                aRad = Math.atan((dY) / (dX));
                aFinal = 360 + ((aRad * 180) / Math.PI);
            }
            else if(dX < 0)
            //gives us 0->90, add 180 and we've got it (180->270)
            {
                aRad = Math.atan((dY) / (dX));
                aFinal =180 + ((aRad * 180) / Math.PI);
            }
            else if(dX == 0)
            //should be able to handle -90 becomes 270 in the new system
            {
                aFinal = 270;
            }
        }
        else if(dY==0)
        {
            if(dX > 0)
            {
                aFinal = 0;
            }
            else if(dX < 0)
            {
                aFinal = 180;
            }
        }
        return aFinal;
    }


    //Angle change and Direction Change; this may represent a substitute for beat cross calculations, however, I have not determined its usefullness and it does appear to be highly dependent on frame rate
    public double[] myAngleDelta(double oldA, double storA, double directionVCL) {
        //aChange is element 0
        //directionVCL is element 1
        //0 is clockwise, 1 is counter-clockwise
        double arrayAngle[] = new double[2];
        //use the angle from above to calculate the change from the previous angle (all angles are relative to the axis...)
        double holdA1 = oldA + 180;
        double holdA2 = oldA - 180;

        if(oldA > 180)
        {
            if(holdA2 < storA && oldA > storA)
            {
                arrayAngle[0] = oldA - storA;
                arrayAngle[1] = 0;
            }
            else if(oldA < storA)
            {
                arrayAngle[0] = storA - oldA;
                arrayAngle[1] = 1;
            }
            else if(oldA > storA)
            {
                arrayAngle[0] = storA - oldA + 360;
                arrayAngle[1] = 1;
            }
            else if(oldA == storA)
            {
                arrayAngle[0] = 0;
                arrayAngle[1] = directionVCL;
            }
            else if(holdA2 == storA)
            {
                arrayAngle[0] = 0;
                arrayAngle[1] = directionVCL;
            }
        }
        else if(oldA <= 180)
        {
            if(holdA1 > storA  && oldA < storA)
            {
                arrayAngle[0] = storA - oldA;
                arrayAngle[1] = 1;
            }
            else if(oldA > storA)
            {
                arrayAngle[0] = oldA - storA;
                arrayAngle[1] = 0;
            }
            else if(oldA < storA)
            {
                arrayAngle[0] = oldA - storA + 360;
                arrayAngle[1] = 0;
            }
            else if(oldA == storA)
            {
                arrayAngle[0] = 0;
                arrayAngle[1] = directionVCL;
            }
            else if(holdA1 == storA)
            {
                arrayAngle[0] = 0;
                arrayAngle[1] = directionVCL;
            }
        }
        return arrayAngle;
    }



    public void track(ImagePlus imp, float minSize, float maxSize, float maxVelocity) {
        int nFrames = imp.getStackSize();
        if (nFrames<2) {
            IJ.showMessage("Tracker", "Stack required");
            return;
        }

        ImageStack stack = imp.getStack();
        int options = 0; // set all PA options false
        int measurements = CENTROID;
        // Initialize results table
        ResultsTable rt = new ResultsTable();
        rt.reset();

        // create storage for particle positions
        List[] theParticles = new ArrayList[nFrames];
        int trackCount=0;

        // record particle positions for each frame in an ArrayList
        for (int iFrame=1; iFrame<=nFrames; iFrame++) {
            theParticles[iFrame-1]=new ArrayList();
            rt.reset();
            ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
            pa.analyze(imp, stack.getProcessor(iFrame));
            float[] sxRes = rt.getColumn(ResultsTable.X_CENTROID);
            float[] syRes = rt.getColumn(ResultsTable.Y_CENTROID);
            if (sxRes==null)
                return;

            for (int iPart=0; iPart<sxRes.length; iPart++) {
                particle aParticle = new particle();
                aParticle.x=sxRes[iPart];
                aParticle.y=syRes[iPart];
                aParticle.z=iFrame-1;
                theParticles[iFrame-1].add(aParticle);
            }
            IJ.showProgress((double)iFrame/nFrames);
        }

        // now assemble tracks out of the particle lists
        // Also record to which track a particle belongs in ArrayLists
        List theTracks = new ArrayList();
        for (int i=0; i<=(nFrames-1); i++) {
            IJ.showProgress((double)i/nFrames);
            for (ListIterator j=theParticles[i].listIterator();j.hasNext();) {
                particle aParticle=(particle) j.next();
                if (!aParticle.inTrack) {
                    // This must be the beginning of a new track
                    List aTrack = new ArrayList();
                    trackCount++;
                    aParticle.inTrack=true;
                    aParticle.trackNr=trackCount;
                    aTrack.add(aParticle);
                    // search in next frames for more particles to be added to track
                    boolean searchOn=true;
                    particle oldParticle=new particle();
                    particle tmpParticle=new particle();
                    oldParticle.copy(aParticle);
                    for (int iF=i+1; iF<=(nFrames-1);iF++) {
                        boolean foundOne=false;
                        particle newParticle=new particle();
                        for (ListIterator jF=theParticles[iF].listIterator();jF.hasNext() && searchOn;) {
                            particle testParticle =(particle) jF.next();
                            float distance = testParticle.distance(oldParticle);
                            // record a particle when it is within the search radius, and when it had not yet been claimed by another track
                            if ( (distance < maxVelocity) && !testParticle.inTrack) {
                                // if we had not found a particle before, it is easy
                                if (!foundOne) {
                                    tmpParticle=testParticle;
                                    testParticle.inTrack=true;
                                    testParticle.trackNr=trackCount;
                                    newParticle.copy(testParticle);
                                    foundOne=true;
                                }
                                else {
                                    // if we had one before, we'll take this one if it is closer.  In any case, flag these particles
                                    testParticle.flag=true;
                                    if (distance < newParticle.distance(oldParticle)) {
                                        testParticle.inTrack=true;
                                        testParticle.trackNr=trackCount;
                                        newParticle.copy(testParticle);
                                        tmpParticle.inTrack=false;
                                        tmpParticle.trackNr=0;
                                        tmpParticle=testParticle;
                                    }
                                    else {
                                        newParticle.flag=true;
                                    }
                                }
                            }
                            else if (distance < maxVelocity) {
                                // this particle is already in another track but could have been part of this one
                                // We have a number of choices here:
                                // 1. Sort out to which track this particle really belongs (but how?)
                                // 2. Stop this track
                                // 3. Stop this track, and also delete the remainder of the other one
                                // 4. Stop this track and flag this particle:
                                testParticle.flag=true;
                            }
                        }
                        if (foundOne)
                            aTrack.add(newParticle);
                        else
                            searchOn=false;
                        oldParticle.copy(newParticle);
                    }
                    theTracks.add(aTrack);
                }
            }
        }


        // this is my movement assessment bit that runs for each sperm
        //this is the number of points to be included in the VAP roaming average
        int vAPpoints = frameRate / 6;

        if(vAPpoints % 2 > 0){
            vAPpoints ++;
        }

        //arrays to hold on to the total movement for a given male
        double[] vCLS = new double[trackCount];
        double[] vAPS = new double[trackCount];
        double[] vSLS = new double[trackCount];
        double[] aLHS = new double[trackCount];
        double[] maxALH = new double[trackCount];
        double[] sumVAPdistsOrigin = new double[trackCount];
        ArrayList indexSlow = new ArrayList();
        ArrayList indexMedium = new ArrayList();
        ArrayList indexRapid = new ArrayList();

        //angle and beat calculations
        double dY = 0;
        double dX = 0;
        double oldA = 0;
        double storA = 0;
        double aChange = 0;
        double directionVCL = 0;
        double oldDirectionVCL = 0;
        double[] dirChangesVCL = new double [trackCount];
        double storTotalRotation = 0;
        double[] framesBeat = new double [trackCount];
        //used for myangleDelta
        double[] arrayAngle = {0, 0};

        //variables used in finding maximum movement from origin
        double firstVAPx = 0;
        double firstVAPy = 0;
        double holdVAPdistance = 0;
        double storMaxVAPVSL = 0;

        //this keeps track of frames on the VCL path
        int[] frames = new int[trackCount];

        //this tracks for each sperm how many times the VCL path crosses VAP, and how many times a direction change in the VCL path occurs
        double[] beatCross = new double[trackCount];
        double vAPxG = 0;
        double vAPyG = 0;
        double vAPxGlast = 0;
        double vAPyGlast = 0;
        //array to hold the VAP points before roaming average
        double[] xCo = new double[vAPpoints];
        double[] yCo = new double[vAPpoints];

        //these two keep track of the number of times movement between frames is zero or low velocity for each sperm
        double[] numVAPzero = new double[trackCount];
        double[] numLowVAP = new double[trackCount];
        //this carries the zero and low vap instances
        double holdZeroVAPs = 0;
        double holdLowVAPs = 0;


        //points used in VAP distance calculations
        double vAPx = 0;
        double vAPy = 0;
        double vAPxOld = 0;
        double vAPyOld = 0;
        double maxALHold = 0;
        //used in vap calculations, holds total of points for average
        double holdX = 0;
        double holdY = 0;

        //for average distance from origin over whole path
        double holdTotalVAPdistance = 0;
        double tempALH = 0;
        double distFirstVAPcurrent = 0;

        //holds values for percent motility calc
        int motileSperm = 0;
        int motileSLOW = 0;
        int motileMEDIUM = 0;
        int motileRAPID = 0;
        int motileULTRARAPID = 0;
        int motilePM = 0;
        int hyperSPZ = 0;
        double totalSperm = 0;

        //for calculations in VCL distance
        double holdDistance = 0;

        double[] framesVAP = new double[trackCount];
        double holdVAPframe = 0;

        //strings to print out all of the data gathered, point by point, also variables for holding points for printing
        String xyPts = " ";
        String xyVAPpts = " ";
        String firstVAPpts = " ";
        double holdMaxX = 0;
        double holdMaxY = 0;

        //initialize variables
        int frame = 0;
        double x1, y1, x2, y2;
        int trackNr=0;
        int displayTrackNr=0;

        //loop through all sperm tracks
        for (ListIterator iT=theTracks.listIterator(); iT.hasNext();) {
            trackNr++;
            List bTrack=(ArrayList) iT.next();
            if (bTrack.size() >= minTrackLength) {
                //keeps track of the current track
                displayTrackNr++;
                ListIterator jT=bTrack.listIterator();
                particle oldParticle=(particle) jT.next();
                particle firstParticle=new particle();
                firstParticle.copy(oldParticle);
                frames[displayTrackNr-1]=bTrack.size();
                frame = 0;
                holdVAPframe = 0;
                vAPx = 0;
                vAPy = 0;

                holdTotalVAPdistance = 0;
                tempALH = 0;
                for (;jT.hasNext();){
                    particle newParticle=(particle) jT.next();

                    //VCL calculations
                    holdDistance= Math.sqrt(sqr(oldParticle.x-newParticle.x)+sqr(oldParticle.y-newParticle.y));
                    vCLS[displayTrackNr-1] += holdDistance;
                    xyPts = " " + newParticle.x + " " + newParticle.y;

                    //hold the x and y co-ordinates for VAP - will only hold the number of co-ordinates in VAPpoints and will cycle them as this goes
                    for (int j=0; j < (vAPpoints - 1); j++){
                        holdX = xCo[j + 1];
                        holdY = yCo[j + 1];
                        xCo[j] = holdX;
                        yCo[j] = holdY;
                    }

                    xCo[vAPpoints - 1] = newParticle.x;
                    yCo[vAPpoints - 1] = newParticle.y;
                    dX = newParticle.x - oldParticle.x;
                    dY = newParticle.y - oldParticle.y;

                    //calculate the angle of sperm between the two points at hand, but only if we have moved
                    if(holdDistance > 0){
                        //call angle function
                        storA = myAngle(dX, dY);


                        //use the angle from above to calculate the change from the previous angle (all angles are relative to the axis...)
                        if(frame>1){
                            //myAngleDelta returns an array
                            //arrayAngle[0] = aDelta
                            //arrayAngle[1] = rotDir
                            arrayAngle = myAngleDelta(oldA, storA, directionVCL);
                            aChange = arrayAngle[0];
                            directionVCL = arrayAngle[1];

                            //note a direction change by comparing the present clockwise or counter clockwise orientation of angular change to the last
                            if(frame>2){
                                if(oldDirectionVCL != directionVCL){
                                    dirChangesVCL[displayTrackNr-1]++;
                                }
                            }else{
                                storTotalRotation += storA;
                            }
                            if(storTotalRotation >=360){
                                dirChangesVCL[displayTrackNr-1]++;
                                storTotalRotation = storTotalRotation - 360;
                            }
                            oldDirectionVCL = directionVCL;
                            framesBeat[displayTrackNr-1] += 1;
                        }
                        oldA = storA;
                    }

                    vAPx = 0;
                    vAPy = 0;

                    //sum the stored vap co-ordinates
                    for(int k=0; k<(vAPpoints); k++) {
                        vAPx += xCo[k];
                        vAPy += yCo[k];
                    }

                    //generate an average to make a vap point
                    vAPx = vAPx / vAPpoints;
                    vAPy = vAPy / vAPpoints;
                    

                    //if current frame is beyond the number of points used in VAP than calculate distance on VAP path
                    if(frame > (vAPpoints + 1)){

                        xyVAPpts =" " + vAPx + " " + vAPy;

                        //calculates VSL
                        if(holdVAPframe ==0){
                            firstVAPy = vAPy;
                            firstVAPx = vAPx;
                            firstVAPpts = " " + firstVAPx + " " + firstVAPy;
                        }
                        else{
                            firstVAPpts = " " + " ";
                        }

                        //calculate dif of vap point from vcl point for beat cross calculation
                        if(vAPx > xCo[vAPpoints / 2]){
                            vAPxG = 1;
                        }else if(vAPx < xCo[vAPpoints / 2]){
                            vAPxG = 0;
                        }
                        if(vAPy > yCo[vAPpoints / 2]){
                            vAPyG = 1;
                        }else if(vAPy < yCo[vAPpoints / 2]){
                            vAPyG = 0;
                        }

                        if(holdVAPframe > 0){
                            //beat cross calculation
                            if((vAPxG != vAPxGlast) || (vAPyG != vAPyGlast)){
                                beatCross[displayTrackNr-1]++;
                            }
                            vAPxGlast = vAPxG;
                            vAPyGlast = vAPyG;

                            distFirstVAPcurrent = Math.sqrt(sqr(firstVAPx - vAPx) + sqr(firstVAPy - vAPy));
                            sumVAPdistsOrigin[displayTrackNr -1] += distFirstVAPcurrent;
                            if(holdVAPframe ==1){
                                storMaxVAPVSL = distFirstVAPcurrent;
                            }else{
                                if(distFirstVAPcurrent > storMaxVAPVSL){
                                    storMaxVAPVSL = distFirstVAPcurrent;
                                    holdMaxX = vAPx;
                                    holdMaxY = vAPy;
                                }
                            }
                        }

                        //keep tracks of frames used in VAP
                        holdVAPframe ++;

                        //calculate the distance between VAP points and add them to the array
                        holdTotalVAPdistance = Math.sqrt(sqr(vAPxOld-vAPx)+sqr(vAPyOld-vAPy));
                        vAPS[displayTrackNr-1] += holdTotalVAPdistance;
                        tempALH = 2*Math.sqrt(sqr((newParticle.x - oldParticle.x)-(vAPxOld-vAPx))+sqr((newParticle.x - oldParticle.x) - (vAPyOld-vAPy)));
                        aLHS[displayTrackNr-1] += tempALH;
                        maxALH[displayTrackNr-1] = Math.max(maxALHold, tempALH);
                        if(holdTotalVAPdistance == 0){
                            numVAPzero[displayTrackNr-1]++;
                        }

                    }
                    
                    vAPxOld = vAPx;
                    vAPyOld = vAPy;
                    maxALHold = maxALH[displayTrackNr-1];
                    frame++;

                    //save the last two VCL points
                    oldParticle=newParticle;
                    if(printXY !=0){
			// Log the tracked sperms IDs and their XY positions
			//IJ.write(trackNr + xyPts); // Eder mod
			//IJ.write(trackNr + " " + xyPts);
                    }
                }
                //put the calculated frames into an array, stored for each sperm and store VSL
                frames[displayTrackNr-1] = frame;
                framesVAP[displayTrackNr-1] = holdVAPframe;

                vSLS[displayTrackNr-1]=storMaxVAPVSL;
            }
        }
        //hold sum of all changes from frame to frame for a given sperm before generating a per second value for that sperm
        double holdVAP =0;
        double holdVCL = 0;
        double holdVSL = 0;
        double holdLIN = 0;
        double holdWOB = 0;
        double holdDistOrigin = 0;
        double holdBeatCross = 0;
        double holdDirChanges = 0;
        //carry over frame rate so that it can be used in calculation
        double frameRateCalc = frameRate;


        //sum the per frame movements
        for(int z=0; z<displayTrackNr; z++){

            holdVAP = vAPS[z];
            holdVCL = vCLS[z];
            holdVSL = vSLS[z];

            holdZeroVAPs = numVAPzero[z];
            holdLowVAPs = numLowVAP[z];
            holdBeatCross = beatCross[z];
            holdDirChanges = dirChangesVCL[z];
            holdDistOrigin = sumVAPdistsOrigin[z];

            //average the per frame to a per second value
            vCLS[z] = holdVCL * (frameRateCalc / (frames[z])) * microPerPixel;
            vAPS[z] = holdVAP * (frameRateCalc / (framesVAP[z])) * microPerPixel;
            vSLS[z] = holdVSL * (frameRateCalc / (framesVAP[z])) * microPerPixel;
            aLHS[z] = aLHS[z] * microPerPixel / (framesVAP[z]);
            maxALH[z] = maxALH[z] * microPerPixel;
            beatCross[z] = holdBeatCross *(frameRateCalc / (framesBeat[z] - 1));
            dirChangesVCL[z] = holdDirChanges * (frameRateCalc / (framesBeat[z] -1));
            //calculate LIN and WOB
            holdLIN = holdVSL / holdVCL;
            holdWOB = holdVSL / holdVAP;
            numVAPzero[z] = holdZeroVAPs / framesVAP[z];
            numLowVAP[z] = holdLowVAPs / framesVAP[z];
            sumVAPdistsOrigin[z] = holdDistOrigin *(frameRateCalc / (framesVAP[z])) *microPerPixel;


        }

        //keep track of velocities for the whole sperm population
        //these are used for description of motility character for all motile sperm
        double totalVAP = 0;
        double totalVCL = 0;
        double totalVSL = 0;
        double totalLIN = 0;
        double totalSTR = 0;
        double totalEffic = 0;
        double totalWOBbeats = 0;
        double totalProg = 0;
        double totalLINeffic = 0;
        double totalBeatCross = 0;
        double totalALH = 0;
        double totalMaxALH = 0;

        //these hold the average (total value in variable above divided by the number of motile sperm)
        double avgVAP = 0;
        double avgVCL = 0;
        double avgVSL = 0;
        double avgEffic = 0;
        double avgWOBbeats = 0;
        double avgProg = 0;
        double avgLINeffic = 0;
        //these are the values for VAP/ VCL and VSL/ VAP to describe path curvature of all motile sperm
        double avgLin = 0;
        double avgSTR = 0;
        double avgALH = 0;
        double avgMaxALH = 0;
        double avgWob = 0;
        double avgBeats = 0;

        //arrays to hold velocity characteristics for all motile sperm and then to calculate median values
        double[] motileVAP = new double [trackCount];
        double[] motileVSL = new double [trackCount];
        double[] motileVCL = new double [trackCount];
        double[] motileWOB = new double [trackCount];
        double[] motileLIN = new double [trackCount];
        double[] motileSTR = new double [trackCount];
        double[] motileALH = new double [trackCount];
        double[] motileMaxALH = new double [trackCount];
        double[] motileEffic = new double [trackCount];
        double[] motileLINeffic = new double [trackCount];
        double[] motileWOBbeats = new double [trackCount];
        double[] motileBeats = new double [trackCount];
        double[] motileProg = new double [trackCount];
        double[] motileTracks = new double [trackCount];

        //these hold the values to be added to the motility array and the total for motile sperm
        double addLIN = 0;
        double addSTR = 0;
        double addWOB = 0;
        double addBeatCross =  0;
        double addVAP =0;
        double addVSL =0;
        double addVCL =0;
        double addALH = 0;
        double addMaxALH = 0;


        //this is where we determine if a sperm is motile
        if(printSperm != 0){
            IJ.write("Parameters for motile sperms:");
            IJ.write("VCL              VAP              VSL              LIN             STR             WOB             BeatCross             ALH");
        }
        for(int m=0;m<displayTrackNr; m++){

            totalSperm++;
            motileTracks[m] = 0;
            //two tiers for motility determination, first check a set of characteristics then if we meet those criteria, check three other sets... have to meet the first tier and one of the second tiers
            if((vSLS[m] > minVSL) && (vCLS[m] > minVCL)  && (vAPS[m] > minVAP)){
                    addVAP =(vAPS[m]);
                    addVSL =(vSLS[m]);
                    addVCL =(vCLS[m]);
                    addALH = (aLHS[m]);
                    addMaxALH = (maxALH[m]);
                    addBeatCross =beatCross[m];

                    addWOB = addVAP / addVCL;
                    addLIN = addVSL / addVCL;
                    addSTR = addVSL / addVAP;
//                    if(printMedian != 0){
                        motileVAP[motileSperm] = addVAP;
                        motileVSL[motileSperm] = addVSL;
                        motileVCL[motileSperm] = addVCL;
                        motileWOB[motileSperm] = addWOB;
                        motileLIN[motileSperm] = addLIN;
                        motileSTR[motileSperm] = addSTR;
                        motileALH[motileSperm] = addALH;
                        motileMaxALH[motileSperm] = addMaxALH;
                        motileEffic[motileSperm] = sumVAPdistsOrigin[m]/ addVAP;
                        motileBeats[motileSperm] = addBeatCross;
                        motileProg[motileSperm] = sumVAPdistsOrigin[m];
                        motileLINeffic[motileSperm] = (addLIN / (sumVAPdistsOrigin[m] / addVAP))*100;
//                    }
                    if(printSperm != 0){
                        //output characteristics for all motile sperm if needed
String sperm = "" + (float)addVCL + "    " + (float)addVAP + "    " + (float)addVSL + "    " + (float)addLIN + "    " + (float)addSTR + "    " + (float)addWOB + "    "  + (float)addBeatCross +"    " + (float)addALH;
                        IJ.write(sperm);
                    }

                if(motileVAP[motileSperm] < maxVAPSlow){
                    indexSlow.add(m);
                    motileSLOW++;
                } else if (motileVAP[motileSperm] >= maxVAPSlow && motileVAP[motileSperm] < maxVAPMedium){
                    indexMedium.add(m);
                    motileMEDIUM++;
                }  else{
                    indexRapid.add(m);
                    motileRAPID++;
                }

                //TODO PARAMETERS con secondo parametro moltiplicato *100
                if(motileVAP[motileSperm] > cutOffMinVAPforMotile && motileSTR[motileSperm] > (cutOffMinMotileSTR/(float)100)){
                    motilePM++;
                }


                //if the sperm is motile we add in its values into the holding variable to determine the average motility characteristics for the sample
                    totalVAP+=(vAPS[m]);
                    totalVSL+=(vSLS[m]);
                    totalVCL+=(vCLS[m]);
                    totalSTR+= addSTR;
                    totalLIN+= addLIN;
                    totalALH+= addALH;
                    totalMaxALH += addMaxALH;
                    totalBeatCross +=beatCross[m];
                    totalEffic += sumVAPdistsOrigin[m] / addVAP;
                    totalWOBbeats += addWOB / addBeatCross;
                    totalProg += sumVAPdistsOrigin[m];
                    totalLINeffic += (addLIN / (sumVAPdistsOrigin[m] / addVAP)) * 100;

                    motileSperm=motileSperm + 1;
                    motileTracks[m] = 1;
            }
        }
        
        
        
//        IJ.write("----->Motile sperms: " + motileSperm);
//
//        if(printMedian != 0){
//        IJ.write("----->Motile sperms SLOW: " + motileSLOW);
//        IJ.write("----->Motile sperms MEDIUM: " + motileMEDIUM);
//        IJ.write("----->Motile sperms RAPID: " + motileRAPID);
//        IJ.write("----->Motile sperms ULTRARAPID: " + motileULTRARAPID);
//        IJ.write("----->Motile sperms PM: " + motilePM);
//        }

        int middleMotArray = 0;
        double medVAP = 0;
        double medVCL = 0;
        double medVSL = 0;
        double medLIN = 0;
        double medSTR = 0;
        double medWOB = 0;
        double medALH = 0;
        double medMaxALH = 0;

        double sigmaVAP = 0;
        double sigmaVCL = 0;
        double sigmaVSL = 0;
        double sigmaLIN = 0;
        double sigmaSTR = 0;
        double sigmaBeats = 0;
        double sigmaWOB = 0;
        double sigmaALH = 0;
        double sigmaMaxALH = 0;


        double medEffic = 0;
        double medWOBbeats = 0;
        double medBeats = 0;
        double medProg = 0;
        double medLINeffic = 0;

        double[] copyVCL = new double[motileSperm];
        double[] copyVAP = new double[motileSperm];
        double[] copyVSL = new double[motileSperm];
        double[] copyLIN = new double[motileSperm];
        double[] copySTR = new double[motileSperm];
        double[] copyWOB = new double[motileSperm];
        double[] copyBeats = new double[motileSperm];
        double[] copyALH = new double[motileSperm];
        double[] copyMaxALH = new double[motileSperm];

        System.arraycopy(motileVCL,0,copyVCL, 0,motileSperm);
        System.arraycopy(motileVAP,0,copyVAP, 0,motileSperm);
        System.arraycopy(motileVSL,0,copyVSL, 0,motileSperm);
        System.arraycopy(motileLIN,0,copyLIN, 0,motileSperm);
        System.arraycopy(motileSTR,0,copySTR, 0,motileSperm);
        System.arraycopy(motileWOB,0,copyWOB, 0,motileSperm);
        System.arraycopy(motileBeats,0,copyBeats, 0,motileSperm);
        System.arraycopy(motileALH,0,copyALH, 0,motileSperm);
        System.arraycopy(motileMaxALH,0,copyMaxALH, 0,motileSperm);

        Arrays.sort(motileVAP);
        Arrays.sort(motileVSL);
        Arrays.sort(motileVCL);
        Arrays.sort(motileBeats);
        Arrays.sort(motileLIN);
        Arrays.sort(motileALH);
        Arrays.sort(motileMaxALH);
        Arrays.sort(motileSTR);
        Arrays.sort(motileWOB);

        Arrays.sort(motileEffic);
        Arrays.sort(motileWOBbeats);
        Arrays.sort(motileBeats);
        Arrays.sort(motileProg);
        Arrays.sort(motileLINeffic);

        middleMotArray = motileSperm /2;

        if (((motileSperm % 2) == 0) && (motileSperm > 2)) {
            medVAP = (motileVAP[trackCount -middleMotArray] + motileVAP[trackCount - middleMotArray -1]) /2;
            medVCL = (motileVCL[trackCount - middleMotArray] + motileVCL[trackCount - middleMotArray-1]) /2;
            medVSL = (motileVSL[trackCount - middleMotArray] + motileVSL[trackCount - middleMotArray-1])/2;
            medLIN = (motileLIN[trackCount - middleMotArray] + motileLIN[trackCount - middleMotArray - 1]) / 2;
            medSTR = (motileSTR[trackCount - middleMotArray] + motileSTR[trackCount - middleMotArray - 1]) / 2;
            medALH = (motileALH[trackCount - middleMotArray] + motileALH[trackCount - middleMotArray - 1]) / 2;
            medMaxALH = (motileMaxALH[trackCount - middleMotArray] + motileMaxALH[trackCount - middleMotArray - 1]) / 2;
            medWOB = (motileWOB[trackCount - middleMotArray] + motileWOB[trackCount - middleMotArray - 1]) / 2;
            medEffic = (motileEffic[trackCount - middleMotArray] + motileEffic[trackCount - middleMotArray - 1]) / 2;
            medBeats = (motileBeats[trackCount - middleMotArray] + motileBeats[trackCount - middleMotArray-1]) / 2;
            medProg = (motileProg[trackCount - middleMotArray] + motileProg[trackCount - middleMotArray - 1]) / 2;
            medLINeffic = (motileLINeffic[trackCount - middleMotArray] + motileLINeffic[trackCount - middleMotArray - 1]) / 2;
        }else if (motileSperm > 2){
            medVAP = motileVAP[trackCount - middleMotArray];
            medVCL = motileVCL[trackCount - middleMotArray];
            medVSL = motileVSL[trackCount - middleMotArray];
            medLIN = motileLIN[trackCount - middleMotArray];
            medSTR = motileSTR[trackCount - middleMotArray];
            medWOB = motileWOB[trackCount - middleMotArray];
            medALH = motileALH[trackCount - middleMotArray];
            medMaxALH = motileMaxALH[trackCount - middleMotArray];
            medEffic = motileEffic[trackCount - middleMotArray];
            medBeats = motileBeats[trackCount - middleMotArray];
            medProg = motileProg[trackCount - middleMotArray];
            medLINeffic = motileLINeffic[trackCount - middleMotArray];
        }


        //average the velocity characteristics for motile sperm using totals from above
        avgBeats = totalBeatCross / motileSperm;
        avgVAP = totalVAP / motileSperm;
        avgVCL = totalVCL / motileSperm;
        avgVSL = totalVSL / motileSperm;
        avgEffic = totalEffic / motileSperm;
        avgProg = totalProg / motileSperm;
        avgLINeffic = totalLINeffic / motileSperm;

        avgLin = totalLIN/motileSperm;
        avgSTR = totalSTR/motileSperm;
        avgALH = totalALH/motileSperm;
        avgMaxALH = totalMaxALH/motileSperm;
        avgWob = avgVAP/avgVCL;
        //percent motility
        double motility = 0;
        motility = motileSperm;


        for (int i = 0; i <motileSperm ; i++) {
            sigmaVCL += (copyVCL[i] - avgVCL)*(copyVCL[i] - avgVCL);
            sigmaVAP += (copyVAP[i] - avgVAP)*(copyVAP[i] - avgVAP);
            sigmaVSL += (copyVSL[i] - avgVSL)*(copyVSL[i] - avgVSL);
            sigmaLIN += (copyLIN[i] - avgLin)*(copyLIN[i] - avgLin);
            sigmaSTR += (copySTR[i] - avgSTR)*(copySTR[i] - avgSTR);
            sigmaWOB += (copyWOB[i] - avgWob)*(copyWOB[i] - avgWob);
            sigmaBeats += (copyBeats[i] - avgBeats)*(copyBeats[i] - avgBeats);
            sigmaALH += (copyALH[i] - avgALH)*(copyALH[i] - avgALH);
            sigmaMaxALH += (copyMaxALH[i] - avgMaxALH)*(copyMaxALH[i] - avgMaxALH);
        }

        int campione = motileSperm - 1;

        sigmaVCL = Math.sqrt(sigmaVCL/campione);
        sigmaVAP = Math.sqrt(sigmaVAP/campione);
        sigmaVSL = Math.sqrt(sigmaVSL/campione);
        sigmaLIN = Math.sqrt(sigmaLIN/campione);
        sigmaSTR = Math.sqrt(sigmaSTR/campione);
        sigmaWOB = Math.sqrt(sigmaWOB/campione);
        sigmaBeats = Math.sqrt(sigmaBeats/campione);
        sigmaALH = Math.sqrt(sigmaALH/campione);
        sigmaMaxALH = Math.sqrt(sigmaMaxALH/campione);


        String sperm2 = "TotalSperm: " + (float)totalSperm;
        IJ.write(sperm2);
        
        
        if(printMedian != 0){
        String sperm = "TM         PM         AvgVCL     AvgVAP     AvgVSL     AvgLIN     AvgSTR     AvgWOB     AvgBeats    AvgALH";
        String _sperm = (float)motility + "       " + (float)motilePM + "       " + (float)avgVCL + "       " + (float)avgVAP + "       " + (float)avgVSL + "       " + (float)avgLin + "       " + (float)avgSTR + "       " + (float)avgWob +"       " + (float)avgBeats+ "       " + (float)avgALH;
        IJ.write(sperm);
        IJ.write(_sperm);
        }
        
        
        if(printMedian != 0){
        String tipology = "SLOW     MEDIUM     RAPID";
        String _tipology = (float)motileSLOW + "   " + (float)motileMEDIUM + "   " + (float)motileRAPID;
        IJ.write(tipology);
            IJ.write(_tipology);
        String sperm1 = "MedianVCL        MedianVAP        MedianVSL        MedianLIN        MedianSTR        MedianWOB        MedianBeats        MedianALH";
            String _sperm1 = (float)medVCL + "   " + (float)medVAP + "   " + (float)medVSL + "   " + (float)medLIN + "   " + (float)medSTR + "   " + (float)medWOB + "   " + (float)medBeats + "   " + (float)medALH;
        IJ.write(sperm1);
            IJ.write(_sperm1);

        String sigma = "SigmaVCL        SigmaVAP        SigmaVSL        SigmaLIN        SigmaSTR        SigmaWOB        SigmaBeats        SigmaALH";
        String _sigma = (float)sigmaVCL + "   " + (float)sigmaVAP + "   " + (float)sigmaVSL + "   " + (float)sigmaLIN + "   " + (float)sigmaSTR + "   " + (float)sigmaWOB + "   " + (float)sigmaBeats + "   " + (float)sigmaALH;
        IJ.write(sigma);
        IJ.write(_sigma);


        }
        
        IJ.write("#############  --- Color LEGEND  --- #############");
        IJ.write("RED: RAPID");
        IJ.write("ORANGE: MEDIUM");
        IJ.write("GREEN: SLOW");
        IJ.write("BLACK: NOT CLASSIFIED");
        
        // 'map' of tracks
        if (imp.getCalibration().scaled()) {
            IJ.showMessage("MultiTracker", "Cannot display paths if image is spatially calibrated");
            return;
        }
        int upRes = 1;
        ImageProcessor ip = new ColorProcessor(imp.getWidth()*upRes, imp.getHeight()*upRes);
        ip.setColor(Color.white);
        ip.fill();
        int trackCount2=0;
        int trackCount3=0;
        int trackCount4=0;
        Color color;
        for (ListIterator iT=theTracks.listIterator();iT.hasNext();) {
            trackCount2++;
            List zTrack=(ArrayList) iT.next();
            if (zTrack.size() >= minTrackLength) {

                ListIterator jT=zTrack.listIterator();
                particle oldParticle=(particle) jT.next();

                color = Color.black;
                if(motileTracks[trackCount3] > 0){
                    if (indexRapid.contains(trackCount3)){
                        color = Color.red;
                        trackCount4++;
                    }
                    else if (indexMedium.contains(trackCount3)){
                        color = Color.orange;
                        trackCount4++;
                    }
                    else if(indexSlow.contains(trackCount3)){
                        color = Color.green;
                        trackCount4++;
                    }
                    else{
                        color = Color.black;
                    }
                }
                trackCount3++;

                for (;jT.hasNext();) {
                    particle newParticle=(particle) jT.next();
                    ip.setColor(color);
                    ip.moveTo((int)oldParticle.x*upRes, (int)oldParticle.y*upRes);
                    ip.lineTo((int)newParticle.x*upRes, (int)newParticle.y*upRes);
			if (color != Color.black){
			String xyCoord = (int)trackCount3 + " " + (int)newParticle.x + " " + (int)newParticle.y; //Eder mod
			IJ.write(xyCoord);} // identifying the trajectory? (oldParticle.trackNr)

                    oldParticle=newParticle;
                }
                if(color != Color.black) {
			ip.drawString(Integer.toString(trackCount3));}

            }
        }
        new ImagePlus("Paths", ip).show();

    }


    // Utility functions
    double sqr(double n) {return n*n;}

    int doOffset (int center, int maxSize, int displacement) {
        if ((center - displacement) < 2*displacement) {
            return (center + 4*displacement);
        }
        else {
            return (center - displacement);
        }
    }

    double s2d(String s) {
        Double d;
        try {d = new Double(s);}
        catch (NumberFormatException e) {d = null;}
        if (d!=null)
            return(d.doubleValue());
        else
            return(0.0);
    }

}
