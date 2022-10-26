/* Define 3D ROI slice by slice*/


w=getWidth();
h=getHeight();
s=nSlices();
newImage("Mask", "8-bit Black", w, h, s);
setForegroundColor(255, 255, 255);
nb=roiManager("count");
for(i=0;i<nb;i++) {
        roiManager("Select",i);
run("Fill", "slice");
} 