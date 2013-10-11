

#include "tracer.h"
#include "objreader.h"

#define EPSILON 0.000001
using namespace std;

#define IM_DEBUGGING 
int recursionDepth = 0;
int maxDepth = 1;
Array2D<Rgba> p(1500,1500) ;
Array2D<Rgba> Sepia(1500,1500) ;
Array2D<Rgba> BandW(1500,1500) ;
Array2D<Rgba> WildBlue(1500,1500) ;
Array2D<Rgba> TVChannel(1500,1500) ;
Array2D<Rgba> Pixelated(1500,1500) ;

double nearVal = 10000000;
vector<items *> objects;
vector<items *> checkedObjects;
vector<items *> tmp;
vector<lights *> lighting;
vector<materials *> material;
vector<string> camera;
vector<string> lights;
string obj;
float eyex, eyey, eyez, wx, wy, wz, cd, iw, ih;
float ux, uy, uz, uMag, vx, vy, vz;
int pw, ph;
items *sph;
items *temp;
class lights *lig;
materials * mat;
int numLights = 0;
float holdPandT[10][2250000];
char type[2250000];
double normCalc[9][2250000];
int lastMaterialLoaded = -1;
int leafCounter = 0;


float BoxIntersection(vector<items *> box){
    //check for intersection with bounding box
    float nwx, nwy, nwz;
    float nwMag = sqrt(pow(wx, 2) + pow(wy, 2) + pow(wz, 2));
    nwx = -1*(wx/nwMag);
    nwy = -1*(wy/nwMag);
    nwz = -1*(wz/nwMag);
    ux = -1* nwz;
    uy = 0;
    uz = 1*(nwx);
    uMag = sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
    ux = ux/uMag;
    uy = uy/uMag;
    uz = uz/uMag;
    
    // v = u x w
    vx = (uy*nwz) - (uz*nwy);
    vy = (uz*nwx) - (ux*nwz);
    vz = (ux*nwy) - (uy*nwx);
    
    //cout << "iw, ih: " << iw << ", " << ih << endl;
    
    float widthE = iw/2;
    float heightE = ih/2;
    float widthB = widthE * -1;
    float heightB = heightE * -1;
    float du = iw/pw;
    float dv = ih/ph;
    int pixelx = -1;
    int pixely = -1;
    int intersects = 0;
    for(double u=widthB; u<widthE; u+=du){
        pixelx +=1;
        pixely =-1;
        //u += ju;
        for(double v=heightB; v<heightE; v+=dv){
            pixely +=1;
            //cout << "Pixels: " << pixelx << ", " << pixely << endl;
            double directionx = cd*wx + (u)*ux + (v)*vx;
            double directiony = cd*wy + (u)*uy + (v)*vy;
            double directionz = cd*wz + (u)*uz + (v)*vz;
            double dMag = sqrt(pow(directionx, 2) + pow(directiony, 2) + pow(directionz, 2));
            
            directionx = directionx/dMag;
            directiony = directiony/dMag;
            directionz = directionz/dMag;
            nearVal = 10000000;
            
            for (int i = 0; i<box.size(); i++) {
                if (box[i]->getType() == 't'){
                    float vert1x, vert1y, vert1z, vert2x, vert2y, vert2z, vert3x, vert3y, vert3z;
                    vert1x = box[i]->getX1();
                    vert1y = box[i]->getY1();
                    vert1z = box[i]->getZ1();
                    vert2x = box[i]->getX2();
                    vert2y = box[i]->getY2();
                    vert2z = box[i]->getZ2();
                    vert3x = box[i]->getX3();
                    vert3y = box[i]->getY3();
                    vert3z = box[i]->getZ3();
                    
                    float magA = ((vert1x - vert2x) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z))))   + ((vert1x - vert3x) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
                    
                    float beta = ((vert1x - eyex) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) * -1 +  (directionx * (((vert1y-eyey)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-eyez)))) ;
                    
                    float gamma = ((vert1x - vert2x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) + ((vert1x - eyex) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) ;
                    
                    float t = ((vert1x - vert2x) * (((vert1y-vert3y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) * -1 +  ((vert1x - eyex) * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
                    
                    beta = beta/magA;
                    gamma = gamma/magA;
                    t = t/magA;
                    double sumCheck = beta + gamma;
                    
                    if (sumCheck <= 1 && beta>=0 && gamma>=0 && t>0) {
                        
                        if (t<nearVal){
                            int pos = pixelx*1500 + pixely;
                            holdPandT[0][pos] = pixelx;
                            holdPandT[1][pos] = pixely;
                            holdPandT[2][pos] = t;
                            holdPandT[3][pos] = directionx;
                            holdPandT[4][pos] = directiony;
                            holdPandT[5][pos] = directionz;
                            holdPandT[6][pos] = box[i]->getMatID();
                            type[pos] = 't';
                            normCalc[0][pos] = vert1x;
                            normCalc[1][pos] = vert1y;
                            normCalc[2][pos] = vert1z;
                            normCalc[3][pos] = vert2x;
                            normCalc[4][pos] = vert2y;
                            normCalc[5][pos] = vert2z;
                            normCalc[6][pos] = vert3x;
                            normCalc[7][pos] = vert3y;
                            normCalc[8][pos] = vert3z;
                            nearVal = t;
                            intersects =t;
                            //cout << "Boom!" << endl;
                        }
                    }
                }//end of triangle intersection checker
            }//end objects loop
        }
    }
    return intersects;
}

//for the 'axisZ' find the largest and smallest z value
//for the 'axisY' find the largest and smallest y value
//for the 'axisX' find the largest and smallest x value
//sort object array on 'axisZ'

void node(vector<items *> bound, unsigned int treedepth) {
	//vector<items *> bound;
    
    // Build the node bounding box based from the triangle list	
    int N = (int)bound.size();
    //vector<items *> temp;
    float maxX= 0, minX =0, maxY=0, minY=0, maxZ=0, minZ = 0;
    //cout << "N: " << N << endl;
    //find the maximum radius 
    float dMax = 0;
    float planeAdder = 0;
    for (int i = 0; i<N; i++) {
        //cout << "Show Radius: " << bound[i]->getD() << endl;
        if (bound[i]->getD() > dMax) {
            dMax = bound[i]->getD();
        }
    }
    for (int i = 0; i<N; i++) {
        //cout << "Show Radius: " << bound[i]->getD() << endl;
        if (bound[i]->getType() == 'p' ) {
            planeAdder = 1000;
        }
    }
    
    for (int i = 0; i<N; i++) {
        if (i == 0) {
            if (bound[0]->getX1()>=bound[0]->getX2() && bound[0]->getX1()>=bound[0]->getX3()){
                maxX = bound[0]->getX1();
            }
            if (bound[0]->getX1()<=bound[0]->getX2() && bound[0]->getX1()<=bound[0]->getX3()){
                minX = bound[0]->getX1();
            } 
            if (bound[0]->getY1()>=bound[0]->getY2() && bound[0]->getY1()>=bound[0]->getY3()){
                maxY = bound[0]->getY1();
            }
            if (bound[0]->getY1()<=bound[0]->getY2() && bound[0]->getY1()<=bound[0]->getY3()){
                minY = bound[0]->getY1();
            } 
            if (bound[0]->getZ1()>=bound[0]->getZ2() && bound[0]->getZ1()>=bound[0]->getZ3()){
                maxZ = bound[0]->getZ1();
            }
            if (bound[0]->getZ1()<=bound[0]->getZ2() && bound[0]->getZ1()<=bound[0]->getZ3()){
                minZ = bound[0]->getZ1();
            } 
            if (bound[0]->getX2()>=bound[0]->getX1() && bound[0]->getX2()>=bound[0]->getX3()){
                maxX = bound[0]->getX2();
            }
            if (bound[0]->getX2()<=bound[0]->getX1() && bound[0]->getX2()<=bound[0]->getX3()){
                minX = bound[0]->getX2();
            } 
            if (bound[0]->getY2()>=bound[0]->getY1() && bound[0]->getY2()>=bound[0]->getY3()){
                maxY = bound[0]->getY2();
            }
            if (bound[0]->getY2()<=bound[0]->getY1() && bound[0]->getY2()<=bound[0]->getY3()){
                minY = bound[0]->getY1();
            } 
            if (bound[0]->getZ2()>=bound[0]->getZ1() && bound[0]->getZ2()>=bound[0]->getZ3()){
                maxZ = bound[0]->getZ2();
            }
            if (bound[0]->getZ2()<=bound[0]->getZ1() && bound[0]->getZ2()<=bound[0]->getZ3()){
                minZ = bound[0]->getZ2();
            } 
            if (bound[0]->getX3()>=bound[0]->getX2() && bound[0]->getX3()>=bound[0]->getX1()){
                maxX = bound[0]->getX3();
            }
            if (bound[0]->getX3()<=bound[0]->getX2() && bound[0]->getX3()<=bound[0]->getX1()){
                minX = bound[0]->getX3();
            } 
            if (bound[0]->getY3()>=bound[0]->getY2() && bound[0]->getY3()>=bound[0]->getY1()){
                maxY = bound[0]->getY3();
            }
            if (bound[0]->getY3()<=bound[0]->getY2() && bound[0]->getY3()<=bound[0]->getY1()){
                minY = bound[0]->getY3();
            } 
            if (bound[0]->getZ3()>=bound[0]->getZ2() && bound[0]->getZ3()>=bound[0]->getZ1()){
                maxZ = bound[0]->getZ3();
            }
            if (bound[0]->getZ3()<=bound[0]->getZ2() && bound[0]->getZ3()<=bound[0]->getZ1()){
                minZ = bound[0]->getZ3();
            } 
        }//end of if (i = 0)
        else{
            //cout << "MaxX, MinX: " << maxX << ", " << minX << "\nMaxY, MinY: " << maxY << ", " << minY << "\nMaxZ, MinZ: " << maxZ << ", " << minZ << endl;
            if (bound[i]->getX1()>maxX){
                //cout << "test" << maxX << ", " << bound[i]->getX1() << endl;
                maxX = bound[i]->getX1();
                //cout << "test" << maxX << ", " << bound[i]->getX1() << endl;
            }
            if (bound[i]->getX1()<minX){
                minX = bound[i]->getX1();
            } 
            if (bound[i]->getY1()>maxY){
                maxY = bound[i]->getY1();
            }
            if (bound[i]->getY1()<minY){
                minY = bound[i]->getY1();
            } 
            if (bound[i]->getZ1()>maxZ){
                maxZ = bound[i]->getZ1();
            }
            if (bound[i]->getZ1()<minZ){
                minZ = bound[i]->getZ1();
            } 
            if (bound[i]->getX2()>maxX){
                maxX = bound[i]->getX2();
            }
            if (bound[i]->getX2()<minX){
                minX = bound[i]->getX2();
            } 
            if (bound[i]->getY2()>maxY){
                maxY = bound[i]->getY2();
            }
            if (bound[i]->getY2()<minY){
                minY = bound[i]->getY2();
            } 
            if (bound[i]->getZ2()>maxZ){
                maxZ = bound[i]->getZ2();
            }
            if (bound[i]->getZ2()<minZ){
                minZ = bound[i]->getZ2();
            }
            if (bound[i]->getX3()>maxX){
                maxX = bound[i]->getX3();
            }
            if (bound[i]->getX3()<minX){
                minX = bound[i]->getX3();
            } 
            if (bound[i]->getY3()>maxY){
                maxY = bound[i]->getY3();
            }
            if (bound[i]->getY3()<minY){
                minY = bound[i]->getY3();
            } 
            if (bound[i]->getZ3()>maxZ){
                maxZ = bound[i]->getZ3();
            }
            if (bound[i]->getZ3()<minZ){
                minZ = bound[i]->getZ3();
            }
        }//end of else ... i != 0 
    }// end of i for loop
    //cout << "Dmax: " << dMax << endl;
    
    
    maxX = maxX +dMax + planeAdder;
    maxY = maxY +dMax + planeAdder;
    maxZ = maxZ +dMax + planeAdder;
    minX = minX -dMax - planeAdder;
    minY = minY -dMax - planeAdder;
    minZ = minZ -dMax - planeAdder;
    //cout << "MaxX, MinX: " << maxX << ", " << minX << "\nMaxY, MinY: " << maxY << ", " << minY << "\nMaxZ, MinZ: " << maxZ << ", " << minZ << endl;
    //for loop to sort the bound on 'axisY' in a new bound array
    //put sorted value in the array called temp
    //bool swapped = true;
    //int j = 0;
    /*
     while (swapped) {
     swapped = false;
     j++;
     for (int i = 0; i < N-j; i++) {
     if (bound[i]->getY1() > bound[i + 1]->getY1()) {
     temp = new items(bound[i]->getType(), bound[i]->getX1(), bound[i]->getY1(), bound[i]->getZ1(), bound[i]->getX2(), bound[i]->getY2(), bound[i]->getZ2(), bound[i]->getX3(), bound[i]->getY3(), bound[i]->getX3(), bound[i]->getD());
     tmp.push_back(temp) ;
     bound[i] = bound[i + 1];
     bound[i + 1] = tmp[0];
     swapped = true;
     }
     }
     }
     for (int i = 0; i < N; i++) {
     cout << "Sorted bound: " << bound[i]->getY1() << endl;
     }*/
    /*
     cout << "Triangle Vertices for bounding box are" << endl;
     1 cout << maxX << ", " << maxY << ", " << maxZ << endl;
     2 cout << maxX << ", " << maxY << ", " << minZ << endl;
     3 cout << maxX << ", " << minY << ", " << maxZ << endl;
     4 cout << maxX << ", " << minY << ", " << minZ << endl;
     5 cout << minX << ", " << maxY << ", " << maxZ << endl;
     6 cout << minX << ", " << maxY << ", " << minZ << endl;
     7 cout << minX << ", " << minY << ", " << maxZ << endl;
     8 cout << minX << ", " << minY << ", " << minZ << endl;
     */
    vector<items *> tribound;
    //six triangle sides of the bounding box
    items *holder  = new items('t', maxX, maxY, maxZ, maxX, maxY, minZ, maxX, minY, maxZ, 0); //123
    tribound.push_back(holder);
    holder  = new items('t', maxX, maxY, maxZ, maxX, minY, minZ, maxX, minY, maxZ, 0); //143
    tribound.push_back(holder);
    holder  = new items('t', maxX, maxY, maxZ, maxX, minY, minZ, minX, maxY, minZ, 0); //146
    tribound.push_back(holder);
    holder  = new items('t', maxX, maxY, maxZ, minX, maxY, maxZ, minX, maxY, minZ, 0); //156
    tribound.push_back(holder);
    holder  = new items('t', maxX, maxY, maxZ, minX, maxY, maxZ, minX, minY, minZ, 0); //158
    tribound.push_back(holder);
    holder  = new items('t', maxX, maxY, maxZ, maxX, maxY, minZ, minX, minY, minZ, 0); //128
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, minX, minY, minZ, minX, maxY, maxZ, 0); //785
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, minX, minY, minZ, maxX, maxY, minZ, 0); //782
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, maxX, minY, maxZ, maxX, maxY, minZ, 0); //732
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, maxX, minY, maxZ, maxX, minY, minZ, 0); //734
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, minX, maxY, minZ, maxX, minY, minZ, 0); //764
    tribound.push_back(holder);
    holder  = new items('t', minX, minY, maxZ, minX, maxY, minZ, minX, maxY, maxZ, 0); //765
    tribound.push_back(holder);
    
    
    if (BoxIntersection(tribound) != 0) {
        vector<items *> left;
        vector<items *> right;
        int sizeLeft;
        if (N == 1) {
            /*
             cout << "\nLeaf Intersection Found!" << endl;
             cout << "Object coordinates: " << endl;
             cout << bound[0]->getX1() << ", " << bound[0]->getY1() << ", " << bound[0]->getZ1() << endl;
             cout << bound[0]->getX2() << ", " << bound[0]->getY2() << ", " << bound[0]->getZ2() << endl;
             cout << bound[0]->getX3() << ", " << bound[0]->getY3() << ", " << bound[0]->getZ3() << endl;
             cout << "Object type: " << (char)bound[0]->getType() << endl;
             cout << "MaxX, MinX: " << maxX << ", " << minX << "\nMaxY, MinY: " << maxY << ", " << minY << "\nMaxZ, MinZ: " << maxZ << ", " << minZ << endl;
             */
            leafCounter++;
            items *tempo  = new items(bound[0]->getType(), bound[0]->getX1(), bound[0]->getY1(), bound[0]->getZ1(), bound[0]->getX2(), bound[0]->getY2(), bound[0]->getZ2(), bound[0]->getX3(), bound[0]->getY3(), bound[0]->getZ3(), bound[0]->getD());
            tempo->StoreMatID(bound[0]->getMatID());
            checkedObjects.push_back(tempo);
        }
        else if (N%2 == 0) {
            sizeLeft = N/2;
            for (int k = 0; k<sizeLeft; k++) {
                holder  =new items(bound[k]->getType(), bound[k]->getX1(), bound[k]->getY1(), bound[k]->getZ1(), bound[k]->getX2(), bound[k]->getY2(), bound[k]->getZ2(), bound[k]->getX3(), bound[k]->getY3(), bound[k]->getZ3(), bound[k]->getD());
                holder->StoreMatID(bound[k]->getMatID());
                left.push_back(holder);
                
            }
            for (int k = sizeLeft; k<N; k++) {
                holder  =new items(bound[k]->getType(), bound[k]->getX1(), bound[k]->getY1(), bound[k]->getZ1(), bound[k]->getX2(), bound[k]->getY2(), bound[k]->getZ2(), bound[k]->getX3(), bound[k]->getY3(), bound[k]->getZ3(), bound[k]->getD());   
                holder->StoreMatID(bound[k]->getMatID());
                right.push_back(holder);
            }
            node(left, 0);
            node(right, 0);
        }
        else{
            sizeLeft = (N+1)/2;
            for (int k = 0; k<sizeLeft; k++) {
                holder  =new items(bound[k]->getType(), bound[k]->getX1(), bound[k]->getY1(), bound[k]->getZ1(), bound[k]->getX2(), bound[k]->getY2(), bound[k]->getZ2(), bound[k]->getX3(), bound[k]->getY3(), bound[k]->getZ3(), bound[k]->getD());
                holder->StoreMatID(bound[k]->getMatID());
                left.push_back(holder);
                
            }
            for (int k = sizeLeft; k<N; k++) {
                holder  =new items(bound[k]->getType(), bound[k]->getX1(), bound[k]->getY1(), bound[k]->getZ1(), bound[k]->getX2(), bound[k]->getY2(), bound[k]->getZ2(), bound[k]->getX3(), bound[k]->getY3(), bound[k]->getZ3(), bound[k]->getD());
                holder->StoreMatID(bound[k]->getMatID());
                right.push_back(holder);
            }
            node(left, 0);
            node(right, 0);
        }
    }//if hit
    else{
        //cout << "Did not hit bounding box!" << endl;
    }
    
}


// Write an RGBA image using class RgbaOutputFile.
void writeRgba (const char fileName[],const Rgba *pixels,int width,int height)
{
    
    RgbaOutputFile file (fileName, width, height, WRITE_RGBA);
    file.setFrameBuffer (pixels, 1, width);
    file.writePixels (height);
}

float getTokenAsFloat (string inString, int whichToken){
    
    float thisFloatVal = 0.;    // the return value
    
    if (whichToken == 0) {
        cerr << "error: the first token on a line is a character!" << endl;
        exit (-1);
    }
    
    // c++ string class has no super-easy way to tokenize, let's use c's:
    char *cstr = new char [inString.size () + 1];
    
    strcpy (cstr, inString.c_str());
    
    char *p = strtok (cstr, " ");
    if (p == 0) {
        cerr << "error: the line has nothing on it!" << endl;
        exit (-1);
    }
    
    for (int i = 0; i < whichToken; i++) {
        p = strtok (0, " ");
        if (p == 0) {
            cerr << "error: the line is not long enough for your token request!" << endl;
            exit (-1);
        }
    }
    
    thisFloatVal = atof (p);
    
    delete[] cstr;
    
    return thisFloatVal;
}

string getTokenAsString (string inString, int whichToken){
    
    float thisFloatVal = 0.;    // the return value
    
    if (whichToken == 0) {
        cerr << "error: the first token on a line is a character!" << endl;
        exit (-1);
    }
    
    // c++ string class has no super-easy way to tokenize, let's use c's:
    char *cstr = new char [inString.size () + 1];
    
    strcpy (cstr, inString.c_str());
    
    char * p = strtok (cstr, " ");
	char * hold;
    if (p == 0) {
        cerr << "error: the line has nothing on it!" << endl;
        exit (-1);
    }
    
    for (int i = 0; i < whichToken; i++) {
        p = strtok (0, " ");
        hold = p;
        //	cout << "test: " << p << endl;
        if (p == 0) {
            cerr << "error: the line is not long enough for your token request!" << endl;
            exit (-1);
        }
    }
    
    //    thisFloatVal = atof (p);
    //cout << "test: " << p << "hold " << hold  << endl;      
    //    delete[] cstr;
    //cout << "test: " << p << "hold " << hold  << endl;  
    return p;
}



void parseSceneFile (char *filnam)
{
    
    // open the file
    ifstream inFile(filnam);    
    string line;
    
    // if it's not open, error out.
    if (! inFile.is_open ()) {  
        cerr << "can't open scene file" << endl;
        exit (-1);
    }
    
    
    
    while (! inFile.eof ()) {   // go through every line in the file until finished
        
        getline (inFile, line); // get the line
        
        switch (line[0])  {     
                
            case 's': //sphere
                
                float x, y, z, r;
                x = getTokenAsFloat (line, 1); 
                y = getTokenAsFloat (line, 2); 
                z = getTokenAsFloat (line, 3); 
                r = getTokenAsFloat (line, 4); 
                
                
                sph = new items('s', x, y, z, 0, 0, 0, 0, 0, 0, r);
                sph->StoreMatID(lastMaterialLoaded);
                objects.push_back(sph);
                
                
#ifdef IM_DEBUGGING
                // if we're debugging, show what we got:
                cout << "" << endl;
                cout << "got a sphere with ";
                cout << "parameters: " << x << " " << y << " " << z << " " << r << endl;
                cout << "" << endl;
#endif
                break;
                
            case 't':   // triangle
                
                float x1, y1, z1, x2, y2, z2,x3, y3, z3;
                x1 = getTokenAsFloat (line, 1); 
                y1 = getTokenAsFloat (line, 2); 
                z1 = getTokenAsFloat (line, 3); 
                x2 = getTokenAsFloat (line, 4); 
                y2 = getTokenAsFloat (line, 5); 
                z2 = getTokenAsFloat (line, 6);
                x3 = getTokenAsFloat (line, 7); 
                y3 = getTokenAsFloat (line, 8); 
                z3 = getTokenAsFloat (line, 9);
                
                sph = new items('t', x1, y1, z1, x2, y2, z2, x3, y3, z3, 0);
                sph->StoreMatID(lastMaterialLoaded);
                objects.push_back(sph);
                
#ifdef IM_DEBUGGING
                // if we're debugging, show what we got:
                cout << "" << endl;
                cout << "got a triangle with " << endl;
                cout << "x1, y1, z1: " << x1 << "," << y1 << "," << z1 << endl;
                cout << "x2, y2, z2: " << x2 << "," << y2 << "," << z2 << endl;
                cout << "x3, y3, z3: " << x3 << "," << y3 << "," << z3 << endl;
                cout << "" << endl;
#endif
                break;
                
            case 'p':   // plane
                
                float  nx, ny, nz, d; 
                nx = getTokenAsFloat (line, 1); 
                ny = getTokenAsFloat (line, 2); 
                nz = getTokenAsFloat (line, 3); 
                d = getTokenAsFloat (line, 4); 
                
                sph = new items('p', nx, ny, nz, 0, 0, 0, 0, 0, 0, d);
                sph->StoreMatID(lastMaterialLoaded);
                objects.push_back(sph);
                
#ifdef IM_DEBUGGING
                // if we're debugging, show what we got:
                cout << "" << endl;
                cout << "got a plane with "<< endl;
                cout << "nx, ny, nz, and d: " << nx << "," << ny << "," << nz << " and " << d << endl;
                cout << "" << endl;
                
#endif
                break;
                
                
                // camera:
            case 'c':  
                // one trick here: the cameras pixel count (width, height) are integers,
                // so cast them.
                
                eyex = getTokenAsFloat (line, 1); 
                eyey = getTokenAsFloat (line, 2); 
                eyez = getTokenAsFloat (line, 3); 
                wx = getTokenAsFloat (line, 4); 
                wy = getTokenAsFloat (line, 5); 
                wz = getTokenAsFloat (line, 6);
                cd = getTokenAsFloat (line, 7); 
                iw = getTokenAsFloat (line, 8); 
                ih = getTokenAsFloat (line, 9);
                pw = (int)getTokenAsFloat (line, 10);
                ph = (int)getTokenAsFloat (line, 11);
                
#ifdef IM_DEBUGGING
                // if we're debugging, show what we got:
                cout << "" << endl;
                cout << "got a camera with " << endl;
                cout << "position: " << eyex << ", " << eyey << ", " << eyez << endl;
                cout << "cam direction: " << wx << ", " << wy << ", " << wz << endl;
                cout << "focal length: " << cd << endl;
                cout << "image plane size: " << iw << " by " << ih << endl;
                cout << "pixels: " << pw << " by " << ph << endl;
                cout << "" << endl;
                
#endif
                break;
                
                
                // lights:
            case 'l':   
                switch (line[2]) {
                    case 'p':   // point light
                        float lx, ly, lz, lr1, lg1, lb1;
                        
                        lx  = getTokenAsFloat (line, 2); 
                        ly  = getTokenAsFloat (line, 3); 
                        lz  = getTokenAsFloat (line, 4); 
                        lr1 = getTokenAsFloat (line, 5); 
                        lg1 = getTokenAsFloat (line, 6); 
                        lb1 = getTokenAsFloat (line, 7);
                        
                        lig = new class lights('p', lx, ly, lz, lr1, lg1, lb1);
                        lighting.push_back(lig);
                        numLights += 1;
                        
#ifdef IM_DEBUGGING
                        // if we're debugging, show what we got:
                        cout << "" << endl;
                        cout << "got a point light with ";
                        cout << "position: " << lx << ", " << ly << ", " << lz << endl;
                        cout << "RGB values: " << lr1 << ", " << lg1 << ", " << lb1 << endl;
                        cout << "" << endl;
#endif
                        
                        break;
                    case 'd':   // directional light
                        float lvx, lvy, lvz, lr2, lg2, lb2;
                        
                        lvx  = getTokenAsFloat (line, 2); 
                        lvy  = getTokenAsFloat (line, 3); 
                        lvz  = getTokenAsFloat (line, 4); 
                        lr2  = getTokenAsFloat (line, 5); 
                        lg2  = getTokenAsFloat (line, 6); 
                        lb2  = getTokenAsFloat (line, 7);
                        
                        lig = new class lights('d', lvx, lvy, lvz, lr2, lg2, lb2);
                        lighting.push_back(lig);
                        numLights += 1;
                        
#ifdef IM_DEBUGGING
                        // if we're debugging, show what we got:
                        cout << "" << endl;
                        cout << "got a directional light with ";
                        cout << "direction: " << lvx << ", " << lvy << ", " << lvz << endl;
                        cout << "RGB values: " << lr2 << ", " << lg2 << ", " << lb2 << endl;
                        cout << "" << endl;
#endif
                        break;
                    case 'a':   // ambient light
                        float lr3, lg3, lb3;
                        
                        lr3  = getTokenAsFloat (line, 2); 
                        lg3  = getTokenAsFloat (line, 3); 
                        lb3  = getTokenAsFloat (line, 4); 
                        
                        lig = new class lights('a', 0, 0, 0, lr3, lg3, lb3);
                        lighting.push_back(lig);
                        numLights += 1;
                        
#ifdef IM_DEBUGGING
                        // if we're debugging, show what we got:
                        cout << "" << endl;
                        cout << "got an ambient light with ";
                        cout << "RGB values: " << lr3 << ", " << lg3 << ", " << lb3 << endl;
                        cout << "" << endl;
#endif
                        break;
                        
                }
                
                
                break;
                
                
                // materials:
            case 'm': 
                
                float dr, dg, db, sr, sg, sb, phong, ir, ig, ib;
                
                dr = getTokenAsFloat (line, 1); 
                dg = getTokenAsFloat (line, 2); 
                db = getTokenAsFloat (line, 3); 
                sr = getTokenAsFloat (line, 4); 
                sg = getTokenAsFloat (line, 5); 
                sb = getTokenAsFloat (line, 6);
                phong = getTokenAsFloat (line, 7); 
                ir = getTokenAsFloat (line, 8); 
                ig = getTokenAsFloat (line, 9);
                ib = getTokenAsFloat (line, 10);
                
                mat = new materials(dr, dg, db, sr, sg, sb, phong, ir, ig, ib);
                material.push_back(mat);
                lastMaterialLoaded += 1;
                //cout << lastMaterialLoaded << endl;
                
#ifdef IM_DEBUGGING
                // if we're debugging, show what we got:
                cout << "" << endl;
                cout << "got materials with ";
                cout << "diffuse components: " << dr << ", " << dg << ", " << db << endl;
                cout << "specular components: " << sr << ", " << sg << ", " << sb << endl;
                cout << "ideal spec. components: " << ir << ", " << ig << ", " << ib << endl;
                cout << "phong exponents: " << phong << endl;
#endif
                
                break;
                
                
            case '/':
                // don't do anything, it's a comment
                break;
                
                
                
                // options
            case 'o':  
                switch (line[1]) {
                    case 'b':   // obj file
                        char * file ;
                        file = (char*)getTokenAsString (line, 1).c_str();
                        //file = getTokenAsString (line, 1);
                        cout << "" << endl;
                        cout << "Boom " << getTokenAsString (line, 1) << endl ;
                        cout << "got obj file with ";
                        cout << "name: " << file << endl;
                        cout << "" << endl;
                        theFile(file);
                        int arrlength = theFile(file)[0][0];
                        cout << "Length: " <<  theFile(file)[0][0] << endl;
                        for (int i=1; i<=arrlength; i++) {
                            x1 = theFile(file)[i][0];
                            y1 = theFile(file)[i][1];
                            z1 = theFile(file)[i][2];
                            x2 = theFile(file)[i][3];
                            y2 = theFile(file)[i][4];
                            z2 = theFile(file)[i][5];
                            x3 = theFile(file)[i][6];
                            y3 = theFile(file)[i][7];
                            z3 = theFile(file)[i][8];
                            sph = new items('t', x1, y1, z1, x2, y2, z2, x3, y3, z3, 0);
                            sph->StoreMatID(lastMaterialLoaded);
                            objects.push_back(sph);
                            
                            cout << "" << endl;
                            cout << "got obj triangle with " << endl;
                            cout << "x1, y1, z1: " << x1 << "," << y1 << "," << z1 << endl;
                            cout << "x2, y2, z2: " << x2 << "," << y2 << "," << z2 << endl;
                            cout << "x3, y3, z3: " << x3 << "," << y3 << "," << z3 << endl;
                            cout << "" << endl;
                        }
                        
                        break;
                        
                }
                break;
        }
        
    }
}

float intersection(float directionx, float directiony, float directionz, float eyex, float eyey, float eyez, float pixelx, float pixely){
    int checker = 0;
    nearVal = 10000000;
    int intersects = 0;
    for (int i = 0; i<checkedObjects.size(); i++) {
        if(checkedObjects[i]->getType() == 's'){
            //e-c
            float hx = eyex - checkedObjects[i]->getX1();
            float hy = eyey - checkedObjects[i]->getY1();
            float hz = eyez - checkedObjects[i]->getZ1();
            
            float A  = (directionx*directionx) + (directiony*directiony) + (directionz*directionz); //v.v
            float B  = 2 * ((directionx*hx) + (directiony*hy) +(directionz*hz));
            float C  = ((hx*hx) + (hy*hy) +(hz*hz)) - pow(checkedObjects[i]->getD(), 2);
            
            float discriminant = pow(B, 2) - (4*A*C);
            
            if (discriminant > 0) {
                float t1  = ((-1*B) + sqrt(discriminant))/2*A;
                float t2  = ((-1*B) - sqrt(discriminant))/2*A;
                if (t1<nearVal){
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos1: " << pos << endl;
                    holdPandT[0][pos] = pixelx;
                    holdPandT[1][pos] = pixely;
                    holdPandT[2][pos] = t1;
                    holdPandT[3][pos] = directionx;
                    holdPandT[4][pos] = directiony;
                    holdPandT[5][pos] = directionz;
                    holdPandT[6][pos] = checkedObjects[i]->getMatID();
                    type[pos] = 's';
                    normCalc[0][pos] = checkedObjects[i]->getX1();
                    normCalc[1][pos] = checkedObjects[i]->getY1();
                    normCalc[2][pos] = checkedObjects[i]->getZ1();
                    normCalc[3][pos] = checkedObjects[i]->getD();
                    checker = 1;
                    nearVal = t1;
                    intersects =t1;
                }
                if (t2<nearVal){
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos1: " << pos << endl;
                    holdPandT[0][pos] = pixelx;
                    holdPandT[1][pos] = pixely;
                    holdPandT[2][pos] = t2;
                    holdPandT[3][pos] = directionx;
                    holdPandT[4][pos] = directiony;
                    holdPandT[5][pos] = directionz;
                    holdPandT[6][pos] = checkedObjects[i]->getMatID(); 
                    type[pos] = 's';
                    normCalc[0][pos] = checkedObjects[i]->getX1();
                    normCalc[1][pos] = checkedObjects[i]->getY1();
                    normCalc[2][pos] = checkedObjects[i]->getZ1();
                    normCalc[3][pos] = checkedObjects[i]->getD();
                    checker = 1;
                    nearVal = t2;
                    intersects =t2;
                }
            }
        }//end of sphere intersection checker
        else if (checkedObjects[i]->getType() == 't'){
            float vert1x, vert1y, vert1z, vert2x, vert2y, vert2z, vert3x, vert3y, vert3z;
            vert1x = checkedObjects[i]->getX1();
            vert1y = checkedObjects[i]->getY1();
            vert1z = checkedObjects[i]->getZ1();
            vert2x = checkedObjects[i]->getX2();
            vert2y = checkedObjects[i]->getY2();
            vert2z = checkedObjects[i]->getZ2();
            vert3x = checkedObjects[i]->getX3();
            vert3y = checkedObjects[i]->getY3();
            vert3z = checkedObjects[i]->getZ3();
            
            float magA = ((vert1x - vert2x) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z))))   + ((vert1x - vert3x) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
            
            float beta = ((vert1x - eyex) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) * -1 +  (directionx * (((vert1y-eyey)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-eyez)))) ;
            
            float gamma = ((vert1x - vert2x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) + ((vert1x - eyex) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) ;
            
            float t = ((vert1x - vert2x) * (((vert1y-vert3y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) * -1 +  ((vert1x - eyex) * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
            
            beta = beta/magA;
            gamma = gamma/magA;
            t = t/magA;
            double sumCheck = beta + gamma;
            
            if (sumCheck <= 1 && beta>=0 && gamma>=0 && t>0) {
                
                if (t<nearVal){
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos1: " << pos << endl;
                    holdPandT[0][pos] = pixelx;
                    holdPandT[1][pos] = pixely;
                    holdPandT[2][pos] = t;
                    holdPandT[3][pos] = directionx;
                    holdPandT[4][pos] = directiony;
                    holdPandT[5][pos] = directionz;
                    holdPandT[6][pos] = checkedObjects[i]->getMatID();
                    type[pos] = 't';
                    normCalc[0][pos] = vert1x;
                    normCalc[1][pos] = vert1y;
                    normCalc[2][pos] = vert1z;
                    normCalc[3][pos] = vert2x;
                    normCalc[4][pos] = vert2y;
                    normCalc[5][pos] = vert2z;
                    normCalc[6][pos] = vert3x;
                    normCalc[7][pos] = vert3y;
                    normCalc[8][pos] = vert3z;
                    nearVal = t;
                    intersects =t;
                }
            }
        }//end of triangle intersection checker
        else if (checkedObjects[i]->getType() == 'p'){
            //(-d-(e.norm)/dir.norm)
            float normx, normy, normz, dist;
            normx = checkedObjects[i]->getX1();
            normy = checkedObjects[i]->getY1();
            normz = checkedObjects[i]->getZ1();
            dist = checkedObjects[i]->getD();
            
            float t = -1*((dist) + (eyex*normx + eyey*normy + eyez*normz))/(directionx*normx + directiony*normy +directionz*normz);
            //cout << "plane t= " << t << endl;
            if (t>0) {
                //cout << "plane t= " << t << endl;
                if (t<nearVal){
                    
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos1: " << pos << endl;
                    holdPandT[0][pos] = pixelx;
                    holdPandT[1][pos] = pixely;
                    holdPandT[2][pos] = t;
                    holdPandT[3][pos] = directionx;
                    holdPandT[4][pos] = directiony;
                    holdPandT[5][pos] = directionz;
                    holdPandT[6][pos] = checkedObjects[i]->getMatID();
                    type[pos] = 'p';
                    normCalc[0][pos] = normx;
                    normCalc[1][pos] = normy;
                    normCalc[2][pos] = normz;
                    nearVal = t; 
                    intersects =t;
                }
            }
            //cout << "plane t value: " << t << endl;
        }//end of plane intersection checker
        
    }//end checkedObjects loop
    
    return intersects;
}



float SpecIntersection(float directionx, float directiony, float directionz, float eyex, float eyey, float eyez, float pixelx, float pixely, float ir, float ig, float ib){
    int checker = 0;
    nearVal = 10000000;
    int intersects = 0;
    for (int i = 0; i<checkedObjects.size(); i++) {
        if(checkedObjects[i]->getType() == 's'){
            //e-c
            float hx = eyex - checkedObjects[i]->getX1();
            float hy = eyey - checkedObjects[i]->getY1();
            float hz = eyez - checkedObjects[i]->getZ1();
            
            float A  = (directionx*directionx) + (directiony*directiony) + (directionz*directionz); //v.v
            float B  = 2 * ((directionx*hx) + (directiony*hy) +(directionz*hz));
            float C  = ((hx*hx) + (hy*hy) +(hz*hz)) - pow(checkedObjects[i]->getD(), 2);
            
            float discriminant = pow(B, 2) - (4*A*C);
            
            if (discriminant > 0) {
                float t1  = ((-1*B) + sqrt(discriminant))/2*A;
                float t2  = ((-1*B) - sqrt(discriminant))/2*A;
                if (t1<nearVal){
                    //cout << "Stored at " << pixelx <<  ", " << pixely << endl;
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos: " << pos << endl;
                    holdPandT[7][pos] += ir;
                    holdPandT[8][pos] += ig;
                    holdPandT[9][pos] += ib;
                    checker = 1;
                    nearVal = t1;
                    intersects =t1;
                }
                if (t2<nearVal){
                    //cout << "Stored at " << pixelx <<  ", " << pixely << endl;
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos: " << pos << endl;
                    holdPandT[7][pos] += ir;
                    holdPandT[8][pos] += ig;
                    holdPandT[9][pos] += ib;
                    checker = 1;
                    nearVal = t2;
                    intersects =t2;
                }
            }
        }//end of sphere intersection checker
        else if (checkedObjects[i]->getType() == 't'){
            float vert1x, vert1y, vert1z, vert2x, vert2y, vert2z, vert3x, vert3y, vert3z;
            vert1x = checkedObjects[i]->getX1();
            vert1y = checkedObjects[i]->getY1();
            vert1z = checkedObjects[i]->getZ1();
            vert2x = checkedObjects[i]->getX2();
            vert2y = checkedObjects[i]->getY2();
            vert2z = checkedObjects[i]->getZ2();
            vert3x = checkedObjects[i]->getX3();
            vert3y = checkedObjects[i]->getY3();
            vert3z = checkedObjects[i]->getZ3();
            
            float magA = ((vert1x - vert2x) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z))))   + ((vert1x - vert3x) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
            
            float beta = ((vert1x - eyex) * (((vert1y-vert3y)*directionz) - (directiony*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) * -1 +  (directionx * (((vert1y-eyey)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-eyez)))) ;
            
            float gamma = ((vert1x - vert2x) * (((vert1y-eyey)*directionz) - (directiony*(vert1z-eyez)))) + ((vert1x - eyex) * (((vert1y-vert2y)*directionz) - (directiony*(vert1z-vert2z)))) * -1 +  (directionx * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) ;
            
            float t = ((vert1x - vert2x) * (((vert1y-vert3y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-vert2y)*(vert1z-eyez)) - ((vert1y-eyey)*(vert1z-vert2z)))) * -1 +  ((vert1x - eyex) * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
            
            beta = beta/magA;
            gamma = gamma/magA;
            t = t/magA;
            double sumCheck = beta + gamma;
            
            if (sumCheck <= 1 && beta>=0 && gamma>=0 && t>0) {
                
                if (t<nearVal){
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos: " << pos << endl;
                    holdPandT[7][pos] += ir;
                    holdPandT[8][pos] += ig;
                    holdPandT[9][pos] += ib;
                    nearVal = t;
                    intersects =t;
                }
            }
        }//end of triangle intersection checker
        else if (checkedObjects[i]->getType() == 'p'){
            float normx, normy, normz, dist;
            normx = checkedObjects[i]->getX1();
            normy = checkedObjects[i]->getY1();
            normz = checkedObjects[i]->getZ1();
            dist = checkedObjects[i]->getD();
            
            float t = -1*((dist) + (eyex*normx + eyey*normy + eyez*normz))/(directionx*normx + directiony*normy +directionz*normz);
            //cout << "plane t= " << t << endl;
            if (t>0) {
                //cout << "plane t= " << t << endl;
                if (t<nearVal){
                    int pos = pixelx*1500 + pixely;
                    //cout << "Pos: " << pos << endl;
                    holdPandT[7][pos] += ir;
                    holdPandT[8][pos] += ig;
                    holdPandT[9][pos] += ib;
                    nearVal = t; 
                    intersects =t;
                }
            }
            //cout << "plane t value: " << t << endl;
            
        }//end of plane intersection checker
        
        float normX;
        float normY;
        float normZ; 
        if (intersects != 0 && recursionDepth<=maxDepth) {
            int pos = pixelx*1500 + pixely;
            
            int pointX = eyex + intersects* directionx;
            int pointY = eyex + intersects* directiony;
            int pointZ = eyex + intersects* directionz;
            if (type[pos] == 's') {
                //calculate normal for sphere
                float centerX = normCalc[0][pos];
                float centerY = normCalc[1][pos];
                float centerZ = normCalc[2][pos];
                float Radius = normCalc[3][pos];
                normX = (pointX - centerX)/Radius;
                normY = (pointY - centerY)/Radius;
                normZ = (pointZ - centerZ)/Radius;
            }
            else if (type[pos] == 't') {
                //calculate normal for triangle
                //cout << "ID: " << ID << endl;
                float X1 = normCalc[0][pos];
                float Y1 = normCalc[1][pos];
                float Z1 = normCalc[2][pos];
                float X2 = normCalc[3][pos];
                float Y2 = normCalc[4][pos];
                float Z2 = normCalc[5][pos];
                float X3 = normCalc[6][pos];
                float Y3 = normCalc[7][pos];
                float Z3 = normCalc[8][pos];
                float side1x = X2 - X1;
                float side1y = Y2 - Y1;
                float side1z = Z2 - Z1;
                float side2x = X3 - X1;
                float side2y = Y3 - Y1;
                float side2z = Z3 - Z1;
                //side1 x side2
                normX = (side1y*side2z) - (side1z*side2y);
                normY = (side1x*side2z) - (side1z*side2x);
                normZ = (side1x*side2y) - (side1y*side2x);
                float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                normX = normX/normMag;
                normY = normY/normMag;
                normZ = normZ/normMag;
            }
            
            else if (type[pos] == 'p') {
                //calculate normal for plane
                normX = normCalc[0][pos];
                normY = normCalc[1][pos];
                normZ = normCalc[2][pos];
                float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                normX = normX/normMag;
                normY = normY/normMag;
                normZ = normZ/normMag;
            }
            
            float dn = directionx*normX + directiony*normY + directionz*normZ;
            float reflectx = directionx - (2*dn*normX);
            float reflecty = directiony - (2*dn*normY);
            float reflectz = directionz - (2*dn*normZ);
            int ID = holdPandT[6][pos] ;   
            recursionDepth++;
            //recursion
            SpecIntersection(reflectx, reflecty, reflectz, pointX, pointY, pointZ, pixelx, pixely, material[ID]->getdr(), material[ID]->getdg(), material[ID]->getdb());
        }
    }//end checkedObjects loop
    
    return intersects;
}
//handle ideal specular reflection
void ideal_reflection(){
    //compute ray direction for individual pixel
    float nwx, nwy, nwz;
    float nwMag = sqrt(pow(wx, 2) + pow(wy, 2) + pow(wz, 2));
    nwx = -1*(wx/nwMag);
    nwy = -1*(wy/nwMag);
    nwz = -1*(wz/nwMag);
    ux = -1* nwz;
    uy = 0;
    uz = 1*(nwx);
    uMag = sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
    ux = ux/uMag;
    uy = uy/uMag;
    uz = uz/uMag;
    
    // v = u x w
    vx = (uy*nwz) - (uz*nwy);
    vy = (uz*nwx) - (ux*nwz);
    vz = (ux*nwy) - (uy*nwx);
    
    
    float widthE = iw/2;
    float heightE = ih/2;
    float widthB = widthE * -1;
    float heightB = heightE * -1;
    float du = iw/(pw*2);
    float dv = ih/(ph*2);
    int pixelx = -1;
    int pixely = -1;
    for(double u=widthB; u<widthE; u+=du){
        pixelx +=1;
        pixely =-1;
        //u += ju;
        for(double v=heightB; v<heightE; v+=dv){
            //adds/subtracts random value to/from the subpixel coodinates to create the jitter sample
            
            pixely +=1;
            double directionx = cd*wx + (u)*ux + (v)*vx;
            double directiony = cd*wy + (u)*uy + (v)*vy;
            double directionz = cd*wz + (u)*uz + (v)*vz;
            double dMag = sqrt(pow(directionx, 2) + pow(directiony, 2) + pow(directionz, 2));
            
            directionx = directionx/dMag;
            directiony = directiony/dMag;
            directionz = directionz/dMag;
            
            
            float intersecter =  intersection(directionx, directiony, directionz, eyex, eyey, eyez, pixelx, pixely);
            
            if (intersecter != 0) {
                int pos = pixelx*1500 + pixely;
                //int pixelx  = holdPandT[0][pos] ;
                //int pixely  = holdPandT[1][pos] ;
                float t     = holdPandT[2][pos] ;
                float dirx  = holdPandT[3][pos] ;
                float diry  = holdPandT[4][pos] ;
                float dirz  = holdPandT[5][pos] ;
                
                float px= eyex + t*dirx;
                float py= eyey + t*diry;
                float pz= eyez + t*dirz;
                float ambR = 0;
                float ambG = 0;
                float ambB = 0;
                int ID = (int)holdPandT[6][pos] ;
                
                Rgba &pix = p[pixelx][pixely]; 
                pix.r = 0;
                pix.g = 0;
                pix.b = 0;
                pix.a = 1;
                
                float normX;
                float normY;
                float normZ;
                float dr, dg, db, sr, sg, sb, phong, ir, ig, ib;
                //cout << pixelx << ", " << pixely << endl; 
                dr = material[ID]->getdr();
                dg = material[ID]->getdg();
                db = material[ID]->getdb();
                sr = material[ID]->getsr();
                sg = material[ID]->getsg();
                sb = material[ID]->getsb();
                ir = material[ID]->getir();
                ig = material[ID]->getig();
                ib = material[ID]->getib();
                phong = material[ID]->getr();
                if (type[pos] == 's') {
                    //calculate normal for sphere
                    float centerX = normCalc[0][pos];
                    float centerY = normCalc[1][pos];
                    float centerZ = normCalc[2][pos];
                    float Radius = normCalc[3][pos];
                    normX = (px - centerX)/Radius;
                    normY = (py - centerY)/Radius;
                    normZ = (pz - centerZ)/Radius;
                }
                else if (type[pos] == 't') {
                    //calculate normal for triangle
                    //cout << "ID: " << ID << endl;
                    float X1 = normCalc[0][pos];
                    float Y1 = normCalc[1][pos];
                    float Z1 = normCalc[2][pos];
                    float X2 = normCalc[3][pos];
                    float Y2 = normCalc[4][pos];
                    float Z2 = normCalc[5][pos];
                    float X3 = normCalc[6][pos];
                    float Y3 = normCalc[7][pos];
                    float Z3 = normCalc[8][pos];
                    float side1x = X2 - X1;
                    float side1y = Y2 - Y1;
                    float side1z = Z2 - Z1;
                    float side2x = X3 - X1;
                    float side2y = Y3 - Y1;
                    float side2z = Z3 - Z1;
                    //side1 x side2
                    normX = (side1y*side2z) - (side1z*side2y);
                    normY = (side1x*side2z) - (side1z*side2x);
                    normZ = (side1x*side2y) - (side1y*side2x);
                    float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                    normX = normX/normMag;
                    normY = normY/normMag;
                    normZ = normZ/normMag;
                }
                else if (type[pos] == 'p') {
                    //calculate normal for plane
                    normX = normCalc[0][pos];
                    normY = normCalc[1][pos];
                    normZ = normCalc[2][pos];
                    float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                    normX = normX/normMag;
                    normY = normY/normMag;
                    normZ = normZ/normMag;
                }
                
                float dn = directionx*normX + directiony*normY + directionz*normZ;
                float reflectx = directionx - (2*dn*normX);
                float reflecty = directiony - (2*dn*normY);
                float reflectz = directionz - (2*dn*normZ);
                float magReflect = sqrt(pow(reflectx, 2) + pow(reflecty, 2) + pow(reflectz, 2));
                reflectx = reflectx/magReflect;
                reflecty = reflecty/magReflect;
                reflectz = reflectz/magReflect;
                
                
                for (int i = 0; i<numLights; i++){
                    float lx;
                    float ly;
                    float lz;
                    if (lighting[i]->getTy() == 'p') {
                        //calculate the direction rays for a point light
                        lx = lighting[i]->getX() - px;
                        ly = lighting[i]->getY() - py;
                        lz = lighting[i]->getZ() - pz;
                        float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        lx = lx/lMag;
                        ly = ly/lMag;
                        lz = lz/lMag;
                    }//if point light
                    if (lighting[i]->getTy() == 'd'){
                        lx = -1*lighting[i]->getX();
                        ly = -1*lighting[i]->getY();
                        lz = -1*lighting[i]->getZ();
                        float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        lx = lx/lMag;
                        ly = ly/lMag;
                        lz = lz/lMag;
                    }//if directional light 
                    if (lighting[i]->getTy() == 'p' || lighting[i]->getTy() == 'd'){ //common code
                        float IR = lighting[i]->getR();
                        float IG = lighting[i]->getG();
                        float IB = lighting[i]->getB();
                        float hx = lx + vx;
                        float hy = lx + vy;
                        float hz = lx + vz;
                        //float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        float hMag = sqrt(pow(hx, 2) + pow(hy, 2) + pow(hz, 2));
                        //float vMag = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
                        //float nMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                        
                        
                        hx = hx/hMag;
                        hy = hy/hMag;
                        hz = hz/hMag;
                        //cout << "Phong:" << phong << endl;
                        //hMag = sqrt(pow(hx, 2) + pow(hy, 2) + pow(hz, 2));
                        //cout << "Magnitude Test: " << lMag << ", " << nMag << ", " << vMag << ", " << hMag <<endl;
                        float nl = normX*lx +normY*ly + normZ*lz;
                        if (0 > nl){
                            nl = 0;
                        }
                        float nh = normX*hx +normY*hy + normZ*hz;
                        if (0 > nh){
                            nh = 0;
                        }
                        float diffuseR = nl*dr;
                        float diffuseG = nl*dg;
                        float diffuseB = nl*db;
                        
                        float specularR = pow(nh, phong)*sr;
                        float specularG = pow(nh, phong)*sg;
                        float specularB = pow(nh, phong)*sb;
                        
                        pix.r += (diffuseR + specularR)*IR;
                        pix.g += (diffuseG + specularG)*IG;
                        pix.b += (diffuseB + specularB)*IB;
                    }
                    
                    
                    if (lighting[i]->getTy() == 'a'){
                        ambR = lighting[i]->getR()*dr;
                        ambG = lighting[i]->getG()*dg;
                        ambB = lighting[i]->getB()*db;
                        
                        pix.r += ambR;
                        pix.g += ambG;
                        pix.b += ambB;
                    }//if ambient light 
                }//for all lights on a pixel
                //cout << pix.r << ", " << pix.b << ", " << pix.g << endl;
                float ff = 0;
                SpecIntersection(reflectx, reflecty, reflectz, px, py, pz, pixelx, pixely, pix.r, pix.g, pix.b);
                
                if (ff != 0) {
                    //cout << "Boom!" << endl;
                }
            }//end of if intersects
            
        }//end of v loop
    }//end of u loop
    
    
}

void draw(){
    
    //compute ray direction for individual pixel
    float nwx, nwy, nwz;
    float nwMag = sqrt(pow(wx, 2) + pow(wy, 2) + pow(wz, 2));
    nwx = -1*(wx/nwMag);
    nwy = -1*(wy/nwMag);
    nwz = -1*(wz/nwMag);
    ux = -1* nwz;
    uy = 0;
    uz = 1*(nwx);
    uMag = sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
    ux = ux/uMag;
    uy = uy/uMag;
    uz = uz/uMag;
    
    // v = u x w
    vx = (uy*nwz) - (uz*nwy);
    vy = (uz*nwx) - (ux*nwz);
    vz = (ux*nwy) - (uy*nwx);
    
    
    float widthE = iw/2;
    float heightE = ih/2;
    float widthB = widthE * -1;
    float heightB = heightE * -1;
    float du = iw/(pw*2);
    float dv = ih/(ph*2);
    int pixelx = -1;
    int pixely = -1;
    for(double u=widthB; u<widthE; u+=du){
        pixelx +=1;
        pixely =-1;
        srand ((int) time(NULL) );
        int jitterU = rand() % 10;
        double ju = jitterU;
        while (ju > (du)/2) {
            ju = ju/2;
        }
        //random way to choose whether to add or subtract jitter value
        int posORneg = rand() % 2;
        if (posORneg) {
            ju = -1*ju;
        }
        //u += ju;
        for(double v=heightB; v<heightE; v+=dv){
            //adds/subtracts random value to/from the subpixel coodinates to create the jitter sample
            
            int jitterV = rand() % 10;
            double jv = jitterV;
            while (jv > (dv)/2) {
                jv = jv/2;
            }
            
            //random way to choose whether to add or subtract jitter value
            int posORneg = rand() % 2;
            //cout << "Pos or Neg " << posORneg << endl;
            if (posORneg) {
                jv = -1*jv;
            }
            //v += jv;
            //cout << ju << ", " << jv << endl;
            pixely +=1;
            double directionx = cd*wx + (u+ju)*ux + (v+jv)*vx;
            double directiony = cd*wy + (u+ju)*uy + (v+jv)*vy;
            double directionz = cd*wz + (u+ju)*uz + (v+jv)*vz;
            double dMag = sqrt(pow(directionx, 2) + pow(directiony, 2) + pow(directionz, 2));
            
            directionx = directionx/dMag;
            directiony = directiony/dMag;
            directionz = directionz/dMag;
            
            intersection(directionx, directiony, directionz, eyex, eyey, eyez, pixelx, pixely);
            
        }//end of v loop
    }//end of u loop
    
    
    
}


void paint(){
    for (int i=0; i<pw*2; i++) {
        for (int j=0; j<ph*2; j++) {
            int pos = i*1500 + j;
            //cout << pos << endl;
            if (holdPandT[2][pos] != 0) {
                int pixelx  = holdPandT[0][pos] ;
                int pixely  = holdPandT[1][pos] ;
                float t     = holdPandT[2][pos] ;
                float dirx  = holdPandT[3][pos] ;
                float diry  = holdPandT[4][pos] ;
                float dirz  = holdPandT[5][pos] ;
                float vx = -1*dirx;
                float vy = -1*diry;
                float vz = -1*dirz;
                int ID      = (int)holdPandT[6][pos] ;
                float reflectedR      = holdPandT[7][pos] ;
                float reflectedG      = holdPandT[8][pos] ;
                float reflectedB      = holdPandT[9][pos] ;
                if (reflectedB+reflectedG+reflectedR != 0) {
                    //cout << "Reflected Components: " <<  reflectedR << ", " << reflectedG << ", " << reflectedB << endl;
                }
                
                float ambR = 0;
                float ambG = 0;
                float ambB = 0;
                float dr, dg, db, sr, sg, sb, phong, ir, ig, ib;
                //cout << pixelx << ", " << pixely << endl; 
                
                dr = material[ID]->getdr();
                dg = material[ID]->getdg();
                db = material[ID]->getdb();
                sr = material[ID]->getsr();
                sg = material[ID]->getsg();
                sb = material[ID]->getsb();
                ir = material[ID]->getir();
                ig = material[ID]->getig();
                ib = material[ID]->getib();
                phong = material[ID]->getr();
                
                float RspecR = reflectedR *ir;
                float RspecG = reflectedG *ig;
                float RspecB = reflectedB *ib;
                /*
                 if (specularG != 0) {
                 //cout << ig << ", "  << specularG << ", " << reflectedG << endl;
                 }
                 */
                float px= eyex + t*dirx;
                float py= eyey + t*diry;
                float pz= eyez + t*dirz;
                
                Rgba &pix = p[pixelx][pixely]; 
                pix.r = RspecR;
                pix.g = RspecG;
                pix.b = RspecB;
                pix.a = 1;
                
                float normX;
                float normY;
                float normZ;
                if (type[pos] == 's') {
                    //calculate normal for sphere
                    float centerX = normCalc[0][pos];
                    float centerY = normCalc[1][pos];
                    float centerZ = normCalc[2][pos];
                    float Radius = normCalc[3][pos];
                    normX = (px - centerX)/Radius;
                    normY = (py - centerY)/Radius;
                    normZ = (pz - centerZ)/Radius;
                }
                else if (type[pos] == 't') {
                    //calculate normal for triangle
                    //cout << "ID: " << ID << endl;
                    float X1 = normCalc[0][pos];
                    float Y1 = normCalc[1][pos];
                    float Z1 = normCalc[2][pos];
                    float X2 = normCalc[3][pos];
                    float Y2 = normCalc[4][pos];
                    float Z2 = normCalc[5][pos];
                    float X3 = normCalc[6][pos];
                    float Y3 = normCalc[7][pos];
                    float Z3 = normCalc[8][pos];
                    float side1x = X2 - X1;
                    float side1y = Y2 - Y1;
                    float side1z = Z2 - Z1;
                    float side2x = X3 - X2;
                    float side2y = Y3 - Y2;
                    float side2z = Z3 - Z2;
                    //side1 x side2
                    normX = (side1y*side2z) - (side1z*side2y);
                    normY = (side1x*side2z) - (side1z*side2x);
                    normZ = (side1x*side2y) - (side1y*side2x);
                    float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                    normX = normX/normMag;
                    normY = normY/normMag;
                    normZ = normZ/normMag;
                }
                else if (type[pos] == 'p') {
                    //calculate normal for plane
                    normX = normCalc[0][pos];
                    normY = normCalc[1][pos];
                    normZ = normCalc[2][pos];
                    float normMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                    normX = normX/normMag;
                    normY = normY/normMag;
                    normZ = normZ/normMag;
                }
                
                for (int k = 0; k<numLights; k++){
                    float lx;
                    float ly;
                    float lz;
                    if (lighting[k]->getTy() == 'p') {
                        //calculate the direction rays for a point light
                        lx = lighting[k]->getX() - px;
                        ly = lighting[k]->getY() - py;
                        lz = lighting[k]->getZ() - pz;
                        float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        lx = lx/lMag;
                        ly = ly/lMag;
                        lz = lz/lMag;
                    }//if point light
                    if (lighting[k]->getTy() == 'd'){
                        lx = -1*lighting[k]->getX();
                        ly = -1*lighting[k]->getY();
                        lz = -1*lighting[k]->getZ();
                        float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        lx = lx/lMag;
                        ly = ly/lMag;
                        lz = lz/lMag;
                    }//if directional light 
                    if (lighting[k]->getTy() == 'p' || lighting[k]->getTy() == 'd'){ //common code
                        float IR = lighting[k]->getR();
                        float IG = lighting[k]->getG();
                        float IB = lighting[k]->getB();
                        if (IR + IG + IB == 0) {
                            //cout << "Pow!" << endl;
                        }
                        float hx = lx + vx;
                        float hy = lx + vy;
                        float hz = lx + vz;
                        //float lMag = sqrt(pow(lx, 2) + pow(ly, 2) + pow(lz, 2));
                        float hMag = sqrt(pow(hx, 2) + pow(hy, 2) + pow(hz, 2));
                        //float vMag = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
                        //float nMag = sqrt(pow(normX, 2) + pow(normY, 2) + pow(normZ, 2));
                        
                        
                        hx = hx/hMag;
                        hy = hy/hMag;
                        hz = hz/hMag;
                        //cout << "Phong:" << phong << endl;
                        //hMag = sqrt(pow(hx, 2) + pow(hy, 2) + pow(hz, 2));
                        //cout << "Magnitude Test: " << lMag << ", " << nMag << ", " << vMag << ", " << hMag <<endl;
                        float nl = normX*lx +normY*ly + normZ*lz;
                        float nh = normX*hx +normY*hy + normZ*hz;
                        if (0 > nl){
                            nl = -1*nl;
                            
                        }
                        
                        if (0 > nh){
                            nh = 0;
                        }
                        float diffuseR = nl*dr;
                        float diffuseG = nl*dg;
                        float diffuseB = nl*db;
                        
                        float specularR = pow(nh, phong)*sr;
                        float specularG = pow(nh, phong)*sg;
                        float specularB = pow(nh, phong)*sb;
                        if (k == 0) {
                            pix.r += (RspecR)*IR;
                            pix.g += (RspecG)*IG;
                            pix.b += (RspecB)*IB;
                        }
                        pix.r += (diffuseR + specularR)*IR;
                        pix.g += (diffuseG + specularG)*IG;
                        pix.b += (diffuseB + specularB)*IB;
                        if (diffuseR + diffuseG + diffuseB == 0) {
                            cout << "Pow!" << endl;
                            cout << "nl, dr " << nl << ", " << dr << endl;
                        }
                    }
                    
                    
                    if (lighting[k]->getTy() == 'a'){
                        ambR = lighting[k]->getR()*dr;
                        ambG = lighting[k]->getG()*dg;
                        ambB = lighting[k]->getB()*db;
                        
                        pix.r += ambR*dr;
                        pix.g += ambG*dg;
                        pix.b += ambB*db;
                    }//if ambient light 
                    
                    for (int i = 0; i<checkedObjects.size(); i++) {
                        if((char)checkedObjects[i]->getType() == 's'){
                            //e-c
                            float hx = px - checkedObjects[i]->getX1();
                            float hy = py - checkedObjects[i]->getY1();
                            float hz = pz - checkedObjects[i]->getZ1();
                            
                            float A  = (lx*lx) + (ly*ly) + (lz*lz); //v.v
                            float B  = 2 * ((lx*hx) + (ly*hy) +(lz*hz));
                            float C  = ((hx*hx) + (hy*hy) +(hz*hz)) - pow(checkedObjects[i]->getD(), 2);
                            
                            float discriminant = pow(B, 2) - (4*A*C);
                            
                            if (discriminant > 0 &&  type[pos] == 'p') {
                                //if there is an intersection, set shadow to ambient light color
                                pix.r = ambR*dr;
                                pix.g = ambG*dg;
                                pix.b = ambB*db;
                                
                                //cout << "Pixel X = " << pixelx  << ", Pixel Y = " << pixely << endl;
                            }
                        }//end of sphere intersection checker
                        else if ((char)checkedObjects[i]->getType() == 't'){
                            float vert1x, vert1y, vert1z, vert2x, vert2y, vert2z, vert3x, vert3y, vert3z;
                            vert1x = checkedObjects[i]->getX1();
                            vert1y = checkedObjects[i]->getY1();
                            vert1z = checkedObjects[i]->getZ1();
                            vert2x = checkedObjects[i]->getX2();
                            vert2y = checkedObjects[i]->getY2();
                            vert2z = checkedObjects[i]->getZ2();
                            vert3x = checkedObjects[i]->getX3();
                            vert3y = checkedObjects[i]->getY3();
                            vert3z = checkedObjects[i]->getZ3();
                            
                            float magA = ((vert1x - vert2x) * (((vert1y-vert3y)*lz) - (ly*(vert1z-vert3z))))   + ((vert1x - vert3x) * (((vert1y-vert2y)*lz) - (ly*(vert1z-vert2z)))) * -1 +  (lx * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
                            
                            float beta = ((vert1x - px) * (((vert1y-vert3y)*lz) - (ly*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-py)*lz) - (ly*(vert1z-pz)))) * -1 +  (lx * (((vert1y-py)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-pz)))) ;
                            
                            float gamma = ((vert1x - vert2x) * (((vert1y-py)*lz) - (ly*(vert1z-pz)))) + ((vert1x - px) * (((vert1y-vert2y)*lz) - (ly*(vert1z-vert2z)))) * -1 +  (lx * (((vert1y-vert2y)*(vert1z-pz)) - ((vert1y-py)*(vert1z-vert2z)))) ;
                            
                            float t = ((vert1x - vert2x) * (((vert1y-vert3y)*(vert1z-pz)) - ((vert1y-py)*(vert1z-vert3z)))) + ((vert1x - vert3x) * (((vert1y-vert2y)*(vert1z-pz)) - ((vert1y-py)*(vert1z-vert2z)))) * -1 +  ((vert1x - px) * (((vert1y-vert2y)*(vert1z-vert3z)) - ((vert1y-vert3y)*(vert1z-vert2z)))) ;
                            
                            
                            beta = beta/magA;
                            gamma = gamma/magA;
                            double sumCheck = beta + gamma;
                            
                            if (sumCheck < 1 && beta>0 && gamma>0 && t>0  && type[pos] == 'p') {
                                //if there is an intersection, set shadow to ambient light color 
                                //cout << "Pixel X = " << pixelx  << ", Pixel Y = " << pixely << endl;
                                pix.r = ambR*dr;
                                pix.g = ambG*dg;
                                pix.b = ambB*db;
                            }
                        }//end of triangle intersection checker
                    }//end checkedObjects loop
                    
                    
                    
                }//for all lights on a pixel
                if (pix.r + pix.g + pix.b == 0) {
                    //cout << "Boom!" << endl;
                }
            }//if t!=0
        }//inner for loop
    }//outer for loop
}// paint function

int main (int argc, char *argv[])
{
    
    if (argc != 2) {
        cerr << "useage: raytra scenefilename" << endl;
        return -1;
    }
    
    parseSceneFile (argv[1]);
    cout << "Performing boundary checking... " << endl;
    node(objects, 0);
    cout << "BVH Structure found " << leafCounter << " leaf intersection" << endl;
    cout << "Please wait ... Rendering image ..." << endl;
    ideal_reflection();
    draw();
    
    cout << "Still rendering ... Ray intersections are quite a pain to render" << endl;
    paint();
    int j = 0;
    
    /*
     image we created above was 2ph by 2pw for anti-aliasing purposes. 
     We now average every 4 pixels into one to perform aliasing 
     and to obtain a sharper/clearer ph by pw image
     
     Here we put the values into in clearer average from which we can then
     easily average our consecutive values
     */
    float averager[4][3];
    int position = 0;
    int pixelx = 0;
    int pixely = 0;
    for (int i = 0; i<pw*2; ) {
        for ( ; j<ph*2; ) {
            if (j == ph*2-1 && i == pw*2-1) {
                //end loop ... average up pixels in our array
                Rgba &pix = p[i][j]; 
                averager[position][0] = pix.r;
                averager[position][1] = pix.g;
                averager[position][2] = pix.b;
                
                Rgba &avg = p[pixelx][pixely]; 
                avg.r = 0;
                avg.g = 0;
                avg.b = 0;
                
                for (int k = 0; k<4 ; k++) {
                    avg.r += averager[k][0];
                    avg.g += averager[k][1];
                    avg.b += averager[k][2];
                }
                avg.r = avg.r/4;
                avg.g = avg.g/4;
                avg.b = avg.b/4;
                
                if (pixely == ph-1) {
                    pixelx++;
                    pixely = 0;
                }else{
                    pixely++;
                }
                //cout << "Position: " << position << endl;
                position = 0;
                j++;
                i++;
            }
            else if (j == (ph*2)-1 && i%2 != 0) {
                //make j=0 and add 1 to i ... average up pixels in our array
                Rgba &pix = p[i][j]; 
                averager[position][0] = pix.r;
                averager[position][1] = pix.g;
                averager[position][2] = pix.b;
                
                Rgba &avg = p[pixelx][pixely]; 
                avg.r = 0;
                avg.g = 0;
                avg.b = 0;
                
                for (int k = 0; k<4 ; k++) {
                    avg.r += averager[k][0];
                    avg.g += averager[k][1];
                    avg.b += averager[k][2];
                }
                avg.r = avg.r/4;
                avg.g = avg.g/4;
                avg.b = avg.b/4;
                //problem
                //cout << "Position: " << position << endl;
                if (pixely == ph-1) {
                    pixelx++;
                    pixely = 0;
                }else{
                    pixely++;
                }
                position = 0;
                
                j = 0;
                i = i +1;
                
            }
            else if (j%2 != 0 && i%2 != 0) {
                //subtract 1 from i - Moving to next pixel quadriple ... average up pixels in our array
                Rgba &pix = p[i][j]; 
                averager[position][0] = pix.r;
                averager[position][1] = pix.g;
                averager[position][2] = pix.b;
                //cout << pixelx << ", " << pixely << endl;
                Rgba &avg = p[pixelx][pixely]; 
                avg.r = 0;
                avg.g = 0;
                avg.b = 0;
                
                for (int k = 0; k<4 ; k++) {
                    avg.r += averager[k][0];
                    avg.g += averager[k][1];
                    avg.b += averager[k][2];
                }
                avg.r = avg.r/4;
                avg.g = avg.g/4;
                avg.b = avg.b/4;
                
                if (pixely == ph-1) {
                    pixelx++;
                    pixely = 0;
                }else{
                    pixely++;
                }
                //cout << "Position: " << position << endl;
                position = 0;
                
                i = i - 1;
                j++;
            }
            else if (j%2 != 0) {
                //subtract 2 from j, add 1 to i (move 2 in sequence)
                Rgba &pix = p[i][j]; 
                averager[position][0] = pix.r;
                averager[position][1] = pix.g;
                averager[position][2] = pix.b;
                position++;
                j = j - 1;
                i++;
            }else{//automatic first and third move, no turning point
                Rgba &pix = p[i][j]; 
                averager[position][0] = pix.r;
                averager[position][1] = pix.g;
                averager[position][2] = pix.b;
                j++;
                position++;
            }
        }
    }
    
    //zero out the alpha channels in other pixels of the formerly enlarge imaged.
    for (int i = 0; i<1500; i++) {
        for (int j = ph; j<1500; j++){
            Rgba &avg = p[i][j]; 
            avg.a = 0;
        }
    }
    for (int i = pw; i<1500; i++) {
        for (int j = 0; j<1500; j++){
            Rgba &avg = p[i][j]; 
            avg.a = 0;
        }
    }
    
    for (int i = 0; i<1500; i++) {
        for (int j = 0; j<1500; j++){
            Rgba &avg = p[i][j]; 
            Rgba &bw = Sepia[i][j];
            Rgba &bw1 = BandW[i][j];
            Rgba &bw2 = WildBlue[i][j]; 
            Rgba &bw3 = Pixelated[i][j]; 
            Rgba &bw4 = TVChannel[i][j]; 
            bw.a = avg.a;
            bw1.a = avg.a;
            bw2.a = avg.a;
            bw3.a = avg.a;
            bw4.a = avg.a;
        }
    }
    
    for (int i = 0; i<pw; i++) {
        for (int j = 1; j<ph; j++) {
            Rgba &bw = Pixelated[i][j]; 
            Rgba &avg = p[i][j]; 
            
            bw.r = avg.r;
            bw.g = avg.g;
            bw.b = avg.b;
            
        }
    }
     
    for (int i = 0; i<pw; i++) {
        for (int j = 1; j<ph; j++) {
            Rgba &bw = Sepia[i][j]; 
            Rgba &avg = p[i][j]; 
            float sum = avg.r*.3 + avg.g*.59 + avg.b*.11;
            
            bw.r = sum * 0.44;
            bw.g = sum * 0.26;
            bw.b = sum * 0.08;
            
        }
    }
    
    for (int i = 0; i<pw; i++) {
        for (int j = 1; j<ph; j++) {
            Rgba &bw = BandW[i][j]; 
            Rgba &avg = p[i][j]; 
            float sum = avg.r*.15 + avg.g*.295 + avg.b*.055;
            
            bw.r = sum;
            bw.g = sum;
            bw.b = sum;
        }
    }
    
    for (int i = 0; i<pw; i+=2) {
        for (int j = 1; j<ph; j+=2) {
            Rgba &bw = TVChannel[i][j]; 
            Rgba &avg = p[i][j]; 
            float sum = avg.r*.3 + avg.g*.59 + avg.b*.11;
            
            bw.r = sum * .54;
            bw.g = sum * .58;
            bw.b = sum * .90;
        }
    }
    
    for (int i = 0; i<pw; i++) {
        for (int j = 1; j<ph; j++) {
            Rgba &bw = WildBlue[i][j]; 
            Rgba &avg = p[i][j]; 
            float sum = avg.r*.3 + avg.g*.59 + avg.b*.11;
            
            bw.r = sum * .54;
            bw.g = sum * .58;
            bw.b = sum * .90;
        }
    }


    for (int i = 0; i<pw; i+=2) {
        for (int j = 1; j<ph; j+=2) {
            Rgba &avg = Pixelated[i][j]; 
            
            avg.a = 0;
        }
    }
    
    writeRgba ("output.exr", &p[0][0], 1500, 1500);
    cout << "The Color image has been stored as color.exr" << endl;
    
    writeRgba ("sepia.exr", &Sepia[0][0], 1500, 1500);
    cout << "The Sephia image has been stored as sepia.exr" << endl;
    
    writeRgba ("blackandwhite.exr", &BandW[0][0], 1500, 1500);
    cout << "The Black and White image has been stored as blackandwhite.exr" << endl;
    
    writeRgba ("wildBlue.exr", &WildBlue[0][0], 1500, 1500);
    cout << "The Wild Blue image has been stored as wildBlue.exr" << endl;
    
    writeRgba ("TVChannel.exr", &TVChannel[0][0], 1500, 1500);
    cout << "The 'TV Channel' image has been stored as TVChannel.exr" << endl;
    
    writeRgba ("Pixelated.exr", &Pixelated[0][0], 1500, 1500);
    cout << "The Pixelated image has been stored as Pixelated.exr" << endl;
    
    for (int pixelx = 0; pixelx < pw*2; pixelx++) {
        for (int pixely=0; pixely < ph*2; pixely++) {
            int pos = pixelx*1500 + pixely;
            
            if (holdPandT[7][pos] != 0) {
                //cout << "Hold P and T: " << holdPandT[0][pos] << ", " << holdPandT[1][pos] << ", " << holdPandT[7][pos] << endl;
            }
        }
    }
    
    return 0;
}
