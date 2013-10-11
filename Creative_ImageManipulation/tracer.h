#include "objects.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <vector>
#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <math.h>
#include <iostream>

using namespace std;
using namespace Imf;
using namespace Imath;


void writeRgba (const char fileName[],
                const Rgba *pixels,
                int width,
                int height);

void
readRgba (const char fileName[],
          Array2D<Rgba> &pixels,
          int &width,
          int &height);

double processPoint(double ** matrix, int x, int y, double ** kernel, int direction);

double* scale(double * pixels,int w1,int h1, double scale);

float getTokenAsFloat (string inString, int whichToken);

string getTokenAsString (string inString, int whichToken);

void draw();

void parseSceneFile (char *filnam);

void node(vector<items *> objects, unsigned int treedepth) ;

float intersection(float directionx, float directiony, float directionz, float eyex, float eyey, float eyez, float pixelx, float pixely);

void paint();

float BoxIntersection(vector<items *> box);