#define MAXNUM 120000

float** theFile(char * fName){

    FILE *fp;
    float v[MAXNUM][3];
    float vn[MAXNUM][3];
    //float triangles[MAXNUM][9];
    float **triangles = 0;
    
    triangles = new float*[MAXNUM];
    for (int i = 0; i < MAXNUM; ++i) {
        triangles[i] = new float[9];
    }
    

    int vertexCount = 0;
    int vertexNcount = 0;
    int triangleCount = 0;
    int pos;
//cout << "test: "  << fName <<endl;

//cout << "test" << endl;        
    if((fp = fopen(fName, "r"))==NULL)
    {
        printf("Cannot open file!");
    }
	else{
//		cout << "Found file " << fName << endl;
	}

//	cout << "test" << endl;    
    char buffer[MAXNUM];
    char temp1[MAXNUM];
    while(!feof(fp))
    {
//	cout << "test: "  << fName <<endl;
        char  hold[100];
        memset(buffer,0, MAXNUM-1);
        fgets(buffer, MAXNUM, fp);
        //printf("%s\n", buffer);
        if( strncmp("vn ", buffer, 3) == 0)
        {
            //normal
            
            sscanf(buffer, "%s %f %f %f", hold, vn[vertexNcount][0], vn[vertexNcount][1], vn[vertexNcount][2]);                  
            vertexCount++;
        }
        else if(strncmp("v ",buffer, 2) == 0)
        {
            //vertex
            float hold1, hold2, hold3;
            //cout << vertexCount << ", Max: " << sizeof(v)/sizeof(float) << "\nBuffer: " << buffer<< endl;
            //sscanf(buffer, "%s %f %f %f", hold, v[vertexCount][0], v[vertexCount][1], v[vertexCount][2]);
            sscanf(buffer, "%s %f %f %f", hold, &hold1, &hold2, &hold3);
            v[vertexCount][0] = hold1;
            v[vertexCount][1] = hold2;
            v[vertexCount][2] = hold3;
            vertexCount++;
        }
        else if(strncmp("f ", buffer, 2) == 0)
        {
            //faces
            int num_words = 0;
            char delims[] = " ";
            char *result = NULL;
            strncpy(temp1, buffer, MAXNUM);
            result = strtok( temp1, " " );
            int check = 0;
            int first_vertex = 0;
            int third_vertex = 0;
            
            int first_norm = 0;
            int third_norm = 0;
            while( result != NULL && strlen(result)!=0) {
                if(num_words==0){
                    //do nothing
                }
                if(num_words==1){
                    
                    //t[triangleCount].v1 = atoi(result);
                    pos = atoi(result)-1;
                    triangles[triangleCount+1][0] = v[pos][0];
                    triangles[triangleCount+1][1] = v[pos][1];
                    triangles[triangleCount+1][2] = v[pos][2];
                }
                if(num_words==2 && check==0){
                    pos = atoi(result)-1;
                    triangles[triangleCount+1][3] = v[pos][0];
                    triangles[triangleCount+1][4] = v[pos][1];
                    triangles[triangleCount+1][5] = v[pos][2];
                }
                if(num_words==3 && check==0){
                    pos = atoi(result)-1;
                    triangles[triangleCount+1][6] = v[pos][0];
                    triangles[triangleCount+1][7] = v[pos][1];
                    triangles[triangleCount+1][8] = v[pos][2];
                }
                
                if(num_words==3 && check==1){
                    pos = atoi(result)-1;
                    triangles[triangleCount+1][0] = triangles[triangleCount][0];
                    triangles[triangleCount+1][1] = triangles[triangleCount][1];
                    triangles[triangleCount+1][2] = triangles[triangleCount][2];
                    triangles[triangleCount+1][3] = triangles[triangleCount][6];
                    triangles[triangleCount+1][4] = triangles[triangleCount][7];
                    triangles[triangleCount+1][5] = triangles[triangleCount][8];
                    triangles[triangleCount+1][6] = v[pos][0];
                    triangles[triangleCount+1][7] = v[pos][1];
                    triangles[triangleCount+1][8] = v[pos][2];
                }
                result = strtok( NULL, " " );
                
                num_words++;
                
                if(num_words==4){
                    
                    check = 1;
                    num_words = 3;
                    triangleCount++;
                }
            }
        }
        
    }
    fclose(fp);
    
    //cout << "Triangle Count : " << triangleCount << endl;
    triangles[0][0] = triangleCount;
    /*
    for (int i = 1; i<= triangleCount; i++) {
        cout << "Triangle Vert1: " << triangles[i][0] << ", " << triangles[i][3] <<  ", " << triangles[i][6] <<  endl;
    } */
    
    return triangles;
 }
