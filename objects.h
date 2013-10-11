

class Points {
public:
    float xval, yval, zval;
    void setPoints(float x, float y, float z){
        xval = x;
        yval = y;
        zval = z;
    }
    
    int getX(){
        return xval;
    }
    int getY(){
        return yval;
    }
    int getZ(){
        return zval;
    }
    
};


class lights {
public:
    float r, g, b, x, y, z;
    char ty;
    
    lights( char type,  float parx, float pary, float parz,
                        float parr, float parg, float parb){
        ty = type;
        r  = parr;
        g  = parg;
        b  = parb;
        x  = parx;
        y  = pary;
        z  = parz;
    }   
    
    char getTy(){
        return ty;
    }
    float getX(){
        return x;
    }
    float getY(){
        return y;
    }
    float getZ(){
        return z;
    }
    float getR(){
        return r;
    }
    float getG(){
        return g;
    }
    float getB(){
        return b;
    }
        
    
};


class items {
public:
    float firstx, firsty, firstz, secondx, secondy, secondz, thirdx, thirdy, thirdz, scalar;
    float dr, dg, db,  sr, sg, sb,  r; 
    float ir, ig ,ib ;
    char ty;
    int MatID;
    items(char type, float x1, float y1, float z1,
                     float x2, float y2, float z2,
                     float x3, float y3, float z3, float d){
        firstx  = x1;
        secondx = x2;
        thirdx  = x3;
        firsty  = y1;
        secondy = y2;
        thirdy  = y3;
        firstz  = z1;
        secondz = z2;
        thirdz  = z3;
        ty      = type;
        scalar  = d;
  	MatID = 0;      
    }
    
    float getX1(){
        return firstx;
    }
    float getY1(){
        return firsty;
    }
    float getZ1(){
        return firstz;
    }
    float getX2(){
        //secondx = 0;
        return secondx;
    }
    float getY2(){
        return secondy;
    }
    float getZ2(){
        return secondz;
    }
    float getX3(){
        return thirdx;
    }
    float getY3(){
        return thirdy;
    }
    float getZ3(){
        return thirdz;
    }
    float getType(){
        return ty;
    }
    float getD(){
        return scalar;
    }
    
    void StoreMatID(int mat){
        MatID = mat;
    }
    int getMatID(){
        return MatID;
    }
    void StoreReflectedLight(float parir, float parig, float parib){
        ir = parir;
        ig = parig;
        ib = parib;
        
    }
    float getIR(){
        return ir;
    }
    float getIG(){
        return ig;
    }
    float getIB(){
        return ib;
    }
};


class materials {
    public:
    float dr, dg, db,  sr, sg, sb,  r, ir,  ig,ib;
    materials(float pardr, float pardg, float pardb, float parsr, float parsg, float parsb, float parr, float parir, float parig, float parib ){
        dr = pardr;
        dg = pardg;
        db = pardb;
        sr = parsr;
        sg = parsg;
        sb = parsb;
        ir = parir;
        ig = parig;
        ib = parib;
        r = parr;
    }
    
    float getdr(){
        return dr;
    }
    float getdg(){
        return dg;
    }
    float getdb(){
        return dg;
    }
    float getsr(){
        return sr;
    }
    float getsg(){
        return sg;
    }
    float getsb(){
        return sb;
    }
    float getir(){
        return ir;
    }
    float getig(){
        return ig;
    }
    float getib(){
        return ib;
    }
    float getr(){
        return r;
    }
     
};
