// la regle de programmation 3  
#include "assertion.hpp" 
#include <cstdlib>
using namespace std;
// definition R
#include "Label.hpp"
#include "R2.hpp"
  

class Vertex2 : public R2,public Label {
  friend  ostream& operator <<(ostream& f, const Vertex2 & v )
    { f << (R2) v << ' ' << (Label &) v   ; return f; }
  friend  istream& operator >> (istream& f,  Vertex2 & v )
      { f >> (R2 &) v >> (Label &) v ; return f; }
   public:
   Vertex2() : R2(),Label(){};
  //   Vertex2(R2 P,int r=0): R2(P),Label(r){}
   private: // pas de copie pour ne pas perdre l'adresse
     Vertex2(const Vertex2 &);
     void operator=(const Vertex2 &); 
};







class Triangle2: public Label {
  // variable prive 
  Vertex2 *vertices[3]; // an array of 3 pointer to vertex
  public:
  R area;
  Triangle2() :area() {vertices[0]=vertices[1]=vertices[2]=0;}; // constructor empty for array
  Vertex2 & operator[](int i) const {
     ASSERTION(i>=0 && i <3);
     return *vertices[i];} // to see traingle as a array of vertex
  Vertex2 * & operator()(int i) {
     ASSERTION(i>=0 && i <3);
     return vertices[i];} // to see traingle as a array of vertex
  void init(Vertex2 * v0,int * iv,int r) 
    { vertices[0]=v0+iv[0];
      vertices[1]=v0+iv[1];
      vertices[2]=v0+iv[2]; 
      R2 AB(*vertices[0],*vertices[1]);
      R2 AC(*vertices[0],*vertices[2]);
      area = (AB^AC)*0.5;
      lab=r;
      ASSERTION(area>0);    
    }
  R2 Edge(int i) const {ASSERTION(i>=0 && i <3);
      return R2(*vertices[(i+1)%3],*vertices[(i+2)%3]);}// opposite edge vertex i
  R2 H(int i) const { ASSERTION(i>=0 && i <3);
      R2 E=Edge(i);return E.perp()/(2*area);} // heigth 

  void Gradlambda(R2 * GradL) const
  {
    GradL[1]= H(1);
    GradL[2]= H(2);
    GradL[0]=-GradL[1]-GradL[2];
  }

  R lenEdge(int i) const {ASSERTION(i>=0 && i <3);
      R2 E=Edge(i);return sqrt((E,E));}

  // Transformation:  $\hat{K} \mapsto K$ 
  R2 operator()(const R2 & Phat) const   {
    const R2 &A =*vertices[0];
    const R2 &B =*vertices[1];
    const R2 &C =*vertices[2];
    return (1-Phat.x- Phat.y)* A + Phat.x *B  +Phat.y*C  ;} 
    
 private:
   Triangle2(const Triangle2 &); // pas de construction par copie
   void operator=(const Triangle2 &);// pas affectation par copy 
  // Ajoute FH  dernier cours ----
public:
  R  mesure() const {return area;}
  static const int nv=3;  
};





// la classe pour le maillage du bord, ici des Segment: 
class Seg: public Label {
public:
  static const int nv=2;     // the nomber of vertices 
private:
  Vertex2 *vertices[nv]; // an array of 2 pointer to vertex
public:
  R l; // lhe lenght of the segment

  Seg() :l() {vertices[0]=vertices[1]=0;} // constructor empty for array

  Vertex2 & operator[](int i) const {
    ASSERTION(i>=0 && i <nv);
    return *vertices[i];} // to see the segment  as a array of vertex

  void init(Vertex2 * v0,int * iv,int r)  // the true constructor
  { 
    vertices[0]=v0+iv[0];
    vertices[1]=v0+iv[1];
    R2 AB(*vertices[0],*vertices[1]);
    l=  AB.norme();
    lab=r;
    ASSERTION(l>0);    
  }
  
  // Transformation:  $[0,1]  \mapsto K$ 
  R2 operator()(const R & Phat) const   {
    const R2 &A =*vertices[0];
    const R2 &B =*vertices[1];
    return (1-Phat)* A + Phat *B  ;} 

  R  mesure() const {return l;}

private:
  Seg(const Seg &); // pas de construction par copie
  void operator=(const Seg &);// pas affectation par copy 
};




class Mesh2 
{ 
public:
  typedef Triangle2 Element;   
  typedef Seg  BorderElement;
  typedef R2 Rd;
  typedef Vertex2  Vertex;
  int nv,nt, nbe;
  R area,peri;
  Vertex2 *vertices;
  Triangle2 *triangles;
  Seg * borderelements;
  Triangle2 & operator[](int i) const {return triangles[CheckT(i)];}
  Vertex2 & operator()(int i) const {return vertices[CheckV(i)];}
  Seg & be(int i) const {return borderelements[CheckBE(i)];} // boundary element ..
  Mesh2(const char * filename); // read on a file
  Mesh2(int nnv,int nnt,int nnbe, Vertex2 *v,Triangle2 *t, Seg *b_e)
    : nv(nnv),nt(nnt),nbe(nnbe),vertices(v),triangles(t),borderelements(b_e),area(0.),peri(0.)  
  {
    for(int i=0;i<nt;++i) area += triangles[i].mesure();
    for(int i=0;i<nbe;++i) peri += borderelements[i].mesure();
  }
  int operator()(const Triangle2 & t) const {return CheckT(&t - triangles);}
  int operator()(const Triangle2 * t) const {return CheckT(t - triangles);}
  int operator()(const Vertex2 & v) const {return CheckV(&v - vertices);}
  int operator()(const Vertex2  * v) const{return CheckV(v - vertices);}
  int operator()(const Seg & v) const {return CheckBE(&v - borderelements);}
  int operator()(const Seg  * v) const{return CheckBE(v - borderelements);}
  int operator()(int it,int j) const {return (*this)(triangles[it][j]);}// Nu vertex j of triangle it
  //  to check the bound 
  int CheckV(int i) const { ASSERTION(i>=0 && i < nv); return i;} 
  int CheckT(int i) const { ASSERTION(i>=0 && i < nt); return i;}
  int CheckBE(int i) const { ASSERTION(i>=0 && i < nbe); return i;}
  ~Mesh2() { 
    delete [] vertices; 
    delete [] triangles;
    delete [] borderelements;}
 private:
   Mesh2(const Mesh2 &); // pas de construction par copie
   void operator=(const Mesh2 &);// pas affectation par copy 
};


 
