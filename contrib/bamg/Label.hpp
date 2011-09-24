#ifndef LABEL_HPP
#define LABEL_HPP
 
class Label {  // reference number for the physics
  public: 
  int lab;
  Label(int r=0):lab(r){}
  bool onGamma() const { return lab;} 
  };
 inline ostream& operator <<(ostream& f,const Label & r  )
    { f <<  r.lab ; return f; }
 inline istream& operator >>(istream& f, Label & r  )
    { f >>  r.lab ; return f; }

#endif
