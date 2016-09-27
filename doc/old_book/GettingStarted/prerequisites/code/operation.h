#ifndef OPERATION_H
#define OPERATION_H
 
#include <iostream>
using namespace std;
 
class operation
{
 public:
  /**
   *Constructors
   */
  operation();
  operation(double x, double y);
 
 /**
  *Accessors and mutators
  */
  void setX(double x);
  void setY(double y);
  double getX() const;
  double getY() const;
 
  
  /**
   *functions
   */
  double addition();
  double substruction();
  double multiplication();
  double division();
 
   /**
    *@param x is a variable
    *@param y is a variable too
    */
 private: 
  double x,y;
};
/**
 *@brief     operation.h
 *@date      june 23, 2015
 *@author:   kyoshe winstone
 *@version   1.0
 */

 
#endif
