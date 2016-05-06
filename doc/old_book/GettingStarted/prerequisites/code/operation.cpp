#include "operation.h"
#include <iostream>
using namespace std;



operation::operation() : x(0), y(0)
{}
  /**
   *Constructs an objects with two arguments x and y
   *@param x is a variable
   *@param y is a variable
   */
operation::operation(double x, double y) : x(x), y(y)
{}

  /**
   *sets the x value
   */
void operation::setX(double x)
{
  this->x = x;
}
 
  /**
   *sets the y value 
   */
void operation::setY(double y)
{
  this->y = y;
}

  /**
   *@return the x value
   */
double operation::getX() const
{
  return this->x;
} 

  /**
   *@return the  y value
   */
double operation::getY() const
{
  return this->y;
} 
 

  /**
   *@returns the sum of two numbers
   */
double operation::addition( )
{
  return (x + y);
}

  /**
   *@returns the diff of two numbers
   */
double operation::substruction( )
{
  return (x - y);
}

  /**
   *@returns the product of two numbers
   */
double operation::multiplication( )
{
  return (x * y);
}

  /**
   *@returns the quotient of two numbers
   */
double operation::division( )
{
 if(y != 0)
    {
      return x / y;

    }
    
  cout << "Error: division by zero.\n";
  return 0;
}
/**
 * @brief     operation.h
 * @date      june 19, 2015
 * @author   kyoshe winstone
 * @version   1.1
 */

