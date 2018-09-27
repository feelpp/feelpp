// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#ifndef BASISFACTORY_H
#define BASISFACTORY_H

#include "MElement.h"
#include "MPyramid.h"
#include "nodalBasis.h"
#include "JacobianBasis.h"

class BasisFactory
{
  typedef std::map<std::pair<int, int>, bezierBasis*> Cont_bezierBasis;

 private:
  static std::map<int, nodalBasis*> fs;
  static std::map<int, JacobianBasis*> js;
  static Cont_bezierBasis bs;
  // store bezier bases by parentType and order (no serendipity..)

 public:
  // Caution: the returned pointer can be NULL
  static const nodalBasis* getNodalBasis(int tag);
  static const JacobianBasis* getJacobianBasis(int tag);
  static const bezierBasis* getBezierBasis(int parentTag, int order);
  static inline const bezierBasis* getBezierBasis(int tag) {
      return getBezierBasis(ElementType::ParentTypeFromTag(tag),
                            ElementType::OrderFromTag(tag) );
    }

  static void clearAll();
};

#endif
