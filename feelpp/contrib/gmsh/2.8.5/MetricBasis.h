// Gmsh - Copyright (C) 1997-2013 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#ifndef _METRIC_BASIS_H_
#define _METRIC_BASIS_H_

#include "MElement.h"
#include "JacobianBasis.h"
#include "fullMatrix.h"
#include <fstream>
#include <cmath>

class MetricBasis {
  friend class MetricCoefficient;
  friend class GMSH_AnalyseCurvedMeshPlugin;
private:
  const JacobianBasis *_jacobian;
  const GradientBasis *_gradients;
  const bezierBasis *_bezier;
  static double _tol;
  static int _which;

  int __maxdepth, __numSubdivision, __TotSubdivision;
  std::vector<int> __numSub;
  MElement *__curElem;

  std::fstream file;

  class IneqData {
  public:
    int i, j, k;
    double val;

  public:
    IneqData(double val, int i, int j, int k = -1) : i(i), j(j), k(k), val(val) {}
  };

  class MetricData {
   public:
    fullMatrix<double> *_metcoeffs;
    fullVector<double> *_jaccoeffs;
    double _RminBez;
    int _depth, _num;

   public:
    MetricData(fullMatrix<double> *m, fullVector<double> *j, double r, int d, int num) :
      _metcoeffs(m), _jaccoeffs(j), _RminBez(r), _depth(d), _num(num) {}
    ~MetricData() {
      delete _metcoeffs;
      delete _jaccoeffs;
    }
  };

  std::map<int, std::vector<IneqData> > _ineqJ2, _ineqP3, _ineqA;

public:
  MetricBasis(int elementTag);

  static void setTol(double tol) {_tol = tol;}
  static double getTol() {return _tol;}
  static void setWhich(int which) {_which = which;}

  double getBoundRmin(MElement*, MetricData*&, fullMatrix<double>&);
  double getMinR(MElement*, MetricData*&, int) const;
  bool notStraight(MElement*, double &metric, int order) const;
  static double boundMinR(MElement *el);
  static double sampleR(MElement *el, int order);
  //double getBoundRmin(int, MElement**, double*);
  //static double boundRmin(int, MElement**, double*, bool sameType = false);

  void interpolate(const MElement*, const MetricData*, const double *uvw, double *minmaxQ, bool write = false) const;

  static int metricOrder(int tag);
  void printTotSubdiv(double n) const {
    Msg::Info("SUBDIV %d, %g", __TotSubdivision, __TotSubdivision/2776.);
  }

private:
  void _fillInequalities(int order);
  void _lightenInequalities(int&, int&, int&); //TODO change

  void _computeRmin(const fullMatrix<double>&, const fullVector<double>&,
                    double &RminLag, double &RminBez, int depth, bool debug = false) const;
  void _computeRmax(const fullMatrix<double>&, const fullVector<double>&,
                    double &RmaxLag) const;
  void _computeTermBeta(double &a, double &K, double &dRda,
                        double &term1, double &phip) const;
  void _getMetricData(MElement*, MetricData*&) const;

  double _subdivideForRmin(MetricData*, double RminLag, double tol, int which) const;

  double _minp(const fullMatrix<double>&) const;
  double _minp2(const fullMatrix<double>&) const;
  double _minq(const fullMatrix<double>&) const;
  double _maxp(const fullMatrix<double>&) const;
  double _maxq(const fullMatrix<double>&) const;
  void _minMaxA(const fullMatrix<double>&, double &min, double &max) const;
  void _minJ2P3(const fullMatrix<double>&, const fullVector<double>&, double &min) const;
  void _maxAstKpos(const fullMatrix<double>&, const fullVector<double>&,
                 double minK, double beta, double &maxa) const;
  void _maxAstKneg(const fullMatrix<double>&, const fullVector<double>&,
                 double minK, double beta, double &maxa) const;
  void _maxKstAfast(const fullMatrix<double>&, const fullVector<double>&,
                 double mina, double beta, double &maxK) const;
  void _maxKstAsharp(const fullMatrix<double>&, const fullVector<double>&,
                 double mina, double beta, double &maxK) const;
  void _minMaxJacobianSqr(const fullVector<double>&, double &min, double &max) const;

  double _Rsafe(double a, double K) const {
    const double x = .5 * (K - a*a*a + 3*a);
    const double phi = std::acos(x) / 3;
    return (a + 2*std::cos(phi + 2*M_PI/3)) / (a + 2*std::cos(phi));
  }
  bool _chknumber(double val) const {
#if defined(_MSC_VER)
    return _isnan(val) || !_finite(val);
#else
    return std::isnan(val) || std::isinf(val);
#endif
  }
  bool _chka(double a) const {return _chknumber(a) || a < 1;}
  bool _chkK(double K) const {return _chknumber(K) || K < 0;}
  int _chkaK(double a, double K) const {
    if (_chka(a)) return 1;
    if (_chkK(K)) return 2;
    if (std::abs(K - a*a*a + 3*a) > 2) {
      Msg::Warning("x = %g", .5 * (K - a*a*a + 3*a));
      return 3;
    }
    return 0;
  }
  bool _chkR(double R) const {return _chknumber(R) || R < 0 || R > 1;}
  int _chkaKR(double a, double K, double R) const {
    const int aK = _chkaK(a, K);
    if (aK) return aK;
    if (_chkR(R)) return 4;
    const double myR = _Rsafe(a, K);
    if (std::abs(myR-R) > 1e-10) return 5;
    return 0;
  }

private:
  class gterIneq {
   public:
    bool operator()(const IneqData &id1, const IneqData &id2) const {
      return id1.val > id2.val;
    }
  };
  struct lessMinB {
    bool operator()(const MetricData *md1, const MetricData *md2) const {
      return md1->_RminBez > md2->_RminBez;
    }
  };
};

#endif
