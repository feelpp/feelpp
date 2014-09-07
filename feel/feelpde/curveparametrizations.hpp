#ifndef _CURVEPARAMETRIZATION_HPP
#define _CURVEPARAMETRIZATION_HPP 1


/*
   This header gives parametrization of few curves
   functions return a tuple containing : x(t), y(t), tStart, tEnd
   needed to make the curve.
   It can be given directly to the method DistToCurve::fromParametrizedCurve
*/


namespace Feel
{

    class CurveParametrization
    {

    public :

        static std::tuple< std::function<double(double)>, std::function<double(double)>, double, double >
        circle( double r, double x0, double y0 )
        {
            auto x = [=] (double t) -> double { return x0 + r * std::cos(t); };
            auto y = [=] (double t) -> double { return y0 + r * std::sin(t); };
            return std::make_tuple( x, y, 0, 2*M_PI );
        }



        static std::tuple< std::function<double(double)>, std::function<double(double)>, double, double >
        ellipse( double a_ell, double b_ell, double x0, double y0 )
        {
            auto x_ell = [=](double t) -> double {return x0 + a_ell * std::cos(t); };
            auto y_ell = [=](double t) -> double {return y0 + b_ell * std::sin(t); };
            return std::make_tuple( x_ell, y_ell, 0, 2*M_PI );
        }



        static std::tuple< std::function<double(double)>, std::function<double(double)>, double, double >
        tiltedEllipse( double a_ell, double b_ell, double teta, double x0, double y0 )
        {
          auto x_ell = [=](double t) -> double {return x0 + (a_ell * std::cos(t)) * std::cos(teta) - (b_ell * std::sin(t)) * std::sin(teta); };
          auto y_ell = [=](double t) -> double {return y0 + (b_ell * std::sin(t)) * std::cos(teta) + (a_ell * std::cos(t)) * std::sin(teta); };
          return std::make_tuple( x_ell, y_ell, 0, 2*M_PI );
        }



        /* don't known the tStart, tEnd .. should be added in the future */
        static std::tuple< std::function<double(double)>, std::function<double(double)> >
        epitrochoid( double a_epi, double b_epi, double scaleFactor, double x0, double y0 )
        {
            auto x_epi = [=](double t) -> double { return ((1+a_epi) * std::cos(a_epi*t) - a_epi*b_epi * std::cos( (1+a_epi) * t )) / scaleFactor + x0; };
            auto y_epi = [=](double t) -> double { return ((1+a_epi) * std::sin(a_epi*t) - a_epi*b_epi * std::sin( (1+a_epi) * t )) / scaleFactor + y0; };
            return std::make_tuple( x_epi, y_epi );
        }



        /* siclke cell shape or "moon" or "croissant/crescent" shape */
        /*  R1 should be the bigger circle */
        static  std::tuple< std::function<double(double)>, std::function<double(double)>, double, double >
        crescent( double R1, double R2, double H, double l, double x0, double y0)
        {
            const double dt1=l/R1;
            const double dt2=l/R2;

            const double t1=std::asin(H/R1);
            const double t2=std::asin(H/R2);

            const double X=R1*std::cos(t1)-R2*std::cos(t2);

            // check for nan numbers
            CHECK( (t1==t1) && (t2==t2) ) << "problem in crescent shape parametrization : an asin function is resulting in Nan\n";

            auto x_sc = [=](double t) -> double {
                if ((-t1+dt1 <= t) && (t <= t1-dt1))
                    return R1*std::cos(t) + x0 - R1*std::cos(t1) ;

                else if ((t1 - dt1 < t) && (t <= t1 + dt2))
                    return R1*std::cos(t1-dt1) + x0 + (t-(t1-dt1))/(t1+dt2-(t1-dt1)) *(R2*std::cos(t2-dt2) + X + x0 - (R1*std::cos(t1-dt1)+x0 )) - R1*std::cos(t1);

                else if ((t1+dt2 < t) && (t <= t1+2*t2-dt2))
                    return R2*std::cos(-(t-t1-t2)) + X + x0 - R1*std::cos(t1) ;

                else if ((t1+2*t2-dt2 < t) && (t <= t1+2*t2+dt1))
                    return R2*std::cos(-(t2-dt2))+X+x0 + (t-(t1+2*t2-dt2))/(t1+2*t2+dt1-(t1+2*t2-dt2))*(R1*std::cos(-t1+dt1)+x0-(R2*std::cos(-(t2-dt2))+X+x0 )) - R1*std::cos(t1) ;

                else
                    CHECK( false ) << "parameter range might be wrong"<<std::endl;
            };

            auto y_sc = [=](double t) -> double {
                if ((-t1+dt1 <= t) && (t <= t1-dt1))
                    return R1*std::sin(t)+y0;

                else if ((t1 - dt1 < t) && (t <= t1 + dt2))
                    return R1*std::sin(t1-dt1)+y0 +  (t-(t1-dt1))/(t1+dt2-(t1-dt1))*(R2*std::sin(t2-dt2)+y0-(R1*std::sin(t1-dt1)+y0 ));

                else if ((t1+dt2 < t) && (t <= t1+2*t2-dt2))
                    return R2*std::sin(-(t-t1-t2))+y0;

                else if ((t1+2*t2-dt2 < t) && (t <= t1+2*t2+dt1))
                    return R2*std::sin(-(t2-dt2)) + y0 + (t-(t1+2*t2-dt2))/(t1+2*t2+dt1-(t1+2*t2-dt2))*(R1*std::sin(-t1+dt1)+y0-(R2*std::sin(-(t2-dt2))+y0));

                else
                    CHECK( false ) << "parameter range might be wrong"<<std::endl;
            };

            return std::make_tuple( x_sc, y_sc, -t1+dt1, t1+2*t2+dt1 );
        }

    }; // CurveParametrization

} // namespace Feel



#endif
