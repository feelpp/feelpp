namespace Feel
{
template <typename T = float>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;
class MultiScaleImage
{
    public :
    int T(holo3_image<float> im, std::pair<double,double> c)
    {
     double x = c.first;
     double y = c.second;
     
     int i = x/dx;
     int j = y/dy;
     /*
     auto Xhc = Pch<1>();
     Hbf2Feelpp h2f(im.cols(),im.rows,Xhc);
    */
    return im(i,j);
    }

    int T(holo3_image<float> im, std::pair<double,double> c, int L)
    {
     double x = c.first;
     double y = c.second;
     
     int i = x/dx;
     int j = y/dy;
     /*
     auto Xhc = Pch<1>();
     Hbf2Feelpp h2f(im.cols(),im.rows,Xhc);
    */
    return im(L*i,L*j);
    }
     
    private :
    double dx =8.9e-3;
    double dy =8.9e-3;
};

}
