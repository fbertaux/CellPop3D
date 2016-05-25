#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "nr3.hpp"

struct RanGen
{
    unsigned long long int u,v,w ;
    RanGen (unsigned long long int j) : v(4101842887655102017LL), w(1)
    {
        u = j ^ v; int64() ;
        v = u; int64() ;
        w = v; int64();
    }
    inline unsigned long long int int64()
    {
        u = u * 2862933555777941757LL + 7046029254386353087LL ;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8 ;
        w = 4294957665U*(w & 0xffffffff) + (w >> 32) ;
        unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4 ;
        return (x + v) ^ w ;
    }
    inline double doub01() { return 5.42101086242752217E-20 * int64() ; }
    inline double gaussian() { return sqrt (-2.0*log (doub01 ())) * cos (2.0 * M_PI *doub01 ()) ; }
    inline unsigned int int32() { return (unsigned int)int64() ; }
    double norm(double mu, double sig)
    {
        return mu+sig*gaussian();
    }
    double logNorm(double mean,double std)
    {
        double cv = std/mean ;
        double mu = log(mean)-0.5*log(1+cv*cv) ;
        double sig = sqrt( log(cv*cv+1) ) ;
        return exp(norm(mu,sig)) ;
    }
    void unit2dVector(double &ex, double &ey)
    {
        ex = gaussian() ;
        ey = gaussian() ;
        double normvec = sqrt(ex*ex+ey*ey) ;
        ex /= normvec ;
        ey /= normvec ;
    }   
    void unit3dVector(double &ex, double &ey, double &ez)
    {
        ex = gaussian() ;
        ey = gaussian() ;
        ez = gaussian() ;
        double normvec = sqrt(ex*ex+ey*ey+ez*ez) ;
        ex /= normvec ;
        ey /= normvec ;
        ez /= normvec ;
    }
} ;

#endif // RANDOM_HPP
