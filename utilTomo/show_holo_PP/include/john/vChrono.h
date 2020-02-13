#include <boost/chrono.hpp>

// requiert boost 1.4.8

// types possibles pour boostClock: 
//
// boost::chrono::system_clock
// boost::chrono::steady_clock
// boost::chrono::high_resolution_clock


#ifndef __VECTRA_CHRONO__
#define __VECTRA_CHRONO__


template< class boostClock >
class vChrono
{
  typename boostClock::time_point start;
  
 public:

  // constructeur
  //vChrono() : start( boostClock::now() ) {}
  // se comprend comme { start = boostClock::now(); }
  vChrono()
    { this -> reset(); }

  // méthode elapsed, retourne une valeur de type boostClock::duration
  // on est obligé de préfixer par typename, car boostClock est un nom-template
  typename boostClock::duration elapsed() const
    {
      return boostClock::now() - start;
    }

  // retourne la durée écoulée depuis la construction de l'instance
  // en secondes
  double seconds() const
  {
    return elapsed().count() * ((double)boostClock::period::num/boostClock::period::den);
  }

  // en millisecondes
  double milliseconds() const
  {
    return 1000 * this -> seconds();
  }

  // remet le compteur à zéro
  void reset()
  {
    start = boostClock::now();
  }
};




#endif
