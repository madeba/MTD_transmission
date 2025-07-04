// *****************************************************************************
// constructors
// *****************************************************************************


template <class T>
vPlotterT<T>::vPlotterT(T (*fun_x)(T), T (*fun_y)(T))
{
  this -> fun_x = fun_x;
  this -> fun_y = fun_y;

  nb_steps = i = 0;
  repeat_counter = t = t_min = t_max = x = y =0;
  repeat_mode_p = regular_p = adaptative_p = start_zero_p = functions_set = range_set = ready = false;
  data_x = data_y = NULL; data_allocated = false;

  repeat_steps = factor_x = factor_y = 1;
}


template <class T>
vPlotterT<T>::vPlotterT()
{
  nb_steps = i = 0;
  repeat_counter = t = t_min = t_max = x = y =0;
  repeat_mode_p = regular_p = adaptative_p = start_zero_p = functions_set = range_set = ready = false;
  data_x = data_y = NULL; data_allocated = false;

  repeat_steps = factor_x = factor_y = 1;
}


// *****************************************************************************
// config normal mode
// *****************************************************************************


template <class T>
void
vPlotterT<T>::setFunctions(T (*fun_x)(T), T (*fun_y)(T))
{
  this -> fun_x = fun_x;
  this -> fun_y = fun_y;

  functions_set = true;
  if ((regular_p || adaptative_p) && range_set)
    ready = true;
}


template <class T>
void
vPlotterT<T>::set_range(T t_min, T t_max)
{
  this -> t_min = t_min; this -> t_max = t_max;
  t = t_min; i = 0;
  range_set = true;
  if ((regular_p || adaptative_p) && functions_set)
    ready = true;
}


template <class T>
void
vPlotterT<T>::set_sampling_regular(size_t steps)
{
  this -> nb_steps = steps;
  regular_p = true;
  ASSERT(! adaptative_p);
  if (functions_set && range_set) ready = true;
}



template <class T>
void
vPlotterT<T>::set_sampling_adaptative(T (*fun_next_t)(T))
{
  this -> fun_updateT = fun_next_t;
  adaptative_p = true;
  ASSERT(! regular_p);
  if (functions_set && range_set) ready = true;
}


// *****************************************************************************
// compute
// *****************************************************************************

template <class T>
void
vPlotterT<T>::compute_next()
{
  ASSERT(ready);

  if (start_zero_p && (i == 0))
    {
      start_zero_p = false;
      x = y = 0;
      return;
    }

  if (data_allocated)
    {
      x = data_x[i] * factor_x;
      y = data_y[i] * factor_y;
    }
  else
    {
      if (regular_p)
	t = ( t_max - t_min ) * (T)i / T(nb_steps);
      if (adaptative_p)
	t += fun_updateT(t);


      x = fun_x(t) * factor_x;
      y = fun_y(t) * factor_y;
    }

  if (! repeat_mode_p)
    i++; // cas général
  else if (regular_p)
    {
      if (repeat_counter == repeat_steps)
	{i++; repeat_counter = 0;}
      else
	repeat_counter++;
    }
}
