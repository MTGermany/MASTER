
#include <stdio.h>
#include <math.h>
#include "general.h"
#include "master.const"
#include "master.h"
#include "Ctrl_ramp.h"

Ctrl_ramp::Ctrl_ramp(int index, 
            double x_center, double length, 
            double max_flow)

{ 
  this->index=index;
  this->x_center=x_center;
  this->length=length;
  this->max_flow=max_flow;
  nveh_wait=0;
  in_flow=desired_flow=rmp_flow=0;

  sprintf(inflow_fname,"%s.ctrl_rmp%i",namepar, index);
  get_array(inflow_fname, n_data, times_data
                          inflow_data);
 
}

