// ***********************************************************************
// On-ramps with flow control
// ***********************************************************************

#ifndef CTRL_RAMP_H
#define CTRL_RAMP_H

const int NLINESMAX = 1000;

class Ctrl_ramp
{
 public:
  Ctrl_ramp(int index, 
            double x_center, double length, 
            double max_flow);

  void update(double time, double dt);
  void print_status();

 private:
  char   inflow_fname[MAXSTR];         // Name of the file with inflow data

  int index;
  double nveh_wait;    // fractional number of queuing vehicles
  double x_center;     // ramp at [x_center-length/2,x_center-length/2]
  double length;

  int n_data;         // how many data lines in inflow file
  double times_data[NLINESMAX];     // from .ctrl_rmp<index> File
  double inflow_data[NLINESMAX];     // from .ctrl_rmp<index> File

  double inflow;    // onramp blocks if inflow > max_flow
  double max_flow;    // onramp blocks if in_flow > max_flow
  double desired_flow;     // in_flow if nveh_wait=0; max_flow otherwise
  double rmp_flow;    // min(desired_flow, max_flow)
};

#endif
