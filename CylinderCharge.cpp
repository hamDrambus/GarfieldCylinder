#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include <Garfield/ViewCell.hh>
#include <Garfield/ViewDrift.hh>
#include <Garfield/ViewSignal.hh>

#include <Garfield/GarfieldConstants.hh>
#include <Garfield/ComponentAnalyticField.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/DriftLineRKF.hh>
#include <Garfield/AvalancheMC.hh>
#include <Garfield/TrackHeed.hh>
#include <Garfield/Random.hh>
#include <Garfield/Shaper.hh>

#include "MaterialCustomDrift.hh"

using namespace Garfield;

std::string strtoken(std::string &in, std::string break_symbs)
{
  std::string out;
  auto erase_marker = in.begin();
  auto copy_from = in.begin(), copy_to = in.end();
  while (erase_marker!=in.end()) {
    bool break_ = false;
    for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
      if (*erase_marker == *h) {
        break_ = true;
        break;
      }
    if ((break_) && copy_from == erase_marker) {
      copy_from = ++erase_marker;
      continue;
    }
    ++erase_marker;
    if (break_)
      break;
    copy_to = erase_marker;
  }
  out = std::string(copy_from, copy_to);
  in.erase(in.begin(), erase_marker);
  return out;
}

std::pair<std::vector<double>, std::vector<double>> read_xy_from_file(std::string fname)
{
  std::vector<double> xs, ys;
  std::ifstream str(fname);
  if (!str.is_open()) {
    std::cerr<<"Failed to open \""<<fname<<"\" file."<<std::endl;
    return std::make_pair(xs, ys);
  }
	std::string line, word;
	int line_n = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		if (line.size() >= 2) // Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		try {
			word = strtoken(line, "\t \r");
			double x = std::stod(word);
			word = strtoken(line, "\t \r");
			double val = std::stod(word);
			xs.push_back(x);
      ys.push_back(val);
		} catch (std::exception &e) {
      continue;
		}
	}
  return std::make_pair(xs, ys);
}

void scale_xy(std::pair<std::vector<double>, std::vector<double>> &xys, double x_scale, double y_scale = 1.0) {
  for (auto & x : xys.first)
    x *= x_scale;
  for (auto & y : xys.second)
    y *= y_scale;
}

bool str_in(const std::string &what, const std::vector<std::string> & set) {
  for (const auto & str : set)
    if (str == what)
      return true;
  return false;
}

bool str_not_in(const std::string &what, const std::vector<std::string> & set) {
  return !str_in(what, set);
}

auto find_str_in(const std::string &what, std::vector<std::pair<std::string, std::string>> & dict) {
  for (std::size_t i = 0, i_end_ = dict.size(); i!=i_end_; ++i)
    if (dict[i].first == what)
      return dict.begin()+i;
  return dict.end();
}

auto find_str_in(const std::string &what, const std::vector<std::pair<std::string, std::string>> & dict) {
  for (std::size_t i = 0, i_end_ = dict.size(); i!=i_end_; ++i)
    if (dict[i].first == what)
      return dict.begin()+i;
  return dict.end();
}

bool parse_argument(double& value, std::string& str_value, const std::string& key, const std::vector<std::pair<std::string, std::string>> & dict)
{
  auto it = find_str_in(key, dict);
  if (it == dict.end()) {
    return true;
  } else {
    try {
      value = std::stod(it->second);
      str_value = it->second;
    } catch (const std::exception &e) {
      std::cerr<<"Ivalid argument \"-"<<it->first<<"="<<it->second<<"\""<<std::endl;
      std::cerr<<e.what()<<std::endl;
      return false;
    }
  }
  return true;
}

bool parse_argument(std::size_t& value, std::string& str_value, const std::string& key, const std::vector<std::pair<std::string, std::string>> & dict)
{
  auto it = find_str_in(key, dict);
  if (it == dict.end()) {
    return true;
  } else {
    try {
      value = std::stoul(it->second);
      str_value = it->second;
    } catch (const std::exception &e) {
      std::cerr<<"Ivalid argument \"-"<<it->first<<"="<<it->second<<"\""<<std::endl;
      std::cerr<<e.what()<<std::endl;
      return false;
    }
  }
  return true;
}

bool parse_argument(bool& value, std::string& str_value, const std::string& key, const std::vector<std::pair<std::string, std::string>> & dict)
{
  std::vector<std::string> bool_args = {"on", "off", "true", "false"};
  auto it = find_str_in(key, dict);
  if (it == dict.end()) {
    return true;
  } else {
    try {
      if (str_not_in(it->second, bool_args))
        throw std::invalid_argument(it->second + " cannot be converted to bool value.");
      if (it->second == "on" || it->second == "true") {
        value = true;
        str_value = "on";
      } else {
        value = false;
        str_value = "off";
      }
    } catch (const std::exception &e) {
      std::cerr<<"Ivalid argument \"-"<<it->first<<"="<<it->second<<"\""<<std::endl;
      std::cerr<<e.what()<<std::endl;
      return false;
    }
  }
  return true;
}

bool parse_argument(std::string& value, std::string& str_value, const std::string& key, const std::vector<std::pair<std::string, std::string>> & dict)
{
  auto it = find_str_in(key, dict);
  if (it == dict.end()) {
    return true;
  } else {
    value = it->second;
    str_value = it->second;
  }
  return true;
}

int main(int argc, char * argv[]) {
  // Define input parameters and set to defaults.
  double V0 = 20; // Wire voltage [kV]. Tube is grounded
  double rWire = 400; // Wire radius [um]
  double rWireOffset = 0; // Wire offset from center [um]
  double rTube = 0.5; // Outer radius of the tube [cm]
  double amplifier_shaping = 0.0; // in [ns]
  // Version of cross sections set which was used to calculate transport parameters (Geant4 NBrS project).
  std::string XS_name = "Var8";
  bool use_diffusion = true;
  bool do_plot = false;
  bool do_plot_tracks = false;
  std::string out_fname = "";
  std::size_t N = 1000;

  // Parameters as string for output files.
  std::string str_V0 = "20.0",
              str_rWire = "400",
              str_rWireOffset = "0",
              str_rTube = "0.5",
              str_XS_name = "Var8",
              str_use_diffusion = "on",
              str_N = "1000",
              str_amplifier_shaping = "0";
  std::string dummy_str;

  // Parse program inputs.
  // Normally one should use boost library for this, but it and ROOT do not get along. 
  bool print_help = false;
  if (argc != 1) {
    std::vector<std::string> pars;
    for (std::size_t i = 1; i!= static_cast<std::size_t>(argc); ++i)
      pars.push_back(argv[i]);
    std::vector<std::string> supported_args = {"V0", "rWire", "rTube", "diff", "XS", "rWireOffset", "N", "plot", "plot_tracks", "out_fname", "sh"};
    std::vector<std::string> help_args = {"-help", "--help", "help", "-h", "--h"};
    std::vector<std::pair<std::string, std::string>> keys_args;
    for (const auto& str : pars) {
      if (str_in(str, help_args)) {
        print_help = true;
        continue;
      }
      std::string arg = str;
      std::string key = strtoken(arg, "=");
      if (arg.empty() || arg.find_first_not_of(" \t\r") == std::string::npos) {
        std::cerr<<"Warning: argument "<<key<<" has invalid value \""<<arg<<"\""<<std::endl;
        continue;
      }
      if (key.empty() || key.size()<2 || key[0]!='-') {
        std::cerr<<"Warning: invalid argument \""<<str<<"\" was provided."<<std::endl;
        print_help = true;
        continue;
      }
      auto temp = std::string(key.begin()+1, key.end());
      key = temp;
      if (str_not_in(key, supported_args)) {
        continue;
      }
      keys_args.push_back(std::make_pair(key, arg));
    }

    if (!parse_argument(V0, str_V0, "V0", keys_args) ||
        !parse_argument(rWire, str_rWire, "rWire", keys_args) ||
        !parse_argument(rWireOffset, str_rWireOffset, "rWireOffset", keys_args) ||
        !parse_argument(rTube, str_rTube, "rTube", keys_args) ||
        !parse_argument(N, str_N, "N", keys_args) ||
        !parse_argument(do_plot, dummy_str, "plot", keys_args) ||
        !parse_argument(do_plot_tracks, dummy_str, "plot_tracks", keys_args) ||
        !parse_argument(XS_name, str_XS_name, "XS", keys_args)  ||
        !parse_argument(amplifier_shaping, str_amplifier_shaping, "sh", keys_args) ||
        !parse_argument(use_diffusion, str_use_diffusion, "diff", keys_args) ||
        !parse_argument(out_fname, dummy_str, "out_fname", keys_args))
      return -2;

    if (rWireOffset + rWire + Small > rTube*1e4) {
      std::cerr<<"Error: Wire intersects or is outside tube (rWireOffset + rWire > rTube).";
      return -3;
    }
  } else {
    std::cout << "Using default settings:" << std::endl;
    std::cout << "-V0=" << str_V0
              << " -rWire=" << str_rWire
              << " -rWireOffset=" << str_rWireOffset
              << " -rTube=" << str_rTube
              << " -N=" << str_N
              << " -sh=" << str_amplifier_shaping
              << " -diff=" << str_use_diffusion
              << " -XS=" << str_XS_name
              << std::endl;
  }
  if (print_help) {
    std::cout<<"Program arguments are:\n";
    std::cout<<"-V0={number}\t\tVoltage at the anode wire in kilovolts (default 20).\n";
    std::cout<<"-rWire={number}\t\tWire radius in micrometers (default 400).\n";
    std::cout<<"-rWireOffset={number}\tWire offset from tube center in micrometers (default 0).\n";
    std::cout<<"-rTube={number}\t\tRadius of cylinder tube in centimeters (default 0.5).\n";
    std::cout<<"-diff={bool}\t\tFlag defining whether to use diffusion in electron drift simulation (default is on).\n";
    std::cout<<"-N={number}\t\tNumber of electrons to simulate (default is 1000).\n";
    std::cout<<"-sh={number}\t\tShaping time of amplifier in ns (default is 0).\n";
    std::cout<<"-plot={bool}\t\tFlag defining whether to plot signal or run program in a batch mode (default is off).\n";
    std::cout<<"-plot_tracks={bool}\tFlag defining whether to plot tracks in addition to signal (default is off).\n"
                "\t\t\tOnly works if -plot is set to on.\n";
    std::cout<<"-XS={string}\t\tVersion of cross sections set which was used to calculate transport parameters in Geant4 NBrS project (default is Var8).\n";
    std::cout<<"-out_fname={string}\t\tOutput filename for signal. If not specified, default filename from input paramters is generated.\n";
    std::cout<<"{bool}=={on|off|true|false}\n";
    std::cout<<std::endl;
  }

  V0 *= 1e3; // kV to V
  rWire *= 1e-4; // um to cm
  rWireOffset *= 1e-4; // um to cm

  bool use_amplifier_shaper = amplifier_shaping > Small;
  double amplifier_gain = 1.0;

  char* cstr = getenv("HOME");
  std::string home_str = cstr == nullptr ? "../../../" : std::string(cstr) + "/";
  std::string data_folder = home_str + ("Documents/Detector_geant4/NBrS_THGEM_LAr_v0/data_cache_LAr_v5.4_exact_")+XS_name+"/";
  std::string out_filename = std::string("r") + str_rWire + "um_R" + str_rTube + "cm_offset" + str_rWireOffset + "um_" +
            "XS_" + str_XS_name + "_V0_" + str_V0 + "kV_diff_" + str_use_diffusion + "_N" + str_N + "_sh" + str_amplifier_shaping + "ns.txt";
  out_filename = out_fname.empty() ? out_filename : out_fname;

  TApplication app("app", &argc, argv);
 
  // Load electron transfer parameters
  // and convert from geant4 to garfield units.
  auto velocity = read_xy_from_file(data_folder + "drift_velocity.dat");
  scale_xy(velocity, 1e7, 1e-1);
  auto diff_long = read_xy_from_file(data_folder + "diffusion_longitudinal.dat");
  scale_xy(diff_long, 1e7, 1e-2);
  auto diff_trans = read_xy_from_file(data_folder + "diffusion_transverse.dat");
  scale_xy(diff_trans, 1e7, 1e-2);
  
  if (velocity.first.empty() || diff_long.first.empty() || diff_trans.first.empty())
    return -1;

  // Crate drift medium
  MaterialCustomDrift gas(velocity.first, velocity.second,
                          diff_long.first, diff_long.second,
                          diff_trans.first, diff_trans.second);

  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.DisableDebugging();
  cmp.SetMedium(&gas);
  // Voltages
  const double vWire = V0;
  const double vTube = 0.;
  // Add the wire in the centre.
  cmp.AddWire(rWireOffset, 0, 2 * rWire, vWire, "s");
  cmp.AddTube(rTube, vTube, 0, "t");
  // Request calculation of the weighting field. 
  cmp.AddReadout("s");

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "s");
  // Set the signal time window.
  const double tstep = 0.1;
  const double tmin = -1.* tstep;
  const unsigned int nbins = 30000;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  Shaper shaper(1, 1, 1, "bipolar");
  if (use_amplifier_shaper) {
    // Set the delta reponse function.
    shaper = Shaper(3, amplifier_shaping, amplifier_gain, "unipolar");
    sensor.SetTransferFunction(shaper);
  }
  sensor.ClearSignal();

  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(3e-4);
  drift.DisableAttachment();
  if (use_diffusion)
    drift.EnableDiffusion();
  else
    drift.DisableDiffusion();
  drift.EnableSignalCalculation(true);
  //drift.UseWeightingPotential(false); // has no effect
  //drift.UseWeightingPotential(true); // has no effect
 
  TCanvas* cD = nullptr;
  ViewCell cellView;
  ViewDrift driftView;
  bool plotDrift = do_plot && do_plot_tracks;
  if (plotDrift) {
    cD = new TCanvas("cD", "", 1000, 800);
    cD->ToggleEventStatus(); cD->ToggleToolBar();
    cellView.SetCanvas(cD);
    cellView.SetComponent(&cmp);
    driftView.SetCanvas(cD);
    drift.EnablePlotting(&driftView);
  }
 
  TCanvas* cS = nullptr;
  ViewSignal signalView;
  bool plotSignal = do_plot;
  if (plotSignal) {
    cS = new TCanvas("cS", "", 1000, 800);
    cS->ToggleEventStatus(); cS->ToggleToolBar();
    signalView.SetCanvas(cS);
    signalView.SetSensor(&sensor);
    signalView.SetLabelY("Charge signal [fC]");
  } 

  const double rTrack = rTube - 1.e-4;
  for (unsigned int j = 0; j < N; ++j) {
    double phi = RndmUniform() * TwoPi;
    double xe = rTrack * cos(phi), ye = rTrack * sin(phi), ze = 0., te = 0.;
    drift.DriftElectron(xe, ye, ze, te);
    sensor.NewSignal();
  }
  if (use_amplifier_shaper)
    sensor.ConvoluteSignal("s");
  if (plotDrift) {
    cellView.Plot2d();
    constexpr bool twod = true;
    constexpr bool drawaxis = false;
    driftView.Plot(twod, drawaxis);
  }
  if (plotSignal)
    signalView.PlotSignal("s");

  sensor.ExportSignal("s", out_filename);

  if (do_plot)
    app.Run(kTRUE);
}
