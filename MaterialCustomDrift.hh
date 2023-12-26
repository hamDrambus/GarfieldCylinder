#ifndef MATERIAL_CUSTOM_DRIFT_HH
#define MATERIAL_CUSTOM_DRIFT_HH

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Garfield/Medium.hh>

namespace Garfield {

class MaterialCustomDrift : public Medium
{
protected:
  double m_v_drift;
  std::vector<double> m_v_drift_e, m_v_drift_v;
  std::vector<double> m_diff_long_e, m_diff_long_v;
  std::vector<double> m_diff_trans_e, m_diff_trans_v;
public:
  /// Default constructor for drift velocity 1e5 cm/s for all fields and zero diffusion
  MaterialCustomDrift(double v_drift = 1e-4);
  /// Constructor for tabulated drift velocity (in cm/ns) vs fields (in V/cm) and zero diffusion
  MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift);
  /// Constructor for tabulated drift velocity (in cm/ns) vs fields (in V/cm)
  /// and tabulated isotropic (D_L = D_T) diffusion (in cm^2/ns) coefficients vs electric field.
  MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift,
              const std::vector<double> & diff_fields, const std::vector<double> & diff_coef);
  /// Constructor for tabulated drift velocity (in cm/ns) vs electric fields (in V/cm)
  /// and tabulated longitudinal (D_L) and transverse diffusion coefficients (in cm^2/ns).
  MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift,
              const std::vector<double> & diff_long_fields, const std::vector<double> & diff_long_coef,
              const std::vector<double> & diff_trans_fields, const std::vector<double> & diff_trans_coef);
  ~MaterialCustomDrift() = default;

  // Transport parameters for electrons
  /// Drift velocity [cm / ns]
  virtual bool ElectronVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz) override;
  /// Longitudinal and transverse diffusion coefficients [cm1/2]
  virtual bool ElectronDiffusion(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& dl,
                                 double& dt) override;
};

}
#endif // MATERIAL_CUSTOM_DRIFT_HH
