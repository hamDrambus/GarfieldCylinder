#include <Garfield/FundamentalConstants.hh>
#include <Garfield/GarfieldConstants.hh>
#include <Garfield/Numerics.hh>
#include "MaterialCustomDrift.hh"

namespace {
double Interpolate_(const std::vector<double> &x, const std::vector<double> &y,
                  const double xx, const unsigned int order) {
  if (xx < x.front())
    return y.front();
  if (xx > x.back())
    return y.back();
  if (order > 1)
    return Garfield::Numerics::Divdif(y, x, x.size(), xx, order);
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend())
    return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < Garfield::Small)
    return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

} // namespace

namespace Garfield {

MaterialCustomDrift::MaterialCustomDrift(double v_drift) :
    Medium(), m_v_drift(v_drift)
{
  m_className = "MaterialCustonDrift";
  m_driftable = true; 
}

MaterialCustomDrift::MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift) :
    MaterialCustomDrift(0)
{
  std::size_t sz = std::min(drift_fields.size(), v_drift.size());
  if (sz == 0) {
    m_v_drift = 1e-4;
    if (m_debug)
      std::cerr<<__PRETTY_FUNCTION__<<": empty input velocity data. Setting drift velocity to constant (1e-4 cm/ns)."<<std::endl;
    return;
  }
  if (m_debug && drift_fields.size() != v_drift.size())
    std::cerr<<__PRETTY_FUNCTION__<<": warning: input velocity data size mismatch."<<std::endl;
  m_v_drift_e.reserve(sz);
  m_v_drift_v.reserve(sz);
  m_v_drift_e = std::vector<double>(drift_fields.begin(), drift_fields.begin()+sz);
  m_v_drift_v = std::vector<double>(v_drift.begin(), v_drift.begin()+sz);
}

MaterialCustomDrift::MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift,
            const std::vector<double> & diff_fields, const std::vector<double> & diff_coef) :
    MaterialCustomDrift(drift_fields, v_drift)
{
  std::size_t sz = std::min(diff_fields.size(), diff_coef.size());
  if (sz == 0) {
    if (m_debug)
      std::cerr<<__PRETTY_FUNCTION__<<": warning: empty diffusion coefficients data. Diffusion is off."<<std::endl;
    return;
  }
  if (m_debug && diff_fields.size() != diff_coef.size())
    std::cerr<<__PRETTY_FUNCTION__<<": warning: input diffusion data size mismatch."<<std::endl;
  m_diff_long_e.reserve(sz);
  m_diff_long_v.reserve(sz);
  m_diff_long_e = std::vector<double>(diff_fields.begin(), diff_fields.begin()+sz);
  m_diff_long_v = std::vector<double>(diff_coef.begin(), diff_coef.begin()+sz);
  for (std::size_t i = 0; i!=sz; ++i) {
    double E = m_diff_long_e[i];
    double vel = m_v_drift_e.empty() ? m_v_drift : Interpolate_(m_v_drift_e, m_v_drift_v, E, 1);
    // Recalculate normal diffusion coefficients (in cm^2/ns) to Garfield++ ones (cm^(1/2))
    m_diff_long_v[i] = sqrt(2.0 * m_diff_long_v[i] / std::fabs(vel));
  }
  m_diff_trans_e = m_diff_long_e;
  m_diff_trans_v = m_diff_long_v;
}

MaterialCustomDrift::MaterialCustomDrift(const std::vector<double> & drift_fields, const std::vector<double> & v_drift,
            const std::vector<double> & diff_long_fields, const std::vector<double> & diff_long_coef,
            const std::vector<double> & diff_trans_fields, const std::vector<double> & diff_trans_coef) :
     MaterialCustomDrift(drift_fields, v_drift, diff_long_fields, diff_long_coef)
{
  std::size_t sz = std::min(diff_trans_fields.size(), diff_trans_coef.size());
  if (sz == 0) {
    if (m_debug)
      std::cerr<<__PRETTY_FUNCTION__<<": warning: empty transverse diffusion coefficients data. Diffusion is off."<<std::endl;
    return;
  }
  if (m_debug && diff_trans_fields.size() != diff_trans_coef.size())
    std::cerr<<__PRETTY_FUNCTION__<<": warning: input transverse diffusion data size mismatch."<<std::endl;
  m_diff_trans_e.reserve(sz);
  m_diff_trans_v.reserve(sz);
  m_diff_trans_e = std::vector<double>(diff_trans_fields.begin(), diff_trans_fields.begin()+sz);
  m_diff_trans_v = std::vector<double>(diff_trans_coef.begin(), diff_trans_coef.begin()+sz);
  for (std::size_t i = 0; i!=sz; ++i) {
    double E = m_diff_trans_e[i];
    double vel = m_v_drift_e.empty() ? m_v_drift : Interpolate_(m_v_drift_e, m_v_drift_v, E, 1);
    // Recalculate normal diffusion coefficients (in cm^2/ns) to Garfield++ ones (cm^(1/2))
    m_diff_trans_v[i] = sqrt(2.0 * m_diff_trans_v[i] / std::fabs(vel));
  }
}

bool MaterialCustomDrift::ElectronVelocity(const double ex, const double ey,
                                          const double ez, const double bx,
                                          const double by, const double bz, double& vx,
                                          double& vy, double& vz)
{
  vx=vy=vz=0.0;
  if (!((0==bx)&&(0==by)&&(0==bz))) {
     std::cerr<<__PRETTY_FUNCTION__<<": error: magnetic field is not supported."<<std::endl;
     return false;
  }
  double E = std::sqrt(ex*ex+ey*ey+ez*ez);
  double V = m_v_drift_e.empty() ? m_v_drift : Interpolate_(m_v_drift_e, m_v_drift_v, E, 1);
  vx = -ex*V/E;
  vy = -ey*V/E;
  vz = -ez*V/E;
  return true;
}

bool MaterialCustomDrift::ElectronDiffusion(const double ex, const double ey,
                                            const double ez, const double bx,
                                            const double by, const double bz, double& dl,
                                            double& dt)
{
  dt = dl = 0.0;
  if (!((0==bx)&&(0==by)&&(0==bz))) {
     std::cerr<<__PRETTY_FUNCTION__<<": error: magnetic field is not supported."<<std::endl;
     return false;
  }
  if (m_diff_long_e.empty() && m_diff_trans_e.empty())
    return true;
  double E = std::sqrt(ex*ex+ey*ey+ez*ez);
  if (m_diff_long_e.empty()) {
    dl = dt = Interpolate_(m_diff_trans_e, m_diff_trans_v, E, 1);
    return true;
  }
  if (m_diff_trans_e.empty()) {
    dt = dl = Interpolate_(m_diff_long_e, m_diff_long_v, E, 1);
    return true;
  }
  dt = Interpolate_(m_diff_trans_e, m_diff_trans_v, E, 1);
  dl = Interpolate_(m_diff_long_e, m_diff_long_v, E, 1);
  return true;
}

}
