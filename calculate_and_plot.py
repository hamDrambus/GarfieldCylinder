#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import os
import subprocess
from scipy.interpolate import CubicSpline, PchipInterpolator, make_interp_spline
import matplotlib.pyplot as plt
from pathlib import Path as path
import numpy as np
import math

"""
All calculations here are done in cm, kV, and ns.
"""

class DriftMedium:
    def __init__(self, transport_pars_folder):
        fname = transport_pars_folder / 'drift_velocity.dat'
        xs, ys = np.loadtxt(fname, comments='//', delimiter='\t', unpack=True)
        xs, ys = xs * 1e4, ys * 0.1 # Scale from Geant4 units to kV/cm, cm/ns
        self.Vdrift = make_interp_spline(xs, ys, k=1)

        fname = transport_pars_folder / 'diffusion_longitudinal.dat'
        xs, ys = np.loadtxt(fname, comments='//', delimiter='\t', unpack=True)
        xs, ys = xs * 1e4, ys * 0.01 # Scale from Geant4 units to kV/cm, cm^2/ns
        self.DiffL = make_interp_spline(xs, ys, k=1)

        fname = transport_pars_folder / 'diffusion_transverse.dat'
        xs, ys = np.loadtxt(fname, comments='//', delimiter='\t', unpack=True)
        xs, ys = xs * 1e4, ys * 0.01 # Scale from Geant4 units to kV/cm, cm^2/ns
        self.DiffT = make_interp_spline(xs, ys, k=1)

    def get_Vdrift(self, field):
        Vdrift = self.Vdrift
        if type(field) in [list, np.ndarray]:
            return [Vdrift(f) if f is not None else None for f in field]
        return Vdrift(field) if field is not None else None
    
    def get_DiffL(self, field):
        DiffL = self.DiffL
        if type(field) in [list, np.ndarray]:
            return [DiffL(f) if f is not None else None for f in field]
        return DiffL(field) if field is not None else None
    
    def get_DiffT(self, field):
        DiffT = self.DiffT
        if type(field) in [list, np.ndarray]:
            return [DiffT(f) if f is not None else None for f in field]
        return DiffT(field) if field is not None else None

class CylinderSetup:
    def __init__(self, Rmin, Rmax, drift_medium):
        self.Rmin = min(Rmin, Rmax)
        self.Rmax = max(Rmax, Rmin)
        self.drift_medium = drift_medium

    def get_field(self, radius, V0):
        Rmin = self.Rmin
        Rmax = self.Rmax
        logR = math.log(Rmax/Rmin)
        if type(radius) in [list, np.ndarray]:
            return [V0/logR/r if r <= Rmax and r >= Rmin else None for r in radius]
        return V0/logR/radius if radius <= Rmax and radius >= Rmin else None
    
    def get_max_field(self, V0):
        return self.get_field(self.Rmin, V0)
    
    def get_min_field(self, V0):
        return self.get_field(self.Rmax, V0)

    def get_Vdrift(self, radius, V0):
        field = self.get_field(radius, V0)
        return self.drift_medium.get_Vdrift(field)

    def get_field_gradient(self, radius, V0):
        Rmin = self.Rmin
        Rmax = self.Rmax
        logR = math.log(Rmax/Rmin)
        if type(radius) in [list, np.ndarray]:
            return [V0/logR/r/r if r <= Rmax and r >= Rmin else None for r in radius]
        return V0/logR/radius/radius if radius <= Rmax and radius >= Rmin else None
    
    def get_max_field_gradient(self, feild):
        Rmin = self.Rmin
        if type(feild) in [list, np.ndarray]:
            return [f/Rmin if f is not None else None for f in feild]
        return field/Rmin if field is not None else None
    
    def get_charge_signal(self, V0, rel_precision = 0.01):
        Rmin = self.Rmin
        Rmax = self.Rmax
        num_split = Rmax * (1.0/Rmin - 1.0/Rmax) / rel_precision
        num_split = int(round(num_split) + 1)
        radii = np.linspace(1.0/Rmax, 1.0/Rmin, num=num_split) # uniform spacing over E, that is 1/r
        radii = [1/r for r in radii] # Sorted from Rmax to Rmin
        fields = self.get_field(radii, V0)
        v_drift = self.get_Vdrift(radii, V0)
        signal = []
        times = []
        time = 0.0
        r_prev = radii[0]
        for r, f, vdr in zip(radii, fields, v_drift):
            time += (r_prev - r) / vdr
            signal.append(f * vdr)
            times.append(time)
            r_prev = r
        return times, signal, radii
    
class ParallelPlateSetup:
    def __init__(self, L, drift_medium):
        self.L = abs(L)
        self.drift_medium = drift_medium

    def get_field(self, x, V0):
        L = self.L
        if type(x) in [list, np.ndarray]:
            return [V0/L if r <= L and r >= 0 else None for r in x]
        return V0/L if x <= L and x >= 0 else None
    
    def get_max_field(self, V0):
        return V0/self.L
    
    def get_min_field(self, V0):
        return V0/self.L

    def get_Vdrift(self, x, V0):
        field = self.get_field(x, V0)
        return self.drift_medium.get_Vdrift(field)

    def get_field_gradient(self, x, V0):
        L = self.L
        if type(x) in [list, np.ndarray]:
            return [0 if r <= L and r >= 0 else None for r in x]
        return 0 if x <= L and x >= 0 else None
    
    def get_max_field_gradient(self, feild):
        L = self.L
        if type(feild) in [list, np.ndarray]:
            return [0 if f is not None else None for f in feild]
        return 0 if field is not None else None
    
    def get_charge_signal(self, V0, rel_precision = 0.01):
        x_min = -0.01
        x_max = self.L + 0.01
        num_split = (x_max - x_min) / rel_precision
        num_split = int(round(num_split) + 1)
        xs = np.linspace(x_min, x_max, num=num_split)
        fields = self.get_field(xs, V0)
        v_drift = self.get_Vdrift(xs, V0)
        signal = []
        times = []
        time = 0.0
        x_prev = xs[0]
        for x, f, vdr in zip(xs, fields, v_drift):
            if f is None:
                time += 0.1
                signal.append(0)
            else:
                time += (x - x_prev) / vdr
                signal.append(f * vdr)
            times.append(time)
            x_prev = x
        return times, signal, radii
    
class CylinderSetupGarfield:
    def __init__(self, Rmin, Rmax, XS_name, WireOffset=0.0):
        self.Rmin = min(Rmin, Rmax)
        self.Rmax = max(Rmax, Rmin)
        self.XS_name = XS_name
        self.WireOffset = WireOffset

    def get_field(self, radius, V0):
        Rmin = self.Rmin
        Rmax = self.Rmax
        logR = math.log(Rmax/Rmin)
        if type(radius) in [list, np.ndarray]:
            return [V0/logR/r if r <= Rmax and r >= Rmin else None for r in radius]
        return V0/logR/radius if radius <= Rmax and radius >= Rmin else None
    
    def get_max_field(self, V0):
        return self.get_field(self.Rmin, V0)
    
    def get_min_field(self, V0):
        return self.get_field(self.Rmax, V0)

    def get_field_gradient(self, radius, V0):
        Rmin = self.Rmin
        Rmax = self.Rmax
        logR = math.log(Rmax/Rmin)
        if type(radius) in [list, np.ndarray]:
            return [V0/logR/r/r if r <= Rmax and r >= Rmin else None for r in radius]
        return V0/logR/radius/radius if radius <= Rmax and radius >= Rmin else None
    
    def get_max_field_gradient(self, feild):
        Rmin = self.Rmin
        if type(feild) in [list, np.ndarray]:
            return [f/Rmin if f is not None else None for f in feild]
        return field/Rmin if field is not None else None
    
    def get_charge_signal(self, V0, rel_precision = 0.001, use_diffusion=True, shaping=0.0):
        N = round(1.0 / abs(rel_precision))
        use_diff = 'on' if use_diffusion else 'off'
        args = f"-V0={V0:.1f} -N={N:d} -rWire={(self.Rmin*1e4):.0f} " + \
               f"-rTube={self.Rmax:.2f} -rWireOffset={(self.WireOffset*1e4):.0f} " + \
               f"-diff={use_diff} -XS={self.XS_name} -sh={shaping:.1f} -out_fname=tmp.txt"
        # print(f"./CylinderCharge {args}")
        p1 = subprocess.Popen(['./CylinderCharge'] + args.split(), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        p1.communicate()
        if p1.returncode is None:
            raise ChildProcessError(f"ERROR while running 'CylinderCharge {args}'.")
        if p1.returncode != 0:
            raise ChildProcessError(f"'CylinderCharge {args}'\nreturned {p1.returncode}.")
        times, signal, _1, _2 = np.loadtxt('tmp.txt', delimiter=',', unpack=True, skiprows=2)
        path('tmp.txt').unlink()
        return times, signal
    
def trimmed_tail(xs, ys, *args, ys_rel_threshold = 1e-3, xs_rel_margin = 0.05):
    max_y = max(abs(np.max(ys)), abs(np.min(ys)))
    th = max_y * ys_rel_threshold
    t0 = xs[0]
    t1 = xs[-1]
    for x, y in zip(reversed(xs), reversed(ys)):
        if abs(y) > th:
            t1 = x
            break
    t1 = (t1 - t0) * (1.0 + xs_rel_margin)
    xys = [(x, y) for x, y in zip(xs, ys) if x <= t1]
    out_xs, out_ys = zip(*xys)
    return out_xs, out_ys, *args

def renormalized(xs, ys, *args, norm = 1.0):
    norm1 = np.trapz(ys, x=xs)
    out_ys = [y * (norm / norm1) for y in ys]
    return xs, out_ys, *args

if __name__ == "__main__":
    r_wire = 400 * 1.e-4 # Radii are in cm
    r_tube = 0.5
    input_folder = path.home() / 'Documents' / 'Detector_geant4' / 'NBrS_THGEM_LAr_v0' / 'data_cache_LAr_v5.4_exact_Var9'
    output_folder = path('LAr_v5.4_exact_Var9')
    LAr = DriftMedium(input_folder)
    setup = CylinderSetup(r_wire, r_tube, LAr)
    plate_setup = ParallelPlateSetup(0.2, LAr)
    garfield_setup = CylinderSetupGarfield(r_wire, r_tube, 'Var9')
    garfield_setup_askew = CylinderSetupGarfield(r_wire, r_tube, 'Var9', WireOffset=r_wire/2.0)

    input_folder1 = path.home() / 'Documents' / 'Detector_geant4' / 'NBrS_THGEM_LAr_v0' / 'data_cache_LAr_v5.4_exact_Var8'
    output_folder1 = path('LAr_v5.4_exact_Var8')
    LAr1 = DriftMedium(input_folder1)
    setup1 = CylinderSetup(r_wire, r_tube, LAr1) 
    plate_setup1 = ParallelPlateSetup(0.2, LAr1)
    garfield_setup1 = CylinderSetupGarfield(r_wire, r_tube, 'Var8')
    garfield_setup1_askew = CylinderSetupGarfield(r_wire, r_tube, 'Var8', WireOffset=r_wire/2.0)

    fugure_sz = (14, 9)
    
    V0s = [8, 14, 20]
    fig, grad = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    grad = grad.flatten()
    grad[0].set_title(fr"Field gradient in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm")
    grad[0].set_xlabel(fr"E (kV/cm)")
    grad[0].set_ylabel(fr"dE/dr (kV/cm$^2$)")

    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(fr"Induced current (charge signal) in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")

    fig, position = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    position = position.flatten()
    position[0].set_title(fr"Electron drift in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm")
    position[0].set_xlabel(fr"Time (ns)")
    position[0].set_ylabel(fr"Position (cm)")

    fig, field = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    field = field.flatten()
    field[0].set_title(fr"Field affecting electrons in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm")
    field[0].set_xlabel(fr"Time (ns)")
    field[0].set_ylabel(fr"Field (kV/cm)")

    for V0 in V0s:
        times, sig, radii = setup.get_charge_signal(V0, rel_precision=0.01)
        fields = setup.get_field(radii, V0)
        gadients = setup.get_field_gradient(radii, V0)
        E_max = setup.get_max_field(V0)
        grad[0].plot(fields, gadients, '-', label=f"Var9, V0={V0}kV")
        signal[0].plot(times, sig, '-', label=f"Var9, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
        position[0].plot(times, radii, '-', label=f"Var9, V0={V0}kV")
        field[0].plot(times, fields, '-', label=f"Var9, V0={V0}kV")
        norm = np.trapz(sig, x=times)
        times1, sig1, radii1 = renormalized(*setup1.get_charge_signal(V0, rel_precision=0.01), norm=norm)
        signal[0].plot(times1, sig1, '--', label=f"Var8, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")

    fields = np.linspace(1, 200, num=300)
    max_gradient = setup.get_max_field_gradient(fields)
    fig, max_grad = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    max_grad = max_grad.flatten()
    max_grad[0].set_title(fr"Maximum field gradient in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm")
    max_grad[0].set_xlabel(fr"E (kV/cm)")
    max_grad[0].set_ylabel(fr"dE/dr (kV/cm$^2$)")
    max_grad[0].plot(fields, max_gradient, '-', label=f"Var9")

    grad[0].legend(loc='best')
    signal[0].legend(loc='best')
    position[0].legend(loc='best')
    field[0].legend(loc='best')

    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(fr"Induced current (charge signal) in parallel plate setup L={plate_setup.L}cm")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = plate_setup.get_max_field(V0)
        times, sig, radii = plate_setup.get_charge_signal(V0, rel_precision=0.001)
        norm = np.trapz(sig, x=times)
        signal[0].plot(times, sig, '-', label=f"Var9, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
        times1, sig1, radii1 = renormalized(*plate_setup1.get_charge_signal(V0, rel_precision=0.001), norm=norm)
        signal[0].plot(times1, sig1, '--', label=f"Var8, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
    signal[0].legend(loc='best')

    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(f"Induced current (charge signal) in cylindrical chamber Rmin={setup.Rmin*1e4:.0f}um, Rmax={setup.Rmax}cm\n \
                        Testing Garfield++ vs Python. Cross sections Var9.")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = setup.get_max_field(V0)
        times, sig, radii = trimmed_tail(*setup.get_charge_signal(V0, rel_precision=0.001))
        norm = np.trapz(sig, x=times)
        signal[0].plot(times, sig, '-.', label=f"Python, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
        times, sig = renormalized(*trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001, use_diffusion=False)), norm=norm)
        signal[0].plot(times, sig, '--', label=f"Garfield, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
        times, sig = renormalized(*trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001, use_diffusion=True)), norm=norm)
        signal[0].plot(times, sig, '-', label=f"Garfield with diffusion, V0 = {V0} kV, Emax = {E_max:.0f} kV/cm")
    signal[0].legend(loc='best')

    V0s = [10]
    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(f"Induced current (charge signal) in cylindrical chamber Rmin={garfield_setup.Rmin*1e4:.0f}um, Rmax={garfield_setup.Rmax}cm\n \
                        Garfield++ with diffusion, V0 = {V0s[0]} kV, E_max = {garfield_setup.get_max_field(V0s[0]):.0f} kV/cm")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = garfield_setup.get_max_field(V0)
        times, sig = trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001))
        norm = -np.trapz(sig, x=times)
        times, sig = renormalized(times, sig, norm=norm)
        signal[0].plot(times, sig, '-', label=f"Var9")
        times, sig = renormalized(*trimmed_tail(*garfield_setup_askew.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, '-.', label=f"Var9, wire offset {garfield_setup_askew.WireOffset * 1e4:.0f} um")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, '--', label=f"Var8")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1_askew.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, ':', label=f"Var8, wire offset {garfield_setup1_askew.WireOffset * 1e4:.0f} um")
    signal[0].legend(loc='best')

    V0s = [20]
    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(f"Induced current (charge signal) in cylindrical chamber Rmin={garfield_setup.Rmin*1e4:.0f}um, Rmax={garfield_setup.Rmax}cm\n \
                        Garfield++ with diffusion, V0 = {V0s[0]} kV, E_max = {garfield_setup.get_max_field(V0s[0]):.0f} kV/cm")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = garfield_setup.get_max_field(V0)
        times, sig = trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001))
        norm = -np.trapz(sig, x=times)
        times, sig = renormalized(times, sig, norm=norm)
        signal[0].plot(times, sig, '-', label=f"Var9")
        times, sig = renormalized(*trimmed_tail(*garfield_setup_askew.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, '-.', label=f"Var9, wire offset {garfield_setup_askew.WireOffset * 1e4:.0f} um")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, '--', label=f"Var8")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1_askew.get_charge_signal(V0, rel_precision=0.001)), norm=norm)
        signal[0].plot(times, sig, ':', label=f"Var8, wire offset {garfield_setup1_askew.WireOffset * 1e4:.0f} um")
    signal[0].legend(loc='best')

    V0s = [20]
    sh = 10 # ns
    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(f"Induced current (charge signal) in cylindrical chamber Rmin={garfield_setup.Rmin*1e4:.0f}um, Rmax={garfield_setup.Rmax}cm\n\
                        Garfield++ with diffusion, V0 = {V0s[0]} kV, E_max = {garfield_setup.get_max_field(V0s[0]):.0f} kV/cm\n\
                        Amplifier shaping = {sh} ns.")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = garfield_setup.get_max_field(V0)
        times, sig = trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001, shaping=sh))
        norm = -np.trapz(sig, x=times)
        times, sig = renormalized(times, sig, norm=norm)
        signal[0].plot(times, sig, '-', label=f"Var9")
        times, sig = renormalized(*trimmed_tail(*garfield_setup_askew.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, '-.', label=f"Var9, wire offset {garfield_setup_askew.WireOffset * 1e4:.0f} um")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, '--', label=f"Var8")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1_askew.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, ':', label=f"Var8, wire offset {garfield_setup1_askew.WireOffset * 1e4:.0f} um")
    signal[0].legend(loc='best')

    V0s = [20]
    sh = 50 # ns
    fig, signal = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    signal = signal.flatten()
    signal[0].set_title(f"Induced current (charge signal) in cylindrical chamber Rmin={garfield_setup.Rmin*1e4:.0f}um, Rmax={garfield_setup.Rmax}cm\n\
                        Garfield++ with diffusion, V0 = {V0s[0]} kV, E_max = {garfield_setup.get_max_field(V0s[0]):.0f} kV/cm\n\
                        Amplifier shaping = {sh} ns.")
    signal[0].set_xlabel(fr"Time (ns)")
    signal[0].set_ylabel(fr"Signal (arb. u.)")
    for V0 in V0s:
        E_max = garfield_setup.get_max_field(V0)
        times, sig = trimmed_tail(*garfield_setup.get_charge_signal(V0, rel_precision=0.001, shaping=sh))
        norm = -np.trapz(sig, x=times)
        times, sig = renormalized(times, sig, norm=norm)
        signal[0].plot(times, sig, '-', label=f"Var9")
        times, sig = renormalized(*trimmed_tail(*garfield_setup_askew.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, '-.', label=f"Var9, wire offset {garfield_setup_askew.WireOffset * 1e4:.0f} um")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, '--', label=f"Var8")
        times, sig = renormalized(*trimmed_tail(*garfield_setup1_askew.get_charge_signal(V0, rel_precision=0.001, shaping=sh)), norm=norm)
        signal[0].plot(times, sig, ':', label=f"Var8, wire offset {garfield_setup1_askew.WireOffset * 1e4:.0f} um")
    signal[0].legend(loc='best')
    
    plt.tight_layout()
    plt.show()
