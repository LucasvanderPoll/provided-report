from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from numpy import arange, min, max
import scipy.io as sio
import csv
import os

# cwd = os.getcwd()
# mat_data = sio.loadmat(os.path.join(cwd,'data/matlab.mat'))
mat_contents = sio.loadmat('../data/matlab.mat')
structure = mat_contents['flightdata']
flightdata = structure[0, 0]

lh_engine_FU = flightdata['lh_engine_FU'][0][0][0] / 2.20462                # convert lbs to kg
rh_engine_FU = flightdata['rh_engine_FU'][0][0][0] / 2.20462                # convert lbs to kg
lh_engine_FMF = flightdata['lh_engine_FMF'][0][0][0] / 2.20462 / 3600       # convert lbs/hr to kg/s
rh_engine_FMF = flightdata['rh_engine_FMF'][0][0][0] / 2.20462 / 3600       # convert lbs/hr to kg/s

t = arange(0, 48220, 1) # arange(0, 48320, 1)

BEW = 9165.0 / 2.20462
W_f = 2700.0 / 2.20462 # 4050.0 / 2.20462
W_p = (80.0 + 102.0 + 101.0 + 76.0 + 84.0 + 81.0 + 76.0 + 85.0 + 86.0) # (68.0 + 78.0 + 75.0 + 66.0 + 61.0 + 86.0 + 74.0 + 92.0 + 95.0)
W_zf = BEW + W_p
W_ramp = BEW + W_f + W_p


def get_fuel_mass(instant):
    fuel_mass = W_f - (lh_engine_FU[instant] + rh_engine_FU[instant])
    return fuel_mass


def get_mass(instant):
    total_mass = W_ramp - (lh_engine_FU[instant] + rh_engine_FU[instant])
    return total_mass


def get_fuel_moment(fuel_mass_input):               # input in kg

    f_mass = fuel_mass_input * 2.20462              # convert back to lbs

    mass = []
    moment = []

    with open('../data/fuelmasses.txt', 'r') as mass_data:
        masses = csv.reader(mass_data)
        for row in masses:
            mass.append(float(row[0]))

    with open('../data/fuelarms.txt', 'r') as moment_data:
        moments = csv.reader(moment_data)
        for row in moments:
            moment.append(float(row[0]))

    spline = interp1d(mass, moment, kind='linear', fill_value='extrapolate')    # create extra/interpolation spline
    fuel_arm_imp = spline.__call__(f_mass)                                      # get moment arm as (pound*inch)/100
    fuel_arm = fuel_arm_imp * 100 * 0.0254 / 2.20462                            # get moment arm as (kg*m)
    return fuel_arm.tolist()


def get_xcg(fuel_mass, moment_payload):
    moment_BEW = 2672953.5 * 0.0254 / 2.20462
    f_moment = get_fuel_moment(fuel_mass)
    xcg = (moment_BEW + moment_payload + f_moment) / (BEW + W_p + fuel_mass)
    return xcg


if __name__ == "__main__":
    W_op_raw = []
    for i in t:
        W_op_raw.append(get_mass(i))

    W_op = []
    for m in W_op_raw:
        nice_m = m.tolist()
        W_op.append(nice_m)
    moment_BEW = 2672953.5 * 0.0254 / 2.20462

    moment_payload = 4163.2378 # 3734.0848

    moment_zf = moment_BEW + moment_payload
    moment_fuel = get_fuel_moment(W_f)
    moment_ramp = moment_zf + moment_fuel

    xgc_BEW = 291.65 * 0.0254
    xgc_zf = moment_zf / W_zf

    print(f'BEW: {BEW}, moment_BEW: {moment_BEW}, xgc_BEW: {xgc_BEW}')
    print(f'W_p: {W_p}, moment_payload: {moment_payload}')
    print(f'W_zf: {W_zf}, moment_zf: {moment_zf}, xgc_zf: {xgc_zf}')
    print(f'W_f: {W_f}, moment_fuel: {moment_fuel}')
    print(f'W_ramp: {W_ramp}, moment_ramp: {moment_ramp}, xcg_ramp: {moment_ramp/W_ramp}')

    plot_all = True

    if plot_all:
        f_masses = []

        for i in t:
            f_masses.append(get_fuel_mass(i))

        print((moment_BEW+4020.2412+get_fuel_moment(f_masses[31570][0]))/(BEW+W_p+f_masses[31570][0]))

        Xcg = []
        for i in t:
            mass = f_masses[i][0]
            f_moment = get_fuel_moment(mass)
            xcg = (moment_BEW + moment_payload + f_moment) / (BEW + W_p + mass)
            Xcg.append(xcg)

        x = [min(Xcg), max(Xcg)]
        y = [min(W_op), max(W_op)]

        instants = [11480, 12880, 14170, 15550, 17780, 19110]

        for instant in instants:
            print('W for instant', instant, '=', W_op[instant], 'Xcg for instant', instant, '=', Xcg[instant])

        fig1 = plt.figure(figsize=(9, 4), dpi=100)
        plt.plot(t, W_op, color='royalblue', linewidth='1.5')
        axes = plt.gca()
        axes.set_xlim([0, t[-1]])
        axes.set_ylim([W_op[t[-1]][0], W_ramp]) # axes.set_ylim([6053.929947859395, W_ramp])
        fig1.suptitle('Gross Aircraft Mass Change', fontsize=14)
        plt.xlabel('Time Index', fontsize=10)
        plt.ylabel('Aircraft Operative Mass [Kg]', fontsize=10)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='grey', linestyle='-')
        plt.grid(b=True, which='minor', color='lightgrey', linestyle='dotted')
        plt.savefig('../plots/mass_change.png', dpi=600, orientation='portrait')
        plt.show()

        fig2 = plt.figure(figsize=(9, 4), dpi=100)
        plt.plot(Xcg, W_op, color='royalblue', linewidth='1.5')
        plt.plot(x, y, color='lightslategrey', linestyle='dotted')
        axes = plt.gca()
        axes.set_xlim(min(Xcg), max(Xcg))
        axes.set_ylim(min(W_op), max(W_op))
        fig2.suptitle('Longitudinal Position of Center of Gravity', fontsize=14)
        plt.ylabel('Aircraft Operative Mass [Kg]', fontsize=10)
        plt.xlabel('Center of Gravity [m]', fontsize=10)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='grey', linestyle='-')
        plt.grid(b=True, which='minor', color='lightgrey', linestyle='dotted')
        plt.savefig('../plots/cg_position.png', dpi=600, orientation='portrait')
        plt.show()

        fig3 = plt.figure(figsize=(9, 4), dpi=100)
        plt.plot(t, Xcg, color='royalblue', linewidth='1.5')
        plt.plot(x, y, color='lightslategrey', linestyle='dotted')
        axes = plt.gca()
        axes.set_xlim([0, t[-1]])
        axes.set_ylim(min(Xcg), max(Xcg))
        fig3.suptitle('Change in Position of Center of Gravity', fontsize=14)
        plt.ylabel('Center of Gravity [m]', fontsize=10)
        plt.xlabel('Time Index', fontsize=10)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='grey', linestyle='-')
        plt.grid(b=True, which='minor', color='lightgrey', linestyle='dotted')
        plt.savefig('../plots/cg_shift.png', dpi=600, orientation='portrait')
        plt.show()
