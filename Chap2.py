import matplotlib.pyplot as plt
import numpy as np

# ######## Example 2.1 ############

Vm = 100
thetav = 0   # Volage amplitude and phase angle
Z = 1.25
gama = 60   # Impedance magnitude and phase angle
def runExample21(Vm,thetav,Z,gama):
    thetai = thetav - gama
    theta = (thetav-thetai)*np.pi/180   # Degree to radian
    Im = Vm/Z                           # Current amplitude
    wt = np.arange(0.0, 2*np.pi, 0.05)  # wt from 0 to 2*pi
    v = Vm*np.cos(wt)                   # Instantaneous voltage
    i = Im*np.cos(wt+thetai*np.pi/180)  # Instantaneous current
    p = v*i                             # Instantaneous power
    V = Vm/np.sqrt(2)                   # rms voltage and current
    I = Im/np.sqrt(2)
    P = V*I*np.cos(theta)               # Average power
    Q = V*I*np.sin(theta)               # Reactive power
    S = P+1j*Q                          # Complex power
    pr = P*(1+np.cos(2*(wt+thetav)))    # Eq. (2.6)
    px = Q*np.sin(2*(wt+thetav))        # Eq. (2.8)
    PP = P*np.ones(len(wt))             # Average power of length w for plot
    xline = np.zeros(len(wt))           # Generates a zero vector
    wt = wt*180/np.pi                   # Converting radians to degree

    # Creating the figure
    fig, axes = plt.subplots(2, 2)
    for line in [v, i, xline]:
        axes[0, 0].plot(wt, line)
    axes[0, 0].set_title("v(t)=Vm cos(wt), i(t)=Im cos(wt+{0})".format(thetai))
    axes[0, 0].set_xlabel("wt, degree")

    for line in [p, xline]:
        axes[0, 1].plot(wt, line)
    axes[0, 1].set_title("p(t)=v(t)*i(t)")
    axes[0, 1].set_xlabel("wt, degree")

    for line in [pr, PP, xline]:
        axes[1, 0].plot(wt, line)
    axes[1, 0].set_title("pr(t)  -  Eq. 2.6")
    axes[1, 0].set_xlabel("wt, degree")

    for line in [px, xline]:
        axes[1, 1].plot(wt, line)
    axes[1, 1].set_title("px(t)  -  Eq. 2.8")
    axes[1, 1].set_xlabel("wt, degree")

    fig.subplots_adjust(
        wspace=0.3, hspace=0.5)  # Legger til mellomrom mellom figurer.
    return fig
fig = runExample21(Vm,thetav,Z,gama)

#plt.show()


# ######### Example 2.2 - Complex power balance ########
print("Example 2.2")
V = 1200 + 1j * 0  # Input Voltage
Z1 = 60 + 1j * 0   # Input impedance 1
Z2 = 6 + 1j * 12   # Input impedance 2
Z3 = 30 - 1j * 30  # Input impedance 3

I1, I2, I3 = [V / Z for Z in [Z1, Z2, Z3]]           # I = V/Z
print("I: ", I1, I2, I3)
S1, S2, S3 = [V * np.conj(I) for I in [I1, I2, I3]]  # S = U*I'
print("S: ", S1, S2, S3)
print("S sum:", S1+S2+S3)


# ######### Example 2.3 - Power factor correction #######
print("Example 2.3")
V = 200 + 1j * 0
Z1 = 100 + 1j * 0
Z2 = 10 + 1j * 20
# Find total real power, total reactive power, PF, total current
S1, S2 = [V**2 / np.conj(Z) for Z in [Z1, Z2]]
Stot = S1 + S2
print(Stot)
Ptot = np.real(Stot)
Qtot = np.imag(Stot)
Itot = np.conj(Stot) / np.conj(V)
print("Total real P: ", Ptot)
print("Total reactive Q: ", Qtot)
print("Total current I:", Itot)
print("PF: ", np.cos(np.angle(Stot)))

# Find the capacitance of the capacitor connected across the loads to
# improve the overall power factor to 0.8 lagging.

# Total power is equal, both for 0,6 and 0,8 power factor.
newTheta = np.arccos(0.8)
newQtot = Ptot * np.tan(newTheta)
Qc = Qtot - newQtot
Zc = V**2 / (1j * Qc)
C = 10**6 / (2 * np.pi * 50 * Zc)           # myF: 10E6/2pi*f*Z
print("New theta", np.arccos(0.8), np.arccos(0.8) * 180 / np.pi)
print("New total Q: ", newQtot)
print("Zc: ", Zc)
print("Installed C:{0} microF".format(C))

Snew = Ptot + 1j * newQtot
Inew = np.conj(Snew) / np.conj(V)
print("New S:", Snew)
print("New I:{0}, ({1} v{2})".format(Inew, np.abs(Inew),
                                     np.angle(Inew) * 180 / np.pi))


# ###### Example 2.4 #########
print("Example 2.4")
V = 1400 + 1j * 0
freq = 60
# Load 1, inductiv, 125 kva at 0.28 pf
# Load 2, capacitiv, 10 kw and 40 kvar
# Load 3, Resistive load of 15 kwar
kv = 1000
theta = np.arccos(0.28)
S1 = (125*np.cos(theta)+1j*125*np.sin(theta))*kv
S2 = (10-1j*40)*kv
S3 = 15 * kv
print("S1:", S1)
print("S2:", S2)
print("S3:", S3)

# Find the total kW, kvar, kVA and the supply PF
Stot = S1 + S2 + S3
print("Total kW: ", np.real(Stot))
print("Total kvar: ", np.imag(Stot))
print("Total kVA: ", np.abs(Stot))
print("PF: ", np.real(Stot) / np.abs(Stot))

# Determine the kvar ratoing of a capacitor and the capacitance in
# microF to improve the total power factor to be 0.8 lagging.
newTheta = np.arccos(0.8)   # Finding the new angle with PF = 0.8
newQtot = np.real(Stot)*np.tan(newTheta)  # Finding the new total Q with PF 0.8
Qc = np.imag(Stot)-newQtot  # Installed Qc = Qold - Qnew
Xc = V**2/np.conj(0-1j*Qc)  # Convert to Z
C = 10**6/(2*np.pi*freq*np.abs(Xc))  # Convert from Z to C

print("newTheta: ", newTheta, newTheta*180/np.pi)
print("newQtot: ", newQtot)
print("Qc kvar: ", Qc)
print("Xc Ohms: ", Xc)
print("C microF: ", C)


# ########## Example 2.5 #################
# Two voltage sources over a line Z      #
def polarToCart(arg, angle=0, deg=False):
    if deg:
        angle = angle * np.pi / 180
    a = arg * (np.cos(angle) + 1j * np.sin(angle))
    return a


def cartToPolar(z, deg=False):
    scale = 1
    if deg:
        scale = 180 / np.pi
    a = [np.abs(z), np.angle(z) * scale]
    return a


V1 = polarToCart(120, -5, deg=True)
V2 = 100
Z = 1 + 1j * 7

I12 = (V1 - V2) / Z
I21 = (V2 - V1) / Z
print("I12: ", cartToPolar(I12, deg=True))
print("I21: ", cartToPolar(I21, deg=True))
S12 = V1 * np.conj(I12)
S21 = V2 * np.conj(I21)
print("S12: {0} W = {1} W".format(S12, cartToPolar(S12, deg=True)))
print("S21: {0} W = {1} W".format(S21, cartToPolar(S21, deg=True)))

# Line losses
Sloss = S12 + S21
print("Sloss: ", Sloss)

# ######### Example 2.6 #########
# Create input from to voltage sources, #1 with variable plus/minus 30
# degree from given value. Magnitude of V1, V2 and angle V2 is kept constant
# Compute complex power for each source and the line loss.
# Tabulate P1, P2, PL against V1 angle.
E1 = 120                    # Source #1 Voltage Mag
a1 = -5                     # Source #1 Phase angle
E2 = 100                    # Source #2 Voltage Mag
a2 = 0                      # Source #2 Phase angle
R = 1                       # Line Resistance
X = 7                       # Line Reactance
Z = R + 1j * X              # Line Impedance
a1 = np.arange(a1 - 30, a1 + 35, 5)  # Finding plus/minus 30 degr
a1r = a1*np.pi/180          # Convert degrees to radians
a2 = np.ones(len(a1))*a2    # Create col. array of same length for a2
a2r = a2*np.pi/180          # Convert degrees to radians
V1 = E1*np.cos(a1r) + 1j * E1*np.sin(a1r)
V2 = E2*np.cos(a2r) + 1j * E2*np.sin(a2r)
I12 = (V1-V2)/Z
I21 = (V2-V1)/Z
S1 = V1*np.conj(I12)        # Finding Apparent power S1
P1 = np.real(S1)
Q1 = np.imag(S1)
S2 = V2*np.conj(I21)        # Finding Apparent power S2
P2 = np.real(S2)
Q2 = np.imag(S2)
Sloss = S1+S2               # Finding Apparent power Sloss
Ploss = np.real(Sloss)
Qloss = np.imag(Sloss)
Result1 = [a1, P1, P2, Ploss]
print("{0} \t {1} \t {2} \t {3}".format("Delta 1", "P-1", "P-2", "Ploss"))
for i in np.arange(len(a1)):
    print("{0} \t {1:10.4f} \t {2:10.4f} \t {3:10.4f}".format(
        Result1[0][i], Result1[1][i], Result1[2][i], Result1[3][i]))
fig, ax = plt.subplots()
for i in Result1[1:]:
    ax.plot(a1, i)
ax.set_xlabel("Source #1 Voltage Phase Angle")
ax.set_ylabel("P, Watts")
ax.text(-26, -550, "P1")
ax.text(-26, 600, "P2")
ax.text(-26, 100, "Ploss")
ax.grid(True)
#plt.show()

# ######### Example 2.7 ############
# A three-phase line has an impedance of 2 + 4j Ohm. Line-to-Line Voltage
# magnitude is |V| = 207.85 V. The load connects to twho balanced three-phase
# loads, one Delta-connected and one Y-connected. The Y-connected has impedance
# of 30 + 40j per phase, and the Delta-connected has 60 - j45 Ohm.
# Find a: Current, real power and reactive power drawn from the supply
# Find b: The line voltage at the combined loads.
# Find c: the current per phase in each load
# Find d: The total real and reactive powers in each load and the line.
print("Example 2.7")
Vsup = 207.85                   # Supply voltage magnitude
Zline = 2+1j*4                  # Line impedance
Zload1_delta = 60 - 1j * 45      # Line-to-line Impedance of Delta-load
Zload2_y = 30 + 1j * 40          # Line-to-neutral Impedance of Y-load
Zload1_y = Zload1_delta/3.0     # Converting the Delta-connected load into Y-eq
print("Z_delta as Y: ", Zload1_y)
# Finding the total phase impedance
Ztot_phase = Zline + (Zload1_y * Zload2_y) / (Zload1_y + Zload2_y)
Vsup_phase = Vsup / np.sqrt(3)
I_phase = Vsup_phase / Ztot_phase
print("Phase current: ", I_phase)
Stot = 3 * Vsup_phase * np.conj(I_phase)
print("Total real power from supply: ", np.real(Stot))
print("Total reactive power from supply: ", np.imag(Stot))

# b)
Vload_phase = Vsup_phase - Zline * I_phase  # V2 = V1 - Z*I
Vload_line = Vload_phase * np.sqrt(3) * polarToCart(
    1, 30, deg=True)  # Convert phase to line
print("Voltage at load, Line-to-Line: ", Vload_line,
      cartToPolar(Vload_line, deg=True))

# c)
I_load_1_y = Vload_phase / Zload1_y
I_load_2_y = Vload_phase / Zload2_y
I_load_1_delta = I_load_1_y / (np.sqrt(3) * polarToCart(1, -30, deg=True))
print("I1, Y: ", I_load_1_y, cartToPolar(I_load_1_y, deg=True))
print("I2, Y: ", I_load_2_y, cartToPolar(I_load_2_y, deg=True))
print("I1ab, delta: ", I_load_1_delta, cartToPolar(I_load_1_delta, deg=True))

# d)
S1 = 3*Vload_phase*I_load_1_y
S2 = 3*Vload_phase*I_load_2_y
Sline = 3*Zline*I_phase**2
print("S1: ", S1, cartToPolar(S1, deg=True))
print("S2: ", S2, cartToPolar(S2, deg=True))
print("Sline: ", Sline, cartToPolar(Sline, deg=True))

# ########## Example 2.8 ##############
# The first load absorbs 560.1 kVA at 0.707 power factor lagging.
# Current lagging after voltage. The second load absorbs 132 kW at unity
# power factor.
# Find a: The magnitude of the line voltage atthe source end of the line.
# Find b: Total real and reactive power loss in the line.
# Find c: Real power and reactive power supplied at the sending end of the line

Z_line_phase = 0.4 + 1j * 2.7   # Phase impedance
Vload_line = 3810.5             # Voltage magnitude at load
Vload_phase = Vload_line/(np.sqrt(3))  # Voltage pr phase at load
theta = np.arccos(0.707)
S1_mag = 560.1*1000
S1 = S1_mag * np.cos(theta) + 1j * S1_mag * np.sin(theta)
S2 = 132*1000
Stot = S1+S2
print(Stot)
I_phase = np.conj(Stot) / (3 * np.conj(Vload_phase))  # S = 3*V*I'
print("I_phase: ", I_phase, cartToPolar(I_phase, deg=True))
Vsup_phase = Vload_phase + I_phase * Z_line_phase
Vsup_line = Vsup_phase * np.sqrt(3)
print("Voltage line magnitude at supply: ", np.abs(Vsup_line))

# b)
Sline = 3 * I_phase**2 * Z_line_phase
print("Real and reactive power loss in the line: ", Sline)

# c)
Ssup = 3 * Vsup_phase * np.conj(I_phase)
print("Real and reactive power supplied at the sending end of the line [kVA]: ",
      Ssup/1000)


# ######## Problems Chap 2, 2.1 #############
# Modify Example 2.1 to take in new values for Vm, thetav, Zmag, gama
# Test an inductive, capacitive and resistive load. Compare.
# Find b: Calculate real and ractive power for each load and compare the
#         results obatined from the curves.
# Find c: Find the sum S of all the loads.

Vm = 100        # Peak voltage Vm,
thetav = 0      # Voltage phase angle,
Zmag = [
    polarToCart(1.25, 60, deg=True),
    polarToCart(2.0, -30, deg=True),
    polarToCart(2.5, 0, deg=True)
]
gama = [60, -30, 0]
for i in np.arange(len(Zmag)):
    print(Zmag[i])
    fig = runExample21(Vm, thetav, np.abs(Zmag[i]), gama[i])
    #plt.show()

# b)
S1 = Vm**2/np.conj(Zmag[0])/2       # Divide by to because = sqrt(2)**2
S2 = Vm**2/np.conj(Zmag[1])/2
S3 = Vm**2/np.conj(Zmag[2])/2
Stot = S1+S2+S3
print("S1: ", S1)
print("S2: ", S2)
print("S3: ", S3)
print("Stot: ", Stot)

# ######### Problem 2.2 ##############
# Find a: the complex power supplied to the load
# Find b: the instantaneous current i(t) and the rms value of the current
#          supplied to the load.
# Find c: The load impedance
# Plot d: v(t), p(t), i(t)
print("Problem 2.2")
v_mag = 200
w = 377

# instant power p(t) = |V||I|*cos(theta)+|V||I|*cos(theta)cos(2(wt+thetav))
thetav = -36.87
# Vm = |V|*sqrt(2)
P = 800
Px = 1000
# 800 = Vm*Im*cos(theta)/(sqrt(2)*sqrt(2))
I_mag = 2*800/(v_mag * np.cos(thetav*np.pi/180))
print("I", I_mag, cartToPolar(I_mag, deg=True))
print("Irms:", I_mag / np.sqrt(2), cartToPolar(I_mag / np.sqrt(2), deg=True))
V = v_mag / np.sqrt(2)
Q = np.sqrt(Px**2 - P**2)
S = P + 1j * Q
print("S: ", S)
P = I_mag * v_mag * np.cos(thetav * np.pi / 180) / 2
Q = I_mag * v_mag * np.sin(thetav * np.pi / 180) / 2
print("S: ", S)

# c)
Z = (V**2) / np.conj(S)
print("Z: ", np.round(Z, 4), np.round(cartToPolar(Z, deg=True), 4))
t = np.arange(0, 16.67+0.1, 0.1)/1000
v_t = v_mag * np.cos(377 * t)
p_t = 800 + 1000 * np.cos((754 * t - 36.87*np.pi/180))
i_t = p_t / v_t

# Figure
fig, axes = plt.subplots(3, 1)
axes[0].plot(t, v_t)
axes[0].set_title("v(t)")
axes[0].set_xlabel("t ")
axes[1].plot(t, p_t)
axes[1].set_title("p(t)")
axes[1].set_xlabel("t")
axes[2].plot(t, i_t)
axes[2].set_title("i(t)=p(t)/v(t)")
axes[2].set_xlabel("t")
axes[2].set_ylim(-15, 15)
plt.show()
